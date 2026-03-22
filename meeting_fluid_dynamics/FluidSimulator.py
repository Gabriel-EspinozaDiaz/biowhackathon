import numpy as np
from scipy.ndimage import map_coordinates
from .Particle import Particle


class FluidSimulator:
    """
    Incompressible 2-D Navier-Stokes solver.

    Grid layout
    -----------
    (vel_y * vel_x) collocated nodes.  Node (j, i) sits at:
        x = -box_size/2 + i * dx
        y = -box_size/2 + j * dy

    Boundary conditions : no-slip (u = v = 0) on all four walls.
    External forcing    : one or more vacuum (suction) point sources.
    """

    def __init__(self, box_size: float, vel_x: int, vel_y: int,
                 viscosity: float = 0.1):
        self.box_size  = box_size
        self.vel_x     = vel_x
        self.vel_y     = vel_y
        self.viscosity = viscosity
        self.dx = box_size / (vel_x - 1)
        self.dy = box_size / (vel_y - 1)

        # Velocity and pressure fields  –  shape: (vel_y, vel_x)
        self.u = np.zeros((vel_y, vel_x))   # x-velocity component
        self.v = np.zeros((vel_y, vel_x))   # y-velocity component
        self.p = np.zeros((vel_y, vel_x))   # pressure

        self.particles: list[Particle] = []
        self.vacuum_sources: list[dict] = []

        # Wave parameters - list to support multiple waves
        self.time = 0.0
        self.waves = []  # list of dicts: {"amplitude", "frequency", "direction", "speed"}

    # ── Public API ────────────────────────────────────────────────────────────

    def add_particle(self, particle: Particle) -> None:
        """Register a Particle so it is updated on every step()."""
        self.particles.append(particle)

    def add_vacuum(self, x: float, y: float,
                   strength: float = 5.0) -> None:
        """
        Place a vacuum (suction) source at physical position (x, y)
        measured from the box centre (cm).  strength > 0 creates inward pull.
        """
        self.vacuum_sources.append(
            {"x": float(x), "y": float(y), "strength": float(strength)}
        )

    def step(self, dt: float = 0.04) -> None:
        """Advance the simulation by one time step *dt* (seconds)."""
        self._apply_body_forces(dt)
        self._apply_wave_forces(dt)
        self.u, self.v = self._advect(dt)
        self._diffuse(dt)
        self._project(iterations=40)
        self._enforce_no_slip()
        for particle in self.particles:
            particle.update(
                self.u, self.v, dt,
                self.box_size, self.vel_x, self.vel_y
            )
        self.time += dt

    # ── Private helpers ───────────────────────────────────────────────────────

    def grid_coords(self) -> tuple[np.ndarray, np.ndarray]:
        """Return 2-D arrays (X, Y) of physical coordinates for every node."""
        half = self.box_size / 2.0
        xs   = np.linspace(-half, half, self.vel_x)
        ys   = np.linspace(-half, half, self.vel_y)
        return np.meshgrid(xs, ys)          # each has shape (vel_y, vel_x)

    def _apply_body_forces(self, dt: float) -> None:
        """Add vacuum-induced body forces to u and v."""
        if not self.vacuum_sources:
            return
        X, Y = self.grid_coords()
        for vac in self.vacuum_sources:
            # Vector from each grid node toward the vacuum source
            rx = vac["x"] - X
            ry = vac["y"] - Y
            r2 = np.maximum(rx**2 + ry**2, (0.5 * self.dx) ** 2)
            r  = np.sqrt(r2)
            # Inverse-square-law body force directed toward the source
            force = vac["strength"] / r2
            self.u += force * (rx / r) * dt
            self.v += force * (ry / r) * dt

    def _apply_wave_forces(self, dt: float) -> None:
        """Add wave-induced velocity perturbations."""
        if not self.waves:
            return
        X, Y = self.grid_coords()
        for wave in self.waves:
            freq = wave["frequency"]
            amp = wave["amplitude"]
            speed = wave["speed"]
            direction = wave["direction"]
            phase_offset = wave.get("phase", 0.0)
            if direction == "left":
                phase = 2 * np.pi * freq * (self.time - X / speed) + phase_offset
                self.u += amp * np.sin(phase)
            elif direction == "right":
                phase = 2 * np.pi * freq * (self.time + X / speed) + phase_offset
                self.u += amp * np.sin(phase)
            elif direction == "top":
                phase = 2 * np.pi * freq * (self.time - Y / speed) + phase_offset
                self.v += amp * np.sin(phase)
            elif direction == "bottom":
                phase = 2 * np.pi * freq * (self.time + Y / speed) + phase_offset
                self.v += amp * np.sin(phase)

    def _advect(self, dt: float) -> tuple[np.ndarray, np.ndarray]:
        """
        Semi-Lagrangian back-trace advection.
        For each grid node, trace a particle backward by dt and
        interpolate the velocity there.
        """
        ny, nx = self.vel_y, self.vel_x
        II, JJ = np.meshgrid(np.arange(nx), np.arange(ny))

        # Back-trace source positions in fractional grid indices
        src_i = np.clip(II - self.u * dt / self.dx, 0, nx - 1)
        src_j = np.clip(JJ - self.v * dt / self.dy, 0, ny - 1)
        coords = [src_j.ravel(), src_i.ravel()]

        u_new = map_coordinates(
            self.u, coords, order=1, mode="nearest"
        ).reshape(ny, nx)
        v_new = map_coordinates(
            self.v, coords, order=1, mode="nearest"
        ).reshape(ny, nx)
        return u_new, v_new

    def _diffuse(self, dt: float) -> None:
        """Explicit central-difference Laplacian diffusion (interior only)."""
        nu   = self.viscosity
        dx2  = self.dx ** 2
        dy2  = self.dy ** 2
        for field in (self.u, self.v):
            lap = np.zeros_like(field)
            lap[1:-1, 1:-1] = (
                (field[1:-1, 2:] - 2*field[1:-1, 1:-1] + field[1:-1, :-2]) / dx2
              + (field[2:, 1:-1] - 2*field[1:-1, 1:-1] + field[:-2, 1:-1]) / dy2
            )
            field += nu * dt * lap      # in-place (field is a view of u/v)

    def _project(self, iterations: int = 40) -> None:
        """
        Pressure-projection: solve ∇²p = ∇·u (Gauss-Seidel), then
        subtract ∇p from the velocity to enforce ∇·u = 0.
        """
        dx, dy = self.dx, self.dy

        # Velocity divergence at interior nodes
        div = np.zeros_like(self.p)
        div[1:-1, 1:-1] = (
            (self.u[1:-1, 2:] - self.u[1:-1, :-2]) / (2 * dx)
          + (self.v[2:, 1:-1] - self.v[:-2, 1:-1]) / (2 * dy)
        )

        # Iterative Poisson solve
        self.p[:] = 0.0
        coeff = 2.0 / dx**2 + 2.0 / dy**2
        for _ in range(iterations):
            p = self.p
            self.p[1:-1, 1:-1] = (
                (p[1:-1, 2:] + p[1:-1, :-2]) / dx**2
              + (p[2:, 1:-1] + p[:-2, 1:-1]) / dy**2
              - div[1:-1, 1:-1]
            ) / coeff
            # Neumann BCs: zero normal pressure gradient at walls
            self.p[ 0, :] = self.p[ 1, :]
            self.p[-1, :] = self.p[-2, :]
            self.p[:,  0] = self.p[:,  1]
            self.p[:, -1] = self.p[:, -2]

        # Correct velocity with pressure gradient
        self.u[1:-1, 1:-1] -= (self.p[1:-1, 2:] - self.p[1:-1, :-2]) / (2 * dx)
        self.v[1:-1, 1:-1] -= (self.p[2:, 1:-1] - self.p[:-2, 1:-1]) / (2 * dy)

    def _enforce_no_slip(self) -> None:
        """Set velocity to zero on all four boundary walls."""
        self.u[ 0, :] = 0;  self.u[-1, :] = 0
        self.u[:,  0] = 0;  self.u[:, -1] = 0
        self.v[ 0, :] = 0;  self.v[-1, :] = 0
        self.v[:,  0] = 0;  self.v[:, -1] = 0
    
    def wave_from_direction(self, amplitude: float, frequency: float, direction: str, speed: float = 1.0, phase: float = 0.0):
        """Add a soundwave coming from the specified direction. Can be called multiple times for multiple waves."""
        self.waves.append({
            "amplitude": amplitude,
            "frequency": frequency,
            "direction": direction,
            "speed": speed,
            "phase": phase
        })

    def clear_waves(self) -> None:
        """Remove all active waves."""
        self.waves = []