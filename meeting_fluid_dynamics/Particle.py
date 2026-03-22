import numpy as np
from scipy.ndimage import map_coordinates


class Particle:

    def __init__(self, x: float, y: float, diameter: float,
                 color: str = "red"):
        self.x        = float(x)
        self.y        = float(y)
        self.diameter = float(diameter)
        self.radius   = diameter / 2.0
        self.color    = color
        self.trail: list[tuple[float, float]] = [(self.x, self.y)]

    # ------------------------------------------------------------------
    def update(self, u: np.ndarray, v: np.ndarray,
               dt: float, box_size: float,
               vel_x: int, vel_y: int) -> None:

        dx   = box_size / (vel_x - 1)
        dy   = box_size / (vel_y - 1)
        half = box_size / 2.0

        # Physical position → fractional grid index  (col=i, row=j)
        gi = (self.x + half) / dx
        gj = (self.y + half) / dy
        gi = float(np.clip(gi, 0, vel_x - 1))
        gj = float(np.clip(gj, 0, vel_y - 1))

        # Bilinear interpolation; map_coordinates uses [row, col] order
        vx = map_coordinates(u, [[gj], [gi]], order=1, mode="nearest").item()
        vy = map_coordinates(v, [[gj], [gi]], order=1, mode="nearest").item()

        # Forward-Euler position update, clamped to box interior
        wall = half - self.radius
        self.x = float(np.clip(self.x + vx * dt, -wall, wall))
        self.y = float(np.clip(self.y + vy * dt, -wall, wall))

        self.trail.append((self.x, self.y))
        if len(self.trail) > 200:       # keep at most 200 history points
            self.trail.pop(0)