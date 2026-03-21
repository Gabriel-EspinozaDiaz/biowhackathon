#!/usr/bin/env python3
"""
Run the 2D Navier-Stokes fluid simulation.

Adjust the parameters below, add vacuum sources and particles, then run:
    python run.py

Grid coordinates are centred on the box, so (0, 0) is the middle.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.animation import FuncAnimation

from FluidSimulator import FluidSimulator
from Particle import Particle

# ─── Simulation parameters ────────────────────────────────────────────────────
box_size  = 10.0    # Side-length of the square box (cm)
vel_x     = 40      # Number of velocity grid points along x
vel_y     = 40      # Number of velocity grid points along y
viscosity = 0.15    # Kinematic viscosity  (cm² / s)
dt        = 0.04    # Time-step  (s)
# ──────────────────────────────────────────────────────────────────────────────


def run_simulation():
    """
    Set up the scene and launch the animated matplotlib window.

    Add vacuum sources with sim.add_vacuum(x, y, strength) before calling
    plt.show() to introduce external forcing.
    """
    sim = FluidSimulator(box_size, vel_x, vel_y, viscosity=viscosity)

    # ── Particles  (x, y from centre, diameter in cm) ─────────────────────────
    particle_specs = [
        (-4.0, -4.0, 0.40, "#e74c3c"),
        ( 4.0, -4.0, 0.40, "#3498db"),
        (-4.0,  4.0, 0.40, "#2ecc71"),
        ( 4.0,  4.0, 0.40, "#f39c12"),
        ( 0.0, -4.5, 0.35, "#9b59b6"),
        ( 0.0,  4.5, 0.35, "#1abc9c"),
        (-4.5,  0.0, 0.35, "#e67e22"),
        ( 4.5,  0.0, 0.35, "#e91e63"),
        (-2.0,  3.0, 0.30, "#00bcd4"),
        ( 2.0, -3.0, 0.30, "#cddc39"),
    ]
    for px, py, diam, col in particle_specs:
        sim.add_particle(Particle(px, py, diameter=diam, color=col))

    # ── Figure / axes ─────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect("equal")
    half = box_size / 2.0
    ax.set_xlim(-half, half)
    ax.set_ylim(-half, half)
    ax.set_facecolor("#ffffff")
    fig.patch.set_facecolor("#9c9c9c")
    ax.set_xlabel("x  (cm)", color="grey")
    ax.set_ylabel("y  (cm)", color="grey")
    ax.tick_params(colors="grey")
    for spine in ax.spines.values():
        spine.set_edgecolor("#334155")

    title = ("Navier–Stokes Fluid Simulator  |  "
             f"{vel_x}×{vel_y} grid  |  ν = {viscosity} cm²/s")
    ax.set_title(title, color="white", fontsize=10, pad=10)

    # Box boundary
    box_rect = mpatches.FancyBboxPatch(
        (-half, -half), box_size, box_size,
        boxstyle="square,pad=0",
        linewidth=1.5, edgecolor="#475569", facecolor="none", zorder=5
    )
    ax.add_patch(box_rect)

    # ── Velocity-magnitude heat-map ────────────────────────────────────────────
    speed = np.sqrt(sim.u**2 + sim.v**2)
    im = ax.imshow(
        speed,
        extent=[-half, half, -half, half],
        origin="lower", cmap="inferno",
        alpha=0.55, vmin=0, vmax=4, zorder=1
    )
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("|velocity|  (cm/s)", color="grey", fontsize=8)
    cbar.ax.yaxis.set_tick_params(color="grey")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="grey", fontsize=7)

    # ── Quiver (velocity direction arrows) ────────────────────────────────────
    X, Y = sim.grid_coords()
    stride = max(1, min(vel_x, vel_y) // 14)    # subsample for readability
    Xq = X[::stride, ::stride]
    Yq = Y[::stride, ::stride]
    quiv = ax.quiver(
        Xq, Yq,
        sim.u[::stride, ::stride],
        sim.v[::stride, ::stride],
        color="#7dd3fc", alpha=0.65,
        scale=25, width=0.003, zorder=4
    )

    # ── Particle trails and circles ────────────────────────────────────────────
    trail_lines = [
        ax.plot([], [], "-", color=p.color, alpha=0.35, lw=1.0)[0]
        for p in sim.particles
    ]
    circles = [
        plt.Circle((p.x, p.y), p.radius, color=p.color, zorder=9)
        for p in sim.particles
    ]
    for c in circles:
        ax.add_patch(c)

    time_text = ax.text(
        0.02, 0.975, "t = 0.00 s",
        transform=ax.transAxes, color="white",
        fontsize=9, va="top", family="monospace"
    )

    # ── Animation update ───────────────────────────────────────────────────────
    frame = [0]

    def update(_):
        sim.step(dt)
        frame[0] += 1

        # Heat-map
        speed = np.sqrt(sim.u**2 + sim.v**2)
        im.set_data(speed)

        # Quiver
        quiv.set_UVC(
            sim.u[::stride, ::stride],
            sim.v[::stride, ::stride]
        )

        # Particles
        for p, circ, trail in zip(sim.particles, circles, trail_lines):
            circ.set_center((p.x, p.y))
            if len(p.trail) > 1:
                tx, ty = zip(*p.trail[-80:])
                trail.set_data(tx, ty)

        time_text.set_text(f"t = {frame[0] * dt:.2f} s")
        return [im, quiv] + circles + trail_lines + [time_text]

    ani = FuncAnimation(
        fig, update,
        interval=30,        # ms between frames
        blit=True,
        cache_frame_data=False
    )

    plt.tight_layout()
    plt.show()
    return sim, ani


if __name__ == "__main__":
    sim, ani = run_simulation()
