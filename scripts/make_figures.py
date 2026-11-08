"""Regenerate the figures in figures/ (run from the repository root):

    python scripts/make_figures.py
"""

import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from double_pendulum import Pendulum, positions, rk4, simulate  # noqa: E402

FIG = ROOT / "figures"
BLUE, ORANGE, INK, MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
plt.rcParams.update({"figure.facecolor": "white", "axes.grid": True, "grid.color": "#e4e3df",
                     "axes.spines.top": False, "axes.spines.right": False, "font.size": 11})

P = Pendulum()
Y0 = np.array([np.radians(120), 0.0, np.radians(-10), 0.0])   # same as the Fortran program
LIM = P.L1 + P.L2 + 0.15


def _pendulum_axes(fig):
    ax = fig.add_subplot()
    ax.set(xlim=(-LIM, LIM), ylim=(-LIM, LIM), aspect="equal", xticks=[], yticks=[])
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.plot(0, 0, "o", color=INK, ms=4)
    return ax


def trajectory_plot():
    """Path of the second mass, coloured by time."""
    tr = simulate([np.radians(91), 0, np.radians(25), 0], 20.0, P, n_out=4000)
    _, _, x2, y2 = positions(tr.y, P)
    pts = np.array([x2, y2]).T.reshape(-1, 1, 2)
    lc = LineCollection(np.concatenate([pts[:-1], pts[1:]], axis=1), cmap="viridis", linewidth=1.5)
    lc.set_array(tr.t[:-1])
    fig, ax = plt.subplots(figsize=(6, 5.4), dpi=130)
    ax.add_collection(lc)
    ax.plot(0, 0, "o", color=INK, ms=5)
    ax.set(xlim=(-LIM, LIM), ylim=(-LIM, LIM), aspect="equal", xlabel="x [m]", ylabel="y [m]",
           title="Path of the second mass, 20 s")
    fig.colorbar(lc, ax=ax, label="time [s]")
    fig.tight_layout()
    fig.savefig(FIG / "trajectory.png")
    plt.close(fig)


def animate(trajs, colors, label, path, t_final=12.0, fps=25):
    n = int(t_final * fps) + 1
    xy = [positions(tr.y, P) for tr in trajs]
    fig = plt.figure(figsize=(4, 4), dpi=90)
    ax = _pendulum_axes(fig)
    trails = [ax.plot([], [], "-", color=c, lw=1, alpha=0.5)[0] for c in colors]
    arms = [ax.plot([], [], "o-", color=c, lw=2, ms=7, mec="white", mew=1.5)[0] for c in colors]
    txt = ax.text(0.02, 0.97, "", transform=ax.transAxes, va="top", color=MUTED, fontsize=10)
    ax.text(0.98, 0.97, label, transform=ax.transAxes, va="top", ha="right", color=MUTED, fontsize=10)

    def update(k):
        for (x1, y1, x2, y2), tr_line, arm in zip(xy, trails, arms):
            arm.set_data([0, x1[k], x2[k]], [0, y1[k], y2[k]])
            tr_line.set_data(x2[max(0, k - 40):k + 1], y2[max(0, k - 40):k + 1])
        txt.set_text(f"t = {trajs[0].t[k]:4.1f} s")
        return (*trails, *arms, txt)

    fig.tight_layout()
    FuncAnimation(fig, update, frames=n, blit=True).save(path, writer=PillowWriter(fps=fps))
    plt.close(fig)


def fortran_comparison():
    """Difference between the Python and Fortran RK4 runs (same scheme, dt = 0.01 s)."""
    data = np.genfromtxt(ROOT / "fortran" / "reference_output.csv", delimiter=",", skip_header=1)
    t, th1_f = data[:, 0], data[:, 1]
    tr = rk4(Y0, 0.01, len(t) - 1, P)
    diff = np.maximum(np.abs(tr.y[0] - th1_f), 1e-12)
    fig, ax = plt.subplots(figsize=(7, 3.6), dpi=130)
    ax.semilogy(t, diff, color=BLUE, lw=1.5)
    ax.axhline(5e-11, color=MUTED, ls="--", lw=1)
    ax.text(0.3, 8e-11, "precision of the Fortran output file", color=MUTED, fontsize=9)
    ax.set(xlabel="time [s]", ylabel="|θ₁ Python − θ₁ Fortran| [rad]", ylim=(1e-12, 1e-5),
           title="Python vs Fortran RK4: agreement until round-off is amplified")
    fig.tight_layout()
    fig.savefig(FIG / "python_vs_fortran.png")
    plt.close(fig)


def angle_plot():
    """theta1(t) for two pendulums released 0.001 deg apart."""
    t_f = 20.0
    a = simulate(Y0, t_f, P, n_out=4001)
    b = simulate(Y0 + [np.radians(0.001), 0, 0, 0], t_f, P, n_out=4001)
    fig, ax = plt.subplots(figsize=(8, 3.6), dpi=130)
    ax.plot(a.t, np.degrees(a.y[0]), color=BLUE, lw=1.5, label="θ₁(0) = 120°")
    ax.plot(b.t, np.degrees(b.y[0]), color=ORANGE, lw=1.5, ls="--", label="θ₁(0) = 120.001°")
    ax.set(xlabel="time [s]", ylabel="θ₁ [deg]", title="Two releases 0.001° apart")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "angles.png")
    plt.close(fig)


if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    trajectory_plot()
    angle_plot()
    fortran_comparison()
    one = simulate(Y0, 12.0, P, n_out=301)
    animate([one], [BLUE], "", FIG / "pendulum.gif")
    two = simulate(Y0 + [np.radians(0.001), 0, 0, 0], 12.0, P, n_out=301)
    animate([one, two], [BLUE, ORANGE], "Δθ₁(0) = 0.001°", FIG / "butterfly.gif")
    print("figures written to", FIG)
