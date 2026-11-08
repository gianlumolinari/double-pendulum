"""Integration with SciPy's solve_ivp and with a fixed-step RK4."""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .dynamics import Pendulum, rhs


@dataclass
class Trajectory:
    t: np.ndarray  # (n,)
    y: np.ndarray  # (4, n)


def simulate(y0, t_final, p=Pendulum(), n_out=2001, rtol=1e-10, atol=1e-10):
    """Adaptive integration (DOP853) sampled at n_out equally spaced times."""
    t_eval = np.linspace(0.0, t_final, n_out)
    sol = solve_ivp(rhs, (0.0, t_final), np.asarray(y0, float), method="DOP853",
                    t_eval=t_eval, args=(p,), rtol=rtol, atol=atol)
    return Trajectory(sol.t, sol.y)


def rk4(y0, dt, n_steps, p=Pendulum()):
    """Classical fixed-step RK4, the same scheme as fortran/double_pendulum.f90."""
    y = np.array(y0, dtype=float)
    ys = [y.copy()]
    for _ in range(n_steps):
        k1 = rhs(0.0, y, p)
        k2 = rhs(0.0, y + 0.5 * dt * k1, p)
        k3 = rhs(0.0, y + 0.5 * dt * k2, p)
        k4 = rhs(0.0, y + dt * k3, p)
        y = y + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        ys.append(y.copy())
    return Trajectory(dt * np.arange(n_steps + 1), np.stack(ys, axis=-1))
