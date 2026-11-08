"""Equations of motion and kinematics of the planar double pendulum.

State vector y = [theta1, omega1, theta2, omega2], angles measured from the
downward vertical. The equations follow from the Lagrangian
L = T - V of two point masses on massless rigid rods.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Pendulum:
    """Physical parameters (SI units)."""

    m1: float = 1.0
    m2: float = 1.0
    L1: float = 1.0
    L2: float = 1.0
    g: float = 9.81


def rhs(t, y, p=Pendulum()):
    """Time derivative of the state."""
    th1, w1, th2, w2 = y
    m1, m2, L1, L2, g = p.m1, p.m2, p.L1, p.L2, p.g

    d = th2 - th1
    den1 = (m1 + m2) * L1 - m2 * L1 * np.cos(d) ** 2
    den2 = (L2 / L1) * den1

    dw1 = (
        m2 * L1 * w1**2 * np.sin(d) * np.cos(d)
        + m2 * g * np.sin(th2) * np.cos(d)
        + m2 * L2 * w2**2 * np.sin(d)
        - (m1 + m2) * g * np.sin(th1)
    ) / den1
    dw2 = (
        -m2 * L2 * w2**2 * np.sin(d) * np.cos(d)
        + (m1 + m2) * (g * np.sin(th1) * np.cos(d) - L1 * w1**2 * np.sin(d) - g * np.sin(th2))
    ) / den2

    return np.array([w1, dw1, w2, dw2])


def positions(y, p=Pendulum()):
    """Cartesian positions (x1, y1, x2, y2) of the two masses [m]."""
    th1, th2 = y[0], y[2]
    x1 = p.L1 * np.sin(th1)
    y1 = -p.L1 * np.cos(th1)
    x2 = x1 + p.L2 * np.sin(th2)
    y2 = y1 - p.L2 * np.cos(th2)
    return x1, y1, x2, y2
