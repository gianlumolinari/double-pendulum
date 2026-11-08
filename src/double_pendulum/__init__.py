"""Planar double pendulum: equations of motion and integration."""

from .dynamics import Pendulum, positions, rhs
from .simulate import Trajectory, rk4, simulate

__all__ = ["Pendulum", "rhs", "positions", "simulate", "rk4", "Trajectory"]
