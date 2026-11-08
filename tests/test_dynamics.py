"""Python RK4 vs the Fortran RK4 in fortran/double_pendulum.f90."""

import pathlib
import shutil
import subprocess

import numpy as np

from double_pendulum import rk4

ROOT = pathlib.Path(__file__).resolve().parents[1]


def _fortran_output(tmp_path):
    """Run the Fortran program if gfortran works, otherwise use the stored output."""
    if shutil.which("gfortran"):
        exe = tmp_path / "dp"
        build = subprocess.run(["gfortran", "-O2", str(ROOT / "fortran" / "double_pendulum.f90"), "-o", str(exe)],
                               capture_output=True, text=True)
        if build.returncode == 0:
            subprocess.run([str(exe)], cwd=tmp_path, check=True, capture_output=True)
            return tmp_path / "double_pendulum_output.csv"
    return ROOT / "fortran" / "reference_output.csv"


def test_python_matches_fortran(tmp_path):
    """Same scheme, same step, same initial state: the angles agree to the
    10 significant digits written by the Fortran program for the first 15 s."""
    data = np.genfromtxt(_fortran_output(tmp_path), delimiter=",", skip_header=1)
    t, th1, th2 = data[:, 0], data[:, 1], data[:, 2]
    tr = rk4([np.radians(120), 0, np.radians(-10), 0], 0.01, len(t) - 1)
    mask = t <= 15.0
    np.testing.assert_allclose(tr.y[0, mask], th1[mask], atol=1e-9)
    np.testing.assert_allclose(tr.y[2, mask], th2[mask], atol=1e-9)
