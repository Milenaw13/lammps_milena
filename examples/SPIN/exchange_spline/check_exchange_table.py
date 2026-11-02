#!/usr/bin/env python3
"""Compare a LAMMPS dimer scan (in.dimer_scan) with the table it used.

usage: check_exchange_table.py TABLEFILE KEYWORD [dimer_KEYWORD.dat]

Checks (parallel spins, offset no):
  1. E_LAMMPS(r) = -J(r): J from the same clamped cubic spline LAMMPS builds
     (slope at r[0] from the first 3 points, dJ/dr = 0 at the cutoff)
  2. Fx on atom 2 = -dE/dr = dJ/dr, compared with a finite difference of E
Optional: pass a Python reference with --ref module:function (J in eV).
"""
import sys
import numpy as np
from scipy.interpolate import CubicSpline


def read_section(fname, key):
    lines = [l.split("#")[0].split() for l in open(fname)]
    lines = [l for l in lines if l]
    i = next(k for k, l in enumerate(lines) if l[0] == key)
    n = int(lines[i + 1][1])
    data = np.array([[float(x) for x in l[1:3]] for l in lines[i + 2:i + 2 + n]])
    return data[:, 0], data[:, 1]


def end_slope(x, y):
    return np.polyfit(x - x[0], y, 2)[1]


def main():
    fname, key = sys.argv[1], sys.argv[2]
    dimer = sys.argv[3] if len(sys.argv) > 3 else f"dimer_{key}.dat"
    r, J = read_section(fname, key)
    spl = CubicSpline(r, J, bc_type=((1, end_slope(r[:3], J[:3])), (1, 0.0)))

    d = np.loadtxt(dimer)
    rs, e, fx = d[:, 0], d[:, 1], d[:, 2]
    inside = (rs >= r[0]) & (rs <= r[-1])
    jmax = np.abs(J).max()
    djmax = np.abs(spl(r, 1)).max()

    de = np.abs(-e[inside] - spl(rs[inside])).max()
    df = np.abs(fx[inside] - spl(rs[inside], 1)).max()
    print(f"table {key}: {len(r)} points, r = {r[0]} .. {r[-1]}")
    print(f"  max |E_lammps + J_spline|      = {de:.3e} eV    ({de / jmax:.1e} of max|J|)")
    print(f"  max |F_lammps - dJ/dr_spline|  = {df:.3e} eV/A  ({df / djmax:.1e} of max|dJ/dr|)")

    # energy/force consistency: central difference of E on the scan grid
    fd = -(e[2:] - e[:-2]) / (rs[2:] - rs[:-2])
    ok = inside[1:-1] & inside[2:] & inside[:-2]
    dfd = np.abs(fd[ok] - fx[1:-1][ok]).max()
    print(f"  max |F_lammps - (-dE/dr)_fd|   = {dfd:.3e} eV/A  (finite difference, O(dr^2))")
    return spl


if __name__ == "__main__":
    main()
