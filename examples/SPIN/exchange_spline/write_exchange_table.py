#!/usr/bin/env python3
"""Write tables for pair_style spin/exchange/spline.

Use as a module:

    from write_exchange_table import sample, write_table
    r, J = sample(my_J, rmin=2.0, rc=7.0, n=1001, scale=2.0e-3)
    write_table("J.table", [("FE_FE", r, J, "Fe-Fe, fit of TB2J data")])

my_J is any function of r (in Angstrom), e.g. the eval() method of an lmfit
ModelResult (lambda r: result.eval(x=r)).  scale converts the values to the
convention of LAMMPS (eV, each pair counted once): for TB2J data in meV use
scale = 2.0e-3 (TB2J sums over i != j, i.e. counts each pair twice).

Run as a script, it writes the table used by in.spin.iron_spline: the
Bethe-Slater J(r) of bcc iron (examples/SPIN/iron), shifted so that
J(rc) = 0 and dJ/dr(rc) = 0.
"""

import numpy as np


def sample(func, rmin, rc, n, scale=1.0):
    """Sample scale*func(r) on n equally spaced points from rmin to rc."""
    r = np.linspace(rmin, rc, n)
    return r, scale * np.asarray(func(r), dtype=float)


def check_cutoff(r, J, key, tol_j=1.0e-4, tol_dj=1.0e-3):
    """Same checks LAMMPS does: J and dJ/dr must vanish at the last point."""
    jmax = np.abs(J).max()
    djmax = np.abs(np.diff(J) / np.diff(r)).max()
    dj_end = np.polyfit(r[-3:] - r[-1], J[-3:], 2)[1]
    ok = True
    if abs(J[-1]) > tol_j * jmax:
        print(f"{key}: J(rc) = {J[-1]:.3e} is not zero (LAMMPS will stop)")
        ok = False
    if abs(dj_end) > tol_dj * djmax:
        print(f"{key}: dJ/dr(rc) = {dj_end:.3e} is not zero (LAMMPS will stop)")
        ok = False
    return ok


def write_table(fname, sections, header="spin/exchange/spline table"):
    """sections: list of (keyword, r, J, comment); r in Angstrom, J in eV."""
    with open(fname, "w") as f:
        f.write(f"# {header}\n# r in Angstrom, J in eV, H = -sum_(i<j) J(r_ij) s_i.s_j\n")
        for key, r, J, comment in sections:
            check_cutoff(r, J, key)
            f.write(f"\n# {comment}\n{key}\nN {len(r)}\n\n")
            for i, (ri, ji) in enumerate(zip(r, J)):
                f.write(f"{i + 1} {ri:.10f} {ji:.16e}\n")


def bethe_slater(r, J1=0.02726, J2=0.2171, J3=1.841):
    """J(r) of pair_style spin/exchange: 4 J1 ra (1 - J2 ra) exp(-ra), ra = (r/J3)^2."""
    ra = (r / J3) ** 2
    return 4.0 * J1 * ra * (1.0 - J2 * ra) * np.exp(-ra)


def bethe_slater_deriv(r, J1=0.02726, J2=0.2171, J3=1.841):
    ra = (r / J3) ** 2
    return 8.0 * J1 * r / J3**2 * (1.0 - ra - J2 * ra * (2.0 - ra)) * np.exp(-ra)


if __name__ == "__main__":
    rc = 5.0

    def shifted(r):
        return bethe_slater(r) - bethe_slater(rc) - (r - rc) * bethe_slater_deriv(rc)

    r, J = sample(shifted, rmin=1.5, rc=rc, n=701)
    write_table("Fe_bethe_slater.table",
                [("FE_BS", r, J, "Bethe-Slater J(r) of bcc Fe, shifted to J=dJ/dr=0 at rc")])
    print("wrote Fe_bethe_slater.table")
