#!/usr/bin/env python3
"""Green-Kubo magnon thermal conductivity, split into exchange channels.

Input: the file of

    compute  qs all spin/heat/flux groupA groupB
    fix      QC all ave/correlate Nevery Nrepeat Nfreq c_qs[1] ... c_qs[M] &
                  type full ave running file Qspin_channels.acf

with all M = N*(N+1)/2 channels of compute spin/heat/flux (N atom types, channel
order (1,1) (1,2) ... (1,N) (2,2) ... (N,N)).  The last block of the file is used.

With q_k = Q_k / A (A = interface area) and C_kl(t) = <q_k(t) q_l(0)>:

    kappa_kl = V / (kB T^2) * integral_0^tmax C_kl(t) dt      (trapezoidal rule)
    kappa    = sum_kl kappa_kl         (= kappa from the total flux c_qs)
    share_k  = sum_l (kappa_kl + kappa_lk) / 2      (sum_k share_k = kappa)

share_k is the contribution of channel k including its cross terms with all other
channels; it can be negative.  The diagonal kappa_kk is channel k on its own.
Units: LAMMPS metal units (eV/ps, Angstrom, K) in, W/(m K) out.
"""

import argparse

import numpy as np

KB = 8.617333262e-5           # Boltzmann constant (eV/K)
CONVERT = 1602.176634         # eV/(ps Angstrom K) -> W/(m K)


def read_last_block(fname):
    """Return the rows (index, time delta, count, values...) of the last block."""
    with open(fname) as f:
        lines = [l.split() for l in f if l.strip() and not l.startswith("#")]
    blocks = []
    i = 0
    while i < len(lines):
        nrows = int(lines[i][1])
        blocks.append(np.array(lines[i + 1:i + 1 + nrows], dtype=float))
        i += 1 + nrows
    if not blocks:
        raise SystemExit(f"no data in {fname}")
    return blocks[-1]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("file", help="output file of fix ave/correlate (type full)")
    p.add_argument("--timestep", type=float, required=True, help="MD time step (ps)")
    p.add_argument("--volume", type=float, required=True, help="volume (Angstrom^3)")
    p.add_argument("--area", type=float, required=True, help="interface area (Angstrom^2)")
    p.add_argument("--temperature", type=float, required=True, help="temperature (K)")
    p.add_argument("--names", help="comma separated names of the atom types, e.g. Fe,Rh,Ur")
    a = p.parse_args()

    data = read_last_block(a.file)
    t = data[:, 1] * a.timestep                 # time delta in ps
    corr = data[:, 3:]
    m = int(round(np.sqrt(corr.shape[1])))
    if m * m != corr.shape[1]:
        raise SystemExit("expected M*M correlation columns (fix ave/correlate type full)")
    ntypes = int(round((np.sqrt(8 * m + 1) - 1) / 2))
    if ntypes * (ntypes + 1) // 2 != m:
        raise SystemExit(f"{m} channels is not N*(N+1)/2 for any number of types N")

    names = a.names.split(",") if a.names else [str(i + 1) for i in range(ntypes)]
    if len(names) != ntypes:
        raise SystemExit(f"--names needs {ntypes} names")
    labels = [f"{names[i]}-{names[j]}" for i in range(ntypes) for j in range(i, ntypes)]

    scale = CONVERT * a.volume / (KB * a.temperature**2 * a.area**2)
    kmat = scale * np.trapezoid(corr, t, axis=0).reshape(m, m)
    kappa = kmat.sum()
    share = 0.5 * (kmat.sum(axis=1) + kmat.sum(axis=0))

    print(f"correlation time {t[-1]:.4g} ps, {m} channels")
    print(f"total magnon thermal conductivity: {kappa:.6g} W/mK\n")
    print(f"{'channel':>12} {'share':>12} {'fraction':>9} {'alone':>12}")
    for k in np.argsort(-np.abs(share)):
        frac = share[k] / kappa if kappa != 0.0 else float("nan")
        print(f"{labels[k]:>12} {share[k]:12.5g} {frac:9.3f} {kmat[k, k]:12.5g}")
    print(f"{'sum':>12} {share.sum():12.5g}\n")
    print("kappa_kl matrix (W/mK), rows k, columns l:")
    print(" " * 12 + "".join(f"{l:>12}" for l in labels))
    for k in range(m):
        print(f"{labels[k]:>12}" + "".join(f"{v:12.4g}" for v in kmat[k]))


if __name__ == "__main__":
    main()
