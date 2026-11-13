"""Free FeRh cluster with spins (MTP + spin/exchange/spline, fix nve/spin).

Exact identity: dG/dt = J_hf + J_prec, G = sum_i x_i e_i,
J_hf = compute heat/flux (ke, pe, centroid/stress/atom), J_prec = compute spin/heat/current.
Per-type energy offsets are removed as in analyse_cluster.py; the residual is
compared with the non-convective part (J_hf - J_conv + J_prec).
"""
import numpy as np

def load(tag):
    d = np.loadtxt(f"clusterspin_{tag}.dat")
    return dict(G=d[:, 1:4], Jc=d[:, 4:7], Js=d[:, 7:10], Jconv=d[:, 10:13],
                Jp=d[:, 13:16], X=d[:, 16:25].reshape(-1, 3, 3), V=d[:, 25:34].reshape(-1, 3, 3))

def analyse(tag, dt):
    r = load(tag)
    A = r["V"].transpose(0, 2, 1).reshape(-1, 3)
    c, *_ = np.linalg.lstsq(A, r["Jconv"].reshape(-1), rcond=None)
    Gp = r["G"] - np.einsum("t,stx->sx", c, r["X"])
    Jconvp = r["Jconv"] - np.einsum("t,stx->sx", c, r["V"])
    s = slice(1, -1)
    ref = (Gp[2:] - Gp[:-2]) / (2 * dt) - Jconvp[s]
    nrm = lambda a: np.sqrt((a**2).sum(1)).mean()
    scale = nrm(ref)
    cases = {"centroid + J_prec": r["Jc"] - r["Jconv"] + r["Jp"],
             "centroid, no J_prec": r["Jc"] - r["Jconv"],
             "stress/atom + J_prec": r["Js"] - r["Jconv"] + r["Jp"]}
    print(f"{tag}: mean|non-conv. part| = {scale:.3e}, mean|J_prec| = {nrm(r['Jp']):.3e} eV*A/ps")
    for k, J in cases.items():
        print(f"    {k:22s} rel. error = {nrm(ref - J[s]) / scale:.2e}")
    return r

a = analyse("mid_dt1_np1", 0.0001)
analyse("mid_dt05_np1", 0.00005)
analyse("seq_dt1_np1", 0.0001)
b = analyse("mid_dt1_np4", 0.0001)
for k in ("Jc", "Jp", "G"):
    print(f"np1 vs np4, {k}: max rel. diff = {np.abs(a[k] - b[k]).max() / np.abs(a[k]).max():.1e}")
