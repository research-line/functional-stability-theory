"""
c2cf_fixed_waterline_audit.py

F??hrt das Fixed-Waterline-Moat-Audit f??r Lemma C2cf.1 durch.
Vermessung der spektralen Masse von k au??erhalb des Fensters V_g = [u0, u0 + g0).

Berechnet f??r g0 in {0.25, 0.5, 1, 2, 5} und lambda in {3.0, 5.0, 7.0}:
  - dim V_g (Anzahl der Eigenwerte u_j - u_0 < g0)
  - ||(I-P_g)k|| (tats??chliche Masse au??erhalb des Fensters)
  - R_0/g_0 (theoretische Schranke mit R_0 = ||(A_lambda - u_0)k||)
  - Anzahl der Slush-Moden (Moden in V_g mit |c_j|^2 < 1e-6)
"""

import os, sys, time, json
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import mpmath as mp
from pathlib import Path

# Pfad zu den Hilfsskripten hinzuf??gen
sys.path.insert(0, str(Path(__file__).parent))
from dirichlet_ccm_fourier_mp import (
    build_QW_mp, project_to_parity, chi_trivial,
)
from c2_approximation_test import (
    k_lambda_value, norm,
)
from c2_poisson_decomposition import (
    project_to_fourier,
)

DPS = int(os.environ.get("DPS", 35))
N_FACTOR = int(os.environ.get("N_FACTOR", 12))

G0_VALUES = [0.25, 0.5, 1.0, 2.0, 5.0]
LAMBDAS = [3.0, 5.0, 7.0]

def to_np(M, dim):
    A = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            A[i, j] = float(M[i, j])
    return A

def main():
    print("=" * 90)
    print("FIXED-WATERLINE SPEKTRAL-AUDIT (Lemma C2cf.1)")
    print(f"DPS={DPS}, N_FACTOR={N_FACTOR}")
    print("=" * 90)

    results = []

    for lam in LAMBDAS:
        mp.mp.dps = DPS
        lam_mp = mp.mpf(lam)
        L_mp = 2 * mp.log(lam_mp)
        L = float(L_mp)
        
        # N-Skalierung
        N = max(40, int(np.ceil(N_FACTOR * L)))
        dim = N + 1
        
        print(f"\nBerechne Operator f??r lam={lam:.1f} (N={N}, dim={dim})...")
        t0 = time.time()
        Mq, _ = build_QW_mp(N, lam, chi_trivial, 1, True, -1, verbose=False)
        Mh, _ = build_QW_mp(N, lam, chi_trivial, 1, False, -1, verbose=False)
        Aq, _ = project_to_parity(Mq, N, parity="even")
        Ah, _ = project_to_parity(Mh, N, parity="even")
        
        A_np = to_np(Aq, dim)
        H_np = to_np(Ah, dim)
        
        # Diagonalisierung
        ws, vs = np.linalg.eigh(A_np)
        
        # Poisson-Vektor k
        def kfull(x):
            u = mp.exp(x) / lam_mp
            return k_lambda_value(float(u), lam)
        cf = project_to_fourier(kfull, L_mp, N)
        nf = norm(cf, dim)
        kn = np.array([float(cf[i, 0] / nf) for i in range(dim)])
        dt = time.time() - t0
        print(f"  Operator und Poisson-Vektor berechnet in {dt:.1f}s")
        
        u0 = ws[0]
        deltas = ws - u0
        overlaps = np.array([float(kn @ vs[:, j])**2 for j in range(dim)])
        
        # Residuum R_0 = ||(A - u0)k||
        res_w0 = A_np @ kn - u0 * kn
        R_0 = float(np.linalg.norm(res_w0))
        
        print(f"  R_0 = {R_0:.6e}")
        
        lam_results = {
            "lam": lam,
            "N": N,
            "R_0": R_0,
            "g0_audits": []
        }
        
        print(f"  {'g0':>6}  {'dim V_g':>8}  {'||(I-P_g)k||':>14}  {'R_0/g_0':>12}  {'Slush-Moden':>12}  {'Massen-Marge':>14}")
        for g0 in G0_VALUES:
            # Maske f??r Moden in V_g (u_j - u_0 < g0)
            in_window = deltas < g0
            dim_Vg = int(np.sum(in_window))
            
            # Masse in V_g
            mass_in_Vg = np.sum(overlaps[in_window])
            mass_outside = max(0.0, 1.0 - mass_in_Vg)
            norm_outside = np.sqrt(mass_outside)
            
            # R_0/g_0 Schranke
            bound = R_0 / g0
            
            # Slush-Moden: Moden in V_g mit |c_j|^2 < 1e-6
            slush_mask = in_window & (overlaps < 1e-6)
            num_slush = int(np.sum(slush_mask))
            
            # ??berpr??fung der Ungleichung
            valid = bool(norm_outside <= bound)
            marker = "OK" if valid else "FAIL!"
            
            print(f"  {g0:6.2f}  {dim_Vg:8d}  {norm_outside:14.6e}  {bound:12.6e}  {num_slush:12d}  {mass_in_Vg:13.6f} ({marker})")
            
            lam_results["g0_audits"].append({
                "g0": g0,
                "dim_Vg": dim_Vg,
                "norm_outside": norm_outside,
                "bound": bound,
                "num_slush": num_slush,
                "mass_in_Vg": mass_in_Vg,
                "valid": valid
            })
        
        results.append(lam_results)
        
    # Ergebnisse in JSON speichern
    output_file = Path(__file__).parent / "c2cf_fixed_waterline_results.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nErgebnisse gespeichert unter: {output_file.name}")

if __name__ == "__main__":
    main()


