"""
milestone_2_chi_v2_rho_wideband.py
==================================

Variante von milestone_2_chi_v2_rho.py mit erweiteter Paley-Wiener-Bandbreite
T_PW = 25, so dass alle L-Nullstellen (gamma <= 20) im Bandbegrenzungs-Bereich
der Prolate-Basis liegen.

**WICHTIG**: das bricht die strikte Paley-Wiener-Einschraenkung T < log(lambda)
aus Milestone 1_chi. Diese Variante ist DIAGNOSTISCH: sie zeigt, ob die
H1/H2-Signatur sichtbar wird, wenn die rho-Summe numerisch wirkt.

Bei lambda = sqrt(14), log(lambda) ~ 1.32, und gamma_1(chi_0) ~ 14.13 -- alle
Nullstellen liegen ausserhalb des strikten Paley-Wiener-Bereichs. Das erklaert
strukturell, warum die strikte Version (milestone_2_chi_v2_rho.py) die
Weil-Positivitaet nicht erreichen kann: die Nullstellen sind nicht im
Prolate-Spektrum.

Test mit T_PW = 25: enthalte gamma bis 20 explizit.

Autor: LG (Claude Opus 4.7 [1M], Session 16 Fortsetzung)
Datum: 2026-04-18
"""

import json
import numpy as np
import time
from pathlib import Path
from scipy.special import digamma

LAMBDA = np.sqrt(14.0)
L = np.log(LAMBDA)
T_WIDE = 25.0
N_GRID = 2400
T_PW = 25.0          # >> L = 1.32, bricht strikte PW-Einschraenkung aber diagnostisch
N_GALERKIN = 60
N_L_TERMS = 400
K_EIGENVALUES = 8
RHO_CUTOFF = 20.0

RESULTS_DIR = Path(__file__).parent.parent / "_results"

RIEMANN_ZEROS = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005150, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
    67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
]
CHI4_ZEROS = [
    6.020948, 10.243766, 12.988096, 16.343297, 18.291996,
    21.428273, 23.265376, 26.068044, 28.106108, 30.296575,
]

def build_prolate_basis(n_grid, t_wide, t_pw, n_galerkin):
    t = np.linspace(-t_wide, t_wide, n_grid)
    dt = t[1] - t[0]
    diff = t[:, None] - t[None, :]
    K = np.where(np.abs(diff) < 1e-14,
                 t_pw / np.pi,
                 np.sin(t_pw * diff) / (np.pi * diff))
    K *= dt
    K = 0.5 * (K + K.T)
    ev, evv = np.linalg.eigh(K)
    idx = np.argsort(ev)[::-1]
    H = evv[:, idx][:, :n_galerkin] / np.sqrt(dt)
    return {'t': t, 'dt': dt, 'H': H, 'lambdas': ev[idx][:n_galerkin]}

def project_diag(H, dt, diag):
    return (H.T * diag[None, :]) @ H * dt

def project_dense(H, dt, K):
    return H.T @ K @ H * dt**2

def interp_H(basis, gammas):
    t = basis['t']; H = basis['H']; dt = basis['dt']
    result = np.zeros((len(gammas), H.shape[1]))
    for i, g in enumerate(gammas):
        if g < t[0] or g > t[-1]: continue
        idx = np.searchsorted(t, g)
        if idx <= 0: result[i, :] = H[0, :]
        elif idx >= len(t): result[i, :] = H[-1, :]
        else:
            alpha = (g - t[idx-1]) / dt
            result[i, :] = (1-alpha)*H[idx-1,:] + alpha*H[idx,:]
    return result

def G_arch(basis, parity):
    t = basis['t']
    arg = 0.25 if parity == +1 else 0.75
    v = 2.0*np.real(digamma(arg + 1j*t/2.0))
    return project_diag(basis['H'], basis['dt'], v)

def G_prim(basis, chi, q):
    t = basis['t']; dt = basis['dt']; H = basis['H']; n = len(t)
    p_max = int(LAMBDA**2)+1
    primes = [p for p in range(2, p_max+1)
              if all(p%q_!=0 for q_ in range(2, int(np.sqrt(p))+1))]
    diff = t[:,None] - t[None,:]
    K = np.zeros((n,n), dtype=complex)
    for p in primes:
        if q>1 and p%q==0: continue
        lp=np.log(p); sp=np.sqrt(p)
        m_max = max(1, int(np.log(LAMBDA**2)/lp))
        for m in range(1, m_max+1):
            cpm = chi(p)**m
            if cpm == 0: continue
            w = -cpm*lp/(sp**m)
            K += w*np.cos(diff*m*lp)
    return project_dense(H, dt, K)

def G_rho(basis, zeros, cutoff):
    valid = [g for g in zeros if abs(g) <= cutoff]
    if not valid:
        N = basis['H'].shape[1]
        return np.zeros((N,N), dtype=complex)
    gammas = np.array(valid + [-g for g in valid])
    Hg = interp_H(basis, gammas)
    return (Hg.T @ Hg).astype(complex)

def make_chi(name, q, parity, vals, zeros):
    arr = np.zeros(q, dtype=complex)
    for n,v in vals.items(): arr[n]=v
    def chi(n): return arr[n%q]
    return {'name':name,'q':q,'parity':parity,'chi':chi,'zeros':zeros,'gamma1':zeros[0] if zeros else 0}

def load_chars():
    zf = Path(__file__).parent.parent.parent/"dirichlet_atlas"/"_results"/"zeros_all_chars.json"
    za = json.load(open(zf))
    c0 = make_chi('chi_0',1,+1,{},RIEMANN_ZEROS)
    def chi_t(n): return 1.0+0j
    c0['chi']=chi_t
    chars=[c0]
    chars.append(make_chi('chi_4',4,-1,{1:1.0+0j,3:-1.0+0j},CHI4_ZEROS))
    chars.append(make_chi('chi_5',5,+1,{1:1.0+0j,2:-1.0+0j,3:-1.0+0j,4:1.0+0j},za['chi_5']))
    chars.append(make_chi('chi_8',8,+1,{1:1.0+0j,3:-1.0+0j,5:-1.0+0j,7:1.0+0j},za['chi_8']))
    arr33=[0,1,1,0,1,-1,0,-1,1,0,-1,0,0,-1,-1,0,1,1,0,-1,-1,0,0,-1,0,1,-1,0,-1,1,0,1,1]
    chars.append(make_chi('chi_33',33,+1,
                          {n:float(arr33[n])+0j for n in range(33) if arr33[n]!=0},
                          za['chi_33']))
    return chars

def L_vals(basis, chi):
    t = basis['t']; s = 0.5 + 1j*t
    r = np.zeros_like(t, dtype=complex)
    for n in range(1, N_L_TERMS+1):
        cn=chi(n)
        if cn!=0: r += cn/(n**s)
    return r

def build_Psi(basis, chi): return project_diag(basis['H'], basis['dt'], L_vals(basis, chi))
def build_PW(basis):      return project_diag(basis['H'], basis['dt'], basis['t']**2+1.0)

def main():
    print("="*95)
    print(f"WIDEBAND Milestone 2_chi v2 (rho-Term + T_PW={T_PW} bricht Milestone-1-Schranke)")
    print(f"lambda={LAMBDA:.4f}, L={L:.4f}, T_PW={T_PW}, T_WIDE={T_WIDE}, N_GALERKIN={N_GALERKIN}")
    print("="*95)

    basis = build_prolate_basis(N_GRID, T_WIDE, T_PW, N_GALERKIN)
    print(f"Prolate-Eigenwerte: first5={basis['lambdas'][:5]}, last={basis['lambdas'][-1]:.3e}")

    PW = build_PW(basis)
    chars = load_chars()

    print("\n{:10s} {:>4s} {:>6s} {:>7s} {:>5s} {:>9s} {:>9s} {:>9s} {:>9s} {:>7s}".format(
        'chi','q','par','gamma1','n_rho','||G_arc||','||G_prim||','||G_rho||','minEW','pos?'))
    print("-"*95)

    eigs = []
    defs = []
    for c in chars:
        Ga = G_arch(basis, c['parity'])
        Gp = G_prim(basis, c['chi'], c['q'])
        Gr = G_rho(basis, c['zeros'], RHO_CUTOFF)
        QW = Gr - Ga - Gp
        QWh = 0.5*(QW + QW.conj().T)
        ew = np.linalg.eigvalsh(QWh)
        n_rho = sum(1 for g in c['zeros'] if abs(g)<=RHO_CUTOFF)
        na = float(np.linalg.norm(Ga,2))
        np_= float(np.linalg.norm(Gp,2))
        nr = float(np.linalg.norm(Gr,2))
        pos = "JA" if ew.min() >= -1e-6 else "nein"
        print(f"{c['name']:10s} {c['q']:4d} {c['parity']:+6d} {c['gamma1']:7.3f} {n_rho:5d} "
              f"{na:9.4f} {np_:9.4f} {nr:9.4f} {ew.min():+9.4f} {pos:>7s}")
        eigs.append({'chi':c['name'],'gamma1':c['gamma1'],'min_ew':float(ew.min()),
                     'n_neg':int(np.sum(ew<-1e-8)),'n_rho':n_rho,
                     'norm_G_arch':na,'norm_G_prim':np_,'norm_G_rho':nr})

        # Defekt
        Psi = build_Psi(basis, c['chi'])
        A = QW @ Psi; B = Psi @ PW
        tAB = np.trace(A.conj().T @ B); tBB = np.trace(B.conj().T @ B)
        mu = tAB/tBB if abs(tBB)>1e-14 else 1.0
        D = A - mu*B
        ds = float(np.linalg.norm(D,2))
        nB = float(np.linalg.norm(B,2))
        nPsi = float(np.linalg.norm(Psi,'fro'))
        defs.append({'chi':c['name'],'gamma1':c['gamma1'],
                     'defect_spec':ds,
                     'rel_a': ds/max(nB,1e-14),
                     'rel_b': ds/max(nPsi,1e-14),
                     'mu_re': float(np.real(mu))})

    print("\nDefekt-Norm:")
    print(f"{'chi':10s} {'gamma1':>7s} {'rel_a':>10s} {'rel_b':>10s} {'mu_re':>10s}")
    for d in defs:
        print(f"{d['chi']:10s} {d['gamma1']:7.3f} {d['rel_a']:10.4f} {d['rel_b']:10.4f} {d['mu_re']:+10.4f}")

    # Log-Log
    even = [d for d in defs if d['chi'] in ('chi_5','chi_8','chi_33')]
    if len(even)==3:
        gs = np.array([d['gamma1'] for d in even])
        for m in ('rel_a','rel_b'):
            v = np.array([d[m] for d in even])
            if np.all(v>0):
                sl,_=np.polyfit(np.log(gs),np.log(v),1)
                print(f"Slope ({m} auf Even-Familie): {sl:+.3f}")

    # Speichern
    out = RESULTS_DIR / "MILESTONE_2_CHI_V2_RHO_WIDEBAND_2026-04-18.json"
    json.dump({'config':{'lambda':float(LAMBDA),'L':float(L),'T_PW':T_PW,'T_WIDE':T_WIDE,
                         'N_grid':N_GRID,'N_galerkin':N_GALERKIN,'rho_cutoff':RHO_CUTOFF},
               'eigenvalue_analysis':eigs,'defect_norms':defs},
              open(out,'w'), indent=2)
    print(f"\nRohdaten: {out}")

if __name__ == "__main__":
    main()
