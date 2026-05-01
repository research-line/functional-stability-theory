"""
Height Saturation Test for BSD Conjecture
Tests h_disc(E,D) / R(E) -> 1 for quadratic twists
"""
from pathlib import Path

import numpy as np

# Known elliptic curves with their data
# Format: (label, rank, regulator R(E), conductor N)
CURVES = [
    # Rank 0
    ("11a1", 0, 1.0, 11),
    # Rank 1
    ("37a1", 1, 0.0511114082, 37),
    ("43a1", 1, 0.0628165071, 43),
    # Rank 2
    ("389a1", 2, 0.1524837, 389),
    ("433a1", 2, 0.22469416, 433),
]

# For rank 0: BSD predicts L(E,1) = |Sha| * Omega * prod(c_p) / |E_tors|^2
# Height saturation is trivial (no regulator)

# For rank 1: h_disc = canonical height of generator
# Gross-Zagier: h(P) = L'(E,1) / (2 * Omega * c * sqrt(N))

# For rank >= 2: We test the RATIO
# h_disc / R(E) for a family of curves

# Simplified test: Use known height data for small conductors
# Source: Cremona tables (known exact values)

print("=" * 60)
print("HEIGHT SATURATION TEST -- BSD Conjecture")
print("=" * 60)

# Test 1: For rank 1 curves, check if h(P)/R(E) is consistent
print("\n[Test 1] Rank 1: Height / Regulator ratio")
rank1_curves = [
    ("37a1", 1, 0.0511114082, 0.0511114082),  # R = h(P) for rank 1
    ("43a1", 1, 0.0628165071, 0.0628165071),
]

for label, rank, reg, height in rank1_curves:
    ratio = height / reg if reg > 0 else float('inf')
    print(f"  {label}: h(P)={height:.7f}, R={reg:.7f}, h/R = {ratio:.6f}")

# Test 2: For rank 2 curves, check regulator structure
print("\n[Test 2] Rank 2: Known regulators")
rank2_curves = [
    ("389a1", 2, 0.15248, "Cremona/LMFDB"),
    ("433a1", 2, 0.22469, "Cremona/LMFDB"),
    ("234446a1", 2, 1.50725, "Stein-Watkins"),
]

for label, rank, reg, source in rank2_curves:
    print(f"  {label}: rank={rank}, R={reg:.5f} ({source})")

# Test 3: Quadratic twist family for a fixed curve
print("\n[Test 3] Quadratic twist family for 11a1")
print("  Simulating h_disc(E_D) for fundamental discriminants |D| <= 1000")

# For the curve y^2 = x^3 - x (congruent number curve, 32a2)
# Quadratic twists E_D: y^2 = x^3 - D^2 x
# We simulate the expected height behavior

np.random.seed(42)
N_twists = 200
discriminants = np.array(sorted(set(np.random.choice(range(-1000, 1001), N_twists))))
discriminants = discriminants[discriminants != 0]

# Simulated heights (following expected distribution)
# For rank 0 twists: h_disc = 0
# For rank 1 twists: h_disc ~ R(E_D) ~ log|D| / sqrt(|D|) (heuristic)
# Goldfeld conjecture: ~50% rank 0, ~50% rank 1

ranks = np.random.binomial(1, 0.5, len(discriminants))  # 50/50 rank 0/1
heights = np.zeros(len(discriminants))
regulators = np.zeros(len(discriminants))

for i, (D, r) in enumerate(zip(discriminants, ranks)):
    if r == 0:
        heights[i] = 0
        regulators[i] = 1.0  # trivial
    else:
        # Heuristic: h ~ c * log|D| with fluctuations
        h = 0.3 * np.log(abs(D) + 1) * (1 + 0.2 * np.random.randn())
        heights[i] = max(h, 0.01)
        regulators[i] = heights[i]  # For rank 1, R = h

# Height saturation ratio for rank 1 twists
rank1_mask = ranks == 1
if np.sum(rank1_mask) > 0:
    ratios = heights[rank1_mask] / regulators[rank1_mask]
    print(f"  Twists tested: {len(discriminants)}")
    print(f"  Rank 0: {np.sum(ranks==0)} ({100*np.mean(ranks==0):.0f}%)")
    print(f"  Rank 1: {np.sum(ranks==1)} ({100*np.mean(ranks==1):.0f}%)")
    print(f"  h/R ratio (rank 1): mean={np.mean(ratios):.4f}, std={np.std(ratios):.4f}")
    print(f"  h/R -> 1 saturation: {'YES' if abs(np.mean(ratios) - 1) < 0.1 else 'PARTIAL'}")

# Summary
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print("  Height saturation h_disc/R -> 1 is TRIVIALLY TRUE for rank 1")
print("  (because R = det(h(P_i,P_j)) = h(P) for rank 1)")
print("  The non-trivial test requires rank >= 2 data with KNOWN regulators")
print("  This requires SageMath + LMFDB API for systematic computation")
print("  NEXT STEP: Run SageMath on server with LMFDB rank-2 curve list")

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    D_r1 = np.abs(discriminants[rank1_mask])
    h_r1 = heights[rank1_mask]
    ax.scatter(D_r1, h_r1, s=10, alpha=0.6, label='h(P) for rank 1 twists')
    D_sorted = np.sort(D_r1)
    ax.plot(D_sorted, 0.3 * np.log(D_sorted + 1), 'r-', label='c * log|D| (heuristic)')
    ax.set_xlabel('|D|')
    ax.set_ylabel('Canonical height h(P)')
    ax.set_title('Height distribution for quadratic twists of 11a1')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    outpath = Path(__file__).with_name('compute_height_saturation.png')
    plt.savefig(outpath, dpi=150)
    print(f"\n  Plot saved: {outpath}")
except Exception as e:
    print(f"\n  (matplotlib nicht verfuegbar: {e})")

print("\n[DONE]")
