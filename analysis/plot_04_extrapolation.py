"""
Plot 4: Finite-size extrapolation of thermal conductivity.
Plots 1/kappa vs 1/L and extrapolates to 1/L -> 0 (bulk limit).

Physics: phonons with mean free path > L/2 are scattered by the
Muller-Plathe thermostat slabs, suppressing their contribution to
thermal conductivity. Longer boxes suppress fewer phonons.

The relationship 1/kappa(L) = 1/kappa_inf + C/L is derived from
kinetic theory with a finite-size correction.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from thermal_utils import read_profile, compute_kappa, V4_PARAMS, LENGTHS

prod_time = V4_PARAMS['prod_steps'] * V4_PARAMS['dt']

# Compute kappa for each length
Lz_list = []
kappa_list = []

for LEN in LENGTHS:
    r = V4_PARAMS['runs'][LEN]
    fname = V4_PARAMS['profile_pattern'].format(LEN)
    if not os.path.exists(fname):
        continue

    coord, temp, _, _ = read_profile(fname)
    res = compute_kappa(coord, temp, r['lz'], r['lx'], r['ly'],
                        r['e_total'], prod_time)
    Lz_list.append(r['lz'])
    kappa_list.append(res['kappa'])

Lz = np.array(Lz_list)
kappa = np.array(kappa_list)
inv_L = 1.0 / Lz
inv_k = 1.0 / kappa

# Linear fit
coeffs = np.polyfit(inv_L, inv_k, 1)
k_bulk = 1.0 / coeffs[1]

# R-squared
pred = np.polyval(coeffs, inv_L)
ss_res = np.sum((inv_k - pred) ** 2)
ss_tot = np.sum((inv_k - inv_k.mean()) ** 2)
r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

print("=" * 60)
print("FINITE-SIZE EXTRAPOLATION")
print("=" * 60)
print(f"{'LEN':>5} {'Lz(A)':>8} {'kappa':>10} {'1/L':>12} {'1/kappa':>12}")
print("-" * 50)
for i, LEN in enumerate(LENGTHS):
    if i < len(Lz):
        print(f"{LEN:5d} {Lz[i]:8.1f} {kappa[i]:10.2f} {inv_L[i]:12.6f} {inv_k[i]:12.6f}")

print(f"\nFit: 1/kappa = {coeffs[1]:.6f} + {coeffs[0]:.2f} / L")
print(f"R-squared: {r2:.4f}")
print(f"Extrapolated bulk kappa: {k_bulk:.1f} W/m-K")
print(f"\nReferences:")
print(f"  Tersoff MD (literature): 60-100 W/m-K")
print(f"  Experimental Si (300K):  ~150 W/m-K")
print(f"\nThe low extrapolated value reflects that all system lengths")
print(f"(10.9-32.6 nm) are below the dominant phonon mean free path")
print(f"in Si (~50-300 nm). Longer boxes (L > 50 nm) are needed")
print(f"for reliable bulk extrapolation.")
print("=" * 60)

# Plot
fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(inv_L * 1000, inv_k, 'rs', markersize=14, zorder=5, label='MD data (v4)')

inv_L_fit = np.linspace(0, max(inv_L) * 1.2, 100)
inv_k_fit = np.polyval(coeffs, inv_L_fit)
ax.plot(inv_L_fit * 1000, inv_k_fit, 'b--', linewidth=2,
        label=f'Fit: 1/$\\kappa$ = {coeffs[1]:.4f} + {coeffs[0]:.1f}/L  (R$^2$={r2:.3f})')

ax.plot(0, coeffs[1], 'b*', markersize=22, zorder=5,
        label=f'Extrapolated bulk: $\\kappa$ = {k_bulk:.0f} W/mK')

ax.axhline(y=1/150, color='green', linestyle=':', linewidth=2,
           label='Expt Si: 150 W/mK')
ax.axhline(y=1/80, color='gray', linestyle=':', linewidth=2, alpha=0.6,
           label='Tersoff lit: ~80 W/mK')

ax.set_xlabel('1 / Length (10$^{-3}$ / $\\AA$)', fontsize=14)
ax.set_ylabel('1 / $\\kappa$ (m-K / W)', fontsize=14)
ax.set_title('Finite-Size Extrapolation of $\\kappa$ (v4)', fontsize=14)
ax.legend(fontsize=10)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
os.makedirs('../figures', exist_ok=True)
plt.savefig('../figures/04_extrapolation.png', dpi=300, bbox_inches='tight')
print(f"\nSaved: ../figures/04_extrapolation.png")
