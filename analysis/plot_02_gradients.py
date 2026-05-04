"""
Plot 2: Linear gradient fits in the bulk region of each temperature profile.
Extracts kappa for each system length using Fourier's law:
    kappa = Q / (2 * A * dT/dz)
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from thermal_utils import read_profile, compute_kappa, V4_PARAMS, LENGTHS

colors = {20: '#1f77b4', 30: '#ff7f0e', 40: '#2ca02c', 60: '#d62728'}
prod_time = V4_PARAMS['prod_steps'] * V4_PARAMS['dt']

fig, ax = plt.subplots(figsize=(10, 6))

print("=" * 65)
print("THERMAL CONDUCTIVITY EXTRACTION (v4)")
print("=" * 65)
print(f"{'LEN':>5} {'Lz(A)':>8} {'dT(K)':>7} {'grad(K/A)':>10} {'Q(eV/ps)':>10} {'kappa':>8}")
print("-" * 65)

for LEN in LENGTHS:
    r = V4_PARAMS['runs'][LEN]
    fname = V4_PARAMS['profile_pattern'].format(LEN)
    if not os.path.exists(fname):
        continue

    coord, temp, temp_err, n_snap = read_profile(fname)
    z = coord * r['lz']

    res = compute_kappa(coord, temp, r['lz'], r['lx'], r['ly'],
                        r['e_total'], prod_time)

    # Plot data points
    ax.plot(coord, temp, 'o', color=colors[LEN], markersize=3, alpha=0.4)

    # Plot gradient fits
    z_left = z[res['left_mask']]
    z_right = z[res['right_mask']]
    ax.plot(coord[res['left_mask']],
            np.polyval(res['l_coeffs'], z_left),
            '-', color=colors[LEN], linewidth=2.5,
            label=f'L={LEN}: $\\kappa$ = {res["kappa"]:.1f} W/mK')
    ax.plot(coord[res['right_mask']],
            np.polyval(res['r_coeffs'], z_right),
            '-', color=colors[LEN], linewidth=2.5)

    print(f"{LEN:5d} {r['lz']:8.1f} {temp.max()-temp.min():7.0f} "
          f"{res['avg_grad']:10.5f} {res['Q_evps']:10.3f} {res['kappa']:8.1f}")

print("=" * 65)

ax.set_xlabel('Reduced Position z/L', fontsize=14)
ax.set_ylabel('Temperature (K)', fontsize=14)
ax.set_title('Gradient Fits in Bulk Region (v4)', fontsize=14)
ax.legend(fontsize=10)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
os.makedirs('../figures', exist_ok=True)
plt.savefig('../figures/02_gradient_fits.png', dpi=300, bbox_inches='tight')
print(f"\nSaved: ../figures/02_gradient_fits.png")
