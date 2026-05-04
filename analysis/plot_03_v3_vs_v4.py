"""
Plot 3: Compare kappa vs system length between v3 and v4.
v3: dt=1fs, swap every 200 steps (200 fs interval), 10 slabs
v4: dt=0.5fs, swap every 100 steps (50 fs interval), 20 slabs
Shows that the 4x more frequent swapping gives ~15% higher kappa
but does not resolve the finite-size limitation.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from thermal_utils import read_profile, compute_kappa, V3_PARAMS, V4_PARAMS, LENGTHS

fig, ax = plt.subplots(figsize=(10, 6))

print("=" * 60)
print("v3 vs v4 COMPARISON")
print("=" * 60)
print(f"{'LEN':>5} {'Lz(A)':>8} {'kappa_v3':>10} {'kappa_v4':>10} {'ratio':>8}")
print("-" * 45)

for version, params, marker, lstyle, label in [
    ('v3', V3_PARAMS, 's', '--', 'v3: swap every 200 fs'),
    ('v4', V4_PARAMS, 'o', '-', 'v4: swap every 50 fs')]:

    prod_time = params['prod_steps'] * params['dt']
    Lz_list = []
    kappa_list = []

    for LEN in LENGTHS:
        r = params['runs'][LEN]
        fname = params['profile_pattern'].format(LEN)
        if not os.path.exists(fname):
            continue

        coord, temp, _, _ = read_profile(fname)
        res = compute_kappa(coord, temp, r['lz'], r['lx'], r['ly'],
                            r['e_total'], prod_time)

        Lz_list.append(r['lz'])
        kappa_list.append(res['kappa'])

    ax.plot(Lz_list, kappa_list, f'{marker}{lstyle}', markersize=10,
            linewidth=2, label=label)

# Print comparison table
v3_prod = V3_PARAMS['prod_steps'] * V3_PARAMS['dt']
v4_prod = V4_PARAMS['prod_steps'] * V4_PARAMS['dt']

for LEN in LENGTHS:
    k_v3 = 0
    k_v4 = 0

    fname_v3 = V3_PARAMS['profile_pattern'].format(LEN)
    if os.path.exists(fname_v3):
        r = V3_PARAMS['runs'][LEN]
        coord, temp, _, _ = read_profile(fname_v3)
        res = compute_kappa(coord, temp, r['lz'], r['lx'], r['ly'],
                            r['e_total'], v3_prod)
        k_v3 = res['kappa']

    fname_v4 = V4_PARAMS['profile_pattern'].format(LEN)
    if os.path.exists(fname_v4):
        r = V4_PARAMS['runs'][LEN]
        coord, temp, _, _ = read_profile(fname_v4)
        res = compute_kappa(coord, temp, r['lz'], r['lx'], r['ly'],
                            r['e_total'], v4_prod)
        k_v4 = res['kappa']

    ratio = k_v4 / k_v3 if k_v3 > 0 else 0
    print(f"{LEN:5d} {V4_PARAMS['runs'][LEN]['lz']:8.1f} {k_v3:10.2f} {k_v4:10.2f} {ratio:8.2f}")

ax.set_xlabel('System Length ($\\AA$)', fontsize=14)
ax.set_ylabel('$\\kappa$ (W/m-K)', fontsize=14)
ax.set_title('Effect of Swap Frequency on $\\kappa$(L)', fontsize=14)
ax.legend(fontsize=12)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
os.makedirs('../figures', exist_ok=True)
plt.savefig('../figures/03_v3_vs_v4.png', dpi=300, bbox_inches='tight')
print(f"\nSaved: ../figures/03_v3_vs_v4.png")
