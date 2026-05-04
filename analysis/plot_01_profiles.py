"""
Plot 1: NEMD temperature profiles for all system lengths (v4).
Shows the steady-state temperature gradient established by
Muller-Plathe velocity swapping.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from thermal_utils import read_profile, V4_PARAMS, LENGTHS

colors = {20: '#1f77b4', 30: '#ff7f0e', 40: '#2ca02c', 60: '#d62728'}

fig, ax = plt.subplots(figsize=(10, 6))

for LEN in LENGTHS:
    r = V4_PARAMS['runs'][LEN]
    fname = V4_PARAMS['profile_pattern'].format(LEN)
    if not os.path.exists(fname):
        print(f"WARNING: {fname} not found")
        continue

    coord, temp, temp_err, n_snap = read_profile(fname)
    z = coord * r['lz']

    ax.errorbar(z, temp, yerr=temp_err, fmt='o-', color=colors[LEN],
                markersize=4, linewidth=1.5, capsize=2,
                label=f'L = {r["lz"]:.0f} $\\AA$ ({LEN} cells, {n_snap} snapshots)')

    print(f"LEN={LEN}: T = {temp.min():.1f} - {temp.max():.1f} K "
          f"(delta = {temp.max()-temp.min():.0f} K)")

ax.set_xlabel('Position z ($\\AA$)', fontsize=14)
ax.set_ylabel('Temperature (K)', fontsize=14)
ax.set_title('NEMD Temperature Profiles (Muller-Plathe, v4)', fontsize=14)
ax.legend(fontsize=10)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
os.makedirs('../figures', exist_ok=True)
plt.savefig('../figures/01_temperature_profiles.png', dpi=300, bbox_inches='tight')
print("\nSaved: ../figures/01_temperature_profiles.png")
