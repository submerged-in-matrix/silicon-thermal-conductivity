import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

def read_temp_profile(filename):
    """Read LAMMPS fix ave/chunk output. Returns all snapshots averaged."""
    snapshots = []
    current = []
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) == 3:
                if current:
                    snapshots.append(np.array(current))
                current = []
            elif len(parts) == 4:
                current.append([float(parts[1]), float(parts[3])])
    if current:
        snapshots.append(np.array(current))
    
    all_data = np.array(snapshots)
    avg = np.mean(all_data, axis=0)
    std = np.std(all_data, axis=0)
    n_snap = len(snapshots)
    return avg[:, 0], avg[:, 1], std[:, 1] / np.sqrt(n_snap), n_snap


def compute_kappa(coord_reduced, temp, lz, lx, ly, energy_swapped,
                  production_steps, dt):
    """
    Compute thermal conductivity from Muller-Plathe NEMD.
    kappa = Q / (2 * A * dT/dz)
    """
    production_time = production_steps * dt
    Q_evps = energy_swapped / production_time
    Q_watts = Q_evps * 1.602e-7
    
    z = coord_reduced * lz
    n = len(coord_reduced)
    mid = n // 2
    
    skip = 2
    left_mask = (np.arange(n) >= skip) & (np.arange(n) <= mid - skip)
    right_mask = (np.arange(n) >= mid + skip) & (np.arange(n) <= n - 1 - skip)
    
    left_coeffs = np.polyfit(z[left_mask], temp[left_mask], 1)
    right_coeffs = np.polyfit(z[right_mask], temp[right_mask], 1)
    
    left_grad = abs(left_coeffs[0])
    right_grad = abs(right_coeffs[0])
    avg_grad = (left_grad + right_grad) / 2.0
    
    grad_Km = avg_grad * 1e10
    A_m2 = lx * ly * 1e-20
    
    kappa = Q_watts / (2.0 * A_m2 * grad_Km)
    
    return kappa, Q_evps, avg_grad, left_grad, right_grad, left_mask, right_mask, left_coeffs, right_coeffs


# === MAIN ===
dt = 0.001
prod_steps = 400000

runs = {
    20: {'lz': 108.533, 'lx': 32.560, 'ly': 32.560,
         'energy_total': 404.605, 'atoms': 5760},
    30: {'lz': 162.815, 'lx': 32.563, 'ly': 32.563,
         'energy_total': 419.070, 'atoms': 8640},
    40: {'lz': 217.100, 'lx': 32.565, 'ly': 32.565,
         'energy_total': 432.931, 'atoms': 11520},
    60: {'lz': 325.548, 'lx': 32.555, 'ly': 32.555,
         'energy_total': 454.101, 'atoms': 17280},
}

profile_files = {
    20: 'temp_profile_L20_v3_prod.dat',
    30: 'temp_profile_L30_v3_prod.dat',
    40: 'temp_profile_L40_v3_prod.dat',
    60: 'temp_profile_L60_v3_prod.dat',
}

print("=" * 70)
print("SILICON THERMAL CONDUCTIVITY - NEMD ANALYSIS (Muller-Plathe)")
print("=" * 70)

results = {}

for LEN in [20, 30, 40, 60]:
    r = runs[LEN]
    fname = profile_files[LEN]
    if not os.path.exists(fname):
        print(f"WARNING: {fname} not found, skipping LEN={LEN}")
        continue
    
    coord, temp, temp_err, n_snap = read_temp_profile(fname)
    
    kappa, Q, avg_grad, l_grad, r_grad, l_mask, r_mask, l_coeffs, r_coeffs = \
        compute_kappa(coord, temp, r['lz'], r['lx'], r['ly'],
                     r['energy_total'], prod_steps, dt)
    
    results[LEN] = {
        'kappa': kappa, 'Q': Q, 'grad': avg_grad,
        'l_grad': l_grad, 'r_grad': r_grad,
        'coord': coord, 'temp': temp, 'temp_err': temp_err,
        'l_mask': l_mask, 'r_mask': r_mask,
        'l_coeffs': l_coeffs, 'r_coeffs': r_coeffs,
        'n_snap': n_snap,
    }
    
    print(f"\nLEN = {LEN} ({r['atoms']} atoms, Lz = {r['lz']:.1f} A)")
    print(f"  Snapshots: {n_snap}")
    print(f"  T range: {temp.min():.1f} - {temp.max():.1f} K (delta = {temp.max()-temp.min():.1f} K)")
    print(f"  Heat flux Q: {Q:.3f} eV/ps")
    print(f"  Gradients: left={l_grad:.5f}, right={r_grad:.5f}, avg={avg_grad:.5f} K/A")
    print(f"  kappa = {kappa:.1f} W/m-K")

# --- Finite-size extrapolation ---
lengths = sorted(results.keys())
Lz_arr = np.array([runs[L]['lz'] for L in lengths])
kappa_arr = np.array([results[L]['kappa'] for L in lengths])

inv_L = 1.0 / Lz_arr
inv_kappa = 1.0 / kappa_arr

coeffs_extrap = np.polyfit(inv_L, inv_kappa, 1)
kappa_bulk = 1.0 / coeffs_extrap[1]

print(f"\n{'='*70}")
print("FINITE-SIZE EXTRAPOLATION")
print(f"{'='*70}")
print(f"{'LEN':>5} {'Lz(A)':>8} {'kappa':>10} {'1/L':>12} {'1/kappa':>12}")
print("-" * 55)
for L in lengths:
    print(f"{L:5d} {runs[L]['lz']:8.1f} {results[L]['kappa']:10.1f} {1/runs[L]['lz']:12.6f} {1/results[L]['kappa']:12.6f}")

print(f"\n1/kappa = {coeffs_extrap[1]:.6f} + {coeffs_extrap[0]:.4f} / L")
print(f"Bulk kappa (1/L -> 0): {kappa_bulk:.1f} W/m-K")
print(f"Experimental Si 300K:  ~150 W/m-K")
print(f"Tersoff literature:    60-100 W/m-K")
print(f"{'='*70}")

# === FIGURES ===
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
clrs = {20: '#1f77b4', 30: '#ff7f0e', 40: '#2ca02c', 60: '#d62728'}

# (a) Temperature profiles
ax = axes[0, 0]
for L in lengths:
    r = results[L]
    z = r['coord'] * runs[L]['lz']
    ax.plot(z, r['temp'], 'o-', color=clrs[L], markersize=4, linewidth=1.5,
            label=f'L = {runs[L]["lz"]:.0f} $\\AA$ ({L} cells)')
ax.set_xlabel('Position z ($\\AA$)', fontsize=14)
ax.set_ylabel('Temperature (K)', fontsize=14)
ax.set_title('(a) NEMD Temperature Profiles', fontsize=14)
ax.legend(fontsize=10)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

# (b) Reduced coords with fits
ax = axes[0, 1]
for L in lengths:
    r = results[L]
    z_abs = r['coord'] * runs[L]['lz']
    ax.plot(r['coord'], r['temp'], 'o', color=clrs[L], markersize=3, alpha=0.5)
    z_left = z_abs[r['l_mask']]
    z_right = z_abs[r['r_mask']]
    ax.plot(r['coord'][r['l_mask']], np.polyval(r['l_coeffs'], z_left),
            '-', color=clrs[L], linewidth=2.5,
            label=f'L={L}: $\\kappa$={r["kappa"]:.0f} W/mK')
    ax.plot(r['coord'][r['r_mask']], np.polyval(r['r_coeffs'], z_right),
            '-', color=clrs[L], linewidth=2.5)
ax.set_xlabel('Reduced Position z/L', fontsize=14)
ax.set_ylabel('Temperature (K)', fontsize=14)
ax.set_title('(b) Linear Fits in Bulk Regions', fontsize=14)
ax.legend(fontsize=9)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

# (c) kappa vs length
ax = axes[1, 0]
ax.plot(Lz_arr, kappa_arr, 'ro-', markersize=12, linewidth=2, label='MD (Tersoff)')
ax.axhline(y=kappa_bulk, color='blue', linestyle='--', linewidth=1.5,
           label=f'Extrapolated: {kappa_bulk:.0f} W/mK')
ax.axhline(y=150, color='green', linestyle=':', linewidth=1.5,
           label='Experimental: 150 W/mK')
ax.set_xlabel('System Length ($\\AA$)', fontsize=14)
ax.set_ylabel('$\\kappa$ (W/m-K)', fontsize=14)
ax.set_title('(c) Size Dependence', fontsize=14)
ax.legend(fontsize=11)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

# (d) 1/kappa vs 1/L
ax = axes[1, 1]
ax.plot(inv_L * 1000, inv_kappa, 'rs', markersize=12, label='MD data')
inv_L_fit = np.linspace(0, max(inv_L) * 1.2, 100)
inv_k_fit = np.polyval(coeffs_extrap, inv_L_fit)
ax.plot(inv_L_fit * 1000, inv_k_fit, 'b--', linewidth=2, label='Linear fit')
ax.plot(0, coeffs_extrap[1], 'b*', markersize=20,
        label=f'$\\kappa_{{bulk}}$ = {kappa_bulk:.0f} W/mK')
ax.axhline(y=1/150, color='green', linestyle=':', linewidth=1.5,
           label='Expt: 150 W/mK')
ax.set_xlabel('1 / Length (10$^{-3}$ / $\\AA$)', fontsize=14)
ax.set_ylabel('1 / $\\kappa$ (m-K/W)', fontsize=14)
ax.set_title('(d) Finite-Size Extrapolation', fontsize=14)
ax.legend(fontsize=10)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('thermal_conductivity_results.png', dpi=300, bbox_inches='tight')
print(f"\nFigure saved: thermal_conductivity_results.png")
