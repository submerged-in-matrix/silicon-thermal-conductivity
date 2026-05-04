"""
Shared utilities for NEMD thermal conductivity analysis.
Used by all analysis scripts in this project.
"""
import numpy as np


def read_profile(filename):
    """
    Read LAMMPS fix ave/chunk temperature profile output.
    Returns: coord (reduced), temp (K), temp_stderr, n_snapshots
    """
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


def compute_kappa(coord_reduced, temp, lz, lx, ly, energy_total, prod_time):
    """
    Compute thermal conductivity from Muller-Plathe NEMD.

    kappa = Q / (2 * A * dT/dz)

    Parameters:
        coord_reduced: fractional z-coordinates of chunks
        temp: averaged temperature per chunk (K)
        lz, lx, ly: box dimensions (Angstrom)
        energy_total: total energy swapped during production (eV)
        prod_time: production run duration (ps)

    Returns dict with kappa and intermediate quantities.
    """
    Q_evps = energy_total / prod_time
    Q_watts = Q_evps * 1.602e-7

    z = coord_reduced * lz
    n = len(coord_reduced)
    mid = n // 2
    skip = 3

    left_mask = (np.arange(n) >= skip) & (np.arange(n) <= mid - skip)
    right_mask = (np.arange(n) >= mid + skip) & (np.arange(n) <= n - 1 - skip)

    l_coeffs = np.polyfit(z[left_mask], temp[left_mask], 1)
    r_coeffs = np.polyfit(z[right_mask], temp[right_mask], 1)

    l_grad = abs(l_coeffs[0])
    r_grad = abs(r_coeffs[0])
    avg_grad = (l_grad + r_grad) / 2.0

    grad_Km = avg_grad * 1e10
    A_m2 = lx * ly * 1e-20

    kappa = Q_watts / (2.0 * A_m2 * grad_Km)

    return {
        'kappa': kappa,
        'Q_evps': Q_evps,
        'Q_watts': Q_watts,
        'avg_grad': avg_grad,
        'l_grad': l_grad,
        'r_grad': r_grad,
        'left_mask': left_mask,
        'right_mask': right_mask,
        'l_coeffs': l_coeffs,
        'r_coeffs': r_coeffs,
    }


# Run data for both versions
V3_PARAMS = {
    'dt': 0.001,
    'prod_steps': 400000,
    'label': 'v3 (dt=1fs, swap/200, 10 slabs)',
    'runs': {
        20: {'lz': 108.533, 'lx': 32.560, 'ly': 32.560, 'e_total': 404.605, 'atoms': 5760},
        30: {'lz': 162.815, 'lx': 32.563, 'ly': 32.563, 'e_total': 419.070, 'atoms': 8640},
        40: {'lz': 217.100, 'lx': 32.565, 'ly': 32.565, 'e_total': 432.931, 'atoms': 11520},
        60: {'lz': 325.548, 'lx': 32.555, 'ly': 32.555, 'e_total': 454.101, 'atoms': 17280},
    },
    'profile_pattern': '../outputs/temp_profile_L{}_v3_prod.dat',
}

V4_PARAMS = {
    'dt': 0.0005,
    'prod_steps': 800000,
    'label': 'v4 (dt=0.5fs, swap/100, 20 slabs)',
    'runs': {
        20: {'lz': 108.543, 'lx': 32.563, 'ly': 32.563, 'e_total': 1312.338, 'atoms': 5760},
        30: {'lz': 162.781, 'lx': 32.556, 'ly': 32.556, 'e_total': 1364.795, 'atoms': 8640},
        40: {'lz': 217.050, 'lx': 32.558, 'ly': 32.558, 'e_total': 1383.792, 'atoms': 11520},
        60: {'lz': 325.548, 'lx': 32.555, 'ly': 32.555, 'e_total': 1453.000, 'atoms': 17280},
    },
    'profile_pattern': '../outputs/temp_profile_v4_L{}_prod.dat',
}

LENGTHS = [20, 30, 40, 60]
