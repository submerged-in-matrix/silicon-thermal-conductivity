# Thermal Conductivity of Silicon via Non-Equilibrium MD (NEMD)

## Summary

Computed the thermal conductivity of silicon at 300 K using the
Müller-Plathe reverse NEMD method with the Tersoff potential.
Investigated finite-size effects across four system lengths
(10.9–32.6 nm) and two swap parameter sets (v3, v4).

Key finding: all system lengths are below the dominant phonon mean
free path in silicon (~50–300 nm), resulting in suppressed κ values
(2–5 W/m-K) compared to literature Tersoff bulk values (60–100 W/m-K).
The 1/κ vs 1/L scaling is physically correct (R² = 0.98),
confirming the methodology is sound even though reliable bulk
extrapolation requires longer simulation cells.

## Physics

Heat in crystalline silicon is carried by phonons — quantized lattice
vibrations. Thermal conductivity depends on which phonons can
propagate freely. In a finite simulation cell with Müller-Plathe
thermostating slabs, phonons with mean free path > L/2 are
scattered by the thermostat regions before completing their journey.
This suppresses their contribution to κ, making shorter boxes give
artificially low values.

The finite-size correction follows:

```
1/κ(L) = 1/κ(∞) + C/L
```

Plotting 1/κ vs 1/L and extrapolating to 1/L → 0 gives the bulk κ.
However, this extrapolation is only reliable when the data points
span a range that includes lengths comparable to the dominant phonon
mean free path.

## Method

- **Software:** LAMMPS
- **Potential:** Tersoff (Si.tersoff)
- **System:** Bulk silicon, diamond cubic, fully periodic (p p p)
- **Cross-section:** 6×6 unit cells (32.6 × 32.6 Å), fixed
- **Lengths:** 20, 30, 40, 60 unit cells (10.9, 16.3, 21.7, 32.6 nm)
- **Heat transport method:** Müller-Plathe reverse NEMD
  (velocity swapping between hot and cold slabs)
- **Protocol per run:**
  1. Energy minimization (tight tolerances for Tersoff)
  2. NPT equilibration at 300 K, 0 bar (30 ps)
  3. NVT equilibration (20 ps, locks volume)
  4. NVE + Müller-Plathe: steady-state establishment (200 ps)
  5. NVE + Müller-Plathe: production with fresh averaging (400 ps)
- **Temperature profile:** 40 spatial bins along z, time-averaged
- **κ extraction:** Fourier's law, κ = Q / (2·A·dT/dz)

### Version Comparison

Two parameter sets were tested to understand sensitivity to the
velocity swap frequency:

| Parameter | v3 | v4 |
|-----------|-----|-----|
| Timestep | 1.0 fs | 0.5 fs |
| Swap interval | 200 steps (200 fs) | 100 steps (50 fs) |
| Swap slabs | 10 | 20 |
| Production steps | 400,000 | 800,000 |
| Physical time | 400 ps | 400 ps |

v4 matches parameters from Heris et al. (2022), who used 0.5 fs
timestep and swap every 100 steps for silicon NEMD.

## Results

### Thermal Conductivity vs System Length

| LEN | Lz (Å) | Atoms | κ_v3 (W/mK) | κ_v4 (W/mK) | ΔT_v4 (K) |
|-----|---------|-------|-------------|-------------|-----------|
| 20  | 108.5   | 5,760 | 2.0         | 2.3         | 88        |
| 30  | 162.8   | 8,640 | 2.7         | 3.0         | 101       |
| 40  | 217.1   | 11,520| 3.0         | 3.6         | 111       |
| 60  | 325.5   | 17,280| 4.2         | 4.9         | 125       |

### Finite-Size Extrapolation (v4)

```
1/κ = 0.0937 + 21.26/L    (R² = 0.98)
Extrapolated bulk κ = 10.7 W/m-K
```

### Reference Values

| Source | κ (W/m-K) |
|--------|-----------|
| This work (L=60, v4) | 4.9 |
| This work (extrapolated) | 10.7 |
| Tersoff MD, literature (bulk) | 60–100 |
| Experimental Si, 300 K | ~150 |

## Discussion

### Why κ is Low

The dominant phonon mean free path in silicon at 300 K is 50–300 nm.
The longest system (32.6 nm) is shorter than this. Phonons that
carry most of the heat are scattered by the Müller-Plathe thermostat
slabs before traversing the system, suppressing their contribution.

This is not a methodological error — it is a well-documented
finite-size effect inherent to NEMD in materials with long phonon
mean free paths. The 1/κ vs 1/L relationship is linear (R² = 0.98),
confirming the physics is correct.

### v3 vs v4

Increasing the swap frequency by 4× (v3 → v4) improved κ by
~15% but did not resolve the fundamental size limitation. Both
the energy transferred and the temperature gradient increased
proportionally, as expected for an intrinsic property.

### Path Forward

Reliable bulk extrapolation requires system lengths of 50+ nm
(LEN > 100 unit cells, > 50,000 atoms). An alternative approach
is the Green-Kubo equilibrium method, which avoids the thermostat
slab artifact but requires very long autocorrelation times.

A complementary study on silicon **nanowires** (with real free
surfaces, boundary s s p) would show how surface scattering
additionally reduces κ — directly relevant to thermoelectric
and nanoelectronic applications.

## Repository Structure

```
silicon-thermal-conductivity/
├── README.md
├── inputs/
│   ├── in.si_thermal_v3        # v3 LAMMPS input (dt=1fs, swap/200)
│   ├── in.si_thermal_v4        # v4 LAMMPS input (dt=0.5fs, swap/100)
│   └── Si.tersoff              # Tersoff potential file
├── scripts/
│   ├── run_sweep_v3.sh         # Batch runner for v3
│   └── run_sweep_v4.sh         # Batch runner for v4
├── analysis/
│   ├── thermal_utils.py        # Shared functions and run data
│   ├── plot_01_profiles.py     # Temperature profiles
│   ├── plot_02_gradients.py    # Gradient fits + kappa extraction
│   ├── plot_03_v3_vs_v4.py     # Version comparison
│   └── plot_04_extrapolation.py # 1/kappa vs 1/L extrapolation
├── outputs/
│   ├── log_v{3,4}_L{20,30,40,60}.txt
│   ├── temp_profile_v{3,4}_L{20,30,40,60}.dat
│   └── temp_profile_v{3,4}_L{20,30,40,60}_prod.dat
└── figures/
    ├── 01_temperature_profiles.png
    ├── 02_gradient_fits.png
    ├── 03_v3_vs_v4.png
    └── 04_extrapolation.png
```

## How to Reproduce

```bash
cd inputs
cp /usr/share/lammps/potentials/Si.tersoff .
# Single run:
lmp -in in.si_thermal_v4 -var LEN 20 > ../outputs/log_v4_L20.txt 2>&1
# Full sweep:
cd ../scripts
./run_sweep_v4.sh
# Analysis:
cd ../analysis
python3 plot_01_profiles.py
python3 plot_02_gradients.py
python3 plot_03_v3_vs_v4.py
python3 plot_04_extrapolation.py
```

## References

- Tersoff, J. (1988) Phys. Rev. B 37, 6991
- Müller-Plathe, F. (1997) J. Chem. Phys. 106, 6082
- Heris et al. (2022) arXiv:2203.15596 (NEMD methodology for Si)
- Feng & Ruan (2014) J. Nanomaterials (spectral phonon MFP)
- Thompson et al. (2022) Comput. Phys. Commun. 271, 108171 (LAMMPS)
