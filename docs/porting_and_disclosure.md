# Synthetic model and implementation parity

This repository includes independently executable MATLAB and Python implementations of the same **reduced baseline** model described in the supplied technical report. It is **not** a literal translation of the full report's frequency sweeps, scenario studies, graphics, or 3-D plots.

## Deliberate public-safe changes
- Frequency: synthetic 2.4 GHz; wavelength and all mechanical displacements are recalculated accordingly.
- Geometry: 20 equally spaced panels around a full notional circumference, 12 vertical rows; a contiguous sector is active. Angular spacing is 18 degrees, and radius is computed as `R=N_panel * d_arc/(2*pi)`, where `d_arc=0.52*lambda`. The geometry is not the source system geometry.
- Vertical spacing: `0.48*lambda`. Both coordinates and beamforming phase are recomputed from these settings.
- Two generic analog elevation strategies: `steered_analog` and `fixed_analog`. No original channel count, operating band, hardware identification, or proprietary beam schedule is reproduced.
- Displacement standard deviation: `0.012*lambda`, synthetic, not a measured wind-load model.
- Generic random seed, correlation coefficients, and scan directions. No performance acceptance threshold is asserted.

## Physics and interpretation
Element position is `[R*cos(phi), R*sin(phi), z]`. Tangential displacement `s` has a panel-shared Gaussian component and an element-local Gaussian residual. Panel covariance is `K_ij=rho_panel^abs(i-j)` on the *open active sector* (not periodic). Displacement components are `[-sin(phi)*s+cos(phi)*r, cos(phi)*s+sin(phi)*r, z_shift]`. The response is `exp(j*2*pi/lambda * xyz dot u)`. Nominal weights remain fixed under deformation. The imaginary part of the difference-to-sum ratio is inverted using a monotonic nominal calibration interval.

## Run
Python (from repository root): `python src/python/run_demo.py`. Install `numpy` and `scipy` from `requirements.txt`.
MATLAB (from repository root): `addpath('src/matlab'); run_demo`.

Each writes its own CSV in `results/`. Because MATLAB uses MT19937 and Python uses PCG64 by default, **individual Monte Carlo realizations and CSV rows are not expected to match**. Compare deterministic nominal results and ensemble statistics, not row-wise samples.

## Remaining verification work
- Execute the MATLAB version in MATLAB and compare nominal geometry, weights, beam curves, and calibration with Python.
- Check sample correlation and displacement statistics with sufficient trials.
- Add parameter-sweep, frequency, histogram, and 3-D visualization parity with the report.
- Verify sign conventions, reference peak definition, and invalid-estimate handling.

No empirical results, acceptance thresholds, or original project measurements have been included.
