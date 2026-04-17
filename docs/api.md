# API Reference — `scs` package

This document describes what the installable `scs` Python package actually ships (`scs/__init__.py`). For the richer companion computations (SCS+CIECAM02 fit, MacAdam validation, ΔE_SCS00 with Ridge, V4 reproduction), see `scripts/` — these are reproducibility scripts, not library API.

All information-theoretic quantities (`S`, `L`, `gft_check`) are in **nats** (natural log), matching the paper's convention. `S + L = log 3 ≈ 1.0986` exactly, for every π on the simplex.

---

## Color space conversion

### `scs.to_scs(xyz, Y_ref=1.0, matrix=None) -> SCSColor`

Convert CIE XYZ to SCS coordinates.

- **xyz**: array-like of 3 floats (CIE tristimulus values)
- **Y_ref**: reference white luminance, defaults to 1.0
- **matrix**: optional XYZ→LMS transform; defaults to HPE
- **Returns**: `SCSColor(ell, S, hue, pi)` dataclass with
  - `ell`: luminance ∈ [0, 1], clipped `xyz[1] / Y_ref`
  - `S`: saturation = `D_KL(π ‖ uniform)` in nats, range `[0, log 3]`
  - `hue`: SCS hue angle in degrees `[0, 360)`, formula
    `atan2(√3·(π₅−π₇), 2π₃−π₅−π₇)` on the simplex
  - `pi`: numpy array `(π₃, π₅, π₇)`, γ-weighted simplex coordinates

```python
from scs import to_scs
c = to_scs([0.95, 1.0, 1.09])
print(f"ℓ={c.ell:.3f}  S={c.S:.3f} nats  θ={c.hue:.1f}°")
# ℓ=1.000  S=0.003 nats  θ=25.0°
```

---

## Color difference

### `scs.delta_e(xyz1, xyz2, Y_ref=1.0, matrix=None) -> float`

SCS color difference between two CIE XYZ colors. Pure geodesic, zero fitted parameters.

Formula:
```
ΔE² = (3/4) · d_lum²  +  (1/4) · d_chrom²
```

where
- `d_lum = 2 |arcsin(√ℓ₁) − arcsin(√ℓ₂)|` is the Fisher–Bernoulli geodesic on the p=2 channel
- `d_chrom = 2 · arccos(Σ √(π̃₁·π̃₂))` is the Bhattacharyya geodesic on the γ-weighted simplex
- The weights 3/4 and 1/4 come from `N/(N+1)` with `N = 3` active primes (Theorem T7)

```python
d = scs.delta_e([0.95, 1.0, 1.09], [0.60, 0.50, 0.30])
```

### `scs.delta_e_lab(L1, a1, b1, L2, a2, b2, white=(0.9505, 1.0, 1.089)) -> float`

Same as `delta_e` but from CIELAB coordinates (convenience wrapper).

### `scs.fisher_luminance(Y1, Y2) -> float`

Fisher–Bernoulli geodesic between two luminance values: `d_lum = 2 |arcsin(√Y₁) − arcsin(√Y₂)|`. Exposed as a standalone function because it is also used as a feature in the hybrid `ΔE_SCS00` metric (see `scripts/delta_e_scs00.py`).

```python
d_lum = scs.fisher_luminance(0.5, 0.45)   # → 0.1002
```

---

## Sum rule (S + L = log 3)

### `scs.saturation(pi) -> float`

`D_KL(π ‖ uniform)` in nats. Equivalent to the `S` field of `SCSColor`, exposed as a standalone function for batch processing.

### `scs.luminance_entropy(pi) -> float`

Shannon entropy `H(π)` in nats. Called "chromatic entropy" `L` in the paper.

### `scs.gft_check(pi) -> (S, L, total, error)`

Verify the sum rule `S + L = log 3`, the information-theoretic identity that holds on any three-outcome probability simplex with uniform reference.

- `S`: saturation in nats (`D_KL(π ‖ uniform)`)
- `L`: chromatic entropy in nats (`H(π)`)
- `total`: `S + L` — should equal `log 3 ≈ 1.0986`
- `error`: `|total − log 3|`, typically `< 1e-15`

```python
S, L, total, err = scs.gft_check([0.5, 0.3, 0.2])
# total = 1.098612, err ≈ 2.2e-16
```

---

## Constants

| Name | Value | Source |
|------|-------|--------|
| `scs.MU_STAR` | 15 | Unique fixed point of the sieve dynamics (T7) |
| `scs.GAMMAS` | `np.array([0.8076…, 0.6963…, 0.5955…])` | Anomalous dimensions at `μ* = 15` |
| `scs.PRIMES` | `(3, 5, 7)` | Active chromatic primes |
| `scs.__version__` | `"0.2.0"` | |

The `GAMMAS` values are computed from the closed-form expression
`γ_p = 4p·q^(p−1)·(1−δ)/(μ·(1−q^p)·(2−δ))` with `q = 1 − 2/μ*` and `δ = (1 − q^p)/p`, so they are not fit parameters.

---

## Companion scripts (reproducibility, not library API)

For the larger computations that back the paper's headline numbers, see `scripts/` in the public repository. These are reproducibility artifacts, not part of the installable `scs` package, and they have their own dependencies (`scipy`, `pandas`, `scikit-learn`, `colour-science`; install via `pip install -e .[full]`).

| Script | What it reproduces |
|--------|---------------------|
| `scripts/macadam_test.py` | MacAdam 18/25 ellipse orientation wins, RMS Δθ = 37.8° |
| `scripts/delta_e_scs00.py` | ΔE_SCS00 on COMBVD: `r = 0.893 vs 0.878`, `p < 0.0001`, 5-fold CV |
| `scripts/scs_cam02_hybrid.py` | SCS + CIECAM02 hybrid on COMBVD: `r = 0.824`, 6 weights |
| `scripts/delta_e_scs.py` | Pure SCS combined metric (Fisher + Fubini–Study + bifurcation) used by `macadam_test.py` |
| `scripts/v4_summary.py` | V4 channel weights from pre-computed CSV: `L−M ≈ 0.37` vs `γ₃/Σγ = 0.385` |
| `scripts/v4_refined_analysis.py` | Full V4 fMRI pipeline (requires OpenNeuro ds005521 download) |
| `scripts/v4_neural_extraction.py` | Extraction of `v4_bold_response.csv` from raw fMRI |
| `scripts/v4_hybrid_model.py` | V4 opponent channels as CAM02 proxy on COMBVD |
| `scripts/scs_companion.py` | Self-test (`--verify`) + paper figures (`--all`) |

All listed scripts use the same canonical SCS pipeline as `scs/__init__.py` (sRGB → XYZ → LMS/HPE → γ-weighted simplex), in nats for S and L.
