# Algebraic Fourier estimates for computer-assisted proofs in state-dependent delay equations

MATLAB/INTLAB code, data and interval certificates for the computer-assisted
proofs (Theorem 2) in

> J.-P. Lessard, *Algebraic Fourier Estimates for Computer-Assisted Proofs in
> State-Dependent Delay Equations*, AMMCS 2026 proceedings (submitted).

The code proves the existence of a 2π-periodic solution of the forced
state-dependent delay equation

    u'(t) = -u(t) - ε₁ u(t - 1 - ε₂ u(t)) + cos t + sin t,   (ε₁, ε₂) = (0.5, 2),

near trigonometric approximations of order N_u = 800 and N_u = 1300, with
certified C¹ error bounds r₀ ≤ 9.655·10⁻¹² and r₀ ≤ 2.451·10⁻¹², respectively.

## Requirements

* MATLAB (the certificates were produced with R2024a)
* [INTLAB](https://www.tuhh.de/ti3/rump/intlab/) (version 13), which is not included here

## Contents

| Path | Description |
|---|---|
| `data/ubar_800.txt`, `data/ubar_1300.txt` | Approximate solutions. The first row is β̄ (an approximation of ū(−π)); the remaining rows are the Fourier coefficients Ū_n, n = −N_u, …, N_u (real and imaginary parts). |
| `run_proof.m` | Runs the full proof for a given N_u and saves a log and an interval report. |
| `verify_paper_bounds.m` | Rechecks the bounds printed in Table 2 of the paper from the saved certificates. |
| `selftest_intlab.m` | Small interval regression tests. |
| `check_intlab.m`, `verify_radii.m` | INTLAB sanity check, and verification of the Newton–Kantorovich inequalities (Theorem 1). |
| `private/` | Implementation of the bounds (see below). |
| `results/` | Log and `.mat` interval report of the validation run for each N_u. Each report also stores the exact source code used. |

## Correspondence with the paper

| File in `private/` | Paper |
|---|---|
| `bell_table.m` | Bell polynomial recurrence (Sect. 3) |
| `Cs_crude_SDDE.m` | A priori constants 𝒦_s of Lemma 3, equation (7) |
| `Cs_grid_SDDE.m` | Midpoint enclosure K̃_s(v) of Lemma 4, with Faà di Bruno's formula (6) |
| `alias_bound.m` | Aliasing bound E_s(n,M) of Lemma 2 |
| `F_SDDE_rig_Cs.m` | Enclosure of the Galerkin part of F(x̄), Proposition 1(a) |
| `DF_SDDE_rig_Cs.m` | Enclosure of the Galerkin matrix of DF(x̄), Proposition 1(c)–(e) |
| `residual_norm_Cs.m` | Pointwise residual bound, second estimate of Lemma 5 |
| `computation_bounds_Cs.m` | Y from (8), Lemma 5 and (9); Z₁, Z₂ following the companion preprint (Sect. 2.6) |
| `comp_mat_norm.m` | C⁰ operator-norm estimate of the companion preprint (Sect. 2.6.6) |
| `kolmogorov_theta_N.m` | Kolmogorov constant θ_N (companion preprint, Sect. 2.6.3) |
| `DF_SDDE.m` | Floating-point Jacobian, used only to build the approximate inverse A |
| `fourier_multiplication_matrix.m`, `complex_error_disc.m`, `upper_decimal.m` | Helpers |

## Usage

Initialize INTLAB in MATLAB (`startintlab`), then change to this directory:

```matlab
selftest_intlab                    % interval regression tests
verified = verify_paper_bounds();  % recheck the printed Table 2 bounds
r800  = run_proof(800);            % full proof, N_u = 800
r1300 = run_proof(1300);           % full proof, N_u = 1300
```

`run_proof` does not clear the workspace or change the working directory, and it
resolves data paths relative to itself. Do not add `private/` to the MATLAB path
separately.

## Citation

If you use this code, please cite the paper above and this archive (see
`CITATION.cff`).

## License

MIT License, see `LICENSE`.
