# Sticky-wage / endogenous-inflation track (`_sw`)

Parallel track adding a sticky-wage union block so inflation becomes endogenous.
The pre-sticky-wage code is untouched and reproducible; see `snapshot_pre_sw/`.

---

## KNOWN ISSUE (pre-existing, deliberately NOT fixed)

### McKay transition matrix rows do not sum to exactly 1

`mckay_process/TransitionMatrix.csv` has row sums in
`[0.99999999952747975, 1.00000000029196001]`, i.e. `max |rowsum - 1| = 4.725e-10`.

`calibration_mckay` (`src/calibration.jl:181-188`) reads the CSV and returns it
as-is. It does **not** normalise rows. Compare `calibration_kmv`
(`src/calibration.jl:83`), which does: `Pi_e ./= sum(Pi_e, dims=2)`.

**Consequence.** `forward_pdf!` propagates mass as `sum_s' P'[s',s]*pdf[i,s]`, so
the distribution leaks mass at ~4.7e-10 per period. Over the 800-period
transition that compounds to `(1+4.725e-10)^800 - 1 = 3.78e-07`.

It shows up as a *constant* relative gap of about `1.566e-07` between two ways of
computing the same aggregate:
  - household side: sum over the (a,z) grid weighted by `pdf_t`
    (z-marginal comes from `forward_pdf!`)
  - government side: `w_t*N_t - tau_n_t*Xi_t` using `z_dist = z_state[t,:]`
    (z-marginal comes from repeated multiplication by `P`)
Both are correct; they disagree only because the two marginals drift apart at the
mass-leak rate.

**Status: NOT FIXED, deliberately.** At 1.6e-07 it is ~5 orders of magnitude
below anything that moves a result. Normalising the rows now would change every
number in the paper for no scientific gain.

**It is PRE-EXISTING**, not introduced by the union block and not related to
`[Z-NORM]`. Verified: the `[Z-NORM]` diff to `calibration.jl` touches only
`z_vals` (the productivity GRID); `BigP` is returned unmodified and
`calibration_mckay` has no diff hunk at all. The row sums were
`1 +/- 4.725e-10` before the `[Z-NORM]` edits too.

One-line fix if ever wanted: `P_mc ./= sum(P_mc, dims=2)` in `calibration_mckay`.

---

## VINTAGE WARNING: pre-`[Z-NORM]` results in the Fourth Year folder

`output/matlab/Plot_Intro.m`, `Plots_series.m` and `Plots_series_long.m` each
contain a hardcoded

```
cd("C:\Users\Leonardo Barreto\Google Drive\GoodNotes\Fourth Year\Aiyagari RPI\AiyagariOpen\src\Matlab")
```

and read `model_time_series_ratio.mat` from **there**, not from this repo.
That file is dated Sep 2025 and is the **pre-`[Z-NORM]` vintage**.

Measured against this repo's post-`[Z-NORM]` Stage 0 output:

| series | median Stage0 / FourthYear |
|---|---|
| y, b, bstar, T, Tu, G | 1.5733 |
| k, X | 1.5763 |
| n | 1.5719 |
| C | 1.5724 |
| tau_n | 1.0857 |
| w | 1.0009 |
| r | 0.9883 |
| **pi, tau_k, un** | **1.0000 (identical)** |

Every LEVEL is scaled by ~1.573; every ratio and exogenous series is identical.
That is exactly the `[Z-NORM]` level factor.

**Implication.** Any figure comparison must use the post-`[Z-NORM]` Stage 0
results as the reference. Comparing `_sw` output against the Fourth Year file
would attribute a 57% level shift from `[Z-NORM]` to sticky wages.
The `_sw` plotting scripts therefore read `output/matlab/` and
`output/matlab_sw/` only, and never `cd` outside this repo.

---

## `[Z-NORM]` (applied, uncommitted)

Rescales the productivity grid by `kappa = 1.2530296932321834` so that
`Z = sum(Pss .* z) = 1` in the INITIAL steady state (it was 0.798). Imposed at
the initial steady state ONLY: `Z_t` still falls ~5.5% on impact and is ~0.45%
low at t=40. Tagged `# [Z-NORM]` on every changed line in
`src/calibration.jl`, `src/transition.jl`, `scripts/main_JMP.jl`.

NOT level-neutral: the economy is ~1.573x larger in levels, `is_pol_valid`
returns false with 5.9e-06 of mass at `a_max = 120`, and a/y drifts to 20.281
against the 20.24 target. Accepted, not fixed.

---

## Track layout

| file | role |
|---|---|
| `snapshot_pre_sw/` | verbatim pre-sticky-wage copy of `src/` + `scripts/`. Reference for reproducing the current paper. Do not edit, do not add to any include path. |
| `scripts/main_JMP_sw.jl` | `_sw` entry point. `main_JMP.jl` is byte-identical to the snapshot. |
| `src/union_block.jl` | `labor_supply_union_`, `hsv_base_union_`, `union_aggregates_` |
| `src/transition_sw.jl` | forked solver + GBC helpers + forked household block |
| `output/jld2_sw/`, `output/matlab_sw/` | `_sw` outputs; the pre-SW `output/jld2/` and `output/matlab/` are never written by this track |

Namespace rule: Julia has no file-level namespacing without modules, so every new
or forked function carries a `_union` / `_sw` suffix. Verified dynamically after
each stage that all originals still resolve to their own `file:line`.

Write rule: `main_JMP_sw.jl` writes no `.jld2`/`.mat` unless the solve converged
AND `SW_WRITE=1`. A stalled point must never sit on disk looking like a result.
`output/jld2_sw/baseline_stage0.jld2` is `chmod 444`.

---

## KNOWN ISSUE: functions that read GLOBALS instead of arguments

Two functions in existing `src/` read globals rather than taking them as
arguments. Both are on the live transition path:

| location | function | reads |
|---|---|---|
| `src/time_series.jl:25` | `make_time_series` | global `g` (for `g.π_bar`) |
| `src/household_problem_and_aggregation.jl:455` | `backwards_once_t0_T!` | global `t` (for `t.δ`) |

Nine more are in `do_plots_slides*` / `do_plots_slides_long*`
(`auxiliary_functions.jl:351, 396, 492-493, 524, 619, 659-660, 698, 723`).
None of those is called by either main, so they are inert today, but they are
landmines if plotting is ever wired into Julia.

**Loud failure** (harmless): call `make_time_series` from inside a function
without a global `g` in scope and you get `UndefVarError: g not defined`.

**Silent failure** (the dangerous one): if two runs coexist in one Julia
process, `g` and `t` resolve to whichever was assigned LAST. `make_time_series`
would then build `π_path` from the wrong `π_bar` with no error at all. This is
the real hazard for Stage 4 counterfactuals.

### Hardening applied

1. `_assert_globals_sw(e_init, t, g)` in `src/transition_sw.jl`, called at the
   top of `solve_transition_h_r_sw`. Checks that `Main.t` / `Main.g`, the
   arguments `t` / `g`, and `e_init.t` / `e_init.g` all agree on
   `δ, A, α, π_bar, τ, ϑ`, and that `g.π_bar == e_init.π_`. Errors otherwise.

2. **One run per Julia process**, enforced in code: `_sw_run_count` is a `Ref`
   incremented by the same guard, and a second call to `solve_transition_h_r_sw`
   in the same process raises. Every counterfactual gets a fresh `julia`.
   This removes the coexistence hazard rather than relying on the check.

---

## CORRECTION to the McKay row-sum note above

The `1.566e-07` relative offset does **not** compound over the 800-period
transition. Measured on `baseline_stage3_kw05.jld2`:

```
    t     sum(pdf_t)-1      sum(z_state[t,:])-1
    1  +1.566403e-07        +1.299338e-11
  800  +1.667685e-07        +1.014120e-08
```

It is already `1.5664e-07` at t=1 and grows only `1.0e-08` across the whole
transition. The offset is inherited from `e_init.ws.pdf`:
`update_stationary_pdf!` (`household_problem_and_aggregation.jl:148-165`)
iterates `forward_pdf!` to convergence and assigns `pdf .= pdf_0` **without
renormalising**. At ~4.7e-10 leaked per iteration, `1.5664e-07 / 4.7e-10 ~ 333`
iterations, which is the right order for a persistent process.

The household side then uses a `pdf` summing to `1 + 1.5664e-07`, while the
government side uses `z_state` built from `ergodic(BigP)`, which IS normalised
(`auxiliary_functions.jl:11`) and sums to `1 + 1.3e-11`. The ratio is the
offset, constant at every date because it is one scale factor rather than an
accumulation.

Alternative one-line fix (leaves `P` untouched): `rdiv!(pdf, sum(pdf))` after
`update_stationary_pdf!`. Still NOT applied, same reasoning as above.

---

## FILE INVENTORY: `output/jld2_sw/` (as of 2026-09-09)

Every `.jld2` in the folder, what produced it, and its headline numbers. Sizes are
~820 MB each except where noted (the whole `path` StructArray is stored, one
`HouseholdWorkspace` per period).

**Stored tuple shape.** `reduced_*` and `cf_*` are **4-tuples**
`(e_init, e_final, path, pi_endog)`. Everything older is a **3-tuple**
`(e_init, e_final, path)` — inflation was fed then, so there was no solved `pi` to
store. Read them with `load_object`, and check `length(o)` before reaching for
`o[4]`.

### The two files that matter

| file | produced by | what it is |
|---|---|---|
| `reduced_epsm.jld2` | `scripts/main_JMP_reduced.jl`, `FEED_EPSM = true` | **THE BASELINE TO REPORT.** kappa_w = 0.05, phi = 1.5, `epsm` fed at its backed-out path so `pi` matches the data exactly. |
| `steady_states.jld2` | `main_JMP_reduced.jl` first run (3.4 MB) | `(e_init, e_final)` cache for the BASELINE. Every `_sw` script loads it instead of re-solving. A counterfactual that moves the terminal steady state must NOT reuse `e_final` from here. |

### Full inventory

| file | produced by | knobs | headline numbers |
|---|---|---|---|
| `baseline_stage0.jld2` | pre-extension code (`snapshot_pre_sw/`) | — | Stage 0 reference. **`chmod 444`, read-only, never overwrite.** |
| `warmstart_stage2.jld2` | Stage 2 work (18 KB) | — | small warm-start vector, not a result |
| `baseline.jld2` | `main_JMP_sw.jl`, Stage 2/3 era (2026-08-25) | — | **provenance not recorded.** Superseded; do not report from it. |
| `baseline_flexcheck.jld2` | Stage 2, `s_mode = :flex` (2026-08-28) | union OFF | the Stage 2 test that `:flex` reproduces Stage 0 |
| `baseline_stage3_kw05.jld2` (+ `.txt`) | `main_JMP_sw.jl`, Stage 3 | kappa_w = 0.05, **pi FED** | 3-tuple. Still load-bearing: it is the source of the backed-out `epsm` and of the Stage 3 warm start. CE mean 1.797. |
| `reduced_epsm.jld2` | `main_JMP_reduced.jl`, `FEED_EPSM = true` | kappa_w = 0.05, phi = 1.5, nkpc_scale = 240 | **BASELINE.** pi reproduces the data (peak +6.78pp vs target). r_f = 1.0778% ann, i_ss = 3.0994% ann. CE mean 1.797, poor 3.004, middle 1.971, rich -0.260. |
| `reduced_run3.jld2` | `main_JMP_reduced.jl`, `FEED_EPSM = false` | same, `epsm = 0` | **COUNTERFACTUAL, not the baseline.** The headline result: peak deviation ABOVE target +0.016pp against the data's +7.06pp, i.e. **0.23% of the surge**. pi range over j=2..40 is 1.30-2.02% (so the largest absolute deviation, -0.70pp, is on the DOWNSIDE — quote the upside number when the claim is "explains none of the surge"). CE mean 0.589. |
| `reduced_run3_kw30.jld2` | **`scripts/main_reduced_kw30_sw.jl`** | kappa_w = 0.30, nkpc_scale = 41.4, `epsm = 0` | The kappa_w bound: a 6x steeper wage curve moves the explained share only 0.23% -> 0.78% (+0.055pp above target). CE mean 0.645. Carries the BASELINE final SS (r_f = 1.0778% ann). Was orphaned — `main_JMP_reduced.jl` only ever writes `reduced_epsm` or `reduced_run3` — and has been reconstructed and verified; see below. |
| `cf_lower_foreign.jld2` | **`scripts/main_CF_lower_foreign_sw.jl`** | Lambda_bstar_f 1.148 -> 0.72, foreign holdings -10pp of GDP over j=2..20 | 64 iterations, inf-norm 8.261e-06, 18.6 min. r_f = 1.1305% ann, i_ss = 3.1531% ann. CE mean **1.055** (paper 1.22). Warm-start residual 0.44, against Run 3's 1.50 and the tau_k run's 161 — **this is the run that established that the Stage 5 machinery works and tau_k is a special case.** |
| `cf_cbo_debt.jld2` | **`scripts/main_CF_cbo_debt_sw.jl`** | terminal debt 120% -> 184% of GDP, linear 2025Q1-2055Q4; foreign held at 28.7% of GDP | 162 iterations, inf-norm 9.919e-06, 47.4 min. r_f = 1.4211% ann, i_ss = 3.4495% ann. CE mean **1.691** (paper 2.51). The gap is the **intercept channel**, not an error: the higher intercept costs -0.6975pp of inflation in closed form, against an observed shortfall of -0.37 to -0.51pp, the r path offsetting +0.19 to +0.33pp. |
| `steady_states_cf_lower_foreign.jld2` | `main_CF_lower_foreign_sw.jl` (first run) | — | that counterfactual's `e_final` only, cached. Created on demand. |
| `steady_states_cf_cbo_debt.jld2` | `main_CF_cbo_debt_sw.jl` (first run) | — | same, for the CBO run. |
| `cf_no_capital_subsidy.jld2` | — | — | **DOES NOT EXIST.** Never converged; blocked on the tau_k calibration question. |

### The counterfactual runners were reconstructed, not recovered

`cf_lower_foreign.jld2` and `cf_cbo_debt.jld2` were produced on 2026-09-08 by a
script that was never saved into the repo. It is not in `scripts/`, not in the VS
Code local history (`AppData/Roaming/Code/User/History` has 3924 entries and none
from that day), and not in the Julia REPL history (whose last JMP entry is
2026-09-03; the afternoon of 2026-09-08 was spent in a different project). Both
runs were driven by a script file that is gone.

The two runners now in `scripts/` were rebuilt by **inverting the solver's own
mapping out of the saved results**:

```
qbstar_path[i] = path.bstar[i] / path.b[i]
qb_path[i]     = path.b[i] * q_lag_i / path.y[i]
q_lag_i        = (1/(1+r_i)) * ((ibar+vartheta)/(1+pi_i) + (1-vartheta)/(1+pi_i)*q_i)
```

which recovers every knob exactly rather than by guesswork:

- **lower foreign demand** — `Lambda_bstar_f = 0.18*4.0 = 0.72` (the saved final
  steady state has `bstar/y_ann = 0.180000` exactly), and along `j = 2..20` foreign
  holdings drop 10pp of annual GDP, which in ratio-to-`b` units is
  `-0.1/(qb_path[j]/4)`. Reproduces e.g. `j=2: 0.2694 - 0.1/1.313375 = 0.19326028`.
  An earlier pre-sticky-wage REPL exploration used `0.287*4 - 0.1*4 = 0.748`
  instead; **that is not what the saved run used.**
- **CBO debt** — linear trend on `qb_path` from the 2024Q4 value (index 19,
  `1.218*4 = 4.872`) to `1.84*4 = 7.36` at index 143 (2055Q4), 124 quarters, step
  `0.020064516`; then flat. Foreign holdings held at `0.287*4 = 1.148` of annual
  GDP from index 20 on, so the foreign SHARE falls as debt rises.

Verify with `scripts/verify_cf_inputs_sw.jl <lower_foreign|cbo_debt>`, which
rebuilds the inputs, reconstructs `x` from the saved path, and evaluates the
reduced residual ONCE (no solve, nothing written). Three checks: the rebuilt knobs
against the saved path elementwise, the residual against the converged value the
run reported, and the `pi` recursion against the stored `pi`.

**Verification result (2026-09-09), both files:**

| check | `cf_lower_foreign` | `cf_cbo_debt` |
|---|---|---|
| rebuilt `qb_path` vs saved path | 8.88e-16 | 1.78e-15 |
| rebuilt `qbstar_path` vs saved path | 0.0 | 5.55e-17 |
| reduced residual at the saved solution | **8.261430e-06** | **9.918897e-06** |
| reported inf-norm at convergence | 8.261e-06 | 9.919e-06 |
| `pi` recursion vs stored `pi` | **0.0** | **0.0** |
| `[SW8]`/`[SW9]` guard counters | 0 / 0 / 0 | 0 / 0 / 0 |

The residual lands on the reported converged value and the `pi` recursion reproduces the stored `pi`
bit for bit, which pins `epsm` and `i_taylor_ss` as well as the fiscal knobs. The rebuilt runners
ARE the runs that produced these files.

**One inconsistency in the CBO run, reproduced deliberately.** Its foreign rule
divides by `qb_path = q*b/y` rather than `b/y`, so the terminal foreign share it
implies is `1.148/7.36 = 0.15597826` while the CBO final steady state actually has
`e_final.bstar/e_final.b = 0.15267816` (because `q_final = 0.9788 < 1`). The
transition ends 2.2% above its own terminal steady-state foreign share. It
converged and the gate passed, so it is not fatal, but it is there. Pass
`foreign_mode = :level_b` for the internally consistent variant — that is a
different run and will not reproduce `cf_cbo_debt.jld2`.

### `reduced_run3_kw30.jld2` is no longer orphaned

`scripts/main_reduced_kw30_sw.jl` reproduces it. It is `main_JMP_reduced.jl` with
three knobs changed and its own output tag:

| knob | baseline | this run |
|---|---|---|
| `KAPPA_W` | 0.05 | **0.30** |
| `NKPC_SCALE` | 240.0 | **41.4** |
| `FEED_EPSM` | true | **false** (`epsm = 0`) |

Confirmed from the saved file: it carries the **baseline** final steady state
(`r_f = 1.0778%` ann, `i_taylor_ss = 3.0994%` ann), so no steady state is re-solved
and `steady_states.jld2` is used as-is.

Verified with `scripts/verify_cf_inputs_sw.jl kw30`: rebuilt paths match the saved
path to 8.9e-16 / 2.8e-17; the reduced residual at the saved solution is
**8.292e-06** against the solver ftol of 1e-5; the `pi` recursion reproduces the
stored `pi` **exactly (0.0)** with `epsm = 0`, which is what pins `FEED_EPSM = false`
and the baseline `i_taylor_ss`; guard counters all 0.

**One caveat, stated because it matters for what the check does and does not prove.**
`nkpc_scale` only ROW-SCALES the NKPC block, so it does not change the solution set:
the residual test **bounds** it rather than identifying it (the unscaled NKPC residual
comes back at 2.003e-07). `kappa_w = 0.30` IS identified — it enters the wage Phillips
curve nonlinearly, so a wrong value would put the unscaled residual orders of
magnitude above ftol. If you need `nkpc_scale` itself pinned down, remeasure the
Jacobian diagonals as the original measurement did; 41.4 is the value implied by the
~linear scaling in `kappa_w` (240 / 5.90).

**Reporting note.** The headline statistic for this run is the peak deviation ABOVE
target, because the claim is about explaining the SURGE: +0.055pp against the data's
+7.06pp, i.e. 0.78%. The largest ABSOLUTE deviation is on the DOWNSIDE (-0.60pp) and
is a different number — the script prints both, labelled, so they do not get mixed up.

### CE welfare gains: the exact convention, and the Table 6 bin labels

**The convention, from `src/statistics.jl:6-42` (`generate_tables`).** Byte-identical
across all four welfare call sites — `generate_tables` (`statistics.jl:6`),
`generate_tables_detailed` (`statistics.jl:76`), `auxiliary_functions.jl:94`, and
`w_gains_path` (`auxiliary_functions.jl:160`) — and identical to git HEAD and to
`snapshot_pre_sw/`, so there is exactly one convention in this project:

1. Quantiles are taken over the **asset grid values**, weighted by the **initial
   steady state's invariant wealth marginal**, `sum(e_init.ws.pdf, dims=2)[:,1]`.
   Not the final steady state, not any transition-date distribution.
2. Cut points are `quantile(a_grid, ProbabilityWeights(pdf_marginal), p)` for
   `p = 0.0:0.1:1.0`, each rounded to 10 digits — **deciles**.
3. Membership is half-open, `qu[k] <= a_grid[i] < qu[k+1]`, with the top edge
   nudged `+10*eps()` so `a_max` is included.
4. Weights inside a bin are the same initial-SS **joint** pdf over (a, z), zeroed
   outside the bin. A group's CE is a pdf-weighted mean over the (a,z) cells whose
   asset value falls in the bin.
5. Poor = bin 1, Middle = plain (not mass-weighted) average of bins 5 and 6,
   Rich = bin 10. So **10 / 20 / 10**, matching the Table 6 labels.
6. The object is `100 * (path.v[1] ./ e_init.ws.v .- 1)`: value at the
   **announcement date** relative to initial-SS value, in percent.

**THE LABELS ARE OFF, and by more than a rounding.** There is an atom at the
borrowing constraint: **7.11% of households sit exactly at `a_min = -0.02`**, and
only 3 grid points lie below the 10th percentile. The weighted quantile cannot
split the atom, so `qu_10[2]` lands at `-0.00828423`, which is the **15.26%** point
of the mass. Every bin above the first is then shifted up by about 5pp:

| bin | realized mass | cumulative | grid points |
|---|---|---|---|
| 1 — labelled "bottom 10%" | **15.26%** | 0 - 15.3% | 3 |
| 2 | 10.16% | 15.3 - 25.4% | 5 |
| 3 | 9.53% | 25.4 - 35.0% | 8 |
| 4 | 9.27% | 35.0 - 44.2% | 11 |
| 5 | 9.08% | 44.2 - 53.3% | 15 |
| 6 | 9.25% | 53.3 - 62.6% | 21 |
| 7 | 9.40% | 62.6 - 72.0% | 30 |
| 8 | 9.31% | 72.0 - 81.3% | 44 |
| 9 | 9.39% | 81.3 - 90.7% | 90 |
| 10 — labelled "top 10%" | 9.33% | 90.7 - 100% | 772 |

So Table 6's **"Poor (bottom 10%)" is really the bottom 15.3%**, **"Middle (40-60%)"
is about 44-63%**, and **"Rich (top 10%)" is the top 9.3%**. This is a labelling fix
for the paper, independent of the sticky-wage extension — the same bins are used by
the published baseline.

### The three columns, computed with that exact convention

Full decile profile of the CE gain (%), announcement date, initial-SS bins:

| bin | mass | `reduced_epsm` (baseline) | `cf_lower_foreign` | `cf_cbo_debt` |
|---|---|---|---|---|
| 1 | 15.26% | 3.0038 | 1.8136 | 3.0544 |
| 2 | 10.16% | 2.6353 | 1.5000 | 2.5836 |
| 3 | 9.53% | 2.5576 | 1.5187 | 2.2839 |
| 4 | 9.27% | 2.4667 | 1.5168 | 2.0128 |
| 5 | 9.08% | 2.1531 | 1.2944 | 1.6509 |
| 6 | 9.25% | 1.7885 | 1.0294 | 1.2642 |
| 7 | 9.40% | 1.4297 | 0.7949 | 0.9242 |
| 8 | 9.31% | 0.9336 | 0.4425 | 0.5411 |
| 9 | 9.39% | 0.4260 | 0.1360 | 0.2930 |
| 10 | 9.33% | **-0.2598** | **-0.0119** | **1.3512** |

Table 6 form:

| | baseline | lower foreign | (paper) | CBO debt | (paper) |
|---|---|---|---|---|---|
| mean | 1.797 | **1.055** | 1.22 | **1.691** | 2.51 |
| min | -2.436 | -4.409 | — | -7.675 | — |
| Poor (bottom 10%) | 3.004 | **1.814** | 2.03 | **3.054** | 4.00 |
| Middle (40-60%) | 1.971 | **1.162** | 1.34 | **1.458** | 2.62 |
| Rich (top 10%) | -0.260 | **-0.012** | 0.06 | **1.351** | 0.50 |

**These supersede the group splits reported on 2026-09-08** (lower foreign
1.688 / 1.216 / 0.224; CBO 2.813 / 1.591 / 0.730). The means are unchanged, 1.055
and 1.691 — the mean is binning-independent. Where the 09-08 group numbers came from
is not known: they are not reproduced by the initial-SS bins (the paper convention
above), nor by final-SS bins, nor by terminal-path bins.

**CONSEQUENCE FOR THE WRITE-UP.** The inference that "the rich move the wrong way in
both counterfactuals" does NOT survive: against the paper the model's top decile is
0.07pp BELOW in the lower-foreign run (-0.012 vs 0.06) but 0.85pp ABOVE in the CBO
run (1.351 vs 0.50) — **opposite signs, so there is no common direction.** With the
09-08 numbers both sat above the paper, which is where that reading came from.

What IS there, and is a better story: the CBO run's decile profile is **non-monotone**
in wealth. It falls 3.05 -> 0.29 through bin 9 and then jumps to **1.351** in bin 10,
while the baseline and the lower-foreign run both fall monotonically to roughly zero
or below. The top decile holds the bonds, so a permanently higher long-run real rate
(r_f 1.0778% -> 1.4211% ann) hands it a higher return path — the same mechanism as
the intercept channel, showing up in the cross-section. That is worth its own
paragraph, and it is robust to the bin-label issue because bin 10 is the one bin
whose realized mass (9.33%) is close to its label.

---

## RULE-2 TRACK: files added 2026-09-09

The **working baseline is now `rule2_epsm.jld2`**, not `reduced_epsm.jld2`. Rule 2
allocates Table E2's "Other" (state/local grants, education, public services, federal
administration, global assistance — grants-in-aid, NOT NIPA federal consumption) at
`j = 1` only, by Table E2's own internal logic: shares checks .475556 / UI .431111 /
firms .093333 of X = $400.1bn. `T_path[1]` 0.28560000 -> 0.32446192,
`Tu_path[1]` 0.27760000 -> 0.31282997, `tau_k[1]` −0.47 -> −0.47520201. Everything
else — `qb_path`, `G_path`, `bstar_path_r`, the lockdown block, `tau_k` for j >= 2 —
is bit-identical to `reduced_epsm`.

**Why:** the debt path is fed from observed federal debt, which financed the whole
package including ~$1.3tn of "Other" that the model has no component for. So the model
borrowed for spending it never did and `tau_n` absorbed the residue — $424bn at j=1,
which is what pushed HSV retention above 1 (`1 - tau_n = -0.0022`). Rule 2 fixes that:
`1 - tau_n = +0.10810933`, +0.60% from steady state.

`epsm` was **re-derived**, not reused: the saved one was backed out of a solve with the
OLD fiscal paths. Two-step, as `baseline_stage3_kw05` -> `reduced_epsm` was built:
Stage 1 solves with `pi` FED (41 iters, 8.1 min) -> `epsm` backed out -> Stage 2 feeds
it to the reduced solver (1 iter). **Regression (b) exact: `max|pi_solved - pi_data| =
0.000000e+00`.** `max|Delta epsm|` vs the old one is 0.0018 ann pp, so holding `epsm`
fixed across counterfactuals stays valid.

**CSI acceptance test PASSED** (`csi_shares_rule2_sw.mat`, same seed 9823, N = 1e6,
both runs in one process): 30.9 / 44.5 / 45.0 / 29.0 / 28.9 / 28.8 against the old
30.9 / 44.0 / 44.7 / 29.1 / 29.0 / 28.9 and data 37 / 43 / 42 / 32 / 32 / 31. Shape
intact; the 2008-19 mean is identical to seven digits because that stretch is
pre-transition.

| file | what it is |
|---|---|
| `rule2_epsm.jld2` | **THE WORKING BASELINE** (4-tuple with `pi_endog`) |
| `rule2_stage1.jld2` | the pi-fed Stage 1 `epsm` was derived from — needed to re-derive it |
| `cf_no_capital_subsidy_rule2.jld2` | 37 iters, 6.80 min. `1 - tau_n = -0.35833` at j=1 |
| `cf_only_checks_rule2.jld2` | 40 iters, 6.62 min. `T_path[1] = 0.97110527` = $4754.5bn |
| `cf_lower_foreign_rule2.jld2` | 64 iters, 14.26 min. own `e_final` (`Lambda_bstar_f = 0.72`) |
| `steady_states_cf_lower_foreign_rule2.jld2` | that run's `e_final` cache |
| `steady_states_cf_cbo_debt_rule2.jld2` | CBO's `e_final` cache (`Lambda_b_f = 7.36`, `r_f = 1.4211%`) |
| `cf_cbo_debt_rule2.jld2` | **DOES NOT EXIST YET** — see below |
| `output/matlab/csi_shares_rule2_sw.mat` | the CSI acceptance test |

`reduced_epsm.jld2`, `csi_shares_sw.mat` and `baseline_stage0.jld2` were **never
touched**, so reverting is just a matter of pointing scripts back at the old tags.

### CBO is the one open column

It moves the terminal steady state (`r_f` 1.0778% -> 1.4211%) AND tightens the period-0
budget enough that `s_flex_1 < Gamma_n` — the lockdown ceiling stops binding, which
`min(Gamma_n, s_flex)` handles correctly but Anderson does not, because the solution
sits ON the kink. Attempt with `s_1` pinned at `Gamma_n` converged to 5.44e-06 in 37
iterations but failed the item (b) gate at j=1 by 2.125e-02 (1.566e-07 elsewhere);
attempt with `s_1` lifted ([SW20]) cycled 1.03e-02 -> 1.56e-01 -> 2.01e-02 -> 7.55e-02
over 336 iterations. `scripts/main_CF_cbo_rule2_outer_sw.jl` runs the manual outer
iteration on `s_1` with the inner 1597 solve that is known to converge.

**If CBO does sit on the s_flex branch, its 2020Q2 output falls MORE than the 8.4%
`Gamma_n` targets** — between 8.4% and 9.14%, since `y_1 ~ s_1^(1-alpha)` with `k[1]`
pinned. That makes 2020Q2 output an OUTCOME in that column rather than a calibration
target, which is a modelling choice to state explicitly.

### New solver machinery, all opt-in and default-off

`src/transition_sw.jl` gained four tagged, opt-in switches; every existing path is
bit-identical with them off:

- **`[SW13]`** `_SW8_ENABLED` (revert the `[SW8]` grid fallback) and `_sw_runs_allowed`
  (permit a staged warm start in one process)
- **`[SW19]`** `sw19` — eliminates `r_k` and `residuals[1]` exactly:
  `r_k* = (alpha*A*dfac)/[mu*(1+tau_k)*(e_init.k/n)^(1-alpha)]`. 1598 -> 1597
- **`[SW20]`** `sw20` — lifts `s_1 = min(Gamma_n, s_flex(w_1, tau_n_1))` into `x`.
  1597 -> 1598. Verified INERT on all four `Gamma_n`-binding runs (margins move by
  <= 5.4e-09, `y[1]` fall 8.4002% throughout, CE reproducing to 4 decimals)
- **`[SW21]`** `s1_pin` — place `s_1` at an arbitrary value so an OUTER iteration can
  put it on the `s_flex` branch while the inner solve stays the well-conditioned 1597

`src/counterfactuals_sw.jl` (`[SW11]`) holds the shared Stage 5 machinery:
`ce_summary_sw`, `gate_lab_tax_sw`, `period0_branch_sw`, `epsm_baseline_sw`,
`warm_start_sw`, the two counterfactual knob builders, `intercept_channel_sw`,
`guard_report_sw`, `print_pi_table_sw`.

### THE THREE EXACT ELIMINATIONS — an appendix paragraph, not a solver note

2396 -> 1597 unknowns through three exact eliminations. **All three are equations that
went trivial after a modelling change and were never revisited**, and together they are
why a three-block heterogeneous-agent transition converges in ~40 Anderson iterations:

| | eliminated | condition inverted |
|---|---|---|
| `[SW4]` | `r_b` | `q_lag = 1`, a pre-announcement initial-SS fact, inverts in closed form |
| `[SW6]` | `pi` | the monetary block inverts as a BACKWARD recursion, contracting at `(1+r)/phi ~ 0.67` **because** `phi > 1` — that IS the Taylor principle, so the bounded solution is unique and the recursion constructs it (no multi-start determinacy check needed) |
| `[SW19]` | `r_k` | period-0 capital is predetermined and period-0 hours are exogenous, so `k[1] = e_init.k` inverts |

**`[SW19]` was available only because the union block made period-0 hours exogenous.**
In the original model `n_0` came from the household labour FOC and depended on `w`,
hence on `r_k` — genuinely simultaneous, which is why `r_k` had to be an unknown.
Rationing at `s*ntilde` with `s_1` pinned makes `n_0` a constant. **So the sticky-wage
extension made the model MORE tractable, not less** — worth saying, because the natural
expectation is the opposite. (Note the exogeneity is conditional on `s_1` being pinned:
with `[SW20]` lifting it, `n_0` depends on `s_1` again, and that second channel is
stabilising because the map is decreasing in `s_1`.)

### TWO TOOLING TRAPS WORTH REMEMBERING

**1. The `[SW8]` counters cannot certify a solve.** They are cumulative over every
residual evaluation, so a non-zero total only means the fallback fired at some
intermediate Anderson iterate — the harmless case. And because `excess_g` is NOT part of
the residual vector, **a contaminated solve reports a perfectly clean inf-norm.** The
only valid test is to reset the counters and evaluate the residual ONCE at the solution
(`stage1_guard_at` / `stage2_guard_at` in `main_JMP_rule2_sw.jl`). Stage 1 fired the
fallback 94 times cumulatively yet was clean at its solution.

**2. Nine loop-locals leak forward into the period-0 block.** The interior loop assigns
`qb, qbstar, tau_k, T, Tu, G, pi_, d, covid`, and they are not reset to their period-0
values until the `qb, qbstar, tau_k, ... = qb_path[1], ...` line several statements into
that block. Code inserted in that window silently inherits **period 2's** values — which
is how the first `[SW19]` attempt computed `r_k*` at `tau_k = +0.1379` instead of
−0.4752, making `k[1]` 3.17x too large. Index the paths directly. Same hazard class as
the globals documented above, leaking forward from a loop rather than in from `Main`.

**And a testing lesson:** the standalone check of `r_k*` PASSED at 2.22e-16 while the
in-solver version was wrong by 3.17x, because the standalone version passed `tau_k`
explicitly and so bypassed the broken path. **A test that reconstructs the computation
independently validates the FORMULA, not the INTEGRATION.** Only the in-solver
comparison at a matched point caught it. Relatedly, `Meta.parseall` does NOT throw on a
syntax error — it returns `Expr(:error, ...)` nodes, so `try Meta.parseall(...) catch`
is a FALSE PASS; walk the tree for `:error`/`:incomplete` heads instead.

---

## [SW22] FORCED PERIOD-0 HOURS: the lockdown as a mandate, and the complete Table 6

### Why the switch exists

On the rule-2 fiscal baseline the CBO-debt counterfactual has **no period-0
equilibrium** under the ceiling specification `s_1 = min(Gamma_n, s_flex(w_1, tau_n_1))`.
Measured on three converged full solves, the map `g(s) = min(Gamma_n, s_flex(w_1(s),
tau_n_1(s)))` has slope **+2.13 to +2.20**, so `h(s) = s - g(s)` is *decreasing* and
`h > 0` at `s = Gamma_n`, `0.8973` and `0.8808` — there is no root at or below
`Gamma_n`. The mechanism is a fiscal-feedback spiral: fewer hours -> smaller labour
tax base -> tighter budget -> lower `tau_n` (a lower *retention level*, i.e. a heavier
tax) -> lower after-tax wage -> fewer hours. Rule 2 and CBO together trigger it;
neither alone does (pre-rule-2 CBO had margin +11.13%, ceiling binding).

`[SW22]` sets hours at period 0 to `Gamma_n * ntilde_z` **exactly**, in both the
household and the government block, rather than `min(Gamma_n, s_flex) * ntilde_z`.
`Gamma_n` is unchanged. The lockdown becomes a **mandated level** of hours rather than
a ceiling. `s_1 = Gamma_n` by construction, so there is nothing for the spiral to act
on. This is a real equilibrium of a slightly different model — households are
*required* to supply `Gamma_n * ntilde_z` rather than *permitted* to supply at most
that — with the budget, the value function and the distribution all solving
consistently around it. The cost is a labour wedge, quantified below.

### Implementation

- `src/counterfactuals_sw.jl` `[SW22]`: `labor_forced_sw`, `labor_income_forced_sw`,
  `disutility_forced_sw`, and `backwards_once_t0_T_forced_sw!` — a verbatim fork of
  `backwards_once_t0_T!` (`src/household_problem_and_aggregation.jl:450`) with the
  three `min()` sites replaced. `t` is passed **explicitly** rather than read from the
  global, closing one instance of the documented global-leak hazard.
- `src/transition_sw.jl`: a `forced_hours = false` keyword threaded through
  `_residuals_h_r_sw!` -> `_residuals_sw_c!` -> `solve_transition_sw_c`, switching the
  period-0 household call. **Default off; nothing in `src/` outside the `_sw` fork was
  modified.**
- `scripts/main_CF_forced_hours_sw.jl` — `ARGS[1] in ("baseline","cbo_debt")`.
  `[SW19]` is deliberately **OFF** for the baseline run (see below) and **ON** for CBO.

### Inertness where the ceiling already binds — CONFIRMED

`[SW19]` is held off in the baseline run on purpose: it pins `k[1] = e_init.k` exactly
rather than to solver tolerance, and `rule2_epsm` was itself solved with `sw19 = false`,
so leaving it on would inject a ~1e-8 difference that has nothing to do with forcing.

The forced-hours baseline **converged in 1 iteration** at inf-norm `7.928507e-06` —
i.e. the rule-2 path already satisfies the forced-hours residual system. Comparison
against `rule2_epsm`:

| series | max abs diff, j>=2 | at j=1 |
|---|---|---|
| `r` | **0.000000e+00** | **0.000000e+00** |
| `pi` | **0.000000e+00** | **0.000000e+00** |
| `w` | 0.000000e+00 | 4.448841e-10 |
| `tau_n` | 5.773160e-15 | 8.111857e-09 |
| `y` | 5.551115e-17 | 2.896263e-10 |

Welfare reproduces to +/-0.000 in every decile. The residual j=1 wobbles (`C` 1.7e-08,
`K` 1.7e-08) are the rule-2 solve's **own** period-0 asset-block tolerance
(`k[1] - e_init.k = -1.611e-08`), not something forcing introduced.

At the baseline solution `s_flex_1 = 1.0356764662`, **+14.78% above `Gamma_n`** — the
ceiling binds with margin, `min()` returns `Gamma_n` identically, and the fork has
nothing to change. So the modification is inert wherever the ceiling already binds,
which makes **CBO the only column where the two specifications differ**.

### The CBO column, forced hours

| | |
|---|---|
| iterations | 35 |
| wall time | 6.07 min |
| final inf-norm | 7.770807e-06 |
| item (b) gate | **PASSED**, worst relative gap 1.566e-07 at j=10 (the mass-leak constant), j=1 included |
| `[SW19]` `k[1] - e_init.k` | -2.220e-16 |
| `[SW8]` bracket / fallback, `[SW9]` floor | 0 / 0 / 0 |
| `s_1` booked by the household block | 0.90230115171930 = `Gamma_n`, diff 1.110e-16 |
| `s_flex_1` at the solution | 0.8924650106, **-1.09% BELOW `Gamma_n`** |

`s_flex_1 < Gamma_n` is the substantive point: under the ceiling specification the two
blocks would book different hours here and the run would be invalid. The mandate is
doing real work in this column, and only in this column.

**The labour wedge at j=1.** `theta * n_z^(1/nu)` against `rho_z * w * z` with
`rho_z = (1-tau) * tau_n * (w z n_z)^(-tau)`:

| | baseline (ceiling binds) | CBO (mandate binds) |
|---|---|---|
| mass-weighted average wedge | **-0.02198677** | **+0.00159929** |
| relative, uniform across z | **-15.03%** | **+1.30%** |
| reading | under-working: the mandate sits *below* free choice | **over-working: the mandate sits *above* free choice** |

The sign is exactly as expected and is a clean diagnostic of which specification
binds. The wedge is uniform in z because GHH preferences make the hours FOC
z-separable up to the HSV curvature, so the ratio is common across types. At +1.30%
the distortion is small: households are pushed 1.3% past their marginal condition.

**`tau_n` (a retention level: lower = heavier tax).**

| | `tau_n` | `1 - tau_n` | vs initial SS |
|---|---|---|---|
| initial SS | 0.88655014 | 0.11344986 | - |
| rule-2 baseline j=1 | 0.89189067 | 0.10810933 | +0.60% |
| **CBO forced j=1** | **0.74813532** | **0.25186468** | **-15.61%** |
| rule-2 baseline j=2 | 0.63830031 | 0.36169969 | -28.00% |
| CBO forced j=2 | 0.64272753 | 0.35727247 | -27.50% |

**2020Q2 output fall = 8.4002%**, matching the baseline as required by
`s_1 = Gamma_n`: `y[1]` differs from the baseline by `-8.552e-12`. `N[1]` is
bit-identical (`+0.000e+00`). The j=2 aggregates then separate (`y` -0.092%,
`C` -0.275%, `N` -0.058%, `K` -0.162%).

**`pi` and `r` against the rule-2 baseline (annualised pp).** The CBO terminal steady
state raises `r_f` from 1.0778% to 1.4211%, hence `i_taylor_ss` from 3.0994% to
3.4495%, and the intercept channel predicts `d(pi) = -d(i_ss)/(phi - 1 - r_f) =
-0.6885 pp`:

| j | solved `d(pi)` | intercept prediction | residual (r-path channel) |
|---|---|---|---|
| 20 | -0.4900 | -0.6885 | +0.1984 |
| 40 | -0.5053 | -0.6885 | +0.1832 |
| 100 | -0.3709 | -0.6885 | +0.3176 |
| 143 | -0.3852 | -0.6885 | +0.3033 |

`r[1]` moves +8.6966 pp annualised (-21.9075 -> -14.4285); `r` is permanently higher,
+0.3311 pp at j=400. `max|d(pi)|` over j=2..60 is 0.6627 pp; terminal `pi` returns to
2.0000% as it must.

### The CBO welfare pattern is INTRINSIC, not an artifact of rule 2 or the mandate

`scripts/check_cbo_welfare_sw.jl` compares the CBO-minus-baseline welfare difference
computed two ways: pre-rule-2 (`cf_cbo_debt` - `reduced_epsm`, ceiling bound
naturally) and rule-2 forced (`cf_cbo_debt_rule2_forced` - `rule2_epsm`, hours
mandated).

| decile | pre-rule-2 (ceiling) | rule-2 (mandate) | gap |
|---|---|---|---|
| 1 | +0.0506 | +0.0453 | -0.0053 |
| 2 | -0.0517 | -0.0560 | -0.0043 |
| 3 | -0.2738 | -0.2760 | -0.0023 |
| 4 | -0.4539 | -0.4562 | -0.0024 |
| 5 | -0.5022 | -0.5042 | -0.0019 |
| 6 | -0.5243 | -0.5260 | -0.0017 |
| 7 | -0.5055 | -0.5069 | -0.0014 |
| 8 | -0.3925 | -0.3936 | -0.0011 |
| 9 | -0.1330 | -0.1345 | -0.0015 |
| 10 | +1.6110 | +1.6085 | -0.0026 |

Agreement to **0.001-0.005 pp in every decile**. Two independent specifications, two
different fiscal baselines, one answer. So the pattern is a property of the CBO
experiment in this model, and the discrepancy against the paper (below) is **not**
caused by rule 2 and **not** caused by the mandate.

### TABLE 6, all four columns on one fiscal baseline

Levels, CE gain in % at the announcement date, paper bin convention (10/20/10 over
initial-SS asset deciles; note the bin-labelling correction documented above):

| | mean | poor | middle | rich | min |
|---|---|---|---|---|---|
| **baseline** model (`rule2_epsm`) | 1.910 | 3.350 | 2.043 | -0.325 | -5.383 |
| baseline paper | 1.82 | 3.03 | 2.00 | -0.20 | - |
| no capital subsidies model | 2.106 | 4.085 | 2.314 | -1.848 | -6.096 |
| no capital subsidies paper | 2.02 | 3.77 | 2.27 | -1.69 | - |
| only stimulus checks model | 4.449 | 8.929 | 4.422 | -2.037 | -7.251 |
| only stimulus checks paper | - | - | - | - | (no like-for-like column) |
| lower foreign demand model | 1.168 | 2.160 | 1.234 | -0.078 | -12.496 |
| lower foreign demand paper | 1.22 | 2.03 | 1.34 | 0.06 | - |
| **CBO debt (mandate)** model | 1.801 | 3.395 | 1.528 | 1.284 | -16.842 |
| CBO debt paper | 2.51 | 4.00 | 2.62 | 0.50 | - |

The comparable object is the **difference from the baseline**, since the levels carry
the known baseline shift:

| column | mean model / paper | poor model / paper | middle model / paper | rich model / paper |
|---|---|---|---|---|
| no capital subsidies | +0.196 / +0.20 | +0.735 / +0.74 | +0.271 / +0.27 | -1.523 / -1.49 |
| only stimulus checks | +2.539 / - | +5.579 / - | +2.378 / - | -1.712 / - |
| lower foreign demand | -0.742 / -0.60 | -1.189 / -1.00 | -0.809 / -0.66 | +0.247 / +0.26 |
| CBO debt (mandate) | -0.109 / +0.69 | +0.045 / +0.97 | -0.515 / +0.62 | +1.608 / +0.70 |

Model-minus-paper gap in the difference (pp):

| column | mean | poor | middle | rich |
|---|---|---|---|---|
| no capital subsidies | -0.004 | -0.005 | +0.001 | **-0.033** |
| lower foreign demand | -0.142 | -0.189 | -0.149 | -0.013 |
| CBO debt (mandate) | **-0.799** | **-0.925** | **-1.135** | **+0.908** |

**`no capital subsidies` reproduces the paper to 0.001-0.033 pp** — that column is
settled. `lower foreign demand` is off by 0.013-0.189 pp, a pre-existing model-vs-paper
gap unrelated to rule 2 (the pre-rule-2 and rule-2 runs agree to 0.0012 pp); it is on
the open-items list.

**`CBO debt` does not reproduce the paper, and the sign pattern differs.** The paper
has all four bins gaining roughly uniformly (+0.62 to +0.97). The model gives a
strongly redistributive profile: only decile 1 (+0.045) and decile 10 (+1.609) gain,
with deciles 3-8 losing about 0.4-0.53. This is not a convergence or specification
artifact — see the intrinsic-pattern check above.

**A hypothesis, not an established result.** The visible mechanism in the model is
that higher terminal debt raises `r_f` (1.0778% -> 1.4211%), which raises
`i_taylor_ss`, which lowers `pi` permanently through the intercept channel (-0.37 to
-0.51 pp solved). Lower inflation on nominal government debt is a transfer *to*
bondholders, and `r` is permanently higher (+0.33 pp at j=400) — both favour decile 10,
while the heavier labour tax needed to service the debt (`tau_n` -15.61% at j=1, and
`s_j` drifting to 0.9466 by j=400 against 0.9893 in the baseline) falls on the middle
deciles, who hold little of the debt. That channel is **absent from the paper's
flexible-wage version**, which would explain a divergence concentrated in this column.
This has not been isolated — testing it needs a run with the intercept channel held
fixed, which is an open item.

### Supporting paths, all four columns

| column | `tau_n` j=1 | `tau_n` j=2 | `r` j=1 (ann%) | `y` j=1 fall % | `pi` j=5 (ann%) | `pi` j=20 (ann%) |
|---|---|---|---|---|---|---|
| baseline | 0.89189067 | 0.63830031 | -21.9075 | 8.4002 | 5.3900 | 2.0000 |
| no capital subsidies | 1.35833061 | 0.63722770 | -21.9529 | 8.4002 | 5.3906 | 2.0033 |
| only stimulus checks | 0.89497590 | 0.63625406 | -22.1263 | 8.4002 | 5.3922 | 2.0148 |
| lower foreign demand | 0.79880636 | 0.65342888 | -23.3526 | 8.4002 | 5.4297 | 2.0406 |
| CBO debt (mandate) | 0.74813532 | 0.64272753 | -14.4285 | 8.4002 | 4.7142 | 1.5035 |

Initial-SS `tau_n = 0.88655014`, `y_ss = 0.0708999383`. All four columns share the
8.4002% 2020Q2 output fall by construction (`s_1 = Gamma_n` in every one).

### Files added by [SW22]

| file | role |
|---|---|
| `src/counterfactuals_sw.jl` (appended) | `[SW22]` fork + wedge diagnostic |
| `src/transition_sw.jl` (`_sw` fork) | `forced_hours` keyword, default off |
| `scripts/main_CF_forced_hours_sw.jl` | `baseline` (inertness) and `cbo_debt` runners |
| `scripts/check_cbo_welfare_sw.jl` | intrinsic-vs-artifact check, read-only |
| `scripts/table6_rule2_sw.jl` | the four-column table, read-only |
| `output/jld2_sw/cf_cbo_debt_rule2_forced.jld2` | the CBO column, converged and gated |

### Open items

1. `lower_foreign`: model -0.742 against paper -0.600 on the mean difference. The two
   model runs agree to 0.0012 pp, so this is model-vs-paper, not rule 2.
2. **`CBO debt`: the -0.80 to -1.14 pp gaps and the rich/middle sign flip.** Shown
   intrinsic to the experiment. The intercept-channel hypothesis above is untested.
3. Bin labelling: Table 6's "bottom 10%" is really the bottom 15.3%.
4. `tau_n = 1.358` in `no capital subsidies` — a retention level above 1, i.e. a net
   labour *subsidy* at j=1. The incidence question there is unsent.
5. 2020Q3's `tau_n` swing, an artifact of quarterly debt matching.

---

## [SW23] INTERCEPT-CHANNEL DECOMPOSITION


### What is held fixed across the six runs, and what the two experiments mean

The rule is `i_j = i_taylor_ss + phi_pi*(pi_j - pi_bar) + epsm_j` throughout.

| | baseline | un-accommodated CF | accommodated CF |
|---|---|---|---|
| `phi_pi` | 1.5 | 1.5 | 1.5 |
| `pi_bar` | 2% ann | 2% ann | 2% ann |
| `i_taylor_ss` | 3.0994% ann | rebuilt from own `r_final` | rebuilt from own `r_final` |
| `epsm` | baseline path | baseline path, untouched | baseline `- d(i_ss)` for `j <= 600` |

The **reaction coefficient and the target are identical in all six runs.** What differs
is the intercept (3.0994% baseline, 3.4495% CBO, 3.1531% lower foreign) and, in the
accommodated runs, a constant level shift to the shock sequence.

**`epsm` is derived once and never re-derived.** From the rule-2 baseline Stage-1 solve,
as the Taylor-rule residual (`epsm_baseline_sw`):

```julia
i_conv[j] = (1 + r[j+1])*(1 + pi[j+1]) - 1                    # Fisher-consistent nominal rate
epsm[j]   = i_conv[j] - i_ss_base - phi_pi*(pi[j] - pi_bar)
```

i.e. whatever the rule fails to explain given the solved real path and the OBSERVED
inflation. Re-deriving it per counterfactual would force each counterfactual to
reproduce observed inflation too, destroying the object of interest — so every
counterfactual inherits the baseline's monetary shocks. The only adjustment ever made
is the constant `- d(i_ss)` over `j <= 600`.

**Fisher is imposed FIRST, globally — not as a terminal-date patch.** The order is:

1. Solve the counterfactual's terminal steady state -> `r_final_cf`.
2. Set `i_taylor_ss = (1 + r_final_cf)*(1 + pi_bar) - 1`. Once, used at every date.
3. Keep `epsm` at the baseline path.

Step 2 is forced by the INTERIOR of the recursion, not by the terminal date as such:
the recursion is pinned at `pi[800] = pi_bar` and the tail `r` converges to
`r_final_cf`, so the bracket `(1+r)(1+pi) - 1 - i_ss - epsm` must vanish there. Skip the
rebuild and the recursion's fixed point leaves `pi_bar` and contradicts the pin — that
is the divergence.

**Steps 2 and 3 together raise the recursion's constant by `d(i_ss)` at every date
relative to the baseline. That IS the intercept channel** — not a modelling choice, but
a consequence of insisting on Fisher consistency at a terminal steady state that has
moved. Switching it off means subtracting it back out over the transition and stopping
before the tail, so step 2's consistency survives. That is what `J_STOP = 600` is for.

So the pair reads:

- **un-accommodated** — same shocks; the rule's intercept **tracks the neutral rate**
  as it rises.
- **accommodated** — same total stance; the central bank **does not let the intercept
  follow `r*`**.

Both are equilibria of the same rule, differing only in the shock sequence. The
difference between them is the intercept channel.
### `i_taylor_ss` and `epsm` are the same instrument

Both enter `solve_transition_sw_c` in **exactly one place**, `src/transition_sw.jl:1092`,
and only as their **sum**:

```julia
pi_vec[j] = pi_bar + ((1 + rj1)*(1 + pi_vec[j+1]) - 1 - i_taylor_ss - epsm_path[j]) / phi_pi
```

Neither is passed to `_residuals_h_r_sw!`. So raising `i_taylor_ss` by `d` and lowering
`epsm` by `d` at the same date is **bit-identical** to doing neither. There is one
degree of freedom here, not two.

**This is why the first specification could not work.** "Rebuild `i_taylor_ss` from the
counterfactual's own `r_final` so Fisher holds at the terminal SS, then shift `epsm` by
`-(i_ss_cf - i_ss_base)` at every date" gives the recursion constant
`i_ss_base + epsm_base[j]` — measured at `+0.000e+00` against the earlier attempt that
held `i_taylor_ss` at the baseline value and diverged. Not similar; identical.

And the sanity check fails in the opposite direction to intuition. Restoring the
constant to baseline is neutral only if the **tail `r`** is also the baseline's, and in
these two counterfactuals it is not. With tail `r -> r_final_cf` and `epsm -> 0` the
recursion's own fixed point becomes

```
pi* = pi_bar + d(i_ss)/(phi - 1 - r_f)
```

| | terminal fixed point | dev from `pi_bar` |
|---|---|---|
| CBO un-accommodated | 2.000000% | +0.0000 pp |
| CBO shift-to-terminal-date | 2.700575% | **+0.6885 pp** |
| lower foreign un-accommodated | 2.000000% | +0.0000 pp |
| lower foreign shift-to-terminal-date | 2.107290% | **+0.1057 pp** |

against the hard pin `pi[800] = pi_bar`. **The un-accommodated runs are the ones with
the consistent terminal condition.**

### The specification used

`i_taylor_ss` stays rebuilt from each counterfactual's own `r_final`. `epsm` is shifted
by `-d(i_ss)` over the transition and returned to baseline for `j > J_STOP = 600`:

```
epsm_accom[j] = epsm_base[j] - d(i_ss)   for j <= 600
              = epsm_base[j]             for j >  600
```

- Backward contraction is `(1+r_f)/phi = 0.669023`, and `0.669^100 = 3.5e-18`, so
  dates `j <= 500` are numerically untouched by the switch point.
- Terminal boundary exactly consistent: measured conflict **+1.313e-09** (CBO) and
  **+1.311e-09** (lower foreign). The floor is the baseline's own
  `epsm_base[799] = -6.519e-10`, which the un-accommodated runs carry too. The
  shift-to-terminal-date variant sits at 0.0017 — rejected by four orders of magnitude.
- `beta^560 = 1.187e-03`, so at the announcement date, where the CE object lives, this
  is welfare-equivalent to a permanent accommodation. Economically: a pre-announced
  return to the original rule 150 years out.
- `d(i_ss)` is applied in **quarterly** units: `8.54511665e-04` (CBO),
  `1.31338607e-04` (lower foreign). Note `400*d` gives `+0.3418` / `+0.0525` pp while
  the compounded annualised differences are `+0.3501` / `+0.0538` pp — same object, two
  annualisations.

`scripts/main_CF_accom_sw.jl` runs the terminal check and `error()`s out **before**
spending a solve.

### Convergence and validity

| | CBO debt | lower foreign |
|---|---|---|
| iterations | 86 | 44 |
| wall time | 16.11 min | 8.22 min |
| final inf-norm | 9.451583e-06 | 8.776202e-06 |
| warm-start inf-norm | 1.374738e-01 | 2.112073e-02 |
| item (b) gate | PASSED, 1.566e-07 | PASSED, 1.566e-07 |
| `[SW19]` `k[1] - e_init.k` | -2.220e-16 | -2.220e-16 |
| `[SW8]`/`[SW9]` counters | 0 / 0 / 0 | 0 / 0 / 0 |
| `s_1` booked | `Gamma_n` to 1.11e-16 | `Gamma_n` to 1.11e-16 |
| `s_flex_1` | 1.0882582558, **+20.61%** | 0.9728407127, **+7.82%** |
| period-0 branch | ceiling binds naturally | **BRANCH OK** |
| forced hours `[SW22]` | on (for comparability) — but **inert**, ceiling binds | off |

The warm-start norms scale with `d(i_ss)` (6.5x), as they should.

**`s_flex_1 = 1.0883` in the accommodated CBO run — the ceiling binds naturally, with
+20.6% of margin.** The fiscal-feedback spiral that left the un-accommodated CBO run
with no period-0 equilibrium was itself substantially the intercept channel's doing.
`[SW22]` was left on for comparability with the un-accommodated run, but since forcing
is inert wherever the ceiling binds (established in the `[SW22]` section), **the
accommodated CBO column does not depend on the mandate at all.**

### The construction removes the intercept channel and nothing else

Inflation, `OFF - unacc` (annualised pp):

| j | CBO | lower foreign |
|---|---|---|
| 2 | +0.6751 | +0.1038 |
| 6 | +0.6769 | +0.1041 |
| 10 | +0.6791 | +0.1044 |
| 20 | +0.6822 | +0.1049 |
| 40 | +0.6863 | +0.1055 |
| 100 | +0.6904 | +0.1061 |
| 143 | +0.6905 | +0.1060 |
| 400 | +0.6900 | +0.1059 |
| 800 | +0.0000 | +0.0000 |
| **mean j=2..500** | **+0.6893** | **+0.1058** |
| **closed form `-d(i_ss)/(phi-1-r_f)`** | **-0.6885** | **-0.1057** |

The measured shift matches the closed form to **0.0008 pp** (CBO) and **0.0001 pp**
(lower foreign), flat across 400 periods, and terminal `pi` is 2.000000% in all six
runs. This is simultaneously a validation of the closed form and of the construction.

**The channel is almost purely nominal.** `r`, `OFF - unacc`, never exceeds
**±0.006 pp** (CBO) or **±0.001 pp** (lower foreign) at any reported date.

### CBO debt: the hypothesis confirmed and quantified

| decile | baseline | un-accom | channel OFF | unacc − base | OFF − base | **channel (OFF − unacc)** |
|---|---|---|---|---|---|---|
| 1 | 3.3498 | 3.3952 | 3.8438 | +0.0453 | +0.4940 | **+0.4487** |
| 2 | 2.9171 | 2.8610 | 3.3726 | -0.0560 | +0.4555 | **+0.5115** |
| 3 | 2.7590 | 2.4829 | 3.1923 | -0.2760 | +0.4334 | **+0.7094** |
| 4 | 2.5988 | 2.1425 | 2.9842 | -0.4562 | +0.3855 | **+0.8417** |
| 5 | 2.2468 | 1.7426 | 2.5433 | -0.5042 | +0.2965 | **+0.8007** |
| 6 | 1.8393 | 1.3133 | 2.0440 | -0.5260 | +0.2046 | **+0.7307** |
| 7 | 1.4451 | 0.9382 | 1.5638 | -0.5069 | +0.1186 | **+0.6255** |
| 8 | 0.9105 | 0.5169 | 0.9676 | -0.3936 | +0.0571 | **+0.4507** |
| 9 | 0.3609 | 0.2265 | 0.4529 | -0.1345 | +0.0920 | **+0.2265** |
| **10** | -0.3249 | **1.2836** | 0.3691 | +1.6085 | +0.6940 | **-0.9144** |
| mean | 1.9100 | 1.8014 | 2.2445 | -0.1087 | +0.3345 | **+0.4431** |
| poor | 3.3498 | 3.3952 | 3.8438 | +0.0453 | +0.4940 | **+0.4487** |
| middle | 2.0431 | 1.5280 | 2.2936 | -0.5151 | +0.2506 | **+0.7657** |
| rich | -0.3249 | 1.2836 | 0.3691 | +1.6085 | +0.6940 | **-0.9144** |
| min | -5.3832 | -16.8420 | -3.0920 | -11.4588 | +2.2911 | **+13.7499** |

**The profile compresses in both directions, exactly as hypothesised.** Decile 10's
+1.6085 falls to +0.6940 — the channel is **57% of the rich's entire gain**. Deciles
3-8's losses of -0.28 to -0.53 become *gains* of +0.06 to +0.43, the channel
contributing +0.45 to +0.84. The mechanism is confirmed and quantified.

`tau_n` at j=1 goes 0.89189067 (baseline) -> 0.74813532 (un-accommodated) ->
**0.94561077** (channel off): with inflation restored, less of the debt has to be
serviced out of labour taxation, so retention rises above even the baseline.

### Against the paper: the channel accounts for the whole sign disagreement

CBO, difference from the baseline:

| | un-accom | channel OFF | paper | residual (OFF − paper) |
|---|---|---|---|---|
| mean | -0.109 | +0.335 | +0.69 | -0.355 |
| poor | +0.045 | +0.494 | +0.97 | -0.476 |
| middle | -0.515 | +0.251 | +0.62 | -0.369 |
| **rich** | **+1.609** | **+0.694** | **+0.70** | **-0.006** |

With the channel off **all four bins gain**, as in the paper, and the rich matches to
**0.006 pp**. What remains is a roughly uniform -0.35 to -0.48 shortfall on the other
three — a different and more tractable residual than the sign flip it replaces. Since
the paper's version is flexible-wage and has no such channel, this is consistent with
the intercept channel being the source of the disagreement.

### Lower foreign demand: the control behaves

| decile | baseline | un-accom | channel OFF | unacc − base | OFF − base | **channel** |
|---|---|---|---|---|---|---|
| 1 | 3.3498 | 2.1605 | 2.2290 | -1.1894 | -1.1208 | **+0.0686** |
| 2 | 2.9171 | 1.7828 | 1.8608 | -1.1343 | -1.0562 | **+0.0780** |
| 3 | 2.7590 | 1.7207 | 1.8287 | -1.0382 | -0.9303 | **+0.1080** |
| 4 | 2.5988 | 1.6490 | 1.7769 | -0.9497 | -0.8219 | **+0.1278** |
| 5 | 2.2468 | 1.3883 | 1.5098 | -0.8586 | -0.7371 | **+0.1215** |
| 6 | 1.8393 | 1.0802 | 1.1909 | -0.7591 | -0.6484 | **+0.1107** |
| 7 | 1.4451 | 0.8101 | 0.9047 | -0.6350 | -0.5404 | **+0.0946** |
| 8 | 0.9105 | 0.4191 | 0.4870 | -0.4914 | -0.4235 | **+0.0679** |
| 9 | 0.3609 | 0.0702 | 0.1042 | -0.2908 | -0.2568 | **+0.0340** |
| 10 | -0.3249 | -0.0781 | -0.2148 | +0.2468 | +0.1101 | **-0.1366** |
| mean | 1.9100 | 1.1683 | 1.2358 | -0.7417 | -0.6743 | **+0.0675** |
| poor | 3.3498 | 2.1605 | 2.2290 | -1.1894 | -1.1208 | **+0.0686** |
| middle | 2.0431 | 1.2343 | 1.3503 | -0.8088 | -0.6927 | **+0.1161** |
| rich | -0.3249 | -0.0781 | -0.2148 | +0.2468 | +0.1101 | **-0.1366** |
| min | -5.3832 | -12.4962 | -10.0633 | -7.1131 | -4.6801 | **+2.4329** |

**The control passes.** In absolute terms the channel moves nothing much: +0.03 to
+0.13 across deciles 1-9 and -0.14 for decile 10, against a counterfactual effect of
-1.19 to +0.25. The reading of the channel survives.

The mean difference goes -0.742 -> -0.674 against the paper's -0.600, so the channel
also accounts for roughly **half** of the pre-existing lower-foreign gap (-0.142 ->
-0.074). Same direction as the CBO finding.

### The channel's welfare effect is linear in `d(i_ss)`

`d(i_ss)` is 6.5065x larger for CBO than for lower foreign. Dividing the CBO channel
contribution by that ratio predicts the lower-foreign one:

| | CBO channel / 6.5065 | lower foreign channel, measured | ratio |
|---|---|---|---|
| mean | +0.0681 | +0.0675 | 0.991 |
| poor | +0.0690 | +0.0686 | 0.995 |
| middle | +0.1177 | +0.1161 | 0.987 |
| rich | -0.1405 | -0.1366 | 0.972 |

Linear to within **1-3%** across two counterfactuals that differ in almost every other
respect. So the channel's welfare contribution can be read off `d(i_ss)` alone, which
makes it cheap to price in any further experiment that moves the terminal steady state.

### Files added by [SW23]

| file | role |
|---|---|
| `scripts/main_CF_accom_sw.jl` | the two runners, `ARGS[1] in ("cbo_debt","lower_foreign")` |
| `scripts/verify_accom_sw.jl` | terminal-consistency verification, read-only, no solve |
| `output/jld2_sw/cf_cbo_debt_rule2_forced_accom.jld2` | CBO, channel off |
| `output/jld2_sw/cf_lower_foreign_rule2_accom.jld2` | lower foreign, channel off |

Saved tuples are `(e_init, e_final, path, pi, epsm)` — the `epsm` actually used is
stored alongside, since it is no longer reconstructible from `rule2_stage1` alone.
Nothing in `src/` was modified by `[SW23]`; both scripts refuse to overwrite an
existing output file.

---

## [SW24] WHY THE RULE AND THE SHOCKS CANNOT BOTH BE HELD FIXED

`scripts/check_terminal_epsm_sw.jl` — read-only, no solve.

### The impossibility

At the terminal steady state `pi = pi_bar` and `r = r_f`, so the rule gives
`i = i_ss + epsm_inf` while Fisher gives `(1+r_f)(1+pi_bar) - 1`. Therefore

```
i_ss + epsm_inf = (1 + r_f)*(1 + pi_bar) - 1
```

The RHS differs across columns whenever `r_f` does:

| column | `r_f` | RHS = `i_ss_cf` |
|---|---|---|
| baseline | 1.0778% ann | 0.00765993239583 |
| CBO debt | 1.4211% ann | 0.00851444406047 |
| lower foreign | 1.1305% ann | 0.00779127100289 |

`i_ss` and `epsm_inf` are the only two objects on the left. **One equation, no free
slot: at least one must differ across columns.** Holding the intercept AND the shock
sequence fixed while terminal inflation returns to target is not a specification the
model admits. This is the sharpest form of the `[SW23]` "same instrument" point.

### A terminal-date-only correction does not work, for two independent reasons

Proposed: `i_taylor_ss = i_ss_base`, `epsm[j] = epsm_base[j]` on 1..799, and
`epsm[800] = (1+r_final_cf)*(1+pi_bar) - 1 - i_ss_base`.

**(a) `epsm[800]` is never read.** The recursion loop is `for j in (n_per-1):-1:2`, so
it reads `epsm_path[2..799]` only. Both `epsm[1]` (pi_1 is fed from data) and
`epsm[800]` are inert. The terminal identity itself closes exactly — `rule i` and
`Fisher i` agree at **0.000e+00** — but it has no effect on the solved path.

**(b) The recursion accumulates rather than decays.** With `x_j = pi_j - pi_bar`,

```
x_j = d(i_ss)/phi + ((1+r_f)/phi) * x_{j+1} ,   x_800 = 0
```

`(1+r_f)/phi = 0.669` is the rate of convergence **TO** the displaced fixed point
`x* = (d/phi)/(1 - 0.669) = d(i_ss)/(phi - 1 - r_f) = +0.6885 pp ann`, **not** decay
toward target. Running backward from `x_800 = 0`, `|x|` GROWS monotonically, reaching
99% of `x*` in **12 periods**. `x_799 = d/phi = +0.2279 pp` predicted, +0.2243 measured
— the one-period prediction is right, the subsequent behaviour is accumulation.

Deviation from target under the proposed spec:

| j | CBO | lower foreign | `beta^j` |
|---|---|---|---|
| 100 | **+0.3170** | +0.1181 | 3.00e-01 |
| 200 | +0.5114 | +0.1077 | 9.02e-02 |
| 300 | +0.6392 | +0.1058 | 2.71e-02 |
| 400 | **+0.6708** | +0.1056 | 8.14e-03 |
| 600 | +0.6822 | +0.1057 | 7.34e-04 |
| 700 | +0.6825 | +0.1057 | 2.20e-04 |
| 795 | +0.5879 | +0.0916 | 7.03e-05 |
| 799 | +0.2243 | +0.0350 | 6.70e-05 |
| 800 | 0.0000 | 0.0000 | 6.62e-05 |

+0.67 pp at j=400 and +0.32 pp at j=100 where `beta^j = 0.30` — inside the
welfare-relevant horizon, so this fails the acceptance criterion.

**(c) Spreading over the last K periods does not help**, because the displacement is
generated in the interior by the constant, not at the tail:

| K | `pi[400]` dev (CBO) | `pi[100]` dev | `pi[770]` dev | `pi[790]` dev |
|---|---|---|---|---|
| 1 | +0.6708 | +0.3170 | +0.6818 | +0.6679 |
| 20 | +0.6708 | +0.3170 | +0.6735 | -0.0082 |
| 30 | +0.6708 | +0.3170 | +0.2212 | -0.0082 |
| 50 | +0.6708 | +0.3170 | -0.0067 | -0.0082 |
| 200 | +0.6708 | +0.3170 | -0.0067 | -0.0082 |
| 799 | -0.0177 | -0.3709 | -0.0067 | -0.0082 |

Only K = 799 — restoring the constant everywhere the recursion reads — works, and that
**is** the un-accommodated run.

### The proposed spec is a system already solved

Recursion constant at j=400 (the recursion only ever sees `i_ss + epsm_j`):

| specification | constant | vs baseline |
|---|---|---|
| baseline | 0.00765984873838 | +0.000e+00 |
| un-accommodated (`i_ss_cf`, `epsm_base`) | 0.00851436040303 | +8.545e-04 |
| PROPOSED (`i_ss_base`, `epsm_base`) | 0.00765984873838 | +0.000e+00 |
| `j<=600` accommodated | 0.00765984873838 | +0.000e+00 |
| **RELABEL (`i_ss_base`, `epsm_base + d`)** | **0.00851436040303** | **+8.545e-04** |

`PROPOSED` minus `j<=600 accommodated` = **+0.000e+00**. It is that run with the tail
treatment stripped out, which is precisely what makes it inconsistent.

### THE RELABEL: a common intercept across all four columns, no new solve

`RELABEL` minus `un-accommodated` = **+0.000e+00**. The pair
`(i_ss = i_ss_base, epsm = epsm_base + d)` and `(i_ss = i_ss_cf, epsm = epsm_base)` are
the **same system**, so the un-accommodated runs already in Table 6 can be *described*
with a common intercept. `d(i_ss)` simply moves from the intercept label to the shock
label. Usable wording:

> The monetary rule is `i_t = i* + 1.5(pi_t - pi_bar) + eps_t` with `i* = 3.0994%`
> annualised and `pi_bar = 2%` in **all four columns**. The baseline monetary-shock
> path `eps_t` is held fixed across columns; the two counterfactuals that shift the
> terminal real rate carry in addition a permanent shock of **+0.34 pp** (CBO debt) and
> **+0.05 pp** (lower foreign demand) annualised, which is exactly the change in the
> neutral nominal rate that Fisher consistency at their terminal steady state requires.

What cannot be claimed without qualification is that the **shock sequence** is also
identical. That is the half the model will not give. The choice is only which label
carries `d(i_ss)`:

- carry it in the **intercept** -> current Table 6 wording, `epsm` identical across
  columns;
- carry it in the **shock** -> common `i*` across columns, `epsm` differs by a
  permanent constant.

Both describe the identical solved paths.

---

## TABLE 6 -- consolidated, post-600 adjustment specification

Generated by `scripts/table6_consolidate_sw.jl` (`[SW25]`), read-only on the saved
`.jld2` files. Plain-text form in `output/table6_sw.txt`; per-run provenance in
`output/jld2_sw/<tag>.txt`.

| # | column | file |
|---|---|---|
| 1 | No capital subsidies | `cf_no_capital_subsidy_rule2.jld2` |
| 2 | Lower foreign demand | `cf_lower_foreign_rule2_accom.jld2` |
| 3 | Only stimulus checks | `cf_only_checks_rule2.jld2` |
| 4 | CBO debt projections | `cf_cbo_debt_rule2_forced_accom.jld2` |
| - | Baseline | `rule2_epsm.jld2` |

### Is the rule common across all four columns?

| column | `i_taylor_ss` | `phi_pi` | `pi_bar` | `r_final` | `epsm` adjustment |
|---|---|---|---|---|---|
| Baseline | 3.0994% ann | 1.50 | 2.00% ann | 1.0778% ann | **none** |
| No capital subsidies | 3.0994% ann | 1.50 | 2.00% ann | 1.0778% ann | **none** |
| Lower foreign demand | 3.1531% ann | 1.50 | 2.00% ann | 1.1305% ann | `-d(i_ss)` over `j <= 600` |
| Only stimulus checks | 3.0994% ann | 1.50 | 2.00% ann | 1.0778% ann | **none** |
| CBO debt projections | 3.4495% ann | 1.50 | 2.00% ann | 1.4211% ann | `-d(i_ss)` over `j <= 600` |

**Columns 1 and 3 carry no adjustment, verified rather than assumed.** Neither moves
the terminal steady state, so `d(i_ss)` is exactly zero and `epsm` is exactly
`epsm_base`, elementwise:

| column | `d(r_final)` | `d(i_ss)` | `max abs(epsm - epsm_base)` |
|---|---|---|---|
| No capital subsidies | +0.0000 pp ann | +0.000000000000 | **0.000e+00** |
| Lower foreign demand | +0.0527 pp ann | +0.000131338607 | (adjusted by construction) |
| Only stimulus checks | +0.0000 pp ann | +0.000000000000 | **0.000e+00** |
| CBO debt projections | +0.3433 pp ann | +0.000854511665 | (adjusted by construction) |

### The rule constant, which is what the model actually uses

`i_taylor_ss` and `epsm` enter the solver only as their sum (`transition_sw.jl:1092`),
so the sum is the rule:

| column | constant at j=100 | constant at j=700 | max abs dev vs baseline, j=2..600 |
|---|---|---|---|
| Baseline | 0.00766440018159 | 0.00765993233657 | **0.000e+00** |
| No capital subsidies | 0.00766440018159 | 0.00765993233657 | **0.000e+00** |
| Lower foreign demand | 0.00766440018159 | 0.00779127094364 | **0.000e+00** |
| Only stimulus checks | 0.00766440018159 | 0.00765993233657 | **0.000e+00** |
| CBO debt projections | 0.00766440018159 | 0.00851444400122 | **0.000e+00** |

**The rule is bit-identical across all five runs at every date `j = 2..600`.** It
differs only for `j > 600`, where columns 2 and 4 revert to the value Fisher
consistency at their own terminal steady state requires. `beta^600 = 7.34e-04`, so
nothing past the switch reaches the announcement-date welfare numbers.

Equivalently, holding `i_taylor_ss` at the baseline 3.0994% in all five columns, the
shock path is exactly `epsm_base` for `j <= 600` everywhere, and columns 2 and 4 carry
a permanent extra shock for `j > 600` of

- Lower foreign demand: **+0.0525 pp** annualised
- CBO debt projections: **+0.3418 pp** annualised

which is exactly the change in the neutral nominal rate. So the rule and the shocks
are held fixed across all four columns over the entire horizon that matters, with the
qualification confined to `j > 600`.

### Validity of each saved solution

| column | iters | wall (min) | inf-norm | gate (rel) | `s_flex_1` margin | forced hours | spec re-verified |
|---|---|---|---|---|---|---|---|
| Baseline | 1 | 0.19 | 7.929e-06 | 1.566e-07 | +14.78% | no | 0.000e+00 |
| No capital subsidies | 37 | 6.80 | 4.704e-06 | 1.566e-07 | +63.90% | no | 0.000e+00 |
| Lower foreign demand | 44 | 8.22 | 8.776e-06 | 1.566e-07 | +7.82% | no | 0.000e+00 |
| Only stimulus checks | 40 | 6.62 | 9.500e-06 | 1.566e-07 | +15.12% | no | 0.000e+00 |
| CBO debt projections | 86 | 16.11 | 9.452e-06 | 1.566e-07 | +20.61% | YES | 0.000e+00 |

Gate tolerance 1e-5; `1.566e-07` is the known mass-leak constant. "Spec re-verified"
re-runs the `pi` recursion on each file's own saved `r` path with the stated
`(i_taylor_ss, epsm)` and compares to the saved `pi` -- so the specification claim is
checked against the files, not against these notes. **A positive `s_flex_1` margin
means the `Gamma_n` branch binds, which every column satisfies** -- including the CBO
column at +20.61%, so `[SW22]` forced hours is inert there and the column does not
depend on it. `s_1 = 0.90230115171930 = Gamma_n` in all five, and the 2020Q2 output
fall is 8.4002% in all five.

### Bin convention and the labelling correction

Quantiles over the asset grid weighted by the initial-SS invariant wealth marginal;
Poor = bin 1, Middle = mean of bins 5 and 6, Rich = bin 10. There is an atom at
`a_min`, so the realized masses are not deciles:

| bin | realized mass | cumulative | grid points |
|---|---|---|---|
| 1 ("bottom 10%") | 15.26% | 15.3% | 3 |
| 2 | 10.16% | 25.4% | 5 |
| 3 | 9.53% | 35.0% | 8 |
| 4 | 9.27% | 44.2% | 11 |
| 5 | 9.08% | 53.3% | 15 |
| 6 | 9.25% | 62.6% | 21 |
| 7 | 9.40% | 72.0% | 30 |
| 8 | 9.31% | 81.3% | 44 |
| 9 | 9.39% | 90.7% | 90 |
| 10 ("top 10%") | 9.33% | 100.0% | 772 |

So **"Poor (bottom 10%)" is really the bottom 15.3%**, "Middle (40-60%)" is about
44.2-62.6%, and "Rich (top 10%)" is the top 9.3%.

### (a) Levels -- CE against the pre-stimulus steady state (%)

| | Baseline | No capital subsidies | Lower foreign demand | Only stimulus checks | CBO debt projections |
|---|---|---|---|---|---|
| **Overall mean** | 1.910 | 2.106 | 1.236 | 4.449 | 2.245 |
| **Overall minimum** | -5.383 | -6.096 | -10.063 | -7.251 | -3.092 |
| **Poor (bottom 10%)** | 3.350 | 4.085 | 2.229 | 8.929 | 3.844 |
| **Middle (40-60%)** | 2.043 | 2.314 | 1.350 | 4.422 | 2.294 |
| **Rich (top 10%)** | -0.325 | -1.848 | -0.215 | -2.037 | 0.369 |

### (b) Differences from the rule-2 baseline (pp)

| | Baseline | No capital subsidies | Lower foreign demand | Only stimulus checks | CBO debt projections |
|---|---|---|---|---|---|
| **Overall mean** | - | +0.196 | -0.674 | +2.539 | +0.334 |
| **Overall minimum** | - | -0.713 | -4.680 | -1.868 | +2.291 |
| **Poor (bottom 10%)** | - | +0.735 | -1.121 | +5.579 | +0.494 |
| **Middle (40-60%)** | - | +0.271 | -0.693 | +2.378 | +0.251 |
| **Rich (top 10%)** | - | -1.523 | +0.110 | -1.712 | +0.694 |

### (c) Published paper numbers, same cells (%)

| | Baseline | No capital subsidies | Lower foreign demand | Only stimulus checks | CBO debt projections |
|---|---|---|---|---|---|
| **Overall mean** | 1.82 | 2.02 | 1.22 | 4.29 | 2.51 |
| **Overall minimum** | -1.77 | -6.18 | -3.72 | -5.80 | -1.27 |
| **Poor (bottom 10%)** | 3.03 | 3.77 | 2.03 | 8.60 | 4.00 |
| **Middle (40-60%)** | 2.00 | 2.27 | 1.34 | 4.28 | 2.62 |
| **Rich (top 10%)** | -0.20 | -1.69 | 0.06 | -1.90 | 0.50 |

### (a) minus (c) -- model level minus published level (pp)

| | Baseline | No capital subsidies | Lower foreign demand | Only stimulus checks | CBO debt projections |
|---|---|---|---|---|---|
| **Overall mean** | +0.090 | +0.086 | +0.016 | +0.159 | -0.265 |
| **Overall minimum** | -3.613 | +0.084 | -6.343 | -1.451 | -1.822 |
| **Poor (bottom 10%)** | +0.320 | +0.315 | +0.199 | +0.329 | -0.156 |
| **Middle (40-60%)** | +0.043 | +0.044 | +0.010 | +0.142 | -0.326 |
| **Rich (top 10%)** | -0.125 | -0.158 | -0.275 | -0.137 | -0.131 |

### (b) minus published differences (pp)

| | Baseline | No capital subsidies | Lower foreign demand | Only stimulus checks | CBO debt projections |
|---|---|---|---|---|---|
| **Overall mean** | - | -0.004 | -0.074 | +0.069 | -0.356 |
| **Overall minimum** | - | +3.697 | -2.730 | +2.162 | +1.791 |
| **Poor (bottom 10%)** | - | -0.005 | -0.121 | +0.009 | -0.476 |
| **Middle (40-60%)** | - | +0.001 | -0.033 | +0.098 | -0.369 |
| **Rich (top 10%)** | - | -0.033 | -0.150 | -0.012 | -0.006 |

### Full decile profile

| bin | mass | Baseline | No capital | Lower foreign | Only stimulus | CBO debt |
|---|---|---|---|---|---|---|
| 1 | 15.26% | 3.3498 | 4.0853 | 2.2290 | 8.9285 | 3.8438 |
| 2 | 10.16% | 2.9171 | 3.6212 | 1.8608 | 7.5301 | 3.3726 |
| 3 | 9.53% | 2.7590 | 3.2713 | 1.8287 | 6.7189 | 3.1923 |
| 4 | 9.27% | 2.5988 | 3.0110 | 1.7769 | 5.8302 | 2.9842 |
| 5 | 9.08% | 2.2468 | 2.5533 | 1.5098 | 4.9122 | 2.5433 |
| 6 | 9.25% | 1.8393 | 2.0743 | 1.1909 | 3.9309 | 2.0440 |
| 7 | 9.40% | 1.4451 | 1.5805 | 0.9047 | 2.9692 | 1.5638 |
| 8 | 9.31% | 0.9105 | 0.9778 | 0.4870 | 1.8761 | 0.9676 |
| 9 | 9.39% | 0.3609 | 0.3425 | 0.1042 | 0.7147 | 0.4529 |
| 10 | 9.33% | -0.3249 | -1.8482 | -0.2148 | -2.0372 | 0.3691 |

Differences from the baseline:

| bin | mass | No capital | Lower foreign | Only stimulus | CBO debt |
|---|---|---|---|---|---|
| 1 | 15.26% | +0.7355 | -1.1208 | +5.5787 | +0.4940 |
| 2 | 10.16% | +0.7042 | -1.0562 | +4.6131 | +0.4555 |
| 3 | 9.53% | +0.5124 | -0.9303 | +3.9600 | +0.4334 |
| 4 | 9.27% | +0.4122 | -0.8219 | +3.2314 | +0.3855 |
| 5 | 9.08% | +0.3064 | -0.7371 | +2.6653 | +0.2965 |
| 6 | 9.25% | +0.2350 | -0.6484 | +2.0915 | +0.2046 |
| 7 | 9.40% | +0.1353 | -0.5404 | +1.5240 | +0.1186 |
| 8 | 9.31% | +0.0674 | -0.4235 | +0.9656 | +0.0571 |
| 9 | 9.39% | -0.0184 | -0.2568 | +0.3538 | +0.0920 |
| 10 | 9.33% | -1.5233 | +0.1101 | -1.7123 | +0.6940 |

### Second panel -- what the extension adds (not in the published table)

| column | `pi` 2022Q2 (ann %) | vs baseline (pp) | `r` at j=20 (ann %) | vs baseline (pp) | peak `pi` dev from target | peak date |
|---|---|---|---|---|---|---|
| Baseline | 9.0600 | - | 0.9862 | - | +6.7823 pp | 2022Q2 |
| No capital subsidies | 9.0624 | +0.0022 | 0.9878 | +0.0016 | +6.7845 pp | 2022Q2 |
| Lower foreign demand | 9.2212 | +0.1510 | 1.0530 | +0.0663 | +6.9333 pp | 2022Q2 |
| Only stimulus checks | 9.0679 | +0.0074 | 0.9928 | +0.0065 | +6.7898 pp | 2022Q2 |
| CBO debt projections | 9.1248 | +0.0607 | 1.1086 | +0.1215 | +6.8430 pp | 2022Q2 |

2022Q2 is `j = 9`. `pi_1` is fed from the data, so the peak is taken over `j >= 2`.

**Read the peak column carefully.** The baseline's `epsm` is backed out so the baseline
reproduces observed inflation, so its +6.78pp peak is the *data's* deviation, not an
endogenous model response. What the extension generates endogenously is the
*difference across columns*, and that is small -- at most **+0.15pp** at 2022Q2 (lower
foreign demand). This is the same finding as the headline result: absent monetary
shocks the fiscal and foreign blocks move inflation very little.

| column | `tau_n` j=1 | `tau_n` j=2 | terminal `pi` (ann %) |
|---|---|---|---|
| Baseline | 0.89189067 | 0.63830031 | 2.000000 |
| No capital subsidies | 1.35833061 | 0.63722770 | 2.000000 |
| Lower foreign demand | 0.82834116 | 0.65306385 | 2.000000 |
| Only stimulus checks | 0.89497590 | 0.63625406 | 2.000000 |
| CBO debt projections | 0.94561077 | 0.64023229 | 2.000000 |

Initial-SS `tau_n = 0.88655014` (a retention *level*: lower means a heavier tax),
`y_ss = 0.0708999383`, `Gamma_n = 0.90230115171930`.

### What matches and what does not

- **Levels agree with the paper closely** in this specification: the mean is within
  0.016-0.265 pp across all five columns, Middle within 0.326 pp, Rich within 0.275 pp.
- **`Overall minimum` is the weak cell.** It is a single `(a,z)` cell with positive
  mass, so it is far more sensitive than any average. The baseline is -5.383 against
  the published -1.77, and lower foreign demand is -10.063 against -3.72; but no
  capital subsidies matches at +0.084. Treat the minimum row as not yet reconciled.
- On **differences**, `no capital subsidies` reproduces the paper to 0.001-0.033 pp and
  `only stimulus checks` to 0.009-0.098 pp. `lower foreign demand` is off by up to
  0.150 pp and `CBO debt` by up to 0.476 pp on Poor -- with Rich now matching to
  0.006 pp, since the intercept channel that produced the old sign flip is switched
  off here (see the `[SW23]` section).

---

## [SW26] MATLAB FIGURES FOR THE RULE-2 BASELINE

### How to run it

```
julia --project -t auto scripts/export_matlab_rule2_sw.jl      # ~1 min, no solve
matlab -batch "run('output/matlab/Plots_series_v0_rule2_sw.m')"    # 7 figures
matlab -batch "run('output/matlab/Plot_new_series_rule2_sw.m')"    # the 6-panel
```

`SW_TAG` and `SW_OUTDIR` are environment overrides, so any other saved run can be
exported the same way, e.g. `SW_TAG=cf_cbo_debt_rule2_forced_accom
SW_OUTDIR=matlab_cbo julia --project scripts/export_matlab_rule2_sw.jl`.

### What was added

| file | role |
|---|---|
| `scripts/export_matlab_rule2_sw.jl` | reads `rule2_epsm.jld2`, writes the three `.mat` files. Read-only, no solve. |
| `output/matlab_rule2/` | **new folder.** `output/matlab_sw/` holds the PRE-rule-2 baseline and is left alone, so the old figures stay reproducible. |
| `output/matlab/Plots_series_v0_rule2_sw.m` | copy of `Plots_series_v0_sw_all.m`; **only paths changed** (model dir and figure names). Column mapping, data series and layouts untouched. |
| `output/matlab/Plot_new_series_rule2_sw.m` | `Plot_new_series_sw.m` extended from 4 to 6 panels. |
| `output/figures_rule2/` | 8 figures. |

`model_new_series_sw.mat` now carries `s_index`, `piw`, `wedge`, `kappa_w`, `phi_pi`,
`pi_bar`, `i_taylor_ss`, `tag`, and the two new series `pi_solved`, `pi_endog`, `epsm`.

### The two new panels

**(f) `epsm`** is the Taylor-rule residual backed out from the rule-2 Stage-1 solve —
exactly the shock sequence the baseline was solved with. It peaks at **-4.61 pp
annualised at 2022Q3**.

**(e) `pi_endog`** re-runs the inflation recursion on the SOLVED baseline `r` path with
`epsm = 0` throughout. It answers *"what inflation would this rule deliver, given this
real allocation, with no monetary shocks?"*

> **This is a decomposition of the inflation equation holding `r` fixed, not a
> re-solved general-equilibrium counterfactual.** With `epsm = 0` the real side would
> itself move. The figure and the script header both say so. A GE version would need a
> full re-solve (~10-20 min) and is not what is plotted.

`j = 1` is fed from the data in both series (`pi_1` is not produced by the recursion),
so the two lines coincide at `j = 1` by construction and separate only from `j = 2`.

**The headline number falls straight out of panel (e).** At the 2022Q2 peak:

| | annualised % |
|---|---|
| `pi` solved, with shocks | **9.0600** |
| `pi` with `epsm = 0` | **2.0326** |
| target `pi_bar` | 2.0000 |

so the model delivers **0.47%** of the peak deviation from target on its own. That is
the "<1% of the 2021-22 inflation without monetary shocks" result, now plotted rather
than asserted.

Decomposition over the first year and a half (annualised %):

| j | date | `pi` solved | `pi`, no shock | `epsm` (pp) |
|---|---|---|---|---|
| 1 | 2020Q2 | 0.65 | 0.65 | +0.72 |
| 4 | 2021Q1 | 2.62 | 1.52 | +2.01 |
| 6 | 2021Q3 | 5.39 | 1.77 | -0.27 |
| 8 | 2022Q1 | 8.54 | 1.98 | -2.71 |
| 9 | 2022Q2 | **9.06** | **2.03** | -4.19 |
| 10 | 2022Q3 | 8.20 | 2.05 | **-4.61** |
| 12 | 2023Q1 | 4.98 | 2.01 | -3.36 |
| 20 | 2025Q1 | 2.00 | 1.87 | -0.08 |

### Checks the exporter runs

- `max|path.pi_ - pi_solved|` = **0.000e+00** — the saved path carries the solved
  inflation, not the fed path, so `piw` built from `path.pi_` and `path.w` is right.
- `s_index` z-spread across all 800 periods = **2.220e-16**, i.e. rationing is exactly
  proportional and `s_j` is well defined.
- Feeding the true `epsm` back through the recursion reproduces the solved `pi` to
  **0.000e+00**, which validates `pi_endog` as the same computation with the shock
  zeroed.

### Figure fixes made while wiring this up

Three were needed, all in the new copy only:

1. `sprintf` with `\%` in a `title` truncated the string — MATLAB's `sprintf` consumes
   the escape before the LaTeX interpreter sees it. Replaced with a **cell-array
   title**, which takes two lines natively and does no escape processing.
2. `\small` is **not** in MATLAB's LaTeX subset and silently kills the whole title
   (`\%` is fine — the `ylabel`s use it). Removed.
3. Panel (b)'s original `'Location','best'` put the legend over the inflation peak.
   Set to `southeast` with `ylim([-2.5 11])`; panel (e) uses `northeast` with
   `ylim([0 14])` so the 9.06 peak clears the box.

### Note on `output/figures_sw/`

The `.mat` inputs in `output/matlab_sw/` are untouched (still dated Sep 1). The PNGs in
`output/figures_sw/` were rewritten at 14:28 and 14:34 on 2026-09-11, which is before
the first MATLAB call in this session (14:41), and both new scripts print only to
`output/figures_rule2/`. Not written by the `[SW26]` work.

### Colour convention (applies to every `_rule2` figure)

Red was doing double duty in the first draft: it marked the data in
`Plots_series_v0_rule2_sw.m` but a *model* series in `Plot_new_series_rule2_sw.m`
(price inflation in panel (b), the baseline run in panel (e)). Fixed so that colour
carries exactly one meaning:

| colour | style | meaning |
|---|---|---|
| **blue** `[0.00 0.20 0.80]` | solid | the baseline model run — **the same run in every panel** |
| **green** `[0.00 0.50 0.25]` | solid | a different *variable* from that same run (wage inflation) |
| **orange** `[0.93 0.49 0.19]` | dash-dot | the counterfactual with the monetary shock switched off |
| **red** `[0.85 0.10 0.10]` | dashed | **data only** |
| grey `[0.25 0.25 0.25]` | dashed | reference lines (target, zero, unity) |

`Plot_new_series_rule2_sw.m` plots no data, so **no red appears in it at all**; the
palette is defined at the top of the file with `c_data` present but commented out, so
the reservation is explicit. `Plots_series_v0_rule2_sw.m` already followed
blue = Model / red dashed = Data and was left untouched.

Consequence worth noting: in panel (b) the two series are both from the baseline run,
so price inflation is blue (matching its colour in `f1_baseline_rule2`) and wage
inflation — a different variable, not a different run — is green.
