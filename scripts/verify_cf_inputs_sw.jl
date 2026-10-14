# ============================================================
# [SW11] VERIFY THE RECONSTRUCTED COUNTERFACTUAL INPUTS
# ============================================================
#   julia --project -t auto scripts/verify_cf_inputs_sw.jl lower_foreign
#   julia --project -t auto scripts/verify_cf_inputs_sw.jl cbo_debt
#
# WHY THIS EXISTS
#   The runners for cf_lower_foreign.jld2 and cf_cbo_debt.jld2 were lost; the knob
#   values in src/counterfactuals_sw.jl were reconstructed by inverting the
#   solver's own mapping out of the saved results.  This script closes the loop
#   WITHOUT paying for a full solve: it rebuilds the inputs exactly as the runner
#   does, reconstructs the unknown vector x from the SAVED path, and evaluates the
#   reduced residual once.
#
#   If the reconstruction is right, the residual comes back at the converged value
#   the run reported (8.261e-06 for lower_foreign, 9.919e-06 for cbo_debt).  If any
#   input is wrong, the residual is orders of magnitude larger.
#
#   e_init and e_final are taken FROM the saved file, so no steady state is
#   re-solved here and nothing is written.  This calls _residuals_sw_c! directly
#   and never solve_transition_sw_c, so it does not touch the one-run-per-process
#   guard.
# ============================================================

using LinearAlgebra
using ProgressMeter
using Roots
using NLsolve
using StructArrays
using StatsBase
using Tullio
using Polyester
using PrettyTables
using CairoMakie
using LaTeXStrings
using DelimitedFiles
using XLSX
using CSV
using DataFrames
using ColorTypes
using Printf
using JLD2
using MAT

const ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(ROOT, "src", "types.jl"))
include(joinpath(ROOT, "src", "preferences.jl"))
include(joinpath(ROOT, "src", "technology_and_fiscal.jl"))
include(joinpath(ROOT, "src", "calibration.jl"))
include(joinpath(ROOT, "src", "auxiliary_functions.jl"))
include(joinpath(ROOT, "src", "household_problem_and_aggregation.jl"))
include(joinpath(ROOT, "src", "stationary_equilibria.jl"))
include(joinpath(ROOT, "src", "transition.jl"))
include(joinpath(ROOT, "src", "union_block.jl"))
include(joinpath(ROOT, "src", "transition_sw.jl"))
include(joinpath(ROOT, "src", "counterfactuals_sw.jl"))
include(joinpath(ROOT, "src", "statistics.jl"))
include(joinpath(ROOT, "src", "time_series.jl"))

const _TOL = 1e-12
const _ZERO_FTOL = 1e-5
const _MAX_ITERS = 20_000
const _ZERO_MAX_ITERS = 450

which = isempty(ARGS) ? "lower_foreign" : ARGS[1]
which in ("lower_foreign", "cbo_debt", "kw30") ||
    error("ARGS[1] must be lower_foreign, cbo_debt or kw30")

PHI_PI = 1.5
if which == "lower_foreign"
    KAPPA_W, NKPC_SCALE, EXPECTED, SAVE_TAG = 0.05, 240.0, 8.261e-06, "cf_lower_foreign"
elseif which == "cbo_debt"
    KAPPA_W, NKPC_SCALE, EXPECTED, SAVE_TAG = 0.05, 240.0, 9.919e-06, "cf_cbo_debt"
else
    # kw30: the kappa_w robustness run.  No inf-norm was recorded for it, so the
    # bar is the solver's own ftol -- still decisive, since a wrong input puts the
    # residual orders of magnitude above it.
    KAPPA_W, NKPC_SCALE, EXPECTED, SAVE_TAG = 0.30, 41.4, 1.0e-05, "reduced_run3_kw30"
end

# ------------------------------------------------------------
# Calib and constructors, verbatim from main_JMP_reduced.jl (only z_scale is used
# here -- e_init/e_final come from the saved file).
# ------------------------------------------------------------
Base.@kwdef struct Calib
    θ::Float64 = 1.0
    ν::Float64 = 1.0
    ies::Float64 = 1.0
    crra::Float64 = 4.0
    β::Float64 = 0.9880431080390226
    a_min::Float64 = -0.02
    a_max::Float64 = 120.0
    n_a::Int = 1000
    z_bar::Int = 12
    Pee::Float64 = 0.975
    Puu::Float64 = 0.4
    Λ_t::Float64 = 0.0
    Λ_tu::Float64 = 0.005
    Λ_a::Float64 = 5.06 * 4.0
    Λ_g::Float64 = 0.05
    Λ_f::Float64 = 0.04
    Λ_b::Float64 = 1.05 * 4.0
    Λ_bstar::Float64 = 0.31 * 4.0
    Λ_rk::Float64 = 0.33
    r_target::Float64 = (1.0 + 0.01)^(1/4) - 1.0
    π_bar::Float64 = (1.0 + 0.02)^(1/4) - 1.0
    mc::Float64 = 5.0 * 4.0
    τ::Float64 = 0.181
    nu_target::Float64 = 0.7
    A::Float64 = 0.035
    μ::Float64 = 1.0
    d_::Float64 = 1.0
    Γ_n::Float64 = 1.0
    Γ_c::Float64 = 1.0
    covid::Float64 = 0.0
    χ::Float64 = 0.44
end
c = Calib()

_, _, z_scale = calibration_unemploy_mckay_(ROOT, c.Pee, c.Puu, 0.0)

# ------------------------------------------------------------
# Baseline final SS (for i_ss_base) and the Stage 3 epsm
# ------------------------------------------------------------
_inp = load_object(joinpath(ROOT, "output", "jld2_sw", "steady_states.jld2"))
e_final_base = _inp[2]

println("loading baseline_stage3_kw05 for epsm ..."); flush(stdout)
_o = load_object(joinpath(ROOT, "output", "jld2_sw", "baseline_stage3_kw05.jld2"))
r_s3 = copy(_o[3].r); π_s3 = copy(_o[3].π_)
_o = nothing; GC.gc()

# ------------------------------------------------------------
# The saved counterfactual
# ------------------------------------------------------------
println("loading ", SAVE_TAG, ".jld2 ..."); flush(stdout)
o = load_object(joinpath(ROOT, "output", "jld2_sw", SAVE_TAG * ".jld2"))
e_init, e_final, path_saved, π_saved = o[1], o[2], o[3], o[4]
h = e_init.h; t = e_init.t; g = e_init.g
π_bar = e_init.π_

i_ss_base = (1 + e_final_base.r) * (1 + π_bar) - 1
i_ss      = (1 + e_final.r)      * (1 + π_bar) - 1
# kw30 sets epsm = 0 (FEED_EPSM = false): pi is the model's own inflation there.
epsm      = which == "kw30" ? zeros(length(r_s3)) :
                              epsm_baseline_sw(r_s3, π_s3, i_ss_base, π_bar, PHI_PI)

# ------------------------------------------------------------
# Rebuild the inputs exactly as the runner does
# ------------------------------------------------------------
b_target = e_init.b / y(e_init)
b_path = let
    T=700; H=100; ρB=0.08; bl = Array{Float64,1}(undef, T+H)
    bl[1]=b_target; bl[2]=b_target; bl[T:end] .= b_target
    for i in 3:T-1; bl[i] = bl[1]*ρB^(i-1) + (1-ρB^(i-1))*b_target; end
    bl
end

if which == "kw30"
    # Baseline paths untouched -- only kappa_w, nkpc_scale and epsm differ.
    Λ_bstar_f = 0.287 * 4.0
    ts = make_time_series_sw(c, e_init, e_final, b_path, Λ_bstar_f, g)
    qb_use  = ts.qb_path
    qbs_use = ts.bstar_path_r
elseif which == "lower_foreign"
    Λ_bstar_f = 0.18 * 4.0
    ts = make_time_series_sw(c, e_init, e_final, b_path, Λ_bstar_f, g)
    qb_use  = ts.qb_path
    qbs_use = cf_qbstar_lower_foreign_sw(ts; drop = 0.10, jlast = 20)
else
    Λ_bstar_f = 0.287 * 4.0
    ts = make_time_series_sw(c, e_init, e_final, b_path, Λ_bstar_f, g)
    qb_use, qbs_use = cf_paths_cbo_debt_sw(ts, e_final;
        start_index = 19, end_index = 143, target = 1.84 * 4.0,
        foreign_to_gdp = Λ_bstar_f, foreign_mode = :as_run)
end

# ------------------------------------------------------------
# CHECK 1 -- do the rebuilt inputs match what the saved path implies?
# Invert the solver's mapping out of the stored path and compare elementwise.
# ------------------------------------------------------------
ibar = (1 + e_init.r) * (1 + g.π_bar) - 1
n_per = length(path_saved.r)
qbs_implied = [ path_saved.bstar[i] / path_saved.b[i] for i in 1:n_per ]
qb_implied  = Vector{Float64}(undef, n_per)
for i in 1:n_per
    qlag = (1/(1+path_saved.r[i])) * ((ibar + g.ϑ)/(1+path_saved.π_[i])
                                      + (1-g.ϑ)/(1+path_saved.π_[i]) * path_saved.q[i])
    qb_implied[i] = path_saved.b[i] * qlag / path_saved.y[i]
end

# j = 1 is excluded: bstar_path_r[1] and qb_path[1] are LEVELS, not ratios, so the
# inversion above does not apply there.
dqb  = maximum(abs.(qb_use[2:end]  .- qb_implied[2:end]))
dqbs = maximum(abs.(qbs_use[2:end] .- qbs_implied[2:end]))
println("\n==================== CHECK 1: rebuilt inputs vs saved path ====================")
@printf("max |qb_rebuilt - qb_implied|         = %.6e\n", dqb)
@printf("max |qbstar_rebuilt - qbstar_implied| = %.6e\n", dqbs)
println(dqb < 1e-10 && dqbs < 1e-10 ?
        "CHECK 1 PASSED -- the knob reconstruction is exact." :
        "*** CHECK 1 FAILED -- a knob is wrong; see the worst offenders below ***")
if !(dqb < 1e-10 && dqbs < 1e-10)
    for j in 2:min(n_per, 160)
        if abs(qb_use[j]-qb_implied[j]) > 1e-10 || abs(qbs_use[j]-qbs_implied[j]) > 1e-10
            @printf("  j=%4d  qb %14.8f vs %14.8f   qbstar %12.8f vs %12.8f\n",
                    j, qb_use[j], qb_implied[j], qbs_use[j], qbs_implied[j])
        end
    end
end

# ------------------------------------------------------------
# CHECK 2 -- evaluate the reduced residual at the SAVED solution
# ------------------------------------------------------------
pars_h = (θ=c.θ, ν=c.ν, ies=c.ies, crra=c.crra, β=c.β, a_min=c.a_min, a_max=c.a_max,
          z_u=0.0, n_a=c.n_a, init_state=e_init.h.Pss, z_scale=z_scale)
h_path, _, _ = h_path_transition_mckay_init(ROOT, ts.Pee_path, ts.Puu_path; pars_h...)

r0, s0 = warm_start_sw(path_saved, e_init, t; method = 1)
x = vcat(r0, s0)
n_r = _n_r_sw(n_per, 1); n_s = n_per - 2

path = _initialize_path(e_init, e_final; qb_path = qb_use)
DY_ss = union_aggregates_(e_init.h, g, e_init.ws; s = 1.0, w = e_init.w, τ_n = e_init.τ_n,
                          n_ss = e_init.ws.n_pol[1,:])
wY_ss = e_init.w * DY_ss.Y
π_vec = Vector{Float64}(undef, n_per)

F = similar(x)
_residuals_sw_c!(F, x; π_vec = π_vec, epsm_path = epsm, i_taylor_ss = i_ss, phi_pi = PHI_PI,
                 path = path, e_init = e_init, e_final = e_final, τ_k_path = ts.τ_k_path,
                 T_path = ts.T_path, Tu_path = ts.Tu_path, G_path = ts.G_path,
                 qb_path = qb_use, qbstar_path = qbs_use, h_path = h_path,
                 π_path = ts.π_path, d_path = ts.d_path, covid_path = ts.covid_path,
                 Γ_n = ts.Γ_n, Γ_c = ts.Γ_c, χ = c.χ, t = t, g = g, method = 1,
                 kappa_w = KAPPA_W, wY_ss = wY_ss, nkpc_scale = NKPC_SCALE, s_mode = :nkpc)

println("\n============ CHECK 2: reduced residual at the SAVED solution ============")
@printf("asset  block max|R| = %.6e\n", maximum(abs, view(F, 1:n_r)))
@printf("nkpc   block max|R| = %.6e   (= %.6e UNSCALED, i.e. / nkpc_scale = %.1f)\n",
        maximum(abs, view(F, (n_r+1):(n_r+n_s))),
        maximum(abs, view(F, (n_r+1):(n_r+n_s))) / NKPC_SCALE, NKPC_SCALE)
@printf("inf-norm            = %.6e\n", maximum(abs, F))
@printf("%s = %.6e\n",
        which == "kw30" ? "solver ftol (no inf-norm was recorded for this run)" :
                          "reported at convergence", EXPECTED)
if which == "kw30"
    println("NOTE nkpc_scale only ROW-SCALES the NKPC block, so it does not change the")
    println("solution set and this check BOUNDS it rather than identifying it.  kappa_w")
    println("IS identified: it enters the wage curve nonlinearly, so a wrong kappa_w")
    println("would put the unscaled residual orders of magnitude above ftol/nkpc_scale.")
end
println(maximum(abs, F) < 10 * EXPECTED ?
        "CHECK 2 PASSED -- these inputs reproduce the saved solution." :
        "*** CHECK 2 FAILED -- the residual is far from converged; an input is wrong ***")

# ------------------------------------------------------------
# CHECK 3 -- does the pi recursion reproduce the stored pi?
# ------------------------------------------------------------
dπ = maximum(abs.(π_vec .- π_saved))
println("\n============ CHECK 3: pi recursion vs stored pi ============")
@printf("max |pi_recursion - pi_saved| = %.6e   (= %.6f ann pp)\n", dπ, 400*dπ)
println(dπ < 1e-9 ? "CHECK 3 PASSED -- epsm and i_taylor_ss are right." :
                    "*** CHECK 3 FAILED -- epsm or i_taylor_ss is wrong ***")

guard_report_sw()
println("\nDONE (nothing written)")
