# ============================================================
# DEBT-FINANCED FISCAL STIMULUS, HETEROGENEITY, AND WELFARE
# STAGE 5 COUNTERFACTUAL 3 -- NO CAPITAL SUBSIDIES
# Author: Leonardo Barreto
# ============================================================
#     julia --project -t auto scripts/main_CF_no_capital_subsidy_sw.jl
#
# tau_k_path[1] is put BACK to its steady-state value (0.13793103) instead of -0.47.
# Everything else is identical to the epsm-fed baseline: qb_path, T_path, Tu_path,
# G_path, bstar_path_r, the lockdown block, epsm and i_taylor_ss.
#
# WHY THIS CAN REUSE THE BASELINE FINAL STEADY STATE
#   tau_k changes ONLY at j = 1.  tau_k_path[j] = e_init.tau_k for every j >= 2
#   already, and stationary_equilibrium_final is called with e_init.tau_k, so the
#   terminal steady state is literally the same object -- hence i_taylor_ss is
#   unchanged and there is NO intercept channel in this counterfactual.  The script
#   ASSERTS both rather than assuming them.
#
# TWO-STAGE WARM START
#   Stage A solves in s_mode = :flex.  The interior s equations become
#   s_j = s_flex(w_j, tau_n_j), so the whole tau_k change is absorbed by r, w and
#   tau_n with no NKPC block in the system.  Stage B then solves s_mode = :nkpc from
#   Stage A's solution.  Measured earlier on the tau_k experiment: this dropped the
#   NKPC warm-start residual from 160.87 to 0.141 (inf-norm 161 -> 0.835), which is
#   what identified the NKPC block as the entire source of the difficulty.
#
# [SW8] IS REVERTED HERE (_SW8_ENABLED[] = false).  The widen-then-grid-argmin
#   fallback makes the residual discontinuous in the unknowns and Anderson cannot
#   handle that; it is the prime suspect for the 5.6e20 blow-up.  With it off, a
#   period budget with no root in the bracket raises instead of silently returning a
#   non-solution.  [SW9] and [SW10] stay ON -- the r_k clamp is a genuine asymmetry
#   fix, not a counterfactual crutch.
#
# IF STAGE A FAILS, STOP.  That would mean the difficulty is NOT the NKPC block and
#   the fiscal calibration has to be fixed first.
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

_φ = 1/Base.MathConstants.φ
ProgressMeter.ijulia_behavior(:clear);

# ============================================================
# KNOBS
# ============================================================
KAPPA_W    = 0.05
PHI_PI     = 1.5
NKPC_SCALE = 240.0
WARM_FROM  = "reduced_epsm"           # converged baseline, for Stage A
EPSM_FROM  = "baseline_stage3_kw05"   # epsm HELD at the baseline path
SAVE_TAG   = "cf_no_capital_subsidy"
RUN_STAGE_B = true                    # false = stop after the :flex stage

# paper numbers for this counterfactual (fed-inflation version)
PAPER_CE = (mean = 2.02, poor = 3.77, middle = 2.27, rich = -1.69)

# [SW13] revert [SW8]; keep [SW9] and [SW10]
_SW8_ENABLED[] = false
# [SW13] two solves in this process (Stage A then Stage B), same e_init/t/g
_sw_runs_allowed[] = 2

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

# ============================================================
# Steady states -- BASELINE, reused unchanged
# ============================================================
println("loading cached steady states ..."); flush(stdout)
_inp = load_object(joinpath(ROOT, "output", "jld2_sw", "steady_states.jld2"))
e_init, e_final = _inp[1], _inp[2]
h = e_init.h; t = e_init.t; g = e_init.g
π_bar = e_init.π_
i_ss  = (1 + e_final.r) * (1 + π_bar) - 1

# ============================================================
# Time series, with the counterfactual tau_k
# ============================================================
b_target = e_init.b / y(e_init)
b_path = let
    T=700; H=100; ρB=0.08; bl = Array{Float64,1}(undef, T+H)
    bl[1]=b_target; bl[2]=b_target; bl[T:end] .= b_target
    for i in 3:T-1; bl[i] = bl[1]*ρB^(i-1) + (1-ρB^(i-1))*b_target; end
    bl
end
ts = make_time_series_sw(c, e_init, e_final, b_path, 0.287*4.0, g)

τ_k_cf = copy(ts.τ_k_path)
τ_k_cf[1] = e_init.τ_k          # <<< THE COUNTERFACTUAL

println("\n---- the counterfactual, and the claim that the final SS is untouched ----")
@printf("tau_k[1]  baseline = %+.8f   this run = %+.8f  (= e_init.tau_k)\n",
        ts.τ_k_path[1], τ_k_cf[1])
@printf("tau_k[2]  = %+.8f   e_init.tau_k = %+.8f   equal? %s\n",
        ts.τ_k_path[2], e_init.τ_k, ts.τ_k_path[2] == e_init.τ_k)
@assert ts.τ_k_path[2] == e_init.τ_k "tau_k deviates beyond j=1; the final SS would move"
@assert all(τ_k_cf[j] == e_init.τ_k for j in 2:length(τ_k_cf)) "tau_k not flat for j>=2"
@assert τ_k_cf == vcat(e_init.τ_k, ts.τ_k_path[2:end]) "tau_k_cf differs from baseline beyond j=1"
@printf("so tau_k_cf is FLAT at e_init.tau_k for every j -- identical to the path the\n")
@printf("final steady state was solved with, hence e_final and i_taylor_ss are unchanged.\n")
@printf("i_taylor_ss = %.10f  (%.4f%% ann)   r_f = %.4f%% ann\n",
        i_ss, 100*((1+i_ss)^4-1), 100*((1+e_final.r)^4-1))
println("EXPECTED, from the baseline: i_ss = 3.0994% ann, r_f = 1.0778% ann  -> NO intercept channel")
# everything else must be bit-identical to the baseline inputs
@printf("qb_path, T_path, Tu_path, G_path, bstar_path_r, d_path, covid_path: taken\n")
@printf("straight from make_time_series_sw, so identical to the baseline by construction.\n")

pars_h = (θ=c.θ, ν=c.ν, ies=c.ies, crra=c.crra, β=c.β, a_min=c.a_min, a_max=c.a_max,
          z_u=0.0, n_a=c.n_a, init_state=e_init.h.Pss, z_scale=z_scale)
h_path, _, _ = h_path_transition_mckay_init(ROOT, ts.Pee_path, ts.Puu_path; pars_h...)

# ============================================================
# epsm HELD at the baseline path, and the Stage A warm start
# ============================================================
println("\nloading ", EPSM_FROM, " for epsm ..."); flush(stdout)
_o = load_object(joinpath(ROOT, "output", "jld2_sw", EPSM_FROM * ".jld2"))
r_s3 = copy(_o[3].r); π_s3 = copy(_o[3].π_)
_o = nothing; GC.gc()
epsm = epsm_baseline_sw(r_s3, π_s3, i_ss, π_bar, PHI_PI)
@printf("max|epsm| = %.8f   epsm at 2022Q2 (j=9) = %+.4f pp ann\n",
        maximum(abs, epsm), 400*epsm[9])

println("\nloading ", WARM_FROM, " for the Stage A warm start ..."); flush(stdout)
_o = load_object(joinpath(ROOT, "output", "jld2_sw", WARM_FROM * ".jld2"))
r0, s0 = warm_start_sw(_o[3], e_init, t; method = 1)
π_base   = copy(_o[4])
r_base   = copy(_o[3].r)
τn_base  = copy(_o[3].τ_n)
path_base_v1 = copy(_o[3].v[1])     # for a like-for-like CE comparison if wanted
_o = nothing; GC.gc()

n_per = length(τ_k_cf)
n_r   = _n_r_sw(n_per, 1)
n_s   = n_per - 2

# shared solver arguments
DY_ss = union_aggregates_(e_init.h, g, e_init.ws; s = 1.0, w = e_init.w, τ_n = e_init.τ_n,
                          n_ss = e_init.ws.n_pol[1,:])
wY_ss = e_init.w * DY_ss.Y

function warm_residual(x, s_mode)
    path_tmp = _initialize_path(e_init, e_final; qb_path = ts.qb_path)
    π_tmp = Vector{Float64}(undef, n_per)
    F = similar(x)
    _residuals_sw_c!(F, x; π_vec = π_tmp, epsm_path = epsm, i_taylor_ss = i_ss, phi_pi = PHI_PI,
        path = path_tmp, e_init = e_init, e_final = e_final, τ_k_path = τ_k_cf,
        T_path = ts.T_path, Tu_path = ts.Tu_path, G_path = ts.G_path,
        qb_path = ts.qb_path, qbstar_path = ts.bstar_path_r, h_path = h_path,
        π_path = ts.π_path, d_path = ts.d_path, covid_path = ts.covid_path,
        Γ_n = ts.Γ_n, Γ_c = ts.Γ_c, χ = c.χ, t = t, g = g, method = 1,
        kappa_w = KAPPA_W, wY_ss = wY_ss, nkpc_scale = NKPC_SCALE, s_mode = s_mode)
    return (asset = maximum(abs, view(F, 1:n_r)),
            sblock = maximum(abs, view(F, (n_r+1):(n_r+n_s))),
            inf = maximum(abs, F))
end

# ============================================================
# STAGE A -- s_mode = :flex
# ============================================================
println("\n", "="^100)
println("STAGE A -- s_mode = :flex   (tau_k change absorbed by r, w, tau_n; no NKPC block)")
println("="^100)
wA = warm_residual(vcat(r0, s0), :flex)
@printf("warm start from %s:  asset max|R| = %.6e   s-block max|R| = %.6e   inf-norm = %.6e\n",
        WARM_FROM, wA.asset, wA.sblock, wA.inf); flush(stdout)

pathA, solA, πA, infoA = solve_transition_sw_c(e_init, e_final;
    τ_k_path=τ_k_cf, T_path=ts.T_path, Tu_path=ts.Tu_path, G_path=ts.G_path,
    qb_path=ts.qb_path, qbstar_path=ts.bstar_path_r, h_path=h_path, π_path=ts.π_path,
    d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
    method=1, kappa_w=KAPPA_W, epsm_path=epsm, phi_pi=PHI_PI, i_taylor_ss=i_ss,
    nkpc_scale=NKPC_SCALE, s_mode=:flex, r_path=r0, s_path=s0)

@printf("\nSTAGE A: iterations = %d   converged = %s   wall = %.1f min\n",
        solA.iterations, solA.f_converged, infoA.wall/60)
guard_report_sw()
if !solA.f_converged
    println("\n*** STAGE A DID NOT CONVERGE ***")
    println("The flexible-wage stage carries NO NKPC block, so the difficulty is NOT the")
    println("wage curve: it is the tau_k shock itself against the fed debt path. The fiscal")
    println("calibration has to be fixed before this counterfactual is attempted again.")
    println("STOPPING as instructed. Nothing written.")
    exit(0)
end
println("\nSTAGE A converged.")
@printf("tau_n[1] = %.8f   tau_n[2] = %.8f   (baseline %.8f, %.8f)\n",
        pathA.τ_n[1], pathA.τ_n[2], τn_base[1], τn_base[2])

# ============================================================
# STAGE B -- s_mode = :nkpc, warm started from Stage A
# ============================================================
if !RUN_STAGE_B
    println("\nRUN_STAGE_B = false; stopping after Stage A. Nothing written.")
    exit(0)
end

rA, sA = warm_start_sw(pathA, e_init, t; method = 1)
println("\n", "="^100)
println("STAGE B -- s_mode = :nkpc, warm started from Stage A")
println("="^100)
wB = warm_residual(vcat(rA, sA), :nkpc)
@printf("warm start from STAGE A:  asset max|R| = %.6e   nkpc max|R| = %.6e   inf-norm = %.6e\n",
        wB.asset, wB.sblock, wB.inf)
wB0 = warm_residual(vcat(r0, s0), :nkpc)
@printf("for comparison, straight from %s: asset %.6e   nkpc %.6e   inf-norm %.6e\n",
        WARM_FROM, wB0.asset, wB0.sblock, wB0.inf)
@printf("=> the two-stage warm start changes the NKPC warm-start residual %.6e -> %.6e\n",
        wB0.sblock, wB.sblock); flush(stdout)

pathB, solB, πB, infoB = solve_transition_sw_c(e_init, e_final;
    τ_k_path=τ_k_cf, T_path=ts.T_path, Tu_path=ts.Tu_path, G_path=ts.G_path,
    qb_path=ts.qb_path, qbstar_path=ts.bstar_path_r, h_path=h_path, π_path=ts.π_path,
    d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
    method=1, kappa_w=KAPPA_W, epsm_path=epsm, phi_pi=PHI_PI, i_taylor_ss=i_ss,
    nkpc_scale=NKPC_SCALE, s_mode=:nkpc, r_path=rA, s_path=sA)

@printf("\nSTAGE B: iterations = %d   converged = %s   wall = %.1f min\n",
        solB.iterations, solB.f_converged, infoB.wall/60)
guard_report_sw()

if !solB.f_converged
    println("\n*** STAGE B DID NOT CONVERGE, but STAGE A DID ***")
    println("That localises the problem to the NKPC block rather than the tau_k shock size.")
    println("Nothing written.")
    exit(0)
end

# ============================================================
# Diagnostics
# ============================================================
println("\n", "="^100)
println("DIAGNOSTICS -- no capital subsidies")
println("="^100)

gate = gate_lab_tax_sw(e_init, pathB, g)          # includes j = 1
gate_ok = print_gate_sw(gate)
branch_ok = period0_branch_sw(e_init, pathB, g, c.ν, ts.Γ_n)

println("\n---- tau_n: the number to watch ----")
@printf("initial SS                 = %.8f\n", e_init.τ_n)
@printf("baseline   j=1 = %.8f   j=2 = %.8f\n", τn_base[1], τn_base[2])
@printf("this run   j=1 = %.8f   j=2 = %.8f\n", pathB.τ_n[1], pathB.τ_n[2])
@printf("change vs baseline  j=1 = %+.8f   j=2 = %+.8f\n",
        pathB.τ_n[1]-τn_base[1], pathB.τ_n[2]-τn_base[2])
@printf("vs the initial SS   j=1 = %+.4f%%   j=2 = %+.4f%%\n",
        100*(pathB.τ_n[1]/e_init.τ_n - 1), 100*(pathB.τ_n[2]/e_init.τ_n - 1))
@printf("retention above 1 at j=1? %s\n", pathB.τ_n[1] > 1 ? "YES" : "no")
println("Removing a subsidy worth 7.70% of annual GDP REDUCES outlays, so with borrowing")
println("fed and unchanged there is MORE spare resource and tau_n should rise FURTHER")
println("above the baseline's 1.00222.")

println("\n---- pi and r against the baseline (annualised %) ----")
@printf("%5s %12s %12s %12s   %12s %12s %12s\n",
        "j", "pi base", "pi cf", "diff pp", "r base", "r cf", "diff pp")
for j in vcat(collect(1:8), [12, 20, 40, 100, 400])
    @printf("%5d %12.4f %12.4f %+12.4f   %12.4f %12.4f %+12.4f\n", j,
            100*((1+π_base[j])^4-1), 100*((1+πB[j])^4-1), 400*(πB[j]-π_base[j]),
            100*((1+r_base[j])^4-1), 100*((1+pathB.r[j])^4-1), 400*(pathB.r[j]-r_base[j]))
end
@printf("\nmax |pi_cf - pi_base| over j=2..60 = %.4f ann pp\n",
        400*maximum(abs.(πB[2:60] .- π_base[2:60])))
@printf("terminal pi = %.6f%% ann (should be pi_bar = %.6f%%)\n",
        100*((1+πB[end])^4-1), 100*((1+π_bar)^4-1))

ce = ce_summary_sw(e_init, pathB)
print_ce_sw(ce; label = "(no capital subsidies)", paper = PAPER_CE)
println("\nfull decile profile (paper convention: deciles of the initial-SS wealth marginal):")
for (k, v) in enumerate(ce.deciles)
    @printf("  bin %2d  %8.4f\n", k, v)
end

# ============================================================
# Save (converged AND gate AND branch only)
# ============================================================
if solB.f_converged && gate_ok && branch_ok
    save_object(joinpath(ROOT, "output", "jld2_sw", SAVE_TAG * ".jld2"),
                (e_init, e_final, pathB, πB))
    println("\nwrote output/jld2_sw/", SAVE_TAG, ".jld2")
else
    println("\nNOT written.  converged = ", solB.f_converged,
            "   gate = ", gate_ok, "   period-0 branch = ", branch_ok)
end
println("\nDONE")
