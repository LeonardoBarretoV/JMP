# ============================================================
# DEBT-FINANCED FISCAL STIMULUS, HETEROGENEITY, AND WELFARE
# STAGE 5 COUNTERFACTUAL 2 -- CBO DEBT PROJECTIONS
# Author: Leonardo Barreto
# ============================================================
# Reproduces output/jld2_sw/cf_cbo_debt.jld2.
#
# RUN IT IN A FRESH JULIA SESSION.  _assert_globals_sw enforces one solve per
# process; a second solve in the same session fails the guard by design.
#
#     julia --project -t auto scripts/main_CF_cbo_debt_sw.jl
#
# Expected, from the 2026-09-08 run:  162 Anderson iterations, inf-norm 9.919e-06,
# 47.4 min, gate passed, guards inert (min x = 2.0e-03, six orders above the floor).
#
# THIS IS THE RUN THAT SHOWS THE INTERCEPT CHANNEL.  Rebuilding i_taylor_ss off the
# new r_final raises the Taylor intercept 3.0994% -> 3.4495% ann, which lowers
# inflation along the whole transition by a predicted -0.70pp, against an observed
# shortfall of -0.37 to -0.51pp; the higher real-rate path offsets +0.19 to +0.33pp.
# The two channels oppose and the intercept wins.  The script prints both.
#
# Everything above the KNOBS block is verbatim from main_JMP_reduced.jl so the
# two share one setup.
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
using Revise
using CairoMakie
using LaTeXStrings
using Polyester
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
include(joinpath(ROOT, "src", "union_block.jl"))          # [SW] sticky-wage track ONLY
include(joinpath(ROOT, "src", "transition_sw.jl"))        # [SW] sticky-wage track ONLY
include(joinpath(ROOT, "src", "counterfactuals_sw.jl"))   # [SW11] Stage 5 machinery
include(joinpath(ROOT, "src", "statistics.jl"))
include(joinpath(ROOT, "src", "time_series.jl"))

const _TOL = 1e-12
const _ZERO_FTOL = 1e-5
const _MAX_ITERS = 20_000
const _ZERO_MAX_ITERS = 450

_φ = 1/Base.MathConstants.φ

ProgressMeter.ijulia_behavior(:clear);

# ============================================================
# Parameters / Calibration
# ============================================================

Base.@kwdef struct Calib
    # Preferences / household
    θ::Float64 = 1.0 # Labor disutility parameter
    ν::Float64 = 1.0 # Frisch elasticity
    ies::Float64 = 1.0 # Intertemporal elasticity of substitution
    crra::Float64 = 4.0 # Risk aversion coefficient
    β::Float64 = 0.9880431080390226  # Discount factor

    # Assets grid
    a_min::Float64 = -0.02
    a_max::Float64 = 120.0
    n_a::Int = 1000

    # Income / stimulus
    z_bar::Int = 12 # Stimulus checks criteria

    # Unemployment discretization
    Pee::Float64 = 0.975 # Probability of remaining employed
    Puu::Float64 = 0.4 # Probability of remaining unemployed

    # Targets
    Λ_t::Float64 = 0.0 # Stimulus checks-to-GDP
    Λ_tu::Float64 = 0.005 # Unemployment insurance-to-GDP
    Λ_a::Float64 = 5.06 * 4.0 # Assets-to-GDP
    Λ_g::Float64 = 0.05 # Government expenditure-to-GDP
    Λ_f::Float64 = 0.04 # Tax revenue from firms-to-GDP
    Λ_b::Float64 = 1.05 * 4.0 # Federal debt-to-GDP
    Λ_bstar::Float64 = 0.31 * 4.0 # External debt-to-GDP
    Λ_rk::Float64 = 0.33 # Capital share

    # Rates
    r_target::Float64 = (1.0 + 0.01)^(1/4) - 1.0 # Real interest rate
    π_bar::Float64 = (1.0 + 0.02)^(1/4) - 1.0 # Inflation rate

    # Markup / pricing
    mc::Float64 = 5.0 * 4.0 # Debt maturity

    # Tax / labor
    τ::Float64 = 0.181 # Progressivity parameter
    nu_target::Float64 = 0.7 # Adjusted Frisch elasticity

    # Technology
    A::Float64 = 0.035 # Total factor productivity
    μ::Float64 = 1.0 #Markup

    # Initial steady-state / COVID block
    d_::Float64 = 1.0 # Productivity parameter
    Γ_n::Float64 = 1.0 # Pandemic labor supply constraint indicator
    Γ_c::Float64 = 1.0 # Pandemic consumption supply constraint indicator
    covid::Float64 = 0.0 # Pandemic indicator
    χ::Float64 = 0.44 # Foreign demand elasticity
end

# Some ratios from the steady state
Λ_bd(c::Calib) = c.Λ_b - c.Λ_bstar
Λ_wn(c::Calib) = 1 - c.Λ_rk
Λ_k(c::Calib)  = c.Λ_a - Λ_bd(c)
i_bar(c::Calib) = (1 + c.r_target) * (1 + c.π_bar) - 1
ϑ(c::Calib) = (1 + i_bar(c)) / c.mc - i_bar(c)

nu_from_tax(c::Calib) = c.nu_target / (1 - c.τ + c.nu_target * c.τ)

# ============================================================
# Constructors for the model
# ============================================================

function make_household(c::Calib)
    v = GHH(θ = c.θ, ν = c.ν)
    P, z_vals, z_scale = calibration_unemploy_mckay_(ROOT, c.Pee, c.Puu, 0.0)   # [Z-NORM]
    u = EZ(ies = c.ies, ra = c.crra, β = c.β)

    a_grid = grid(start = c.a_min, stop = c.a_max, length = c.n_a, scale = :log)
    return Household(u = u, a_grid = a_grid, v = v, P = P, z_grid = z_vals), z_scale   # [Z-NORM]
end

function make_technology(c::Calib)
    δ = (c.Λ_rk - c.Λ_f) / Λ_k(c) - c.r_target
    α = c.Λ_rk
    t = CobbDouglasTechnologyMut(α, c.A^((1 - α)), δ, c.μ)
    return t, δ
end

function make_taxsys(c::Calib, h)
    ind_T  = zeros(length(h.z_grid), 1)
    ind_T[1:c.z_bar] .= 1

    ind_Tu = zeros(length(h.z_grid), 1)
    ind_Tu[1] = 1

    g = TaxSys(c.τ, ϑ(c), c.π_bar, ind_T, ind_Tu)
    return g
end

function τ_k_ss(c::Calib, δ::Float64)
    return c.Λ_f / ((c.r_target + δ) * Λ_k(c))
end

# ============================================================
# Initial steady state
# ============================================================
c = Calib()

h, z_scale = make_household(c)
t, δ = make_technology(c)
g = make_taxsys(c, h)
τk_ss = τ_k_ss(c, δ)
println("[Z-NORM] scale factor = $z_scale;  Z in initial SS = $(sum(h.Pss .* h.z_grid))")


# ============================================================
# KNOBS
# ============================================================
KAPPA_W      = 0.05        # wage Phillips slope.  NKPC_SCALE must be REMEASURED if this changes
PHI_PI       = 1.5         # Taylor coefficient
NKPC_SCALE   = 240.0       # measured at KAPPA_W = 0.05 (asset/NKPC diagonal ratio 1.01).
                           # It scales ~linearly in kappa_w: 41.4 at kappa_w = 0.30.

# --- the counterfactual itself -------------------------------------------------
LAMBDA_B_CBO   = 1.84 * 4.0    # terminal debt-to-GDP: 184% vs the baseline 120% (= c.Lambda_b + 0.15*4).
                               # RECOVERED from cf_cbo_debt.jld2: q*b/y(e_final) = 7.36 exactly.
                               # Equivalently c.Lambda_b + 0.79*4.0, which is how the pre-sticky-wage
                               # REPL version wrote it.
LAMBDA_BSTAR_F = 0.287 * 4.0   # terminal external debt-to-GDP: UNCHANGED from the baseline
TREND_START    = 19            # 2024Q4, the last observed quarter; its value (1.218*4) is left in place
TREND_END      = 143           # 2055Q4.  124 quarters, step 0.020064516
FOREIGN_MODE   = :as_run       # :as_run reproduces the saved file (divides by q*b/y);
                               # :level_b is the internally consistent variant (divides by b/y)
                               # and is a DIFFERENT run -- see src/counterfactuals_sw.jl

# --- warm start / monetary residual -------------------------------------------
WARM_FROM  = "reduced_epsm"           # converged baseline
EPSM_FROM  = "baseline_stage3_kw05"   # epsm is HELD at the baseline path (that is Stage 5)
USE_CACHE  = true                     # reuse the cached counterfactual steady states
SAVE_TAG   = "cf_cbo_debt"             # output/jld2_sw/<SAVE_TAG>.jld2

# paper numbers for the side-by-side (fed-inflation version of this experiment)
PAPER_CE = (mean = 2.51, poor = 4.00, middle = 2.62, rich = 0.50)


# ============================================================
# Steady states
# ============================================================
# e_init is the BASELINE initial steady state (unchanged by this counterfactual).
# e_final MUST be re-solved: terminal debt-to-GDP is what the counterfactual moves,
# and the resulting r_final is what drives the intercept channel.
ss_base = joinpath(ROOT, "output", "jld2_sw", "steady_states.jld2")
ss_cf   = joinpath(ROOT, "output", "jld2_sw", "steady_states_" * SAVE_TAG * ".jld2")

println("loading baseline steady states ..."); flush(stdout)
_inp = load_object(ss_base); e_init, e_final_base = _inp[1], _inp[2]
h = e_init.h; t = e_init.t; g = e_init.g      # required by _assert_globals_sw

if USE_CACHE && isfile(ss_cf)
    println("loading cached counterfactual final steady state ..."); flush(stdout)
    e_final = load_object(ss_cf)
else
    ī_ = (1 + e_init.r) * (1 + c.π_bar) - 1
    println("solving the counterfactual FINAL steady state at Lambda_b = ",
            LAMBDA_B_CBO, " ..."); flush(stdout)
    @time e_final = stationary_equilibrium_final(h, t, g, e_init.τ_k, c.Λ_t, c.Λ_tu, c.Λ_g,
                                                 LAMBDA_B_CBO, LAMBDA_BSTAR_F, ī_;
                                                 r_range = (0.0, 0.008), root_max = 1.0, verbose = true)
    save_object(ss_cf, e_final)
    println("wrote ", ss_cf)
end

@printf("\ninitial SS : r = %.10f (%.4f%% ann)   b/y_ann = %.6f\n",
        e_init.r, 100*((1+e_init.r)^4-1), e_init.b/y(e_init)/4)
@printf("baseline final SS : r = %.10f (%.4f%% ann)   q*b/y_ann = %.6f   bstar/b = %.8f\n",
        e_final_base.r, 100*((1+e_final_base.r)^4-1),
        e_final_base.q*e_final_base.b/y(e_final_base)/4, e_final_base.bstar/e_final_base.b)
@printf("THIS RUN final SS : r = %.10f (%.4f%% ann)   q*b/y_ann = %.6f   bstar/b = %.8f\n",
        e_final.r, 100*((1+e_final.r)^4-1),
        e_final.q*e_final.b/y(e_final)/4, e_final.bstar/e_final.b)
println("EXPECTED (cf_cbo_debt.jld2): r_f = 1.4211% ann, q*b/y_ann = 1.840000, bstar/b = 0.15267816")

# ============================================================
# Time series and the counterfactual paths
# ============================================================
b_target = e_init.b / y(e_init)
b_path = let
    T=700; H=100; ρB=0.08; bl = Array{Float64,1}(undef, T+H)
    bl[1]=b_target; bl[2]=b_target; bl[T:end] .= b_target
    for i in 3:T-1; bl[i] = bl[1]*ρB^(i-1) + (1-ρB^(i-1))*b_target; end
    bl
end

# make_time_series_sw is handed THIS run's e_final, so the qb_path fill is already
# 7.36; cf_paths_cbo_debt_sw then lays the 2025Q1-2055Q4 trend over j = 20..143 and
# rewrites the foreign share from j = 20 on.
ts = make_time_series_sw(c, e_init, e_final, b_path, LAMBDA_BSTAR_F, g)

qb_cf, qbstar_cf = cf_paths_cbo_debt_sw(ts, e_final;
    start_index = TREND_START, end_index = TREND_END,
    target = LAMBDA_B_CBO, foreign_to_gdp = LAMBDA_BSTAR_F, foreign_mode = FOREIGN_MODE)

println("\n---- debt and foreign paths: baseline vs this run ----")
@printf("%5s %14s %14s %14s %14s\n", "j", "qb base", "qb cf", "qbstar base", "qbstar cf")
for j in vcat(collect(17:22), [40, 60, 100, 140, 143, 144, 200, 800])
    @printf("%5d %14.8f %14.8f %14.8f %14.8f\n", j,
            ts.qb_path[j], qb_cf[j], ts.bstar_path_r[j], qbstar_cf[j])
end
println("EXPECTED: qb[20] = 4.89206452, qb[100] = 6.49722581, qb[143] = 7.36,")
println("          qbstar[20] = 0.23466575, qbstar[200] = 0.15597826")

pars_h = (θ=c.θ, ν=c.ν, ies=c.ies, crra=c.crra, β=c.β, a_min=c.a_min, a_max=c.a_max,
          z_u=0.0, n_a=c.n_a, init_state=e_init.h.Pss, z_scale=z_scale)
h_path, _, _ = h_path_transition_mckay_init(ROOT, ts.Pee_path, ts.Puu_path; pars_h...)

# ============================================================
# Warm start and the monetary residual
# ============================================================
n_per = length(qb_cf)
π_bar = e_init.π_

i_ss_base = (1 + e_final_base.r) * (1 + π_bar) - 1
i_ss      = (1 + e_final.r)      * (1 + π_bar) - 1

println("\nloading ", EPSM_FROM, " for the monetary residual ..."); flush(stdout)
_o = load_object(joinpath(ROOT, "output", "jld2_sw", EPSM_FROM * ".jld2"))
r_s3 = copy(_o[3].r); π_s3 = copy(_o[3].π_)
_o = nothing; GC.gc()
epsm = epsm_baseline_sw(r_s3, π_s3, i_ss_base, π_bar, PHI_PI)
@printf("max|epsm| = %.8f   epsm at 2022Q2 (j=9) = %+.4f pp ann\n",
        maximum(abs, epsm), 400*epsm[9])

println("\nloading ", WARM_FROM, " for the warm start ..."); flush(stdout)
_o = load_object(joinpath(ROOT, "output", "jld2_sw", WARM_FROM * ".jld2"))
r0, s0 = warm_start_sw(_o[3], e_init, t; method = 1)
π_base = length(_o) >= 4 ? copy(_o[4]) : copy(_o[3].π_)
_o = nothing; GC.gc()

ic = intercept_channel_sw(i_ss, i_ss_base, PHI_PI, e_final.r)
println("EXPECTED: d(i_ss) = 0.000854511, predicted d(pi) = -0.6975 pp ann")

# ============================================================
# Solve
# ============================================================
path_c, sol_c, π_endog, info_c = solve_transition_sw_c(e_init, e_final;
    τ_k_path=ts.τ_k_path, T_path=ts.T_path, Tu_path=ts.Tu_path, G_path=ts.G_path,
    qb_path=qb_cf, qbstar_path=qbstar_cf, h_path=h_path, π_path=ts.π_path,
    d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
    method=1, kappa_w=KAPPA_W, epsm_path=epsm, phi_pi=PHI_PI, i_taylor_ss=i_ss,
    nkpc_scale=NKPC_SCALE, s_mode=:nkpc, r_path=r0, s_path=s0)

# ============================================================
# Report
# ============================================================
print_pi_table_sw(ts.π_path, π_endog, π_bar; π_ref = π_base, label = "cbo debt")

# The intercept channel against the solved gap, j by j.  The residual is the
# r-path channel; independently, dpi/dr cumulates to about 2 in the recursion.
println("\n---- intercept channel vs the solved inflation gap (ann pp) ----")
@printf("%6s %16s %16s %16s %14s\n", "j", "solved d(pi)", "intercept", "residual = r path", "d(r) ann pp")
for j in (20, 40, 100, 143)
    solved = 400 * (π_endog[j] - π_base[j])
    icp    = 400 * ic.d_pi
    @printf("%6d %16.4f %16.4f %16.4f %14.4f\n", j, solved, icp, solved - icp,
            100*((1+path_c.r[j])^4 - 1) - 100*((1+e_final_base.r)^4 - 1))
end
println("EXPECTED: solved -0.4900 / -0.5052 / -0.3708 / -0.3852 against intercept -0.6975")

ce = ce_summary_sw(e_init, path_c)
print_ce_sw(ce; label = "(CBO debt projections)", paper = PAPER_CE)

gate = gate_lab_tax_sw(e_init, path_c, g)
gate_ok = print_gate_sw(gate)
branch_ok = period0_branch_sw(e_init, path_c, g, c.ν, ts.Γ_n)
guard_report_sw()

@printf("\n---- run record ----\niterations = %d   converged = %s   wall = %.1f min\n",
        sol_c.iterations, sol_c.f_converged, info_c.wall/60)
println("EXPECTED (cf_cbo_debt.jld2): 162 iterations, inf-norm 9.919e-06, 47.4 min,")
println("gate passed, guards inert, CE mean 1.691")

# ============================================================
# Save (converged only)
# ============================================================
if sol_c.f_converged && gate_ok && branch_ok
    save_object(joinpath(ROOT, "output", "jld2_sw", SAVE_TAG * ".jld2"),
                (e_init, e_final, path_c, π_endog))
    println("\nwrote output/jld2_sw/", SAVE_TAG, ".jld2")
else
    println("\nNOT written.  converged = ", sol_c.f_converged,
            "   gate = ", gate_ok, "   period-0 branch = ", branch_ok)
end
