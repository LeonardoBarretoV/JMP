# ============================================================
# DEBT-FINANCED FISCAL STIMULUS, HETEROGENEITY, AND WELFARE
# THE KAPPA_W ROBUSTNESS RUN  --  kappa_w = 0.30, epsm = 0
# Author: Leonardo Barreto
# ============================================================
# Reproduces output/jld2_sw/reduced_run3_kw30.jld2.
#
# WHY THIS FILE EXISTS
#   reduced_run3_kw30.jld2 had no producer on disk.  main_JMP_reduced.jl only ever
#   writes the tags "reduced_epsm" or "reduced_run3" (see its save block), so that
#   file came from a hand-edited copy which was not kept.  It is the file behind the
#   claim that NO wage-stickiness value generates the 2021-22 surge, so it should
#   not be orphaned.
#
#   This is main_JMP_reduced.jl with three knobs changed and its own output tag.
#   Confirmed against the saved file: it carries the BASELINE final steady state
#   (r_f = 1.0778% ann, i_taylor_ss = 3.0994% ann), so no steady state is re-solved
#   and steady_states.jld2 is used as-is.
#
# RUN IT IN A FRESH JULIA SESSION.  _assert_globals_sw enforces one solve per
# process; a second solve in the same session fails the guard by design.
#
#     julia --project -t auto scripts/main_reduced_kw30_sw.jl
#
# THE RESULT THIS PRODUCES
#   A 6x steeper wage curve moves the explained share of the surge only
#   0.23% -> 0.78% (peak deviation ABOVE target +0.016pp -> +0.055pp against the
#   data's +7.06pp).  The chain kappa_w -> s -> r -> pi is throttled at the r link:
#   hours respond (max|ds| = 1.53%) but max|dr| = 0.062 ann pp, because r is pinned
#   by asset-market clearing and then phi = 1.5 pins pi.  As kappa_w -> infinity the
#   NKPC forces the wage wedge to zero, i.e. s -> s_flex, the flexible-wage limit
#   (= the original paper), so the family is BOUNDED and no kappa_w generates the
#   surge.
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
include(joinpath(ROOT, "src", "counterfactuals_sw.jl"))   # [SW11] gate / CE / guard reporting
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
# ------------------------------------------------------------
# The three that differ from main_JMP_reduced.jl are marked <<< CHANGED.
# ============================================================
KAPPA_W    = 0.30     # <<< CHANGED (0.05).  6x steeper wage Phillips curve.
PHI_PI     = 1.5      # Taylor coefficient
NKPC_SCALE = 41.4     # <<< CHANGED (240.0).  nkpc_scale MUST be remeasured when kappa_w
                      # changes: it equalises the asset and NKPC Jacobian diagonals and
                      # scales ~linearly in kappa_w (0.43 at kappa_w=0.05 -> 2.54 at 0.30,
                      # a factor 5.90 against the 6x change), so 240 -> 240/5.90 = 41.4.
                      # NOTE this knob only ROW-SCALES the NKPC block; it moves the solver
                      # path but not the solution set, so it is the one knob the residual
                      # check below cannot pin down (see verify_cf_inputs_sw.jl kw30).
FEED_EPSM  = false    # <<< CHANGED (true).  false => epsm = 0, pi is the MODEL's own.
USE_CACHE  = true     # load output/jld2_sw/steady_states.jld2 instead of re-solving
SAVE_TAG   = "reduced_run3_kw30"

# ============================================================
# Steady states  (BASELINE -- kappa_w does not touch either steady state)
# ============================================================
ss_file = joinpath(ROOT, "output", "jld2_sw", "steady_states.jld2")
if USE_CACHE && isfile(ss_file)
    println("loading cached steady states ..."); flush(stdout)
    _inp = load_object(ss_file); e_init, e_final = _inp[1], _inp[2]
else
    @time e_init = stationary_equilibrium_init(h, t, g, τk_ss, c.Λ_t, c.Λ_tu, c.Λ_g, c.Λ_b, c.Λ_bstar;
                                               r_range = (-0.007*0, 0.01), root_max = 1, verbose = true)
    ī_ = (1 + e_init.r)*(1 + c.π_bar) - 1
    @time e_final = stationary_equilibrium_final(h, t, g, e_init.τ_k, c.Λ_t, c.Λ_tu, c.Λ_g,
                                                 c.Λ_b + 0.15*4.0, 0.287*4.0, ī_;
                                                 r_range = (0.0, 0.005), root_max = 1.0, verbose = true)
    save_object(ss_file, (e_init, e_final))
end
h = e_init.h; t = e_init.t; g = e_init.g
println("Debt-to-GDP = ", e_init.b / y(e_init))
@printf("final SS: r_f = %.4f%% ann   (EXPECTED 1.0778%%, the BASELINE final SS)\n",
        100*((1+e_final.r)^4-1))

# ============================================================
# Time series -- identical to the baseline; only kappa_w changes
# ============================================================
b_target = e_init.b / y(e_init)
b_path = let
    T=700; H=100; ρB=0.08; bl = Array{Float64,1}(undef, T+H)
    bl[1]=b_target; bl[2]=b_target; bl[T:end] .= b_target
    for i in 3:T-1; bl[i] = bl[1]*ρB^(i-1) + (1-ρB^(i-1))*b_target; end
    bl
end
ts = make_time_series_sw(c, e_init, e_final, b_path, 0.287*4.0, g)
pars_h = (θ=c.θ, ν=c.ν, ies=c.ies, crra=c.crra, β=c.β, a_min=c.a_min, a_max=c.a_max,
          z_u=0.0, n_a=c.n_a, init_state=e_init.h.Pss, z_scale=z_scale)
h_path, _, _ = h_path_transition_mckay_init(ROOT, ts.Pee_path, ts.Puu_path; pars_h...)

# ============================================================
# Warm start from the Stage 3 solution, and the monetary residual
# ============================================================
println("loading baseline_stage3_kw05 for the warm start ..."); flush(stdout)
_o = load_object(joinpath(ROOT, "output", "jld2_sw", "baseline_stage3_kw05.jld2")); pc = _o[3]
n_per = length(pc.r)
π_bar = e_init.π_
i_ss  = (1 + e_final.r)*(1 + π_bar) - 1          # Taylor intercept, built from r_FINAL
r0, s0 = warm_start_sw(pc, e_init, t; method = 1)

# epsm backed out of the Stage 3 path (only used when FEED_EPSM)
epsm     = epsm_baseline_sw(copy(pc.r), copy(pc.π_), i_ss, π_bar, PHI_PI)
epsm_use = FEED_EPSM ? epsm : zeros(n_per)
_o = nothing; pc = nothing; GC.gc()

@printf("i_taylor_ss = %.10f  (%.4f%% ann)   EXPECTED 3.0994%% ann\n",
        i_ss, 100*((1+i_ss)^4-1))
@printf("max|epsm| in use = %.10f   (FEED_EPSM = %s)\n", maximum(abs, epsm_use), FEED_EPSM)

# ============================================================
# Solve
# ============================================================
path_c, sol_c, π_endog, info_c = solve_transition_sw_c(e_init, e_final;
    τ_k_path=ts.τ_k_path, T_path=ts.T_path, Tu_path=ts.Tu_path, G_path=ts.G_path,
    qb_path=ts.qb_path, qbstar_path=ts.bstar_path_r, h_path=h_path, π_path=ts.π_path,
    d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
    method=1, kappa_w=KAPPA_W, epsm_path=epsm_use, phi_pi=PHI_PI, i_taylor_ss=i_ss,
    nkpc_scale=NKPC_SCALE, s_mode=:nkpc, r_path=r0, s_path=s0)

# ============================================================
# Report
# ============================================================
print_pi_table_sw(ts.π_path, π_endog, π_bar; label = "kappa_w=0.30")

# The headline statistic for this run is the peak deviation ABOVE target, because the
# claim is about explaining the SURGE.  The largest ABSOLUTE deviation is on the
# downside and is a different number -- keep them apart.
up_model = maximum(π_endog[2:60] .- π_bar)
up_data  = maximum(ts.π_path[2:60] .- π_bar)
println("\n---- the kappa_w bound ----")
@printf("peak deviation ABOVE target: model %+.4f pp   data %+.4f pp   -> share %.2f%%\n",
        400*up_model, 400*up_data, 100*up_model/up_data)
@printf("peak ABSOLUTE deviation:     model %+.4f pp  (downside; do NOT quote as the share)\n",
        400*maximum(abs.(π_endog[2:60] .- π_bar)))
@printf("pi range over j = 2..40:     %.2f%% to %.2f%% ann\n",
        100*(minimum(1 .+ π_endog[2:40])^4-1), 100*(maximum(1 .+ π_endog[2:40])^4-1))
println("EXPECTED (reduced_run3_kw30.jld2): +0.055pp above target, 0.78% of the surge,")
println("pi range 1.39-2.05%, CE mean 0.645")

ce = ce_summary_sw(e_init, path_c)
print_ce_sw(ce; label = "(kappa_w = 0.30, epsm = 0)")

gate = gate_lab_tax_sw(e_init, path_c, g)
gate_ok = print_gate_sw(gate)
branch_ok = period0_branch_sw(e_init, path_c, g, c.ν, ts.Γ_n)
guard_report_sw()

@printf("\n---- run record ----\niterations = %d   converged = %s   wall = %.1f min\n",
        sol_c.iterations, sol_c.f_converged, info_c.wall/60)

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
