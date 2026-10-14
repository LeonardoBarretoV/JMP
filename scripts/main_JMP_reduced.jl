# ============================================================
# DEBT-FINANCED FISCAL STIMULUS, HETEROGENEITY, AND WELFARE
# Author: Leonardo Barreto 
# Created on: 10/16/2024
# Edits history:
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
include(joinpath(ROOT, "src", "union_block.jl"))    # [SW] sticky-wage track ONLY
include(joinpath(ROOT, "src", "transition_sw.jl"))  # [SW] sticky-wage track ONLY
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
    P, z_vals, z_scale = calibration_unemploy_mckay_(ROOT, c.Pee, c.Puu, 0.0)   # [Z-NORM] z_vals now normalized so Z = 1 in the initial SS (was: P, z_vals = ...)
    u = EZ(ies = c.ies, ra = c.crra, β = c.β)

    a_grid = grid(start = c.a_min, stop = c.a_max, length = c.n_a, scale = :log)
    return Household(u = u, a_grid = a_grid, v = v, P = P, z_grid = z_vals), z_scale   # [Z-NORM] also return the scale factor (was: return Household(...))
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

h, z_scale = make_household(c)   # [Z-NORM] z_scale is reused for the transition households (was: h = make_household(c))
t, δ = make_technology(c)
g = make_taxsys(c, h)
τk_ss = τ_k_ss(c, δ)
println("[Z-NORM] scale factor = $z_scale;  Z in initial SS = $(sum(h.Pss .* h.z_grid))")   # [Z-NORM] sanity check: should print 1.0

# ============================================================
# [SW6] REDUCED (pi-eliminated) DRIVER  --  Stage 4 Run 3
# ------------------------------------------------------------
# This is NOT main_JMP_sw.jl.  That script is the Stage 3 solver: 2x2 on (r,s)
# with pi FED from data.  Here pi is ENDOGENOUS, eliminated analytically by the
# backward monetary recursion, so the unknown vector is still [r ; s] (1598) and
# the monetary block holds identically at every residual evaluation.
#
# RUN IT IN A FRESH JULIA SESSION.  _assert_globals_sw enforces one solve per
# process; a second solve in the same session will fail the guard by design.
#
# Everything above this line is verbatim from main_JMP_sw.jl (packages,
# includes, Calib, make_household/technology/taxsys) so the two share one setup.
# ============================================================

# ---- knobs ----------------------------------------------------------------
KAPPA_W    = 0.05     # wage Phillips slope
PHI_PI     = 1.5      # Taylor coefficient
NKPC_SCALE = 240.0    # row scaling; measured, equalises asset/NKPC diagonals (ratio 1.01)
FEED_EPSM  = true     # true = BASELINE (epsm fed at its backed-out path; pi matches data)
                      # true  = regression (b): feed the backed-out epsm, must return pi = pi_data
USE_CACHE  = true     # load output/jld2_sw/steady_states.jld2 instead of re-solving the SS

# ---- steady states --------------------------------------------------------
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

# ---- time series (forked: retargeted period-0 debt, g passed explicitly) ---
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

# ---- warm start from the Stage 3 solution ---------------------------------
println("loading baseline_stage3_kw05 for the warm start ..."); flush(stdout)
_o = load_object(joinpath(ROOT, "output", "jld2_sw", "baseline_stage3_kw05.jld2")); pc = _o[3]
n_per = length(pc.r); n_r = _n_r_sw(n_per, 1); n_s = n_per - 2
π_bar = e_init.π_
i_ss  = (1 + e_final.r)*(1 + π_bar) - 1          # Taylor intercept, built from r_FINAL
r0 = Vector{Float64}(undef, n_r); r0[1] = pc.r_k[1] - t.δ      # r_k_hat; r_b is substituted out
for i in 2:n_per; r0[i] = pc.r[i]; end                          # no +1 shift
ñ   = e_init.ws.n_pol[1,:]
emp = [iz for iz in eachindex(h.z_grid) if h.z_grid[iz] != 0]
s0  = [sum(pc[j].n_pol[1,iz]/ñ[iz] for iz in emp)/length(emp) for j in 1:n_per][2:(n_per-1)]

# epsm backed out of the Stage 3 path (only used when FEED_EPSM)
i_conv = [ j < n_per ? (1+pc.r[j+1])*(1+pc.π_[j+1]) - 1 : i_ss for j in 1:n_per ]
epsm   = [ i_conv[j] - i_ss - PHI_PI*(pc.π_[j] - π_bar) for j in 1:n_per ]
epsm_use = FEED_EPSM ? epsm : zeros(n_per)

# ---- solve ----------------------------------------------------------------
path_c, sol_c, π_endog, info_c = solve_transition_sw_c(e_init, e_final;
    τ_k_path=ts.τ_k_path, T_path=ts.T_path, Tu_path=ts.Tu_path, G_path=ts.G_path,
    qb_path=ts.qb_path, qbstar_path=ts.bstar_path_r, h_path=h_path, π_path=ts.π_path,
    d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
    method=1, kappa_w=KAPPA_W, epsm_path=epsm_use, phi_pi=PHI_PI, i_taylor_ss=i_ss,
    nkpc_scale=NKPC_SCALE, s_mode=:nkpc, r_path=r0, s_path=s0)

# ---- report ---------------------------------------------------------------
println("\n   j  quarter     pi_data(ann%)  pi_model(ann%)     diff(ann pp)")
ql = ["2020Q2","2020Q3","2020Q4","2021Q1","2021Q2","2021Q3","2021Q4","2022Q1","2022Q2","2022Q3","2022Q4","2023Q1"]
for j in 1:12
    @printf("%4d  %s   %+10.4f     %+10.4f     %+12.3e\n", j, ql[j],
            100*((1+ts.π_path[j])^4-1), 100*((1+π_endog[j])^4-1), 400*(π_endog[j]-ts.π_path[j]))
end
for j in (16, 20, 40, 100, 800)
    @printf("%4d           %+10.4f     %+10.4f     %+12.3e\n", j,
            100*((1+ts.π_path[j])^4-1), 100*((1+π_endog[j])^4-1), 400*(π_endog[j]-ts.π_path[j]))
end
@printf("\npeak deviation from target: model %+.4f pp   data %+.4f pp\n",
        400*maximum(abs.(π_endog[2:60] .- π_bar)), 400*maximum(abs.(ts.π_path[2:60] .- π_bar)))

# ---- save (converged only) ------------------------------------------------
if sol_c.f_converged
    tag = FEED_EPSM ? "reduced_epsm" : "reduced_run3"
    save_object(joinpath(ROOT, "output", "jld2_sw", tag * ".jld2"), (e_init, e_final, path_c, π_endog))
    println("\nwrote output/jld2_sw/", tag, ".jld2")
else
    println("\nNOT converged -- nothing written")
end
