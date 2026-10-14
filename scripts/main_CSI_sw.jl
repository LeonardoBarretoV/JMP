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
# [CSI] Consumer-Sentiment "better off" shares for the sticky-wage baselines
# ------------------------------------------------------------
# Run in a FRESH Julia session.  Produces output/matlab/csi_shares_sw.mat, which
# Plot_CSI_answers_sw.m reads -- no hard-coded numbers anywhere.
#
# See src/simulation_sw.jl for the provenance of the estimator.  Short version:
# the numbers in the original CSI_answers_3.png came from Plots_histograms.m in
# the Fourth Year project, NOT from any share_better_* series in Julia, using a
# ~0.02% band that has since been widened to 1% in that script.
# ============================================================
include(joinpath(ROOT, "src", "simulation_sw.jl"))

N_SIM   = 1_000_000
T0_SIM  = 10
TT_SIM  = 40
U_HIGH  = 1.0002        # the band that reproduces CSI_answers_3 exactly
SEED    = 9823
YEARS   = ["mean 2008-2019","2020","2021","2022","2023","2024"]
DATA    = [37 28 35; 43 30 27; 42 30 28; 32 22 47; 32 24 45; 31 27 43]   # Plot_CSI_answers.m
FOURTH  = [30, 43, 44, 29, 27, 28]   # reproduced EXACTLY from sim_data_2.mat, u_high in [1e-7,5e-4]

_inp = load_object(joinpath(ROOT,"output","jld2_sw","steady_states.jld2"))
e_init, e_final = _inp[1], _inp[2]
h = e_init.h; t = e_init.t; g = e_init.g
b_target = e_init.b / y(e_init)
b_path = let
    T=700; H=100; rB=0.08; bl=Array{Float64,1}(undef,T+H)
    bl[1]=b_target; bl[2]=b_target; bl[T:end].=b_target
    for i in 3:T-1; bl[i]=bl[1]*rB^(i-1)+(1-rB^(i-1))*b_target; end
    bl
end
ts = make_time_series_sw(c, e_init, e_final, b_path, 0.287*4.0, g)
pars_h = (θ=c.θ, ν=c.ν, ies=c.ies, crra=c.crra, β=c.β, a_min=c.a_min, a_max=c.a_max,
          z_u=0.0, n_a=c.n_a, init_state=e_init.h.Pss, z_scale=z_scale)
h_path, _, _ = h_path_transition_mckay_init(ROOT, ts.Pee_path, ts.Puu_path; pars_h...)

RUNS = [("BASELINE  (epsm fed, pi = data)", "reduced_epsm.jld2"),
        ("Run 3     (pi endog, ~2%)",  "reduced_run3.jld2")]
res = Dict{String,Any}()
for (lab, file) in RUNS
    println("\n", "="^70, "\n", lab, "   <- ", file); flush(stdout)
    o = load_object(joinpath(ROOT,"output","jld2_sw",file))
    pth = o[3]
    @time sim = simulate_panel_sw(e_init, pth, h_path; N=N_SIM, T0=T0_SIM, T_trans=TT_SIM, seed=SEED)
    for (tag, uu) in (("v", sim.u_v), ("x", sim.u_x))
        s = csi_shares(uu; u_high = U_HIGH)
        res[file*"_"*tag] = s
        @printf("  u = %s :  better %s\n", tag == "v" ? "stored value function" : "flow composite x", string(s.better))
        @printf("            same   %s\n            worse  %s\n", string(s.same), string(s.worse))
        flush(stdout)
    end
    o = nothing; pth = nothing; sim = nothing; GC.gc(); GC.gc()
end

println("\n", "="^70, "\n  BETTER-OFF SHARE (%)  -- all variants vs data\n", "="^70)
@printf("  %-18s", "period"); for y in YEARS; @printf("%16s", y); end; println()
@printf("  %-18s", "DATA"); for i in 1:6; @printf("%16.1f", Float64(DATA[i,1])); end; println()
@printf("  %-18s", "Fourth Year"); for x in FOURTH; @printf("%16.1f", Float64(x)); end; println()
for (lab, file) in RUNS, tag in ("v","x")
    @printf("  %-18s", (tag=="v" ? "v: " : "x: ") * split(lab)[1])
    for x in res[file*"_"*tag].better; @printf("%16.1f", x); end; println()
end

matfolder = joinpath(ROOT, "output", "matlab")
matwrite(joinpath(matfolder, "csi_shares_sw.mat"), Dict(
    "years"          => YEARS,
    "data"           => Float64.(DATA),
    "fourth_year"    => Float64.(FOURTH),
    "epsm_fed_v"     => Float64.(res["reduced_epsm.jld2_v"].better),
    "epsm_fed_x"     => Float64.(res["reduced_epsm.jld2_x"].better),
    "run3_v"         => Float64.(res["reduced_run3.jld2_v"].better),
    "run3_x"         => Float64.(res["reduced_run3.jld2_x"].better),
    "epsm_fed_v_same"=> Float64.(res["reduced_epsm.jld2_v"].same),
    "epsm_fed_v_worse"=>Float64.(res["reduced_epsm.jld2_v"].worse),
    "run3_v_same"    => Float64.(res["reduced_run3.jld2_v"].same),
    "run3_v_worse"   => Float64.(res["reduced_run3.jld2_v"].worse),
    "u_high"         => U_HIGH, "N" => Float64(N_SIM), "T0" => Float64(T0_SIM), "seed" => Float64(SEED)))
println("\nwrote output/matlab/csi_shares_sw.mat")
println("=== CSI done ===")
