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
    P, z_vals = calibration_unemploy_mckay_(ROOT, c.Pee, c.Puu, 0.0)
    u = EZ(ies = c.ies, ra = c.crra, β = c.β)

    a_grid = grid(start = c.a_min, stop = c.a_max, length = c.n_a, scale = :log)
    return Household(u = u, a_grid = a_grid, v = v, P = P, z_grid = z_vals)
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

h = make_household(c)
t, δ = make_technology(c)
g = make_taxsys(c, h)
τk_ss = τ_k_ss(c, δ)

@time e_init = stationary_equilibrium_init(h, t, g, τk_ss, c.Λ_t, c.Λ_tu, c.Λ_g, c.Λ_b, c.Λ_bstar; r_range = (-0.007 * 0, 0.01), root_max = 1, verbose = true)
b_y = e_init.b / y(e_init);
println("Debt-to-GDP = $b_y;")

# ============================================================
# Debt path in the initial steady state
# ============================================================
b_target = e_init.b / y(e_init)
b_path = let   # Smooth path of constant debt-to-GDP 
    T = 700  # period of adjustment of debt
    H = 100   # debt level no longer moving
    ρB = 0.08
    b_list = Array{Float64,1}(undef, T + H)
    b_list[1] = e_init.b / y(e_init) #2019Q4
    b_list[2] = e_init.b / y(e_init) #2020Q1
    b_list[T:end] .= b_target
    for i in 3:T-1
        b_list[i] = b_list[1] * ρB^(i-1) + (1 - ρB^(i-1)) * b_target
    end
    b_list
end;

# ============================================================
# Final steady state
# ============================================================
ī = (1 + e_init.r)*(1+c.π_bar) - 1
Λ_b_f = c.Λ_b + 0.15*4.0 # debt-to-GDP in the final steady state
Λ_bstar_f = 0.287 * 4.0 # external debt-to-GDP in the final steady state
@time e_final = stationary_equilibrium_final(h, t, g, e_init.τ_k, c.Λ_t, c.Λ_tu, c.Λ_g, Λ_b_f, Λ_bstar_f, ī; r_range= (0.0, 0.005), root_max = 1.0, verbose = true)

# ============================================================
# Transition path
# ============================================================
ts = make_time_series(c, e_init, e_final, b_path, Λ_bstar_f)

τ_k_path    = ts.τ_k_path
T_path      = ts.T_path
Tu_path     = ts.Tu_path
G_path      = ts.G_path
b_path      = ts.qb_path
bstar_path_r = ts.bstar_path_r
π_path      = ts.π_path
d_path      = ts.d_path
covid_path  = ts.covid_path
Γ_n         = ts.Γ_n
Γ_c         = ts.Γ_c

Pee_path = ts.Pee_path
Puu_path = ts.Puu_path
pars_h_path = (θ = c.θ, ν = c.ν, ies = c.ies, crra = c.crra, β = c.β, a_min = c.a_min, a_max = c.a_max, z_u = 0.0, n_a = c.n_a, init_state = e_init.h.Pss)
h_path, chain, chain_long = h_path_transition_mckay_init(ROOT, Pee_path, Puu_path; pars_h_path...)

@time path_r = solve_transition_h_r(e_init, e_final;
    τ_k_path = τ_k_path, T_path = T_path, Tu_path = Tu_path, G_path = G_path,
    qb_path = b_path, qbstar_path = bstar_path_r, h_path = h_path,
    π_path = π_path, d_path = d_path, covid_path = covid_path,
    Γ_n = Γ_n, Γ_c = Γ_c, χ = c.χ, t, g,
    method = 1, r_path = [e_init.r for _ in 1:1:length(b_path)+1]
)

# ============================================================
# Generate tables and save 
# ============================================================
# Print tables
tab_detailed = generate_tables(e_init, e_final, path_r, io=String);
println(tab_detailed...)
tab_detailed = generate_tables_detailed(e_init, e_final, path_r, io=String);
println(tab_detailed...)

# Save results
jldfolder = joinpath(ROOT, "output", "jld2")
mkpath(jldfolder)
save_object(joinpath(jldfolder, "baseline.jld2"), (e_init, e_final, path_r))
(e_init, e_final, path_r) = load_object(joinpath(jldfolder, "baseline.jld2"))

# ============================================================
# Welfare gains
# ============================================================
mat_all = w_gains_path(e_init, path_r)
data_to_save = Dict(
    "welfare_mat" => mat_all
    )
outfolder = joinpath(ROOT, "output", "matlab")
mkpath(outfolder)
file_path = joinpath(outfolder, "welfare_time_paths.mat")
matopen(file_path, "w") do f
    for (name, value) in data_to_save
        write(f, name, value)
    end
end

# ============================================================
# Export time paths
# ============================================================
data_to_save = Dict(
    "time_series" => [path_r.y path_r.π_ path_r.un path_r.q path_r.b path_r.bd path_r.bstar path_r.l path_r.lump_sum_T path_r.lump_sum_Tu path_r.r path_r.r_k path_r.k path_r.n path_r.w path_r.C path_r.lab_tax path_r.τ_k path_r.τ_n path_r.X path_r.G],
    
    "ss_init" => [y(e_init) e_init.π_ e_init.h.Pss[1] e_init.q e_init.b e_init.bd e_init.bstar (e_init.b - ((1-g.ϑ) / (1+e_init.π_))*e_init.b) e_init.lump_sum_T e_init.lump_sum_Tu e_init.r e_init.r e_init.k e_init.n e_init.w (y(e_init) - t.δ * e_init.k - e_init.G) e_init.lab_tax e_init.τ_k e_init.τ_n (t.δ * e_init.k) e_init.G],

    "ss_final" => [y(e_final) e_final.π_ e_final.h.Pss[1] e_final.q e_final.b e_final.bd e_final.bstar (e_final.b - ((1-g.ϑ) / (1+e_final.π_))*e_final.b) e_final.lump_sum_T e_final.lump_sum_Tu e_final.r e_final.r e_final.k e_final.n e_final.w (y(e_final) - t.δ * e_final.k - e_final.G) e_final.lab_tax e_final.τ_k e_final.τ_n (t.δ * e_final.k) e_final.G]
)
outfolder = joinpath(ROOT, "output", "matlab")
mkpath(outfolder)
file_path = joinpath(outfolder, "model_time_series_ratio.mat")
matopen(file_path, "w") do f
    for (name, value) in data_to_save
        write(f, name, value)
    end
end