# ============================================================
# [SW19] GATES 1 AND 2 -- regression against rule2_epsm, then block diagonals
# ============================================================
#     julia --project -t auto scripts/check_sw19_sw.jl
#
# GATE 1  Re-solve the rule-2 baseline with [SW19] ON and confirm it reproduces
#         rule2_epsm.jld2:
#           - r_k must come out at the value the old solver FOUND
#           - residuals[1] must be GONE, not merely small: the system is 1597 long
#             and the period-0 capital condition holds identically
#           - paths must match to solver tolerance for j >= 2, with the period-0
#             difference equal to however far the old solve missed k[1] = e_init.k
# GATE 2  Remeasure the asset and NKPC Jacobian block diagonals after the
#         elimination.  nkpc_scale = 240 was measured on a system that still had
#         residuals[1] in it; report whether the ratio moved BEFORE changing it.
#
# Nothing is written.
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
for f in ("types","preferences","technology_and_fiscal","calibration","auxiliary_functions",
          "household_problem_and_aggregation","stationary_equilibria","transition",
          "union_block","transition_sw","counterfactuals_sw","statistics","time_series")
    include(joinpath(ROOT, "src", f * ".jl"))
end
const _TOL = 1e-12; const _ZERO_FTOL = 1e-5
const _MAX_ITERS = 20_000; const _ZERO_MAX_ITERS = 450

KAPPA_W = 0.05; PHI_PI = 1.5; NKPC_SCALE = 240.0
TK_RULE2 = -0.47520201
T1_RULE2 = 0.32446192
TU1_RULE2 = 0.31282997

Base.@kwdef struct Calib
    θ::Float64 = 1.0; ν::Float64 = 1.0; ies::Float64 = 1.0; crra::Float64 = 4.0
    β::Float64 = 0.9880431080390226
    a_min::Float64 = -0.02; a_max::Float64 = 120.0; n_a::Int = 1000; z_bar::Int = 12
    Pee::Float64 = 0.975; Puu::Float64 = 0.4
    Λ_t::Float64 = 0.0; Λ_tu::Float64 = 0.005; Λ_a::Float64 = 5.06*4.0
    Λ_g::Float64 = 0.05; Λ_f::Float64 = 0.04; Λ_b::Float64 = 1.05*4.0
    Λ_bstar::Float64 = 0.31*4.0; Λ_rk::Float64 = 0.33
    r_target::Float64 = (1.0+0.01)^(1/4)-1.0; π_bar::Float64 = (1.0+0.02)^(1/4)-1.0
    mc::Float64 = 5.0*4.0; τ::Float64 = 0.181; nu_target::Float64 = 0.7
    A::Float64 = 0.035; μ::Float64 = 1.0
    d_::Float64 = 1.0; Γ_n::Float64 = 1.0; Γ_c::Float64 = 1.0; covid::Float64 = 0.0
    χ::Float64 = 0.44
end
c = Calib()
_, _, z_scale = calibration_unemploy_mckay_(ROOT, c.Pee, c.Puu, 0.0)

println("loading steady states and the rule-2 baseline ..."); flush(stdout)
_inp = load_object(joinpath(ROOT,"output","jld2_sw","steady_states.jld2"))
e_init, e_final = _inp[1], _inp[2]
h = e_init.h; t = e_init.t; g = e_init.g
δ = t.δ; π_bar = e_init.π_
i_ss = (1 + e_final.r)*(1 + π_bar) - 1

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

T_r2  = copy(ts.T_path);   T_r2[1]  = T1_RULE2
Tu_r2 = copy(ts.Tu_path);  Tu_r2[1] = TU1_RULE2
τk_r2 = copy(ts.τ_k_path); τk_r2[1] = TK_RULE2

_o = load_object(joinpath(ROOT,"output","jld2_sw","rule2_epsm.jld2"))
pOLD = _o[3]; πOLD = copy(_o[4])
rk_old   = pOLD.r_k[1]
k1_old   = pOLD.k[1]
r_old    = copy(pOLD.r); τn_old = copy(pOLD.τ_n); w_old = copy(pOLD.w); y_old = copy(pOLD.y)
_o = nothing; GC.gc()
_o = load_object(joinpath(ROOT,"output","jld2_sw","rule2_stage1.jld2"))
epsm = epsm_baseline_sw(copy(_o[3].r), copy(_o[3].π_), i_ss, π_bar, PHI_PI)
_o = nothing; GC.gc()

n_per = length(τk_r2); n_s = n_per - 2
n_r_old = _n_r_sw(n_per, 1, false)
n_r_19  = _n_r_sw(n_per, 1, true)
@printf("\nunknowns: without [SW19] %d + %d = %d    with [SW19] %d + %d = %d\n",
        n_r_old, n_s, n_r_old+n_s, n_r_19, n_s, n_r_19+n_s)

# what the OLD solve actually achieved on the period-0 capital row
@printf("\nold solve, period-0 capital row: e_init.k = %.12f, path.k[1] = %.12f\n",
        e_init.k, k1_old)
@printf("  residuals[1] at the old solution = %.6e   <- the miss [SW19] removes\n",
        e_init.k - k1_old)
@printf("old solve r_k[1] = %.12f  (r_k_hat = %.12f)\n", rk_old, rk_old - δ)

# closed-form r_k* at the rule-2 tau_k
n0 = labor_supply_union_(h_path[2]; s = ts.Γ_n, n_ss = e_init.ws.n_pol[1,:],
                         z_dist = (let zs = zeros(length(h_path), length(e_init.h.z_grid))
                                       zs[1,:] = transpose(e_init.h.Pss)*h_path[1].P
                                       zs[1,:] end))
r_min = -δ + eps(); r_max = 1/h_path[1].u.β
rk_star = _rk_star_sw19(t, e_init; τ_k = TK_RULE2, n0 = n0, d = ts.d_path[1],
                        covid = ts.covid_path[1], r_min = r_min, r_max = r_max)
@printf("closed-form r_k*  = %.12f  (r_k_hat* = %.12f)\n", rk_star, rk_star - δ)
@printf("  difference from what the old solver found: %.6e  (%.4f%%)\n",
        rk_star - rk_old, 100*(rk_star/rk_old - 1))

# ============================================================
# GATE 1 -- re-solve with [SW19] on
# ============================================================
println("\n", "="^100)
println("GATE 1 -- re-solve the rule-2 baseline with [SW19] ON")
println("="^100)
r0_19, s0_19 = warm_start_sw(pOLD, e_init, t; method = 1, sw19 = true)
@printf("warm start length %d (expected %d)\n", length(r0_19), n_r_19)
path19, sol19, π19, info19 = solve_transition_sw_c(e_init, e_final;
    τ_k_path=τk_r2, T_path=T_r2, Tu_path=Tu_r2, G_path=ts.G_path,
    qb_path=ts.qb_path, qbstar_path=ts.bstar_path_r, h_path=h_path, π_path=ts.π_path,
    d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
    method=1, kappa_w=KAPPA_W, epsm_path=epsm, phi_pi=PHI_PI, i_taylor_ss=i_ss,
    nkpc_scale=NKPC_SCALE, s_mode=:nkpc, r_path=r0_19, s_path=s0_19, sw19=true)
@printf("\n[SW19] iterations %d   converged %s   wall %.2f min\n",
        sol19.iterations, sol19.f_converged, info19.wall/60)
guard_report_sw()

println("\n---- period-0 capital: is residuals[1] GONE, not merely small? ----")
@printf("e_init.k       = %.14f\n", e_init.k)
@printf("path.k[1] [SW19] = %.14f   difference = %.6e\n", path19.k[1], e_init.k - path19.k[1])
@printf("path.k[1] old    = %.14f   difference = %.6e\n", k1_old, e_init.k - k1_old)
@printf("r_k[1] [SW19] = %.12f   old = %.12f   diff = %.3e\n",
        path19.r_k[1], rk_old, path19.r_k[1] - rk_old)

println("\n---- path agreement, [SW19] vs rule2_epsm ----")
@printf("%-16s %14s %14s %14s\n", "series", "max|diff| j>=2", "at j=1", "max rel j>=2")
for (nm, a, b) in (("r",     path19.r,   r_old),
                   ("tau_n", path19.τ_n, τn_old),
                   ("w",     path19.w,   w_old),
                   ("y",     path19.y,   y_old),
                   ("pi",    π19,        πOLD))
    d2 = maximum(abs.(a[2:end] .- b[2:end]))
    rel = maximum(abs.((a[2:end] .- b[2:end]) ./ max.(abs.(b[2:end]), eps())))
    @printf("%-16s %14.6e %14.6e %14.6e\n", nm, d2, abs(a[1]-b[1]), rel)
end

# ============================================================
# GATE 2 -- block diagonals after the elimination
# ============================================================
println("\n", "="^100)
println("GATE 2 -- Jacobian block diagonals, with and without [SW19]")
println("="^100)
DY_ss = union_aggregates_(e_init.h, g, e_init.ws; s=1.0, w=e_init.w, τ_n=e_init.τ_n,
                          n_ss=e_init.ws.n_pol[1,:])
wY_ss = e_init.w*DY_ss.Y

function diag_probe(x, sw19flag, nr; h = 1e-6, idx_asset, idx_nkpc)
    path_t = _initialize_path(e_init, e_final; qb_path = ts.qb_path)
    πv = Vector{Float64}(undef, n_per)
    kw = (; π_vec=πv, epsm_path=epsm, i_taylor_ss=i_ss, phi_pi=PHI_PI, path=path_t,
          e_init=e_init, e_final=e_final, τ_k_path=τk_r2, T_path=T_r2, Tu_path=Tu_r2,
          G_path=ts.G_path, qb_path=ts.qb_path, qbstar_path=ts.bstar_path_r,
          h_path=h_path, π_path=ts.π_path, d_path=ts.d_path, covid_path=ts.covid_path,
          Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g, method=1, kappa_w=KAPPA_W,
          wY_ss=wY_ss, nkpc_scale=NKPC_SCALE, s_mode=:nkpc, sw19=sw19flag)
    F0 = similar(x); _residuals_sw_c!(F0, x; kw...)
    ds = Float64[]
    for k in idx_asset
        xp = copy(x); xp[k] += h
        F1 = similar(x); _residuals_sw_c!(F1, xp; kw...)
        push!(ds, (F1[k] - F0[k])/h)
    end
    dn = Float64[]
    for k in idx_nkpc
        xp = copy(x); xp[k] += h
        F1 = similar(x); _residuals_sw_c!(F1, xp; kw...)
        push!(dn, (F1[k] - F0[k])/h)
    end
    return (asset = ds, nkpc = dn)
end

# probe a handful of interior rows in each block, same PERIODS in both systems
per_asset = [50, 200, 400]      # asset row for period j sits at index j-1 (SW19) or j (old)
per_nkpc  = [50, 200, 400]      # nkpc row for period j sits at n_r + (j-1)
x19 = vcat(collect(path19.r[2:n_per]), collect(sol19.zero[(n_r_19+1):end]))
xold_r = Vector{Float64}(undef, n_r_old); xold_r[1] = rk_old - δ
for i in 2:n_per; xold_r[i] = r_old[i]; end
xold = vcat(xold_r, collect(sol19.zero[(n_r_19+1):end]))

d19 = diag_probe(x19, true,  n_r_19,
                 idx_asset = [j-1 for j in per_asset],
                 idx_nkpc  = [n_r_19 + (j-1) for j in per_nkpc])
dol = diag_probe(xold, false, n_r_old,
                 idx_asset = [j   for j in per_asset],
                 idx_nkpc  = [n_r_old + (j-1) for j in per_nkpc])
@printf("%-8s %16s %16s %16s %16s\n", "period", "asset [SW19]", "asset old", "nkpc [SW19]", "nkpc old")
for (i, j) in enumerate(per_asset)
    @printf("%-8d %16.6e %16.6e %16.6e %16.6e\n", j, d19.asset[i], dol.asset[i],
            d19.nkpc[i], dol.nkpc[i])
end
ma19 = sum(abs, d19.asset)/length(d19.asset); mn19 = sum(abs, d19.nkpc)/length(d19.nkpc)
mao  = sum(abs, dol.asset)/length(dol.asset); mno  = sum(abs, dol.nkpc)/length(dol.nkpc)
@printf("\nmean |asset diagonal|: [SW19] %.6e   old %.6e   ratio %.4f\n", ma19, mao, ma19/mao)
@printf("mean |nkpc  diagonal|: [SW19] %.6e   old %.6e   ratio %.4f\n", mn19, mno, mn19/mno)
@printf("\nasset/nkpc diagonal ratio: [SW19] %.4f   old %.4f\n", ma19/mn19, mao/mno)
println("nkpc_scale = 240 was measured to bring this ratio to ~1. If the [SW19] ratio is")
println("still ~1, 240 stands; if it moved materially, nkpc_scale needs remeasuring.")
println("\nDONE (nothing written)")
