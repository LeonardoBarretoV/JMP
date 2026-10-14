# ============================================================
# [SW15] ALLOCATION RULES FOR THE MISSING "OTHER"  --  nothing re-solved
# ============================================================
#     julia --project -t auto scripts/fiscal_allocation_sw.jl
#
# Option (c): distribute Table E2's "Other" across the three modeled instruments.
# Each rule is sized to the ACTUAL hole in the period-0 budget, not to Other's full
# $1.3tn, because the full size overshoots by ~3.16x (see fiscal_other_sw.jl).
#
# PART A  is the hole a NUMBER or a PATH?  per-period labour-tax residual, j = 1..20
# PART B  the four allocation rules: implied tau_k, tau_n, and the tau_k removal shock
# PART C  per-decile incidence of each channel, from the initial-SS distribution
#
# All of PART B is the PERIOD-0 GBC with r_k, the debt path, pi, hours and the
# household side held at baseline.  PART C is pure incidence arithmetic off e_init:
# it says WHO receives each channel, which is the sign and relative size of the
# welfare effect by decile, but not the CE numbers -- those need a re-solve.
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

GDP_2019 = 21380.0
E2 = (checks = 21.4, ui = 19.4, firms = 30.8, other = 22.9)
TAX_CREDITS_PCT = 3.1        # Table E2 tax credits, % of the FULL package

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

println("loading reduced_epsm.jld2 ..."); flush(stdout)
o = load_object(joinpath(ROOT, "output", "jld2_sw", "reduced_epsm.jld2"))
e_init, e_final, p = o[1], o[2], o[3]
h = e_init.h; t = e_init.t; g = e_init.g
δ = t.δ; π_bar = e_init.π_
ī = (1 + e_init.r)*(1 + π_bar) - 1
y_ss = y(e_init); Y_ANN = 4*y_ss
pct(x) = 100*x/Y_ANN
dol(x) = GDP_2019*x/Y_ANN
labtax_ss = e_init.lab_tax
caprev_ss = e_init.τ_k*(e_init.r + δ)*e_init.k
T_ss = e_init.lump_sum_T; Tu_ss = e_init.lump_sum_Tu; G_ss = e_init.G

# ============================================================
println("\n", "="^100)
println("PART A -- is the hole a NUMBER (j=1 only) or a PATH?")
println("="^100)
println("The hole in period j is the LABOUR-TAX RESIDUAL: how much the fed borrowing")
println("over- or under-covers that period's outlays. Sign convention here:")
println("  hole > 0  => borrowing EXCEEDS outlays, so the labour tax is CUT (tau_n up)")
println("  hole < 0  => borrowing FALLS SHORT, so the labour tax is RAISED (tau_n down)")
println()
@printf("%5s %12s %14s %16s %12s\n", "j", "tau_n", "vs SS (%)", "hole (%ann GDP)", "hole (\$bn)")
for j in 1:20
    hole = -(p.lab_tax[j] - labtax_ss)
    @printf("%5d %12.6f %+14.2f %16.4f %12.0f\n", j, p.τ_n[j],
            100*(p.τ_n[j]/e_init.τ_n - 1), pct(hole), dol(hole))
end
cum20 = sum(-(p.lab_tax[j] - labtax_ss) for j in 1:20)
@printf("\ncumulative j=1..20 : %.4f%% ann GDP = \$%.0f bn\n", pct(cum20), dol(cum20))
cum19 = sum(-(p.lab_tax[j] - labtax_ss) for j in 1:19)
@printf("cumulative j=1..19 : %.4f%% ann GDP = \$%.0f bn\n", pct(cum19), dol(cum19))
holes = [-(p.lab_tax[j] - labtax_ss) for j in 1:19]
@printf("\nsign pattern j=1..19: %s\n", join([h > 0 ? "+" : "-" for h in holes], ""))
@printf("largest surplus  = %+.4f%% at j=%d\n", pct(maximum(holes)), argmax(holes))
@printf("largest shortfall= %+.4f%% at j=%d\n", pct(minimum(holes)), argmin(holes))
@printf("mean |hole| over j=1..19 = %.4f%% ann GDP = \$%.0f bn\n",
        pct(sum(abs, holes)/19), dol(sum(abs, holes)/19))

# ============================================================
# period-0 GBC with channel top-ups
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
z_state = zeros(length(h_path), length(e_init.h.z_grid))
z_state[1,:] = transpose(e_init.h.Pss)*h_path[1].P
for i = 2:lastindex(h_path); z_state[i,:] = transpose(z_state[i-1,:])*h_path[i].P; end
h0 = h_path[2]; z_dist0 = z_state[1,:]
n_ss = e_init.ws.n_pol[1,:]; s_0 = ts.Γ_n
r_k0 = p.r_k[1]; r_b0 = p.r[1]
qb0 = ts.qb_path[1]; qbp0 = ts.qb_path[2]; yp0 = p.y[2]
T0 = ts.T_path[1]; Tu0 = ts.Tu_path[1]; G0 = ts.G_path[1]
d0 = ts.d_path[1]; cov0 = ts.covid_path[1]

function period0(; τ_k = -0.47, dTlev = 0.0, dTulev = 0.0)
    # dTlev / dTulev are LEVEL top-ups; T enters the GBC as T*y so convert by /y.
    pr = ProgressUnknown()
    f = (τ_n) -> _τ_n_excess_t0_port_sw_!(h0, t, g; τ_n, τ_k,
        T = T0 + dTlev/p.y[1], Tu = Tu0 + dTulev/p.y[1], G = G0,
        qbprime = qbp0, qb = qb0, yprime = yp0, d = d0, s = s_0, n_ss, covid = cov0,
        z_dist = z_dist0, r_b = r_b0, r_k = r_k0, verbose = false, p = pr).excess_g
    τ_n = _solve_tau_n_sw(f, 0.0, 8.0)
    out = _τ_n_excess_t0_port_sw_!(h0, t, g; τ_n, τ_k,
        T = T0 + dTlev/p.y[1], Tu = Tu0 + dTulev/p.y[1], G = G0,
        qbprime = qbp0, qb = qb0, yprime = yp0, d = d0, s = s_0, n_ss, covid = cov0,
        z_dist = z_dist0, r_b = r_b0, r_k = r_k0, verbose = false, p = pr)
    fsub = -(τ_k*r_k0*out.k - caprev_ss)
    dT1 = out.lump_sum_T - T_ss; dTu1 = out.lump_sum_Tu - Tu_ss
    mod3 = dT1 + dTu1 + fsub
    hole = -(out.lab_tax - labtax_ss)
    return (; τ_n = out.τ_n, y = out.y, k = out.k, fsub, dT1, dTu1, mod3, hole,
            share_firms = 100*fsub/mod3, share_checks = 100*dT1/mod3, share_ui = 100*dTu1/mod3)
end

base = period0()
@assert abs(base.τ_n - p.τ_n[1]) < 1e-8 "period-0 GBC not reproduced"
@printf("\nperiod-0 GBC reproduced (tau_n = %.8f). base subsidy = %.4f%% ann GDP.\n",
        base.τ_n, pct(base.fsub))

# tau_k that delivers a target subsidy LEVEL
function τ_k_for_subsidy(target; dTlev = 0.0, dTulev = 0.0)
    f = (τ_k) -> period0(; τ_k, dTlev, dTulev).fsub - target
    return find_zero(f, (-0.95, 0.10), Bisection(); xatol = 1e-10)
end

# Solve, for allocation shares (pc, pu, pf), the total X that closes the period-0
# budget with tau_n at a TARGET value, delivering pf*X through tau_k.
function solve_rule(pc, pu, pf; τ_n_target = e_init.τ_n)
    function resid(X)
        dTl  = pc*X; dTul = pu*X
        τ_k  = τ_k_for_subsidy(base.fsub + pf*X; dTlev = dTl, dTulev = dTul)
        a    = period0(; τ_k, dTlev = dTl, dTulev = dTul)
        return a.τ_n - τ_n_target, τ_k, a
    end
    fX = (X) -> resid(X)[1]
    X = find_zero(fX, (0.0, 0.05), Bisection(); xatol = 1e-14)
    _, τ_k, a = resid(X)
    return (; X, τ_k, a, pc, pu, pf)
end

println("\n", "="^100)
println("PART B -- the four allocation rules, each sized to close the j=1 budget")
println("="^100)
den3 = E2.checks + E2.ui + E2.firms
# 1: proportional to the three existing shares, renormalised to 100
R1 = (E2.checks/den3, E2.ui/den3, E2.firms/den3)
# 2: state+local 13.1, education 5.5, public services 1.8 -> labour (split checks:UI
#    in their existing 21.4:19.4 ratio); fed admin 1.8 + global 0.3 -> firms
lab_share = (13.1 + 5.5 + 1.8)/22.5; firm_share = (1.8 + 0.3)/22.5
R2 = (lab_share*E2.checks/(E2.checks+E2.ui), lab_share*E2.ui/(E2.checks+E2.ui), firm_share)
# 4: nothing to firms (minimises the tau_k shock), checks:UI in their existing ratio
R4 = (E2.checks/(E2.checks+E2.ui), E2.ui/(E2.checks+E2.ui), 0.0)

rules = [("1  proportional to existing shares", R1),
         ("2  E2 internal logic (labour-weighted)", R2),
         ("4  nothing to firms (min tau_k shock)", R4)]

results = Any[]
for (lbl, R) in rules
    r = solve_rule(R...)
    push!(results, (lbl, r))
end

# 3: route the tax-credit slice through the labour tax, rest by rule 2.
#    Target tau_n so the retained labour-tax cut equals the tax-credit amount.
P_full = base.mod3/0.771
L_credit = (TAX_CREDITS_PCT/100)*P_full
R3 = R2
function solve_rule3(pc, pu, pf, Lkeep)
    function resid(X)
        dTl = pc*X; dTul = pu*X
        τ_k = τ_k_for_subsidy(base.fsub + pf*X; dTlev = dTl, dTulev = dTul)
        a   = period0(; τ_k, dTlev = dTl, dTulev = dTul)
        return a.hole - Lkeep, τ_k, a
    end
    fX = (X) -> resid(X)[1]
    X = find_zero(fX, (0.0, 0.05), Bisection(); xatol = 1e-14)
    _, τ_k, a = resid(X)
    return (; X, τ_k, a, pc, pu, pf)
end
r3 = solve_rule3(R3..., L_credit)
push!(results, ("3  rule 2 + tax credits via the labour tax", r3))

@printf("\nthe hole to be filled at j=1 (base run) = %.4f%% ann GDP = \$%.0f bn\n",
        pct(base.hole), dol(base.hole))
@printf("Table E2 tax credits = %.1f%% of the full package = %.4f%% ann GDP = \$%.0f bn\n\n",
        TAX_CREDITS_PCT, pct(L_credit), dol(L_credit))

@printf("%-42s %8s %8s %8s %10s %10s %12s %12s\n",
        "rule", "checks", "UI", "firms", "X (\$bn)", "tau_k", "tau_n", "shock %GDP")
for (lbl, r) in results
    @printf("%-42s %8.3f %8.3f %8.3f %10.0f %10.4f %12.6f %12.4f\n",
            lbl, r.pc, r.pu, r.pf, dol(r.X), r.τ_k, r.a.τ_n, pct(r.a.fsub))
end
println("\n'shock %GDP' is the subsidy that the no-capital-subsidy counterfactual has to")
println("remove, i.e. the size of the tau_k shock. Baseline (no Other) is 7.7015%.")

println("\n---- the cost of each rule, relative to the min-shock rule ----")
minshock = minimum(pct(r.a.fsub) for (_, r) in results)
for (lbl, r) in results
    @printf("%-42s shock %.4f%%  = %+.4f pp vs min  (%+.1f%% larger)\n",
            lbl, pct(r.a.fsub), pct(r.a.fsub) - minshock,
            100*(pct(r.a.fsub)/minshock - 1))
end

println("\n---- resulting composition vs Table E2, as shares of the FULL package ----")
@printf("%-42s %9s %9s %9s\n", "rule", "checks", "UI", "firms")
@printf("%-42s %9.2f %9.2f %9.2f\n", "Table E2 target", E2.checks, E2.ui, E2.firms)
for (lbl, r) in results
    # With Other distributed across the three instruments, the modeled package IS the
    # full package, so mod3 is the right denominator for a Table E2 comparison.
    Pf = r.a.mod3
    @printf("%-42s %9.2f %9.2f %9.2f\n", lbl,
            100*r.a.dT1/Pf, 100*r.a.dTu1/Pf, 100*r.a.fsub/Pf)
end
println("(with Other distributed, the modeled package IS the full package, so these are")
println(" directly comparable to Table E2's own shares.)")

# ============================================================
println("\n", "="^100)
println("PART C -- per-decile incidence of each channel (initial-SS distribution)")
println("="^100)
pdfm = sum(e_init.ws.pdf, dims=2)[:,1]
a_grid = round.(e_init.h.a_grid, digits=10)
qu = round.([quantile(a_grid, StatsBase.ProbabilityWeights(pdfm), x) for x in 0.0:0.1:1.0], digits=10)
zg = e_init.h.z_grid
ind_T  = e_init.g.ind_T[:]  .== 1.0
ind_Tu = e_init.g.ind_Tu[:] .== 1.0
nt = e_init.ws.n_pol[1,:]
w0 = e_init.w; τ = g.τ

# incidence weights per (a,z) cell for a $1 aggregate increase in each channel
wT  = zeros(size(e_init.ws.pdf)); wTu = zeros(size(e_init.ws.pdf))
wLT = zeros(size(e_init.ws.pdf))
for i in axes(wT,1), s in axes(wT,2)
    wT[i,s]  = ind_T[s]  ? 1.0 : 0.0                       # equal per recipient
    wTu[i,s] = ind_Tu[s] ? 1.0 : 0.0
    wLT[i,s] = zg[s] == 0 ? 0.0 : (w0*zg[s]*nt[s])^(1-τ)   # HSV: d(net)/d(tau_n)
end
function decile_shares(wt)
    sh = Float64[]
    for (q1,q2) in zip(qu, qu[2:end])
        q2_ = (q2 == qu[end]) ? q2 + 10*eps() : q2
        acc = 0.0
        for i in axes(wt,1), s in axes(wt,2)
            (q1 <= a_grid[i] < q2_) && (acc += e_init.ws.pdf[i,s]*wt[i,s])
        end
        push!(sh, acc)
    end
    return 100 .* sh ./ sum(sh)
end
sT = decile_shares(wT); sTu = decile_shares(wTu); sLT = decile_shares(wLT)
println("share of each channel's total that accrues to each initial-SS wealth decile (%):")
@printf("%6s %12s %12s %14s\n", "decile", "checks (T)", "UI (Tu)", "labour tax")
for k in 1:10
    @printf("%6d %12.2f %12.2f %14.2f\n", k, sT[k], sTu[k], sLT[k])
end
@printf("\nbottom 3 deciles: checks %.1f%%   UI %.1f%%   labour tax %.1f%%\n",
        sum(sT[1:3]), sum(sTu[1:3]), sum(sLT[1:3]))
@printf("top 3 deciles   : checks %.1f%%   UI %.1f%%   labour tax %.1f%%\n",
        sum(sT[8:10]), sum(sTu[8:10]), sum(sLT[8:10]))
@printf("\nunemployed (z=0) mass = %.4f; they receive UI but NOTHING from the labour tax.\n",
        sum(e_init.ws.pdf[:, zg .== 0]))
println("The capital subsidy is not a transfer: it works through w and k, so its incidence")
println("is on labour and capital income rather than on a recipient list.")
println("\nDONE (nothing written)")
