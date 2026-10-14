# ============================================================
# [SW14] THE MISSING "OTHER" COMPONENT  --  period-0 budget, nothing re-solved
# ============================================================
#     julia --project -t auto scripts/fiscal_other_sw.jl
#
# CONTEXT (2026-09-09, corrects the earlier "G is broken" diagnosis).
#   Table E2 (Appendix E.5, from PRAC) classifies the package as Households /
#   Firms / Other.  Its "Other" 22.9% = state and local governments 13.1%,
#   education 5.5%, public services 1.8%, federal programme administration 1.8%,
#   global assistance 0.3%.  That is grants-in-aid and transfers, NOT NIPA federal
#   consumption of goods and services, which is what G is in this model.
#
#   G is CORRECT and is not the problem: G_ss/y_ss = 0.050000 exactly, and the
#   ratio path is 6.00% at 2020Q2 then 5.33, 5.26, 5.48, 5.23, 5.03 ... i.e. 5%
#   rising to ~6% and decaying, exactly as Section 4.2 describes and as Figure
#   D1(a) plots.  An earlier note here claimed G was ~40x too small; that was a
#   LEVEL deviation from the initial steady state, a different object from the
#   ratio the paper targets, and the apparent 22.70/22.9 match against Table E2
#   was a coincidence between two unrelated quantities.  Disregard it.
#
#   THE ACTUAL PROBLEM: "Other" has no counterpart in the model at all.  The debt
#   path is fed from observed federal debt held by the public, which financed the
#   WHOLE package including Other, so the model borrows for spending it never
#   does, and tau_n -- the only free instrument in the GBC -- absorbs the residue.
#
# WHAT THIS SCRIPT DOES
#   PART 1  component shares against the MODELED package (checks + UI + firms)
#   PART 2  the tau_k that hits the firms target under that denominator
#   PART 3  an added "Other" outlay X: tau_n at j=1, and whether it still crosses 1
#   PART 4  the four options (a)-(d), with the tau_n each delivers
#
#   Everything is the PERIOD-0 GBC, holding r_k, the debt path, pi, hours and the
#   household side at their baseline values (w and k respond to tau_k through the
#   firm's FOC).  A reproduction check against the saved tau_n[1] runs first.
#   Welfare cannot be read off this -- it needs a re-solve -- so PART 4 reports
#   the budget consequence and the SIGN of the welfare channel, not CE numbers.
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

# ---- knobs ----------------------------------------------------------------
GDP_2019 = 21380.0
# Table E2, Appendix E.5 (PRAC).  Note the four shares sum to 94.5%, not 100%.
E2 = (checks = 21.4, ui = 19.4, firms = 30.8, other = 22.9)
E2_OTHER_DETAIL = (state_local = 13.1, education = 5.5, public_services = 1.8,
                   fed_admin = 1.8, global_assistance = 0.3)

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

caprev_ss = e_init.τ_k * (e_init.r + δ) * e_init.k
T_ss = e_init.lump_sum_T; Tu_ss = e_init.lump_sum_Tu; G_ss = e_init.G
labtax_ss = e_init.lab_tax
coupon_ss = (ī + g.ϑ)/(1 + π_bar) * e_init.b

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
z_state[1,:] = transpose(e_init.h.Pss) * h_path[1].P
for i = 2:lastindex(h_path); z_state[i,:] = transpose(z_state[i-1,:]) * h_path[i].P; end

h0 = h_path[2]; z_dist0 = z_state[1,:]
n_ss = e_init.ws.n_pol[1,:]; s_0 = ts.Γ_n
r_k0 = p.r_k[1]; r_b0 = p.r[1]
qb0 = ts.qb_path[1]; qbp0 = ts.qb_path[2]; yp0 = p.y[2]
T0 = ts.T_path[1]; Tu0 = ts.Tu_path[1]; G0 = ts.G_path[1]
d0 = ts.d_path[1]; cov0 = ts.covid_path[1]

# ---- the period-0 GBC, with an extra pure outlay X and optional T/Tu top-ups ----
function period0(; τ_k = -0.47, X = 0.0, dT_extra = 0.0, dTu_extra = 0.0)
    pr = ProgressUnknown()
    Tuse  = T0  + dT_extra
    Tuuse = Tu0 + dTu_extra
    Guse  = G0 + X                      # period 0 books G as a LEVEL, so X adds here
    f = (τ_n) -> _τ_n_excess_t0_port_sw_!(h0, t, g; τ_n, τ_k, T = Tuse, Tu = Tuuse,
        G = Guse, qbprime = qbp0, qb = qb0, yprime = yp0, d = d0, s = s_0, n_ss,
        covid = cov0, z_dist = z_dist0, r_b = r_b0, r_k = r_k0, verbose = false, p = pr).excess_g
    τ_n = _solve_tau_n_sw(f, 0.0, 8.0)
    out = _τ_n_excess_t0_port_sw_!(h0, t, g; τ_n, τ_k, T = Tuse, Tu = Tuuse, G = Guse,
        qbprime = qbp0, qb = qb0, yprime = yp0, d = d0, s = s_0, n_ss, covid = cov0,
        z_dist = z_dist0, r_b = r_b0, r_k = r_k0, verbose = false, p = pr)
    fsub = -(τ_k * r_k0 * out.k - caprev_ss)
    dT1  = out.lump_sum_T  - T_ss
    dTu1 = out.lump_sum_Tu - Tu_ss
    mod3 = dT1 + dTu1 + fsub
    return (; τ_n = out.τ_n, y = out.y, k = out.k, w = out.w, fsub, dT1, dTu1, mod3, X,
            share_firms = 100*fsub/mod3, share_checks = 100*dT1/mod3, share_ui = 100*dTu1/mod3,
            lab_tax = out.lab_tax)
end

base = period0()
println("\nreproduction check: tau_n = ", base.τ_n, "   saved = ", p.τ_n[1],
        "   diff = ", base.τ_n - p.τ_n[1])
@assert abs(base.τ_n - p.τ_n[1]) < 1e-8 "period-0 GBC not reproduced"
println("period-0 GBC reproduced exactly.\n")

# ============================================================
println("="^100)
println("PART 1 -- shares against the MODELED package (checks + UI + firms)")
println("="^100)
tgt771 = (checks = E2.checks/0.771, ui = E2.ui/0.771, firms = E2.firms/0.771)
den716 = E2.checks + E2.ui + E2.firms
tgt716 = (checks = 100*E2.checks/den716, ui = 100*E2.ui/den716, firms = 100*E2.firms/den716)
@printf("%-12s %10s %12s %10s %12s %10s\n", "", "model %", "tgt /77.1", "diff", "tgt /71.6", "diff")
@printf("%-12s %10.2f %12.2f %+10.2f %12.2f %+10.2f\n", "checks",
        base.share_checks, tgt771.checks, base.share_checks-tgt771.checks,
        tgt716.checks, base.share_checks-tgt716.checks)
@printf("%-12s %10.2f %12.2f %+10.2f %12.2f %+10.2f\n", "UI",
        base.share_ui, tgt771.ui, base.share_ui-tgt771.ui, tgt716.ui, base.share_ui-tgt716.ui)
@printf("%-12s %10.2f %12.2f %+10.2f %12.2f %+10.2f\n", "firms",
        base.share_firms, tgt771.firms, base.share_firms-tgt771.firms,
        tgt716.firms, base.share_firms-tgt716.firms)
@printf("\n/77.1 targets sum to %.2f%% (Table E2 four shares add to 94.5%%, so 5.5pp is\n",
        tgt771.checks+tgt771.ui+tgt771.firms)
println("unallocated and that gap survives the renormalization). /71.6 sums to 100% by construction.")

# ============================================================
println("\n", "="^100)
println("PART 2 -- the tau_k that hits the firms target, modeled-package denominator")
println("="^100)
@printf("%9s %12s %12s %14s %12s\n", "tau_k", "tau_n", "firms %", "package %GDP", "subsidy %GDP")
for τ_k in -0.62:0.02:-0.36
    a = period0(; τ_k)
    @printf("%9.4f %12.6f %12.2f %14.4f %12.4f\n",
            τ_k, a.τ_n, a.share_firms, pct(a.mod3), pct(a.fsub))
end
for (lbl, tgt) in (("/77.1  (39.95%)", tgt771.firms), ("/71.6  (43.02%)", tgt716.firms))
    fk = (τ_k) -> period0(; τ_k).share_firms - tgt
    z = try find_zero(fk, (-0.65, -0.30), Bisection(); xatol = 1e-8) catch e; NaN end
    if isnan(z)
        println("no tau_k in [-0.65,-0.30] hits the ", lbl, " target")
    else
        a = period0(; τ_k = z)
        @printf("\ntarget %s -> tau_k = %+.6f   (paper -0.47, diff %+.4f)\n", lbl, z, z + 0.47)
        @printf("   tau_n = %.6f (%+.2f%% vs SS)   subsidy = %.4f%% GDP = \$%.0f bn\n",
                a.τ_n, 100*(a.τ_n/e_init.τ_n - 1), pct(a.fsub), dol(a.fsub))
    end
end

# ============================================================
println("\n", "="^100)
println("PART 3 -- adding an 'Other' outlay X: what happens to tau_n at j=1")
println("="^100)
P771 = base.mod3/0.771
X_E2 = 0.229*P771                       # Other at its Table E2 size
@printf("if the three modeled components are 77.1%% of the package:\n")
@printf("  implied package = %.4f%% ann GDP (\$%.0f bn);  Other = %.4f%% (\$%.0f bn)\n",
        pct(P771), dol(P771), pct(X_E2), dol(X_E2))
@printf("\n%14s %14s %12s %10s %14s\n", "X (%ann GDP)", "X (\$bn)", "tau_n", ">1?", "vs SS (%)")
Xgrid = [0.0, 0.005, 0.010, 0.015, 0.0198, 0.025, 0.030, 0.040, X_E2]
for X in sort(unique(Xgrid))
    a = period0(; X)
    @printf("%14.4f %14.0f %12.6f %10s %+14.2f\n",
            pct(X), dol(X), a.τ_n, a.τ_n > 1 ? "YES" : "no", 100*(a.τ_n/e_init.τ_n - 1))
end
fX = (X) -> period0(; X).τ_n - e_init.τ_n
X_flat = try find_zero(fX, (0.0, 0.06), Bisection(); xatol = 1e-14) catch e; NaN end
fX1 = (X) -> period0(; X).τ_n - 1.0
X_one = try find_zero(fX1, (0.0, 0.06), Bisection(); xatol = 1e-14) catch e; NaN end
if !isnan(X_flat)
    @printf("\nX that holds tau_n at its SS value  = %.4f%% ann GDP = \$%.0f bn\n",
            pct(X_flat), dol(X_flat))
end
if !isnan(X_one)
    @printf("X that just brings tau_n down to 1  = %.4f%% ann GDP = \$%.0f bn\n",
            pct(X_one), dol(X_one))
end
aE2 = period0(; X = X_E2)
@printf("\nAT THE FULL TABLE E2 SIZE, X = %.4f%% (\$%.0f bn):  tau_n = %.6f  (%+.2f%% vs SS)\n",
        pct(X_E2), dol(X_E2), aE2.τ_n, 100*(aE2.τ_n/e_init.τ_n - 1))
println(aE2.τ_n < e_init.τ_n ?
    "   -> OVERSHOOTS: tau_n falls BELOW its steady-state value, i.e. the labour tax RISES." :
    "   -> still above the steady-state value.")
@printf("   Other is %.2fx the X that would hold tau_n flat.\n", X_E2/max(X_flat, eps()))

# ============================================================
println("\n", "="^100)
println("PART 4 -- the four options")
println("="^100)
println("Options (a), (b) and (c) are IDENTICAL in the period-0 budget for a given total")
println("outlay X: the GBC only sees the total. They differ in the RESOURCE CONSTRAINT and")
println("in household budgets, which is where welfare comes from:")
println("  path.C = y - X_inv - G - (foreign block), so an outlay booked like G reduces C")
println("  one-for-one, while T and Tu enter household budgets and are redistributed.")
println()
for (lbl, X) in (("(a)/(b) at the tau_n-flat size", isnan(X_flat) ? 0.0 : X_flat),
                 ("(a)/(b) at the full E2 size", X_E2))
    a = period0(; X)
    @printf("%-34s X = %7.4f%% GDP (\$%5.0f bn) -> tau_n = %.6f (%s1)\n",
            lbl, pct(X), dol(X), a.τ_n, a.τ_n > 1 ? ">" : "<")
end
# (c) distribute X across the three modeled components in Table E2 proportions
for (lbl, X) in (("(c) distributed, tau_n-flat size", isnan(X_flat) ? 0.0 : X_flat),
                 ("(c) distributed, full E2 size", X_E2))
    # checks and UI are ratios-to-GDP in T_path/Tu_path, so convert the level share
    shr_c = E2.checks/den716; shr_u = E2.ui/den716; shr_f = E2.firms/den716
    dTx  = shr_c * X / p.y[1]      # T enters as T*y, y ~ p.y[1] at j=1
    dTux = shr_u * X / p.y[1]
    # the firms slice is delivered by a more negative tau_k; solve for it
    want = base.fsub + shr_f * X
    fk = (τ_k) -> period0(; τ_k, dT_extra = dTx, dTu_extra = dTux).fsub - want
    zk = try find_zero(fk, (-0.95, -0.30), Bisection(); xatol = 1e-8) catch e; NaN end
    if isnan(zk)
        @printf("%-34s no tau_k in [-0.95,-0.30] delivers the firms slice\n", lbl)
    else
        a = period0(; τ_k = zk, dT_extra = dTx, dTu_extra = dTux)
        @printf("%-34s X = %7.4f%% GDP -> tau_k = %+.4f, tau_n = %.6f (%s1)\n",
                lbl, pct(X), zk, a.τ_n, a.τ_n > 1 ? ">" : "<")
    end
end
a_d = period0()
@printf("%-34s X = 0                        -> tau_n = %.6f (>1)   [status quo]\n",
        "(d) leave it unmodeled", a_d.τ_n)

println("\n---- what 'Other' actually is, which bears on (a) vs (c) ----")
tot = E2_OTHER_DETAIL.state_local + E2_OTHER_DETAIL.education +
      E2_OTHER_DETAIL.public_services + E2_OTHER_DETAIL.fed_admin +
      E2_OTHER_DETAIL.global_assistance
@printf("state and local governments %.1f, education %.1f, public services %.1f,\n",
        E2_OTHER_DETAIL.state_local, E2_OTHER_DETAIL.education, E2_OTHER_DETAIL.public_services)
@printf("federal administration %.1f, global assistance %.1f   (sum %.1f of the package)\n",
        E2_OTHER_DETAIL.fed_admin, E2_OTHER_DETAIL.global_assistance, tot)
@printf("state+local+education = %.1f of %.1f = %.0f%% of Other -- spending that DOES reach\n",
        E2_OTHER_DETAIL.state_local + E2_OTHER_DETAIL.education, tot,
        100*(E2_OTHER_DETAIL.state_local + E2_OTHER_DETAIL.education)/tot)
@printf("households and firms, just not through the federal budget. Only fed admin + global\n")
@printf("assistance (%.1f, %.0f%% of Other) is plausibly a pure outlay with no household counterpart.\n",
        E2_OTHER_DETAIL.fed_admin + E2_OTHER_DETAIL.global_assistance,
        100*(E2_OTHER_DETAIL.fed_admin + E2_OTHER_DETAIL.global_assistance)/tot)
println("\nDONE (nothing written)")
