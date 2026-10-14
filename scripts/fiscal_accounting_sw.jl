# ============================================================
# [SW12] FISCAL ACCOUNTING OFF THE SAVED BASELINE  --  nothing re-solved
# ============================================================
#     julia --project -t auto scripts/fiscal_accounting_sw.jl
#
# Answers, purely from output/jld2_sw/reduced_epsm.jld2 plus the initial SS:
#   PART A  the government budget identity, verified against the model's own
#           bookkeeping so the accounting below is the model's, not a new one
#   PART B  new borrowing q_t*l_t against each spending component, at 2020Q2 and
#           cumulated over the stimulus window; shares vs Table E2
#   PART C  how far tau_n actually moves in 2020Q2
#   PART D  the tau_k question: what hits the 31% target, what -0.47 delivers, and
#           how tau_n moves under each.  PART D evaluates the PERIOD-0 GBC only,
#           holding r_k, the debt path, pi, hours and the household side at their
#           baseline values.  It is not a re-solve and is labelled as such.
#
# THE IDENTITY (derived from _tau_n_excess_port_sw_! and its period-0 twin).
# Interior j:
#   0 = lab_tax + tau_k*(r+delta)*k + q*b' - (1+r)*q_lag*b - T - Tu - G*y
# Period 0 (qb and G enter as LEVELS, not ratios):
#   0 = lab_tax + tau_k*r_k*k     + q*b' - (1+r_b)*qb    - T - Tu - G
# Using (1+r)*q_lag = (ibar+vartheta)/(1+pi) + (1-vartheta)/(1+pi)*q and
# l = b' - (1-vartheta)/(1+pi)*b, both rearrange to the SAME statement:
#
#   q*l  =  coupon service  +  T + Tu + G  -  lab_tax  -  tau_k*(base)*k
#
# so NEW BORROWING finances the coupon bill plus primary outlays, net of labour
# tax revenue and net of capital tax revenue (tau_k < 0 makes that last term an
# OUTLAY, the capital subsidy).
#
# WHAT THAT MEANS FOR "DO THE COMPONENTS ADD UP".  In deviations from the initial
# steady state the identity is
#
#   d(q*l) = d(coupon) + dT + dTu + dG - d(lab_tax) - d(tau_k*base*k)
#
# The model SOLVES the GBC for tau_n every period, so this holds exactly by
# construction and the components cannot fail to add up.  The meaningful question
# is therefore not whether there is an accounting error but HOW THE GAP between
# new borrowing and the four named components is split between the two unnamed
# terms: d(coupon), which is mechanical given the debt path and inflation, and
# d(lab_tax), which is the residual tau_n absorbs.  That split is what PART B
# reports.
#
# ONE THING THIS SETTLES.  G_path[1] = e_init.G * 1.1 was flagged as a suspected
# level/ratio bug.  It is NOT a bug: the period-0 helper subtracts `- G`, a LEVEL,
# while interior periods subtract `- G*y`, and path.G[1] = G confirms it.  So
# 2020Q2 government consumption is 10% above its steady-state LEVEL, exactly as
# intended, and the earlier "-1.44% of annual GDP" reading of the G component came
# from treating that entry as a ratio.
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

# ---- knobs ----------------------------------------------------------------
BASE_TAG   = "reduced_epsm"     # the epsm-fed baseline
GDP_2019   = 21380.0            # nominal US GDP 2019, $bn (BEA).  Only a scale for the
                                # dollar columns; every share is scale-free.
WINDOW     = 19                 # stimulus window j = 1..19 = 2020Q2..2024Q4, the span
                                # over which make_time_series_sw sets G by hand
TABLE_E2   = (checks = 21.4, ui = 19.4, firms = 30.8, other = 22.9)   # paper Table E2, %

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

println("loading ", BASE_TAG, ".jld2 ..."); flush(stdout)
o = load_object(joinpath(ROOT, "output", "jld2_sw", BASE_TAG * ".jld2"))
e_init, e_final, p, π_endog = o[1], o[2], o[3], o[4]
h = e_init.h; t = e_init.t; g = e_init.g
δ = t.δ
π_bar = e_init.π_
ī = (1 + e_init.r) * (1 + π_bar) - 1
ϑ = g.ϑ
n_per = length(p.r)

y_ss     = y(e_init)
Y_ANN    = 4 * y_ss                      # annual 2019 GDP, model units
pct(x)   = 100 * x / Y_ANN               # share of annual 2019 GDP
dol(x)   = GDP_2019 * x / Y_ANN          # $bn

# ============================================================
# PART A -- the identity, verified against the model
# ============================================================
coupon(j) = (ī + ϑ) / (1 + p.π_[j]) * p.b[j]
capbase(j) = j == 1 ? p.r_k[1] : (p.r[j] + δ)          # period 0 uses gross r_k
caprev(j)  = p.τ_k[j] * capbase(j) * p.k[j]            # REVENUE (tau_k<0 => subsidy)
qell(j)    = p.q[j] * p.l[j]

# initial steady state, same objects
coupon_ss  = (ī + ϑ) / (1 + π_bar) * e_init.b
caprev_ss  = e_init.τ_k * (e_init.r + δ) * e_init.k
ell_ss     = e_init.b - (1 - ϑ) / (1 + π_bar) * e_init.b
qell_ss    = e_init.q * ell_ss
T_ss       = e_init.lump_sum_T
Tu_ss      = e_init.lump_sum_Tu
G_ss       = e_init.G
labtax_ss  = e_init.lab_tax

println("\n", "="^100)
println("PART A -- government budget identity:  q*l = coupon + T + Tu + G - lab_tax - tau_k*base*k")
println("="^100)
@printf("%5s %14s %14s %14s\n", "j", "q*l (model)", "q*l (identity)", "residual")
for j in vcat(collect(1:6), [12, 19, 40, 100, 400])
    rhs = coupon(j) + p.lump_sum_T[j] + p.lump_sum_Tu[j] + p.G[j] - p.lab_tax[j] - caprev(j)
    @printf("%5d %14.10f %14.10f %14.3e\n", j, qell(j), rhs, qell(j) - rhs)
end
worst = maximum(abs(qell(j) - (coupon(j) + p.lump_sum_T[j] + p.lump_sum_Tu[j] + p.G[j]
                               - p.lab_tax[j] - caprev(j))) for j in 1:(n_per-1))
@printf("worst |residual| over j = 1..%d : %.3e\n", n_per-1, worst)
rhs_ss = coupon_ss + T_ss + Tu_ss + G_ss - labtax_ss - caprev_ss
@printf("initial SS: q*l = %.10f   identity = %.10f   residual = %.3e\n",
        qell_ss, rhs_ss, qell_ss - rhs_ss)
println(worst < 1e-9 ? "IDENTITY HOLDS -- the accounting below is the model's own." :
        "*** IDENTITY FAILS -- do not trust the decomposition ***")

println("\ninitial steady-state levels (model units, and % of annual 2019 GDP):")
for (nm, v) in (("coupon service", coupon_ss), ("T (checks)", T_ss), ("Tu (UI)", Tu_ss),
                ("G", G_ss), ("labour tax revenue", labtax_ss),
                ("capital tax revenue", caprev_ss), ("q*l new borrowing", qell_ss))
    @printf("  %-22s %14.8f %10.4f%%\n", nm, v, pct(v))
end
@printf("  check: Lambda_f target 4%% of quarterly GDP -> caprev_ss/y_ss = %.6f\n",
        caprev_ss / y_ss)
@printf("  check: Lambda_g target 5%% of quarterly GDP -> G_ss/y_ss      = %.6f\n",
        G_ss / y_ss)

# ============================================================
# PART B -- components against new borrowing
# ============================================================
# deviations from the initial steady state
dT(j)      = p.lump_sum_T[j]  - T_ss
dTu(j)     = p.lump_sum_Tu[j] - Tu_ss
dG(j)      = p.G[j]           - G_ss
dsub(j)    = -(caprev(j) - caprev_ss)      # OUTLAY: positive = net subsidy to firms
dcoupon(j) = coupon(j)        - coupon_ss
dlabtax(j) = p.lab_tax[j]     - labtax_ss
dqell(j)   = qell(j)          - qell_ss

function report_window(js, label)
    sT  = sum(dT, js);      sTu = sum(dTu, js)
    sS  = sum(dsub, js);    sG  = sum(dG, js)
    sC  = sum(dcoupon, js); sL  = sum(dlabtax, js)
    sQ  = sum(dqell, js)
    pkg = sT + sTu + sS + sG
    gap = sQ - pkg
    println("\n", "-"^100)
    println(label)
    println("-"^100)
    @printf("%-34s %14s %12s %14s\n", "", "model units", "% ann GDP", "\$bn")
    @printf("%-34s %14.8f %12.4f %14.1f\n", "NEW BORROWING  d(q*l)", sQ, pct(sQ), dol(sQ))
    println()
    @printf("%-34s %14.8f %12.4f %14.1f\n", "stimulus checks   dT", sT, pct(sT), dol(sT))
    @printf("%-34s %14.8f %12.4f %14.1f\n", "unemployment ins. dTu", sTu, pct(sTu), dol(sTu))
    @printf("%-34s %14.8f %12.4f %14.1f\n", "capital subsidy   -d(tau_k*base*k)", sS, pct(sS), dol(sS))
    @printf("%-34s %14.8f %12.4f %14.1f\n", "govt consumption  dG", sG, pct(sG), dol(sG))
    @printf("%-34s %14.8f %12.4f %14.1f\n", "  SUM OF COMPONENTS", pkg, pct(pkg), dol(pkg))
    println()
    @printf("%-34s %14.8f %12.4f %14.1f\n", "GAP = borrowing - components", gap, pct(gap), dol(gap))
    @printf("%-34s %14.8f %12.4f %14.1f\n", "  of which d(coupon service)", sC, pct(sC), dol(sC))
    @printf("%-34s %14.8f %12.4f %14.1f\n", "  of which -d(lab_tax)", -sL, pct(-sL), dol(-sL))
    @printf("%-34s %14.3e\n", "  closure check (should be 0)", gap - (sC - sL))
    println("\nshares of the SUM OF COMPONENTS, against Table E2:")
    @printf("%-22s %10s %10s %10s\n", "", "model %", "Table E2", "diff pp")
    for (nm, v, e2) in (("stimulus checks", sT, TABLE_E2.checks), ("unemployment ins.", sTu, TABLE_E2.ui),
                        ("firms (capital)", sS, TABLE_E2.firms), ("other (G)", sG, TABLE_E2.other))
        @printf("%-22s %10.2f %10.2f %+10.2f\n", nm, 100*v/pkg, e2, 100*v/pkg - e2)
    end
    # If a single package total could reconcile the model with Table E2, then
    # component / E2_share would give the SAME implied total for all four.  How far
    # apart these four numbers are says whether the problem is the total's SIZE or
    # its COMPOSITION.
    println("\npackage total IMPLIED by each component and its Table E2 share (% ann GDP):")
    @printf("  from checks %8.4f   from UI %8.4f   from firms %8.4f   from other/G %8.4f\n",
            pct(sT)/(TABLE_E2.checks/100), pct(sTu)/(TABLE_E2.ui/100),
            pct(sS)/(TABLE_E2.firms/100), pct(sG)/(TABLE_E2.other/100))
    println("  (equal across the four => a pure scaling problem; spread => composition is wrong)")
    return (; sT, sTu, sS, sG, sC, sL, sQ, pkg, gap)
end

println("\n", "="^100)
println("PART B -- new borrowing against the spending components")
println("="^100)
@printf("annual 2019 GDP = %.8f model units = \$%.0f bn (assumed)\n", Y_ANN, GDP_2019)
B1  = report_window(1:1,       "2020Q2 ONLY  (j = 1)")
B12 = report_window(1:12,      "CUMULATIVE j = 1..12   (2020Q2 - 2023Q1)")
BW  = report_window(1:WINDOW,  "CUMULATIVE j = 1..$WINDOW   (2020Q2 - 2024Q4, the hand-set G window)")

# ============================================================
# PART C -- how far tau_n moves in 2020Q2
# ============================================================
println("\n", "="^100)
println("PART C -- the tau_n movement, which was the original selection criterion")
println("="^100)
println("tau_n is the HSV RETENTION LEVEL lambda, not a tax rate: y_net = tau_n * y_gross^(1-tau).")
println("tau_n > 1 means labour income is NET SUBSIDISED at the bottom of the earnings")
println("distribution, since retention exceeds gross income there.\n")
@printf("initial SS tau_n            = %.8f\n", e_init.τ_n)
@printf("2020Q2 (j=1) tau_n          = %.8f\n", p.τ_n[1])
@printf("change                      = %+.8f  (%+.4f%%)\n",
        p.τ_n[1] - e_init.τ_n, 100*(p.τ_n[1]/e_init.τ_n - 1))
@printf("does it exceed 1?           %s\n", p.τ_n[1] > 1 ? "YES -- labour income net subsidised" : "no")
println()
@printf("%5s %14s %14s %14s\n", "j", "tau_n", "vs SS (%)", "lab_tax")
for j in vcat(collect(1:8), [12, 19, 40])
    @printf("%5d %14.8f %+14.4f %14.8f\n", j, p.τ_n[j], 100*(p.τ_n[j]/e_init.τ_n - 1), p.lab_tax[j])
end

# ============================================================
# PART D -- the tau_k question, PERIOD-0 GBC ONLY (no re-solve)
# ============================================================
# Rebuild the exact period-0 arguments.  Verified below by reproducing the
# baseline tau_n from the saved r_k.
b_target = e_init.b / y(e_init)
b_path0 = let
    T=700; H=100; ρB=0.08; bl = Array{Float64,1}(undef, T+H)
    bl[1]=b_target; bl[2]=b_target; bl[T:end] .= b_target
    for i in 3:T-1; bl[i] = bl[1]*ρB^(i-1) + (1-ρB^(i-1))*b_target; end
    bl
end
ts = make_time_series_sw(c, e_init, e_final, b_path0, 0.287*4.0, g)
pars_h = (θ=c.θ, ν=c.ν, ies=c.ies, crra=c.crra, β=c.β, a_min=c.a_min, a_max=c.a_max,
          z_u=0.0, n_a=c.n_a, init_state=e_init.h.Pss, z_scale=z_scale)
h_path, _, _ = h_path_transition_mckay_init(ROOT, ts.Pee_path, ts.Puu_path; pars_h...)

z_state = zeros(length(h_path), length(e_init.h.z_grid))
z_state[1,:] = transpose(e_init.h.Pss) * h_path[1].P
for i = 2:lastindex(h_path)
    z_state[i,:] = transpose(z_state[i-1,:]) * h_path[i].P
end
h0      = h_path[2]
z_dist0 = z_state[1,:]
n_ss    = e_init.ws.n_pol[1,:]
s_0     = ts.Γ_n
r_k0    = p.r_k[1]
r_b0    = p.r[1]
qb0     = ts.qb_path[1]
qbp0    = ts.qb_path[2]
yp0     = p.y[2]
T0      = ts.T_path[1]; Tu0 = ts.Tu_path[1]; G0 = ts.G_path[1]
d0      = ts.d_path[1]; cov0 = ts.covid_path[1]

function period0_at(τ_k)
    pr = ProgressUnknown()
    f = (τ_n) -> _τ_n_excess_t0_port_sw_!(h0, t, g; τ_n, τ_k, T = T0, Tu = Tu0, G = G0,
        qbprime = qbp0, qb = qb0, yprime = yp0, d = d0, s = s_0, n_ss, covid = cov0,
        z_dist = z_dist0, r_b = r_b0, r_k = r_k0, verbose = false, p = pr).excess_g
    τ_n = _solve_tau_n_sw(f, 0.0, 5.0)
    out = _τ_n_excess_t0_port_sw_!(h0, t, g; τ_n, τ_k, T = T0, Tu = Tu0, G = G0,
        qbprime = qbp0, qb = qb0, yprime = yp0, d = d0, s = s_0, n_ss, covid = cov0,
        z_dist = z_dist0, r_b = r_b0, r_k = r_k0, verbose = false, p = pr)
    sub = -(τ_k * r_k0 * out.k - caprev_ss)          # capital subsidy OUTLAY at j=1
    dT1 = out.lump_sum_T - T_ss
    dTu1 = out.lump_sum_Tu - Tu_ss
    dG1 = G0 - G_ss                                   # period 0 books G as a LEVEL
    pkg = dT1 + dTu1 + sub + dG1
    return (; τ_k, τ_n = out.τ_n, k = out.k, y = out.y, w = out.w,
            lab_tax = out.lab_tax, sub, pkg, share = 100 * sub / pkg,
            excess = out.excess_g)
end

println("\n", "="^100)
println("PART D -- tau_k, the 31% target, and the tau_n criterion")
println("="^100)
println("PERIOD-0 GBC ONLY.  r_k, r_b, the debt path, pi, hours (s_0 = Gamma_n) and the")
println("household side are HELD at their baseline values; only tau_k moves, and w and k")
println("respond through the firm's first-order conditions at fixed r_k.  A full re-solve")
println("would move r_k and the whole path, so read these as the period-0 budget answer,")
println("which is exactly the object the old tau_n criterion was about.\n")

chk = period0_at(-0.47)
@printf("reproduction check at tau_k = -0.47:  tau_n = %.8f  (saved %.8f, diff %.3e)\n",
        chk.τ_n, p.τ_n[1], chk.τ_n - p.τ_n[1])
println(abs(chk.τ_n - p.τ_n[1]) < 1e-6 ?
        "period-0 GBC reproduced -- the sweep below is on the model's own budget.\n" :
        "*** period-0 GBC NOT reproduced; treat the sweep with suspicion ***\n")

@printf("%9s %14s %12s %14s %14s %12s\n",
        "tau_k", "tau_n", "vs SS (%)", "subsidy %GDP", "package %GDP", "firms share")
for τ_k in -0.60:0.05:0.15
    a = period0_at(τ_k)
    @printf("%9.4f %14.8f %+12.4f %14.4f %14.4f %12.2f\n",
            τ_k, a.τ_n, 100*(a.τ_n/e_init.τ_n - 1), pct(a.sub), pct(a.pkg), a.share)
end

# the two criteria
f_share = (τ_k) -> period0_at(τ_k).share - 31.0
f_taun  = (τ_k) -> period0_at(τ_k).τ_n - e_init.τ_n
println("\n---- where each criterion lands ----")
τk_31 = try find_zero(f_share, (-0.60, 0.15), Bisection(); xatol = 1e-8) catch e; NaN end
τk_fl = try find_zero(f_taun,  (-0.60, 0.15), Bisection(); xatol = 1e-8) catch e; NaN end
if !isnan(τk_31)
    a = period0_at(τk_31)
    @printf("tau_k hitting the 31%% firms share : %+.6f   -> tau_n = %.8f (%+.4f%% vs SS)\n",
            τk_31, a.τ_n, 100*(a.τ_n/e_init.τ_n - 1))
end
if !isnan(τk_fl)
    a = period0_at(τk_fl)
    @printf("tau_k holding tau_n flat at SS    : %+.6f   -> firms share = %.2f%%\n",
            τk_fl, a.share)
end
a47 = period0_at(-0.47)
@printf("tau_k = -0.47 (the paper)         : firms share = %.2f%%   tau_n = %.8f (%+.4f%% vs SS)\n",
        a47.share, a47.τ_n, 100*(a47.τ_n/e_init.τ_n - 1))
if !isnan(τk_31) && !isnan(τk_fl)
    @printf("\nDO THE TWO CRITERIA COINCIDE?  gap = %.6f in tau_k units -> %s\n",
            τk_31 - τk_fl, abs(τk_31 - τk_fl) < 0.01 ? "YES, to within 0.01" : "NO")
end

# ============================================================
# PART E -- if G is the broken component, how big should it be?
# ============================================================
# PART B identifies government consumption as the component that is wrong: its
# share of the package is 0.6% at j=1 and NEGATIVE cumulatively, against Table E2's
# 22.9%.  So run the test the other way round: hold tau_k at the value that hits the
# 31% firms share, and ask what period-0 G makes tau_n sit at its steady-state
# value.  If the G that clears the budget is ALSO the G that delivers roughly 22.9%
# of the package, then all three conditions are mutually consistent once G is fixed,
# and the tau_k tuning was compensating for a mis-specified G.  If not, something
# else is off too.
#
# Period 0 books G as a LEVEL, so G is a free scalar here.
function period0_at_G(τ_k, G)
    pr = ProgressUnknown()
    f = (τ_n) -> _τ_n_excess_t0_port_sw_!(h0, t, g; τ_n, τ_k, T = T0, Tu = Tu0, G = G,
        qbprime = qbp0, qb = qb0, yprime = yp0, d = d0, s = s_0, n_ss, covid = cov0,
        z_dist = z_dist0, r_b = r_b0, r_k = r_k0, verbose = false, p = pr).excess_g
    τ_n = _solve_tau_n_sw(f, 0.0, 5.0)
    out = _τ_n_excess_t0_port_sw_!(h0, t, g; τ_n, τ_k, T = T0, Tu = Tu0, G = G,
        qbprime = qbp0, qb = qb0, yprime = yp0, d = d0, s = s_0, n_ss, covid = cov0,
        z_dist = z_dist0, r_b = r_b0, r_k = r_k0, verbose = false, p = pr)
    sub  = -(τ_k * r_k0 * out.k - caprev_ss)
    dT1  = out.lump_sum_T  - T_ss
    dTu1 = out.lump_sum_Tu - Tu_ss
    dG1  = G - G_ss
    pkg  = dT1 + dTu1 + sub + dG1
    return (; τ_n = out.τ_n, y = out.y, sub, dG1, pkg,
            share_firms = 100*sub/pkg, share_other = 100*dG1/pkg)
end

println("\n", "="^100)
println("PART E -- if G is the broken component, what G clears the period-0 budget?")
println("="^100)
for (lbl, τ_k) in (("tau_k = -0.47 (paper)", -0.47),
                   ("tau_k at the 31% share", isnan(τk_31) ? -0.3844 : τk_31))
    fG = (G) -> period0_at_G(τ_k, G).τ_n - e_init.τ_n
    Gstar = try find_zero(fG, (G_ss, 60*G_ss), Bisection(); xatol = 1e-12) catch e; NaN end
    println("\n", lbl)
    if isnan(Gstar)
        println("  no G in [G_ss, 60*G_ss] sets tau_n to its steady-state value")
        continue
    end
    a = period0_at_G(τ_k, Gstar)
    @printf("  G that holds tau_n at SS      = %.8f   (= %.2f x G_ss, ratio to y_1 = %.4f)\n",
            Gstar, Gstar/G_ss, Gstar/a.y)
    @printf("  baseline G in the model       = %.8f   (= %.2f x G_ss, ratio to y_1 = %.4f)\n",
            G0, G0/G_ss, G0/p.y[1])
    @printf("  dG at that G                  = %.8f  = %.4f%% ann GDP = \$%.1f bn\n",
            a.dG1, pct(a.dG1), dol(a.dG1))
    @printf("  package total                 = %.4f%% ann GDP\n", pct(a.pkg))
    @printf("  firms share                   = %.2f%%   (target 31.0%%)\n", a.share_firms)
    @printf("  other/G share                 = %.2f%%   (Table E2: %.1f%%)\n",
            a.share_other, TABLE_E2.other)
    @printf("  tau_n                         = %.8f  (SS %.8f)\n", a.τ_n, e_init.τ_n)
end

# ============================================================
# PART F -- solve (tau_k, G) JOINTLY for tau_n flat AND firms share = 30.8%
# ============================================================
# PART E shows that the G which clears the budget at the 31%-consistent tau_k also
# reproduces Table E2's "other" share almost exactly (22.70% vs 22.9%).  But adding
# that G enlarges the package, which pushes the firms share DOWN, so tau_k and G have
# to be chosen together.  Two unknowns, two conditions:
#
#   (i)  tau_n = tau_n_ss              the budget closes without moving the labour tax
#   (ii) firms share = 30.8%           Table E2 composition for the capital subsidy
#
# Table E2 has FOUR shares, so imposing two of them leaves the other two as OUTPUT.
# Whether those come out near 21.4% and 22.9% is the real test of whether the
# decomposition can be made consistent at all.
# Bracket G tightly.  A very small G leaves so much spare borrowing that the required
# tau_n exceeds the [0,5] bracket of _solve_tau_n_sw, the [SW8] grid fallback fires and
# the returned tau_n no longer satisfies the GBC -- which silently breaks any outer
# root-find.  Keeping G in [0.5, 12] x G_ss stays well inside the well-posed region.
function G_clearing(τ_k)
    fG = (G) -> period0_at_G(τ_k, G).τ_n - e_init.τ_n
    lo, hi = 0.5 * G_ss, 12.0 * G_ss
    (sign(fG(lo)) == sign(fG(hi))) && return NaN
    return find_zero(fG, (lo, hi), Bisection(); xatol = 1e-14)
end

println("\n", "="^100)
println("PART F -- (tau_k, G) chosen jointly: tau_n flat AND firms share 30.8%")
println("="^100)
println("tau_k grid, with G set at each point to close the budget with tau_n flat:")
@printf("%9s %14s %10s %14s %12s %12s %12s\n",
        "tau_k", "G (clearing)", "G/G_ss", "package %GDP", "firms %", "other/G %", "checks+UI %")
dT1 = p.lump_sum_T[1] - T_ss; dTu1 = p.lump_sum_Tu[1] - Tu_ss
grid_out = NamedTuple[]
for τ_k in -0.550:0.025:-0.300
    Gc = G_clearing(τ_k)
    if isnan(Gc)
        @printf("%9.4f %14s\n", τ_k, "no clearing G in [0.5,12]xG_ss")
        continue
    end
    a = period0_at_G(τ_k, Gc)
    push!(grid_out, (; τ_k, Gc, share_firms = a.share_firms, share_other = a.share_other,
                     pkg = a.pkg, τ_n = a.τ_n))
    @printf("%9.4f %14.8f %10.2f %14.4f %12.2f %12.2f %12.2f\n",
            τ_k, Gc, Gc/G_ss, pct(a.pkg), a.share_firms, a.share_other,
            100*(dT1+dTu1)/a.pkg)
end

τk_joint = NaN
if length(grid_out) >= 2
    lo = findlast(r -> r.share_firms >= TABLE_E2.firms, grid_out)
    hi = findfirst(r -> r.share_firms <= TABLE_E2.firms, grid_out)
    if lo !== nothing && hi !== nothing
        f_joint = (τ_k) -> (a = period0_at_G(τ_k, G_clearing(τ_k)); a.share_firms - TABLE_E2.firms)
        a_, b_ = min(grid_out[lo].τ_k, grid_out[hi].τ_k), max(grid_out[lo].τ_k, grid_out[hi].τ_k)
        τk_joint = try find_zero(f_joint, (a_, b_), Bisection(); xatol = 1e-7) catch e; NaN end
    end
end

if isnan(τk_joint)
    println("\nno tau_k on this grid satisfies both conditions")
else
    Gj = G_clearing(τk_joint)
    a  = period0_at_G(τk_joint, Gj)
    println()
    dT1 = p.lump_sum_T[1] - T_ss; dTu1 = p.lump_sum_Tu[1] - Tu_ss
    @printf("tau_k   = %+.6f      (paper: -0.47)\n", τk_joint)
    @printf("G       = %.8f  = %.2f x G_ss,  ratio to y_1 = %.4f  (model has %.4f)\n",
            Gj, Gj/G_ss, Gj/a.y, G0/p.y[1])
    @printf("tau_n   = %.8f  (SS %.8f)   -> budget closes with NO labour-tax move\n",
            a.τ_n, e_init.τ_n)
    @printf("package = %.4f%% ann GDP  (\$%.0f bn)\n", pct(a.pkg), dol(a.pkg))
    println("\nresulting shares -- the two IMPOSED and the two that come out:")
    @printf("%-22s %10s %10s %10s   %s\n", "", "model %", "Table E2", "diff pp", "")
    @printf("%-22s %10.2f %10.2f %+10.2f   %s\n", "firms (capital)",
            a.share_firms, TABLE_E2.firms, a.share_firms - TABLE_E2.firms, "IMPOSED")
    @printf("%-22s %10.2f %10.2f %+10.2f   %s\n", "other (G)",
            a.share_other, TABLE_E2.other, a.share_other - TABLE_E2.other, "output")
    @printf("%-22s %10.2f %10.2f %+10.2f   %s\n", "stimulus checks",
            100*dT1/a.pkg, TABLE_E2.checks, 100*dT1/a.pkg - TABLE_E2.checks, "output")
    @printf("%-22s %10.2f %10.2f %+10.2f   %s\n", "unemployment ins.",
            100*dTu1/a.pkg, TABLE_E2.ui, 100*dTu1/a.pkg - TABLE_E2.ui, "output")
    @printf("\nchecks + UI = %.2f%% of the package against Table E2's %.2f%%\n",
            100*(dT1+dTu1)/a.pkg, TABLE_E2.checks + TABLE_E2.ui)
    @printf("in levels: model %.4f%% ann GDP vs target %.4f%%  -> excess %.4f%% = \$%.0f bn\n",
            pct(dT1+dTu1), (TABLE_E2.checks + TABLE_E2.ui)/100 * pct(a.pkg),
            pct(dT1+dTu1) - (TABLE_E2.checks + TABLE_E2.ui)/100 * pct(a.pkg),
            dol(dT1+dTu1) - (TABLE_E2.checks + TABLE_E2.ui)/100 * dol(a.pkg))
    println("\nT and Tu are EXOGENOUS in T_path/Tu_path -- they are not free to absorb this,")
    println("so any residual error lands on whichever component is left endogenous.")
end
println("\nDONE (nothing written)")
