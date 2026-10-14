# ============================================================
# STAGE 5 COUNTERFACTUAL -- ONLY STIMULUS CHECKS, off the RULE-2 baseline, [SW19] ON
# ============================================================
#     julia --project -t auto scripts/main_CF_only_checks_rule2_sw.jl
#
# Section 9.1: all additional resources from the debt expansion go exclusively to
# eligible households as stimulus checks, with the DEBT PATH UNCHANGED.
#
#   tau_k[1]  -> e_init.tau_k   (PPP removed)
#   Tu_path[1]-> its STEADY-STATE ratio  (the UI expansion removed IN FULL)
#   T_path[1] -> absorbs the freed UI and PPP resources, PLUS all of rule 2's
#                "Other" allocation, so the Tu instrument returns fully to steady
#                state and the experiment matches its name
#   qb_path, G_path, bstar_path_r, lockdown, and tau_k for j >= 2: UNCHANGED
#
# WHY RULE 2's "OTHER" STAYS IN TOTAL BUT MOVES TO CHECKS
#   Rule 2's Other allocation is a CORRECTION to the baseline's budget coherence --
#   the model was borrowing for ~$400bn of spending it never did -- not a policy
#   being counterfactualled.  Dissolving it here would compare a corrected baseline
#   against an uncorrected counterfactual and would make this column incomparable
#   with the other three.  So the $400bn total is preserved; only its routing
#   changes, entirely into checks.  Under the alternative of keeping rule 2's
#   Other-to-UI slice in UI, this run would still deliver ~$172bn through Tu while
#   being labelled "only stimulus checks".
#
# ############################################################################
# # DOCUMENT THIS WHEN WRITING THE SECTION:                                  #
# # This counterfactual differs from the paper's version in TWO ways, not    #
# # one: (a) the composition change that is the experiment, and (b) rule 2's #
# # "Other" routing, which the published baseline did not have.  The         #
# # published 2.02 / 3.77 / 2.27 / -1.69 came from a baseline WITHOUT the    #
# # Other correction, so the LEVELS are NOT like-for-like and must not be    #
# # presented as though they are.  The counterfactual DIFFERENCE (cf minus   #
# # rule-2 baseline) IS like-for-like, because the Other correction is       #
# # common to both and cancels -- that is the comparison to report.          #
# ############################################################################
#
# SIZING CHECK FIRST.  Only the SPLIT of a fixed total moves, so the period-0 GBC
# must close with tau_n at (or very near) its steady-state value.  The script
# computes the reallocation from the resource amounts, then CHECKS what tau_n comes
# out -- it does not solve for T to force tau_n = SS, which would assume the answer.
# If tau_n drifts more than TOL it STOPS without solving.
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
const _TOL = 1e-12
const _ZERO_FTOL = 1e-5
const _MAX_ITERS = 20_000
const _ZERO_MAX_ITERS = 450
_φ = 1/Base.MathConstants.φ
ProgressMeter.ijulia_behavior(:clear);

KAPPA_W    = 0.05
PHI_PI     = 1.5
NKPC_SCALE = 240.0
TK_RULE2   = -0.47520201
T1_RULE2   = 0.32446192
TU1_RULE2  = 0.31282997
GDP_2019   = 21380.0
SAVE_TAG   = "cf_only_checks_rule2"
TOL_TAUN   = 0.02            # stop if |tau_n/tau_n_ss - 1| exceeds this
PAPER_CE   = (mean = 2.02, poor = 3.77, middle = 2.27, rich = -1.69)   # NOT like-for-like, see above

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

println("loading steady states ..."); flush(stdout)
_inp = load_object(joinpath(ROOT,"output","jld2_sw","steady_states.jld2"))
e_init, e_final = _inp[1], _inp[2]
h = e_init.h; t = e_init.t; g = e_init.g
δ = t.δ; π_bar = e_init.π_
i_ss = (1 + e_final.r)*(1 + π_bar) - 1
y_ss = y(e_init); Y_ANN = 4*y_ss
pct(x) = 100*x/Y_ANN
dol(x) = GDP_2019*x/Y_ANN
caprev_ss = e_init.τ_k*(e_init.r + δ)*e_init.k

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
h0 = h_path[2]; z_dist0 = z_state[1,:]; n_ss = e_init.ws.n_pol[1,:]
r_min = -δ + eps(); r_max = 1/h_path[1].u.β
n0 = labor_supply_union_(h0; s = ts.Γ_n, n_ss = n_ss, z_dist = z_dist0)

# reference runs
println("loading rule2_epsm (baseline) ..."); flush(stdout)
_o = load_object(joinpath(ROOT,"output","jld2_sw","rule2_epsm.jld2"))
pB = _o[3]; π_base = copy(_o[4]); r_base = copy(pB.r); τn_base = copy(pB.τ_n)
y1_B = pB.y[1]; rb_B = pB.r[1]; yp_B = pB.y[2]
r0, s0 = warm_start_sw(pB, e_init, t; method = 1, sw19 = true)
_o = nothing; GC.gc()
println("loading cf_no_capital_subsidy_rule2 (same tau_k[1], for r_b) ..."); flush(stdout)
_o = load_object(joinpath(ROOT,"output","jld2_sw","cf_no_capital_subsidy_rule2.jld2"))
pN = _o[3]; rb_N = pN.r[1]; yp_N = pN.y[2]
rN0, sN0 = warm_start_sw(pN, e_init, t; method = 1, sw19 = true)
_o = nothing; GC.gc()
println("loading rule2_stage1 for epsm ..."); flush(stdout)
_o = load_object(joinpath(ROOT,"output","jld2_sw","rule2_stage1.jld2"))
epsm = epsm_baseline_sw(copy(_o[3].r), copy(_o[3].π_), i_ss, π_bar, PHI_PI)
_o = nothing; GC.gc()

# ============================================================
# SIZING: compute the reallocation from the resource amounts
# ============================================================
Tu_ss_ratio = ts.Tu_path[2]          # make_time_series_sw fills j>=2 with the SS ratio
τk_cf = e_init.τ_k
rk_r2 = _rk_star_sw19(t, e_init; τ_k = TK_RULE2, n0, d = ts.d_path[1],
                      covid = ts.covid_path[1], r_min, r_max)
rk_cf = _rk_star_sw19(t, e_init; τ_k = τk_cf, n0, d = ts.d_path[1],
                      covid = ts.covid_path[1], r_min, r_max)
# [SW19] pins k[1] = e_init.k and n is fixed, so y[1] is the SAME in every one of
# these runs -- the level/ratio conversion below uses one number for all of them.
y1 = t.A * (ts.d_path[1]*e_init.k)^t.α * n0^(1 - t.α)
@printf("\ny[1] under [SW19] = %.10f   (rule-2 baseline path had %.10f, diff %.2e)\n",
        y1, y1_B, y1 - y1_B)

caprev_r2 = TK_RULE2 * rk_r2 * e_init.k     # negative = subsidy outlay
caprev_cf = τk_cf    * rk_cf * e_init.k     # positive = revenue
freed_PPP = caprev_cf - caprev_r2           # resources released by removing the subsidy
freed_UI  = (TU1_RULE2 - Tu_ss_ratio) * y1  # resources released by removing the UI expansion
T1_level_r2 = T1_RULE2 * y1
T1_level_cf = T1_level_r2 + freed_UI + freed_PPP
T1_cf = T1_level_cf / y1

println("\n", "="^100)
println("SIZING -- the reallocation, computed from the resource amounts")
println("="^100)
@printf("%-42s %14s %12s %12s\n", "", "model units", "% ann GDP", "\$bn")
@printf("%-42s %14.8f %12.4f %12.1f\n", "rule-2 checks T[1]",          T1_level_r2, pct(T1_level_r2), dol(T1_level_r2))
@printf("%-42s %14.8f %12.4f %12.1f\n", "  + freed UI expansion",       freed_UI,    pct(freed_UI),    dol(freed_UI))
@printf("%-42s %14.8f %12.4f %12.1f\n", "  + freed PPP subsidy",        freed_PPP,   pct(freed_PPP),   dol(freed_PPP))
@printf("%-42s %14.8f %12.4f %12.1f\n", "  = checks T[1] this run",     T1_level_cf, pct(T1_level_cf), dol(T1_level_cf))
@printf("\nT_path[1] ratio : %.8f -> %.8f\n", T1_RULE2, T1_cf)
@printf("Tu_path[1] ratio: %.8f -> %.8f  (steady-state ratio)\n", TU1_RULE2, Tu_ss_ratio)
@printf("tau_k[1]        : %+.8f -> %+.8f\n", TK_RULE2, τk_cf)
@printf("capital tax     : %.8f (subsidy) -> %+.8f (revenue)\n", caprev_r2, caprev_cf)

# period-0 GBC check: what tau_n does this imply?
function taun_at(; τ_k, T1, Tu1, r_k, r_b, yprime)
    pr = ProgressUnknown()
    f = (τ_n) -> _τ_n_excess_t0_port_sw_!(h0, t, g; τ_n, τ_k, T = T1, Tu = Tu1,
        G = ts.G_path[1], qbprime = ts.qb_path[2], qb = ts.qb_path[1], yprime = yprime,
        d = ts.d_path[1], s = ts.Γ_n, n_ss, covid = ts.covid_path[1], z_dist = z_dist0,
        r_b = r_b, r_k = r_k, verbose = false, p = pr).excess_g
    return _solve_tau_n_sw(f, 0.0, 8.0)
end
τn_chk_base = taun_at(; τ_k = TK_RULE2, T1 = T1_RULE2, Tu1 = TU1_RULE2, r_k = rk_r2,
                      r_b = rb_B, yprime = yp_B)
τn_chk_cf   = taun_at(; τ_k = τk_cf, T1 = T1_cf, Tu1 = Tu_ss_ratio, r_k = rk_cf,
                      r_b = rb_N, yprime = yp_N)
println("\n---- period-0 GBC check (r_b and y' taken from the matching solved run) ----")
@printf("rule-2 baseline reproduces tau_n = %.8f  (saved %.8f, diff %.2e)\n",
        τn_chk_base, τn_base[1], τn_chk_base - τn_base[1])
@printf("this run implies      tau_n = %.8f   vs SS %.8f   -> %+.4f%%\n",
        τn_chk_cf, e_init.τ_n, 100*(τn_chk_cf/e_init.τ_n - 1))
drift = abs(τn_chk_cf/e_init.τ_n - 1)
if drift > TOL_TAUN
    @printf("\n*** tau_n DRIFTS %.4f%% from steady state, above the %.2f%% tolerance.\n",
            100*drift, 100*TOL_TAUN)
    println("*** The reallocation is not sized right -- only the SPLIT of a fixed total")
    println("*** should be moving. STOPPING without solving. Nothing written.")
    exit(0)
end
println("SIZING OK: only the split moved, and the period-0 budget still closes at ~SS.")
flush(stdout)

# ============================================================
# SOLVE
# ============================================================
T_cf  = copy(ts.T_path);   T_cf[1]  = T1_cf
Tu_cf = copy(ts.Tu_path);  Tu_cf[1] = Tu_ss_ratio
τk_cf_path = copy(ts.τ_k_path); τk_cf_path[1] = τk_cf
@assert τk_cf_path == fill(e_init.τ_k, length(τk_cf_path)) "tau_k not flat at e_init.tau_k"
@assert all(T_cf[j]  == ts.T_path[j]  for j in 2:length(T_cf))
@assert all(Tu_cf[j] == ts.Tu_path[j] for j in 2:length(Tu_cf))
@printf("\ntau_k flat at e_init.tau_k for every j -> e_final unchanged, i_taylor_ss = %.4f%% ann\n",
        100*((1+i_ss)^4-1))

n_per = length(τk_cf_path); n_r = _n_r_sw(n_per, 1, true); n_s = n_per - 2
DY_ss = union_aggregates_(e_init.h, g, e_init.ws; s=1.0, w=e_init.w, τ_n=e_init.τ_n, n_ss=n_ss)
wY_ss = e_init.w*DY_ss.Y
function warm_resid(x, s_mode)
    path_t = _initialize_path(e_init, e_final; qb_path = ts.qb_path)
    πv = Vector{Float64}(undef, n_per); F = similar(x)
    _residuals_sw_c!(F, x; π_vec=πv, epsm_path=epsm, i_taylor_ss=i_ss, phi_pi=PHI_PI,
        path=path_t, e_init=e_init, e_final=e_final, τ_k_path=τk_cf_path, T_path=T_cf,
        Tu_path=Tu_cf, G_path=ts.G_path, qb_path=ts.qb_path, qbstar_path=ts.bstar_path_r,
        h_path=h_path, π_path=ts.π_path, d_path=ts.d_path, covid_path=ts.covid_path,
        Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g, method=1, kappa_w=KAPPA_W,
        wY_ss=wY_ss, nkpc_scale=NKPC_SCALE, s_mode=s_mode, sw19=true)
    return (asset=maximum(abs, view(F,1:n_r)), nkpc=maximum(abs, view(F,(n_r+1):(n_r+n_s))),
            inf=maximum(abs, F))
end
wB = warm_resid(vcat(r0, s0), :nkpc)
wN = warm_resid(vcat(rN0, sN0), :nkpc)
@printf("\nwarm start from rule2_epsm            : asset %.6e  nkpc %.6e  inf %.6e\n", wB.asset, wB.nkpc, wB.inf)
@printf("warm start from cf_no_capital_subsidy : asset %.6e  nkpc %.6e  inf %.6e\n", wN.asset, wN.nkpc, wN.inf)
rS, sS = wN.inf < wB.inf ? (rN0, sN0) : (r0, s0)
println("using the closer of the two."); flush(stdout)

println("\n", "="^100); println("SOLVE -- :nkpc, [SW19] ON"); println("="^100)
path_c, sol_c, π_c, info_c = solve_transition_sw_c(e_init, e_final;
    τ_k_path=τk_cf_path, T_path=T_cf, Tu_path=Tu_cf, G_path=ts.G_path,
    qb_path=ts.qb_path, qbstar_path=ts.bstar_path_r, h_path=h_path, π_path=ts.π_path,
    d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
    method=1, kappa_w=KAPPA_W, epsm_path=epsm, phi_pi=PHI_PI, i_taylor_ss=i_ss,
    nkpc_scale=NKPC_SCALE, s_mode=:nkpc, r_path=rS, s_path=sS, sw19=true)

println("\n", "="^100); println("1. CONVERGENCE"); println("="^100)
@printf("iterations %d   converged %s   wall %.2f min\n",
        sol_c.iterations, sol_c.f_converged, info_c.wall/60)
guard_report_sw()
_SW8_WIDEN[] = 0; _SW8_FALLBACK[] = 0
warm_resid(sol_c.zero, :nkpc)
fb_sol = _SW8_FALLBACK[]
@printf("[SW17] at the SOLUTION: [SW8] fallback %d\n", fb_sol)
@printf("\n[SW19] residuals[1] = e_init.k - path.k[1] = %.6e   r_k_hat = %.12f (clamp %s)\n",
        e_init.k - path_c.k[1], path_c.r_k[1] - δ,
        (r_min <= path_c.r_k[1]-δ <= r_max) ? "OK" : "VIOLATED")
sol_c.f_converged || (println("\n*** NOT CONVERGED -- nothing written ***"); exit(0))

println("\n", "="^100); println("2. ITEM (b) GATE and the PERIOD-0 BRANCH"); println("="^100)
gate = gate_lab_tax_sw(e_init, path_c, g); gate_ok = print_gate_sw(gate)
s_flex_1 = s_flex_sw(g, c.ν; w = path_c.w[1], τ_n = path_c.τ_n[1],
                     w_init = e_init.w, τ_n_init = e_init.τ_n)
branch_ok = ts.Γ_n <= s_flex_1
@printf("\nGamma_n = %.10f   s_flex = %.10f   -> %s\n", ts.Γ_n, s_flex_1,
        branch_ok ? "Gamma_n (CORRECT)" : "s_flex (*** INVALID ***)")

println("\n", "="^100); println("3. tau_n AND 1 - tau_n"); println("="^100)
@printf("%-24s %14s %14s %14s\n", "", "tau_n", "1 - tau_n", "vs SS (%)")
@printf("%-24s %14.8f %14.8f %14s\n", "initial SS", e_init.τ_n, 1-e_init.τ_n, "-")
for j in (1,2)
    @printf("%-24s %14.8f %14.8f %+14.2f\n", "rule-2 baseline j=$j", τn_base[j], 1-τn_base[j],
            100*(τn_base[j]/e_init.τ_n - 1))
    @printf("%-24s %14.8f %14.8f %+14.2f\n", "only checks     j=$j", path_c.τ_n[j], 1-path_c.τ_n[j],
            100*(path_c.τ_n[j]/e_init.τ_n - 1))
end

println("\n", "="^100); println("pi AND r vs THE RULE-2 BASELINE (annualised %)"); println("="^100)
@printf("%5s %12s %12s %12s   %12s %12s %12s\n","j","pi base","pi cf","diff pp","r base","r cf","diff pp")
for j in vcat(collect(1:8),[12,20,40,100,400])
    @printf("%5d %12.4f %12.4f %+12.4f   %12.4f %12.4f %+12.4f\n", j,
            100*((1+π_base[j])^4-1), 100*((1+π_c[j])^4-1), 400*(π_c[j]-π_base[j]),
            100*((1+r_base[j])^4-1), 100*((1+path_c.r[j])^4-1), 400*(path_c.r[j]-r_base[j]))
end
@printf("\nmax|pi_cf - pi_base| j=2..60 = %.4f ann pp   terminal pi = %.6f%% ann\n",
        400*maximum(abs.(π_c[2:60] .- π_base[2:60])), 100*((1+π_c[end])^4-1))

println("\n", "="^100); println("5. WELFARE (CE, paper bin convention)"); println("="^100)
ce = ce_summary_sw(e_init, path_c)
ceB = (mean = 1.910, poor = 3.350, middle = 2.043, rich = -0.325)   # rule-2 baseline
println("!"^100)
println("NOT LIKE-FOR-LIKE WITH THE PUBLISHED NUMBERS. This run differs from the paper's")
println("version in TWO ways: the composition change that IS the experiment, and rule 2's")
println("'Other' routing, which the published baseline did not have. The published")
println("2.02 / 3.77 / 2.27 / -1.69 came from a baseline WITHOUT the Other correction.")
println("Report the DIFFERENCE (cf minus rule-2 baseline): the correction is common to")
println("both and cancels, so the difference IS like-for-like. The levels are not.")
println("!"^100)
@printf("\n%-8s %12s %12s %14s   %12s %12s\n",
        "group","rule-2 base","only checks","DIFFERENCE","paper level","paper diff")
pdiff = (mean = 2.02-1.82, poor = 3.77-3.03, middle = 2.27-2.00, rich = -1.69-(-0.20))
for (nm, bv, cv, pv, pd) in (("mean", ceB.mean, ce.mean, PAPER_CE.mean, pdiff.mean),
                             ("poor", ceB.poor, ce.poor, PAPER_CE.poor, pdiff.poor),
                             ("middle", ceB.middle, ce.middle, PAPER_CE.middle, pdiff.middle),
                             ("rich", ceB.rich, ce.rich, PAPER_CE.rich, pdiff.rich))
    @printf("%-8s %12.3f %12.3f %+14.3f   %12.2f %+12.2f\n", nm, bv, cv, cv-bv, pv, pd)
end
@printf("%-8s %12s %12.3f\n", "min", "-2.436/-5.383", ce.min)
println("\nfull decile profile:")
for (k,v) in enumerate(ce.deciles); @printf("  bin %2d  %8.4f\n", k, v); end

ok = sol_c.f_converged && gate_ok && branch_ok && fb_sol == 0 &&
     abs(e_init.k - path_c.k[1]) < 1e-12
println("\n", "="^100)
@printf("converged %s  gate %s  branch %s  [SW8]@sol %d  |res[1]| %.3e\n",
        sol_c.f_converged, gate_ok, branch_ok, fb_sol, abs(e_init.k - path_c.k[1]))
if ok
    save_object(joinpath(ROOT,"output","jld2_sw",SAVE_TAG*".jld2"), (e_init,e_final,path_c,π_c))
    println("wrote output/jld2_sw/", SAVE_TAG, ".jld2")
else
    println("NOT WRITTEN.")
end
println("\n=== done ===")
