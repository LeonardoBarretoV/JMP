# ============================================================
# RE-RUN lower_foreign AND cbo_debt OFF THE RULE-2 BASELINE, [SW19] ON
# ============================================================
#     julia --project -t auto scripts/main_CF_rule2_rerun_sw.jl lower_foreign
#     julia --project -t auto scripts/main_CF_rule2_rerun_sw.jl cbo_debt
#
# Both were solved on the PRE-rule-2 baseline (2026-09-08).  Re-running them off
# rule 2 puts all four columns of Table 6 on ONE fiscal baseline, which is the
# precondition for re-reporting.
#
# BOTH MOVE THE TERMINAL STEADY STATE, so each needs its own e_final AND its own
# i_taylor_ss rebuilt from that new r_final.  That is the INTERCEPT CHANNEL and it
# is the whole reason the CBO column differs from the paper: i_ss enters the Taylor
# rule additively, the backward recursion turns a shift in it into
#     d(pi) = -d(i_ss)/(phi - 1 - r_f)
# so a higher long-run real rate mechanically lowers inflation along the whole
# transition.  Rebuilding i_ss is NOT optional bookkeeping -- it IS the mechanism.
#
#   lower_foreign : Lambda_bstar_f 1.148 -> 0.72, foreign holdings -10pp of annual
#                   GDP over j = 2..20 (in ratio-to-b units, -0.1/(qb_path[j]/4))
#   cbo_debt      : terminal debt 1.84*4 = 7.36, linear on qb_path from index 19's
#                   value to index 143 (2055Q4), foreign held at 0.287*4 of annual
#                   GDP from index 20 on (foreign_mode = :as_run, reproducing the
#                   original run's q*b/y divisor)
#
# Rule-2 fiscal settings are applied in BOTH: T[1] = 0.32446192,
# Tu[1] = 0.31282997, tau_k[1] = -0.47520201.  epsm is HELD at the rule-2 path.
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

WHICH = isempty(ARGS) ? "lower_foreign" : ARGS[1]
WHICH in ("lower_foreign","cbo_debt") || error("ARGS[1] must be lower_foreign or cbo_debt")

KAPPA_W    = 0.05
PHI_PI     = 1.5
NKPC_SCALE = 240.0
TK_RULE2   = -0.47520201
T1_RULE2   = 0.32446192
TU1_RULE2  = 0.31282997
SAVE_TAG   = "cf_" * WHICH * "_rule2"
SS_TAG     = "steady_states_" * SAVE_TAG
PAPER = WHICH == "lower_foreign" ? (mean=1.22, poor=2.03, middle=1.34, rich=0.06) :
                                   (mean=2.51, poor=4.00, middle=2.62, rich=0.50)

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

println("=== RE-RUN: ", WHICH, " off the rule-2 baseline, [SW19] ON ===")
_inp = load_object(joinpath(ROOT,"output","jld2_sw","steady_states.jld2"))
e_init, e_final_base = _inp[1], _inp[2]
h = e_init.h; t = e_init.t; g = e_init.g
δ = t.δ; π_bar = e_init.π_
i_ss_base = (1 + e_final_base.r)*(1 + π_bar) - 1

# ---- the counterfactual FINAL steady state, and its OWN i_taylor_ss ----------
ss_file = joinpath(ROOT,"output","jld2_sw",SS_TAG*".jld2")
Λ_bstar_f = WHICH == "lower_foreign" ? 0.18*4.0 : 0.287*4.0
Λ_b_f     = WHICH == "lower_foreign" ? c.Λ_b + 0.15*4.0 : 1.84*4.0
if isfile(ss_file)
    println("loading cached counterfactual final steady state ..."); flush(stdout)
    e_final = load_object(ss_file)
else
    ī_ = (1 + e_init.r)*(1 + c.π_bar) - 1
    rr = WHICH == "lower_foreign" ? (0.0, 0.005) : (0.0, 0.008)
    @printf("solving the counterfactual FINAL steady state: Lambda_b_f = %.4f, Lambda_bstar_f = %.4f\n",
            Λ_b_f, Λ_bstar_f); flush(stdout)
    @time e_final = stationary_equilibrium_final(h, t, g, e_init.τ_k, c.Λ_t, c.Λ_tu, c.Λ_g,
                                                 Λ_b_f, Λ_bstar_f, ī_;
                                                 r_range = rr, root_max = 1.0, verbose = true)
    save_object(ss_file, e_final)
    println("wrote ", ss_file)
end
i_ss = (1 + e_final.r)*(1 + π_bar) - 1
println("\n---- the INTERCEPT CHANNEL: i_taylor_ss is rebuilt from the NEW r_final ----")
@printf("baseline r_f = %.4f%% ann   i_ss = %.10f (%.4f%% ann)\n",
        100*((1+e_final_base.r)^4-1), i_ss_base, 100*((1+i_ss_base)^4-1))
@printf("THIS RUN r_f = %.4f%% ann   i_ss = %.10f (%.4f%% ann)\n",
        100*((1+e_final.r)^4-1), i_ss, 100*((1+i_ss)^4-1))
ic = intercept_channel_sw(i_ss, i_ss_base, PHI_PI, e_final.r)
@printf("q*b/y(e_final) = %.6f   bstar/y(e_final) = %.6f   bstar/b = %.8f\n",
        e_final.q*e_final.b/y(e_final), e_final.bstar/y(e_final), e_final.bstar/e_final.b)

# ---- time series, built with THIS run's e_final, then rule-2 fiscal settings --
b_target = e_init.b / y(e_init)
b_path = let
    T=700; H=100; ρB=0.08; bl = Array{Float64,1}(undef, T+H)
    bl[1]=b_target; bl[2]=b_target; bl[T:end] .= b_target
    for i in 3:T-1; bl[i] = bl[1]*ρB^(i-1) + (1-ρB^(i-1))*b_target; end
    bl
end
ts = make_time_series_sw(c, e_init, e_final, b_path, Λ_bstar_f, g)
T_r2  = copy(ts.T_path);   T_r2[1]  = T1_RULE2
Tu_r2 = copy(ts.Tu_path);  Tu_r2[1] = TU1_RULE2
τk_r2 = copy(ts.τ_k_path); τk_r2[1] = TK_RULE2
@assert all(τk_r2[j] == e_init.τ_k for j in 2:length(τk_r2)) "tau_k not flat for j>=2"

if WHICH == "lower_foreign"
    qb_use  = ts.qb_path
    qbs_use = cf_qbstar_lower_foreign_sw(ts; drop = 0.10, jlast = 20)
    @printf("\nqbstar: j=2 %.8f  j=20 %.8f  j>=21 %.8f  (expected 0.19326028 / 0.15321854 / %.8f)\n",
            qbs_use[2], qbs_use[20], qbs_use[21], e_final.bstar/e_final.b)
else
    qb_use, qbs_use = cf_paths_cbo_debt_sw(ts, e_final; start_index = 19, end_index = 143,
        target = 1.84*4.0, foreign_to_gdp = 0.287*4.0, foreign_mode = :as_run)
    @printf("\nqb: j=19 %.8f  j=20 %.8f  j=143 %.8f  j=200 %.8f\n",
            qb_use[19], qb_use[20], qb_use[143], qb_use[200])
    @printf("qbstar: j=20 %.8f  j=200 %.8f  (expected 0.23466575 / 0.15597826)\n",
            qbs_use[20], qbs_use[200])
end

pars_h = (θ=c.θ, ν=c.ν, ies=c.ies, crra=c.crra, β=c.β, a_min=c.a_min, a_max=c.a_max,
          z_u=0.0, n_a=c.n_a, init_state=e_init.h.Pss, z_scale=z_scale)
h_path, _, _ = h_path_transition_mckay_init(ROOT, ts.Pee_path, ts.Puu_path; pars_h...)

# ---- epsm HELD at the rule-2 path; warm start from the rule-2 baseline -------
println("\nloading rule2_stage1 for epsm ..."); flush(stdout)
_o = load_object(joinpath(ROOT,"output","jld2_sw","rule2_stage1.jld2"))
epsm = epsm_baseline_sw(copy(_o[3].r), copy(_o[3].π_), i_ss_base, π_bar, PHI_PI)
_o = nothing; GC.gc()
@printf("epsm held at the rule-2 path: max|epsm| = %.8f, j=9 = %+.4f pp ann\n",
        maximum(abs, epsm), 400*epsm[9])
println("loading rule2_epsm for the BASELINE COMPARISON ..."); flush(stdout)
_o = load_object(joinpath(ROOT,"output","jld2_sw","rule2_epsm.jld2"))
pB = _o[3]; π_base = copy(_o[4]); r_base = copy(pB.r); τn_base = copy(pB.τ_n)
r0, s0 = warm_start_sw(pB, e_init, t; method = 1, sw19 = true)
_o = nothing; pB = nothing; GC.gc()

# WARM START SOURCE.  cbo_debt moves the TERMINAL steady state hard (r_f 1.0778% ->
# 1.4211%), so the rule-2 baseline's r path is wrong by ~0.34pp annualised across all
# 800 periods -- that put the warm-start nkpc residual at 6.65 and Anderson took a
# 2.85e+08 step on iteration 2 and blew up.  The PRE-rule-2 CBO solution already has
# the correct terminal r path and identical qb/qbstar paths; it differs from the
# target only by rule 2's period-0 fiscal change.  This is a starting-point choice,
# not a model or solver change: the solution set is untouched.
WARM_SRC = WHICH == "cbo_debt" ? "cf_cbo_debt" : "rule2_epsm"
if WARM_SRC != "rule2_epsm"
    println("warm-starting from ", WARM_SRC, " (its terminal r path is already right) ..."); flush(stdout)
    _o = load_object(joinpath(ROOT,"output","jld2_sw",WARM_SRC*".jld2"))
    global r0, s0 = warm_start_sw(_o[3], e_init, t; method = 1, sw19 = true)
    _o = nothing; GC.gc()
end

n_per = length(τk_r2); n_r = _n_r_sw(n_per, 1, true); n_s = n_per - 2
DY_ss = union_aggregates_(e_init.h, g, e_init.ws; s=1.0, w=e_init.w, τ_n=e_init.τ_n,
                          n_ss=e_init.ws.n_pol[1,:])
wY_ss = e_init.w*DY_ss.Y
function warm_resid(x)
    path_t = _initialize_path(e_init, e_final; qb_path = qb_use)
    πv = Vector{Float64}(undef, n_per); F = similar(x)
    _residuals_sw_c!(F, x; π_vec=πv, epsm_path=epsm, i_taylor_ss=i_ss, phi_pi=PHI_PI,
        path=path_t, e_init=e_init, e_final=e_final, τ_k_path=τk_r2, T_path=T_r2,
        Tu_path=Tu_r2, G_path=ts.G_path, qb_path=qb_use, qbstar_path=qbs_use,
        h_path=h_path, π_path=ts.π_path, d_path=ts.d_path, covid_path=ts.covid_path,
        Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g, method=1, kappa_w=KAPPA_W,
        wY_ss=wY_ss, nkpc_scale=NKPC_SCALE, s_mode=:nkpc, sw19=true)
    return (asset=maximum(abs, view(F,1:n_r)), nkpc=maximum(abs, view(F,(n_r+1):(n_r+n_s))),
            inf=maximum(abs, F))
end
w0 = warm_resid(vcat(r0, s0))
@printf("\nwarm start from rule2_epsm: asset %.6e   nkpc %.6e   inf-norm %.6e\n",
        w0.asset, w0.nkpc, w0.inf); flush(stdout)

println("\n", "="^100); println("SOLVE -- :nkpc, [SW19] ON"); println("="^100)
path_c, sol_c, π_c, info_c = solve_transition_sw_c(e_init, e_final;
    τ_k_path=τk_r2, T_path=T_r2, Tu_path=Tu_r2, G_path=ts.G_path,
    qb_path=qb_use, qbstar_path=qbs_use, h_path=h_path, π_path=ts.π_path,
    d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
    method=1, kappa_w=KAPPA_W, epsm_path=epsm, phi_pi=PHI_PI, i_taylor_ss=i_ss,
    nkpc_scale=NKPC_SCALE, s_mode=:nkpc, r_path=r0, s_path=s0, sw19=true)

println("\n", "="^100); println("REPORT -- ", WHICH, " off rule 2"); println("="^100)
@printf("iterations %d   converged %s   wall %.2f min\n",
        sol_c.iterations, sol_c.f_converged, info_c.wall/60)
guard_report_sw()
_SW8_WIDEN[] = 0; _SW8_FALLBACK[] = 0
warm_resid(sol_c.zero); fb_sol = _SW8_FALLBACK[]
@printf("[SW17] at the SOLUTION: [SW8] fallback %d\n", fb_sol)
r_min = -δ + eps(); r_max = 1/h_path[1].u.β
@printf("[SW19] residuals[1] = %.6e   r_k_hat = %.12f (clamp %s)\n",
        e_init.k - path_c.k[1], path_c.r_k[1]-δ,
        (r_min <= path_c.r_k[1]-δ <= r_max) ? "OK" : "VIOLATED")
sol_c.f_converged || (println("\n*** NOT CONVERGED -- nothing written ***"); exit(0))

gate = gate_lab_tax_sw(e_init, path_c, g); gate_ok = print_gate_sw(gate)
s_flex_1 = s_flex_sw(g, c.ν; w = path_c.w[1], τ_n = path_c.τ_n[1],
                     w_init = e_init.w, τ_n_init = e_init.τ_n)
branch_ok = ts.Γ_n <= s_flex_1
@printf("\nGamma_n = %.10f   s_flex = %.10f  -> %s\n", ts.Γ_n, s_flex_1,
        branch_ok ? "Gamma_n (CORRECT)" : "s_flex (*** INVALID ***)")
@printf("\ntau_n j=1 = %.8f (1-tau_n = %+.8f)   j=2 = %.8f\n",
        path_c.τ_n[1], 1-path_c.τ_n[1], path_c.τ_n[2])

println("\n---- pi and r vs the rule-2 baseline, and the intercept channel ----")
@printf("%6s %12s %12s %14s %14s %14s\n","j","pi base","pi cf","solved d(pi)","intercept","residual")
for j in (20,40,100,143)
    solved = 400*(π_c[j]-π_base[j]); icp = 400*ic.d_pi
    @printf("%6d %12.4f %12.4f %+14.4f %+14.4f %+14.4f\n", j,
            100*((1+π_base[j])^4-1), 100*((1+π_c[j])^4-1), solved, icp, solved-icp)
end
@printf("terminal pi = %.6f%% ann (pi_bar = %.6f%%)\n",
        100*((1+π_c[end])^4-1), 100*((1+π_bar)^4-1))

ce = ce_summary_sw(e_init, path_c)
ceB = (mean = 1.910, poor = 3.350, middle = 2.043, rich = -0.325)
pB_paper = (mean = 1.82, poor = 3.03, middle = 2.00, rich = -0.20)
println("\n---- WELFARE (paper bin convention). LEVELS are NOT like-for-like with the")
println("     published numbers (rule 2's Other correction is absent there); the")
println("     DIFFERENCE from the rule-2 baseline IS, because the correction cancels. ----")
@printf("%-8s %12s %12s %14s   %12s %14s\n",
        "group","rule-2 base","this cf","DIFFERENCE","paper level","paper diff")
for (nm, bv, cv, pl, pbv) in (("mean", ceB.mean, ce.mean, PAPER.mean, pB_paper.mean),
                              ("poor", ceB.poor, ce.poor, PAPER.poor, pB_paper.poor),
                              ("middle", ceB.middle, ce.middle, PAPER.middle, pB_paper.middle),
                              ("rich", ceB.rich, ce.rich, PAPER.rich, pB_paper.rich))
    @printf("%-8s %12.3f %12.3f %+14.3f   %12.2f %+14.2f\n", nm, bv, cv, cv-bv, pl, pl-pbv)
end
@printf("%-8s %12.3f\n", "min", ce.min)
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
println("\n=== ", WHICH, " done ===")
