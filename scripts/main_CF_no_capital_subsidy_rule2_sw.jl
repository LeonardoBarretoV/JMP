# ============================================================
# STAGE 5 COUNTERFACTUAL -- NO CAPITAL SUBSIDIES, off the RULE-2 baseline, [SW19] ON
# ============================================================
#     julia --project -t auto scripts/main_CF_no_capital_subsidy_rule2_sw.jl
#
# tau_k[1] = e_init.tau_k = +0.13793103 instead of the rule-2 −0.47520201.
# EVERYTHING else identical to rule2_epsm: T_path[1] = 0.32446192,
# Tu_path[1] = 0.31282997, qb_path, G_path, bstar_path_r, lockdown, and tau_k for
# all j >= 2.  epsm HELD at the rule-2 path (re-derived there, not the old one).
#
# [SW19] is ON: r_k and residuals[1] are eliminated, unknowns 1597.  Gates 1 and 2
# passed -- the elimination reproduces rule2_epsm (r and pi identical to 0.0,
# path.k[1] = e_init.k to 2.2e-16 against the old solve's -4.9e-10 miss) and the
# Jacobian block diagonals are UNCHANGED (asset/nkpc ratio 1.0184 either way), so
# nkpc_scale = 240 stands.
#
# WHY THIS SHOULD NOW SOLVE.  The previous attempt blew up in the ASSET block at
# 5.6e20.  residuals[1] was carrying essentially the whole warm-start asset residual
# (0.841 of it, against 7.4e-03 for rows 2..n), and [SW19] removes that row exactly
# rather than asking Anderson to solve it.
#
# WRITES output/jld2_sw/cf_no_capital_subsidy_rule2.jld2 only if converged AND the
# item (b) gate passes AND the period-0 branch is Gamma_n AND the [SW8] fallback
# fired 0 times AT THE SOLUTION.
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

# ---- knobs ----------------------------------------------------------------
KAPPA_W    = 0.05
PHI_PI     = 1.5
NKPC_SCALE = 240.0
TK_RULE2   = -0.47520201     # the rule-2 baseline value
T1_RULE2   = 0.32446192
TU1_RULE2  = 0.31282997
SAVE_TAG   = "cf_no_capital_subsidy_rule2"
TRY_FLEX_FIRST = false       # [SW19] should make the direct :nkpc solve viable
PAPER_CE   = (mean = 2.02, poor = 3.77, middle = 2.27, rich = -1.69)
_sw_runs_allowed[] = 2       # in case the :flex fallback is needed

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

# rule-2 fiscal paths, then the counterfactual tau_k
T_r2  = copy(ts.T_path);   T_r2[1]  = T1_RULE2
Tu_r2 = copy(ts.Tu_path);  Tu_r2[1] = TU1_RULE2
τk_base = copy(ts.τ_k_path); τk_base[1] = TK_RULE2
τk_cf   = copy(ts.τ_k_path); τk_cf[1]   = e_init.τ_k

println("\n---- the counterfactual, and the final-steady-state claim CHECKED ----")
@printf("tau_k[1]  rule-2 baseline = %+.8f   this run = %+.8f  (= e_init.tau_k)\n",
        τk_base[1], τk_cf[1])
@assert all(τk_cf[j] == e_init.τ_k for j in eachindex(τk_cf)) "tau_k is not flat at e_init.tau_k"
@assert τk_cf == fill(e_init.τ_k, length(τk_cf)) "tau_k_cf differs from a flat e_init.tau_k path"
println("tau_k_cf is FLAT at e_init.tau_k for EVERY j -- identical to the path")
println("stationary_equilibrium_final was called with, so e_final is the same object.")
@printf("i_taylor_ss = %.10f (%.4f%% ann)   r_f = %.4f%% ann   -- unchanged, NO intercept channel\n",
        i_ss, 100*((1+i_ss)^4-1), 100*((1+e_final.r)^4-1))
@assert T_r2[1] == T1_RULE2 && Tu_r2[1] == TU1_RULE2
@assert all(T_r2[j] == ts.T_path[j] for j in 2:length(T_r2))
@assert all(Tu_r2[j] == ts.Tu_path[j] for j in 2:length(Tu_r2))
println("T, Tu identical to the rule-2 baseline at every j; qb/G/bstar/lockdown untouched.")

# epsm from the rule-2 stage-1 path (the one rule2_epsm was built with)
println("\nloading rule2_stage1 for epsm ..."); flush(stdout)
_o = load_object(joinpath(ROOT,"output","jld2_sw","rule2_stage1.jld2"))
epsm = epsm_baseline_sw(copy(_o[3].r), copy(_o[3].π_), i_ss, π_bar, PHI_PI)
_o = nothing; GC.gc()
@printf("max|epsm| = %.8f   at 2022Q2 (j=9) = %+.4f pp ann\n", maximum(abs, epsm), 400*epsm[9])

println("\nloading rule2_epsm for the warm start and the baseline comparison ..."); flush(stdout)
_o = load_object(joinpath(ROOT,"output","jld2_sw","rule2_epsm.jld2"))
pB = _o[3]
r0, s0 = warm_start_sw(pB, e_init, t; method = 1, sw19 = true)
π_base = copy(_o[4]); r_base = copy(pB.r); τn_base = copy(pB.τ_n); rk_base = pB.r_k[1]
_o = nothing; GC.gc()

n_per = length(τk_cf); n_r = _n_r_sw(n_per, 1, true); n_s = n_per - 2
DY_ss = union_aggregates_(e_init.h, g, e_init.ws; s=1.0, w=e_init.w, τ_n=e_init.τ_n,
                          n_ss=e_init.ws.n_pol[1,:])
wY_ss = e_init.w*DY_ss.Y

# [SW19] r_k* at the counterfactual tau_k, and the clamp check
z_state = zeros(length(h_path), length(e_init.h.z_grid))
z_state[1,:] = transpose(e_init.h.Pss)*h_path[1].P
n0 = labor_supply_union_(h_path[2]; s = ts.Γ_n, n_ss = e_init.ws.n_pol[1,:], z_dist = z_state[1,:])
r_min = -δ + eps(); r_max = 1/h_path[1].u.β
rk_cf = _rk_star_sw19(t, e_init; τ_k = e_init.τ_k, n0 = n0, d = ts.d_path[1],
                      covid = ts.covid_path[1], r_min = r_min, r_max = r_max)
@printf("\n[SW19] r_k* at the counterfactual tau_k = %.12f   r_k_hat* = %.12f\n", rk_cf, rk_cf - δ)
@printf("       inside the [SW10] clamp [%.6f, %.6f]? %s\n",
        r_min, r_max, (r_min <= rk_cf - δ <= r_max) ? "YES" : "NO")
@printf("       rule-2 baseline r_k was %.12f  -> ratio %.6f\n", rk_base, rk_cf/rk_base)

function warm_resid(x, s_mode)
    path_t = _initialize_path(e_init, e_final; qb_path = ts.qb_path)
    πv = Vector{Float64}(undef, n_per); F = similar(x)
    _residuals_sw_c!(F, x; π_vec=πv, epsm_path=epsm, i_taylor_ss=i_ss, phi_pi=PHI_PI,
        path=path_t, e_init=e_init, e_final=e_final, τ_k_path=τk_cf, T_path=T_r2,
        Tu_path=Tu_r2, G_path=ts.G_path, qb_path=ts.qb_path, qbstar_path=ts.bstar_path_r,
        h_path=h_path, π_path=ts.π_path, d_path=ts.d_path, covid_path=ts.covid_path,
        Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g, method=1, kappa_w=KAPPA_W,
        wY_ss=wY_ss, nkpc_scale=NKPC_SCALE, s_mode=s_mode, sw19=true)
    return (asset=maximum(abs, view(F,1:n_r)), nkpc=maximum(abs, view(F,(n_r+1):(n_r+n_s))),
            inf=maximum(abs, F))
end
wn = warm_resid(vcat(r0, s0), :nkpc)
wf = warm_resid(vcat(r0, s0), :flex)
println("\n---- warm start off the rule-2 baseline, with [SW19] ----")
@printf(":nkpc   asset %.6e   nkpc   %.6e   inf-norm %.6e\n", wn.asset, wn.nkpc, wn.inf)
@printf(":flex   asset %.6e   s-blk  %.6e   inf-norm %.6e\n", wf.asset, wf.nkpc, wf.inf)
println("(the pre-[SW19] :flex asset residual at the counterfactual was 8.41e-01)")
flush(stdout)

# ============================================================
# SOLVE
# ============================================================
rS, sS = r0, s0
if TRY_FLEX_FIRST
    println("\n", "="^100); println("STAGE A -- s_mode = :flex"); println("="^100)
    pA, sA, πA, iA = solve_transition_sw_c(e_init, e_final;
        τ_k_path=τk_cf, T_path=T_r2, Tu_path=Tu_r2, G_path=ts.G_path,
        qb_path=ts.qb_path, qbstar_path=ts.bstar_path_r, h_path=h_path, π_path=ts.π_path,
        d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
        method=1, kappa_w=KAPPA_W, epsm_path=epsm, phi_pi=PHI_PI, i_taylor_ss=i_ss,
        nkpc_scale=NKPC_SCALE, s_mode=:flex, r_path=r0, s_path=s0, sw19=true)
    @printf("STAGE A: iterations %d  converged %s  wall %.2f min\n",
            sA.iterations, sA.f_converged, iA.wall/60)
    sA.f_converged || (println("*** STAGE A failed; stopping ***"); exit(0))
    global rS, sS = warm_start_sw(pA, e_init, t; method = 1, sw19 = true)
end

println("\n", "="^100); println("SOLVE -- s_mode = :nkpc, [SW19] ON"); println("="^100)
path_c, sol_c, π_c, info_c = solve_transition_sw_c(e_init, e_final;
    τ_k_path=τk_cf, T_path=T_r2, Tu_path=Tu_r2, G_path=ts.G_path,
    qb_path=ts.qb_path, qbstar_path=ts.bstar_path_r, h_path=h_path, π_path=ts.π_path,
    d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
    method=1, kappa_w=KAPPA_W, epsm_path=epsm, phi_pi=PHI_PI, i_taylor_ss=i_ss,
    nkpc_scale=NKPC_SCALE, s_mode=:nkpc, r_path=rS, s_path=sS, sw19=true)

# ============================================================
# REPORT
# ============================================================
println("\n", "="^100); println("1. CONVERGENCE"); println("="^100)
@printf("warm start  asset %.6e   nkpc %.6e   inf-norm %.6e\n", wn.asset, wn.nkpc, wn.inf)
@printf("iterations %d   converged %s   wall %.2f min\n",
        sol_c.iterations, sol_c.f_converged, info_c.wall/60)
println("cumulative guards:"); cumg = guard_report_sw()

# solution-specific guard + [SW19] verification AT THE SOLUTION
_SW8_WIDEN[] = 0; _SW8_FALLBACK[] = 0
wsol = warm_resid(sol_c.zero, :nkpc)
@printf("\n[SW17] at the SOLUTION: [SW8] fallback %d, widened %d\n", _SW8_FALLBACK[], _SW8_WIDEN[])
fb_sol = _SW8_FALLBACK[]

println("\n", "="^100); println("4. [SW19] VERIFICATION AT THE CONVERGED SOLUTION"); println("="^100)
@printf("e_init.k     = %.14f\n", e_init.k)
@printf("path.k[1]    = %.14f\n", path_c.k[1])
@printf("residuals[1] = e_init.k - path.k[1] = %.6e  <- must be machine zero\n",
        e_init.k - path_c.k[1])
@printf("r_k[1]       = %.12f   r_k_hat = %.12f\n", path_c.r_k[1], path_c.r_k[1] - δ)
@printf("inside the [SW10] clamp [%.6f, %.6f]? %s\n",
        r_min, r_max, (r_min <= path_c.r_k[1] - δ <= r_max) ? "YES" : "NO")

if !sol_c.f_converged
    println("\n*** NOT CONVERGED -- nothing written ***"); exit(0)
end

println("\n", "="^100); println("2. ITEM (b) GATE and the PERIOD-0 BRANCH"); println("="^100)
gate = gate_lab_tax_sw(e_init, path_c, g)
gate_ok = print_gate_sw(gate)
s_flex_1 = s_flex_sw(g, c.ν; w = path_c.w[1], τ_n = path_c.τ_n[1],
                     w_init = e_init.w, τ_n_init = e_init.τ_n)
branch_ok = ts.Γ_n <= s_flex_1
@printf("\nGamma_n = %.10f   s_flex(w_1, tau_n_1) = %.10f\n", ts.Γ_n, s_flex_1)
println("min(Gamma_n, s_flex) selects: ", branch_ok ? "Gamma_n  (CORRECT)" : "s_flex  (*** INVALID ***)")

println("\n", "="^100); println("3. tau_n AND 1 - tau_n"); println("="^100)
@printf("%-22s %14s %14s %14s\n", "", "tau_n", "1 - tau_n", "vs SS (%)")
@printf("%-22s %14.8f %14.8f %14s\n", "initial SS", e_init.τ_n, 1-e_init.τ_n, "-")
for j in (1,2)
    @printf("%-22s %14.8f %14.8f %+14.2f\n", "rule-2 baseline j=$j", τn_base[j], 1-τn_base[j],
            100*(τn_base[j]/e_init.τ_n - 1))
    @printf("%-22s %14.8f %14.8f %+14.2f\n", "NO SUBSIDY      j=$j", path_c.τ_n[j], 1-path_c.τ_n[j],
            100*(path_c.τ_n[j]/e_init.τ_n - 1))
end

println("\n", "="^100); println("pi AND r vs THE RULE-2 BASELINE (annualised %)"); println("="^100)
@printf("%5s %12s %12s %12s   %12s %12s %12s\n",
        "j", "pi base", "pi cf", "diff pp", "r base", "r cf", "diff pp")
for j in vcat(collect(1:8), [12,20,40,100,400])
    @printf("%5d %12.4f %12.4f %+12.4f   %12.4f %12.4f %+12.4f\n", j,
            100*((1+π_base[j])^4-1), 100*((1+π_c[j])^4-1), 400*(π_c[j]-π_base[j]),
            100*((1+r_base[j])^4-1), 100*((1+path_c.r[j])^4-1), 400*(path_c.r[j]-r_base[j]))
end
@printf("\nmax|pi_cf - pi_base| over j=2..60 = %.4f ann pp   terminal pi = %.6f%% ann\n",
        400*maximum(abs.(π_c[2:60] .- π_base[2:60])), 100*((1+π_c[end])^4-1))

println("\n", "="^100); println("5. WELFARE (CE, paper bin convention)"); println("="^100)
ce = ce_summary_sw(e_init, path_c)
print_ce_sw(ce; label = "(no capital subsidies, rule-2 baseline)", paper = PAPER_CE)
println("\nfull decile profile:")
for (k,v) in enumerate(ce.deciles); @printf("  bin %2d  %8.4f\n", k, v); end

ok = sol_c.f_converged && gate_ok && branch_ok && fb_sol == 0 &&
     abs(e_init.k - path_c.k[1]) < 1e-12
println("\n", "="^100)
@printf("converged %s   gate %s   branch %s   [SW8]@solution %d   |res[1]| %.3e\n",
        sol_c.f_converged, gate_ok, branch_ok, fb_sol, abs(e_init.k - path_c.k[1]))
if ok
    save_object(joinpath(ROOT,"output","jld2_sw",SAVE_TAG*".jld2"),
                (e_init, e_final, path_c, π_c))
    println("wrote output/jld2_sw/", SAVE_TAG, ".jld2")
else
    println("NOT WRITTEN.")
end
println("\n=== done ===")
