# ============================================================
# [SW21] CBO OFF RULE 2 -- MANUAL OUTER ITERATION ON s_1
# ============================================================
#     julia --project -t auto scripts/main_CF_cbo_rule2_outer_sw.jl
#
# WHY.  CBO is the only counterfactual whose solution sits ON the Gamma_n kink.
# Three routes failed:
#   1. rule-2 warm start, s_1 pinned at Gamma_n -> blew up at iteration 2 (2.85e+08 step)
#   2. pre-rule-2 warm start, s_1 pinned at Gamma_n -> CONVERGED to 5.44e-06 in 37
#      iterations, but the pinned value was the WRONG BRANCH: s_flex_1 = 0.8914921
#      < Gamma_n, so the blocks booked different hours and the item (b) gate failed
#      at j=1 by 2.125e-02 (1.566e-07 everywhere else)
#   3. [SW20], s_1 lifted into x -> no blow-up but a sustained limit cycle across the
#      kink: 1.03e-02 -> 1.56e-01 -> 2.01e-02 -> 7.55e-02 over 336 iterations
#
# Route 2 is the informative one: a WELL-CONDITIONED solution exists just below
# Gamma_n; the only problem is finding which s_1.  So make s_1 an OUTER unknown and
# keep the inner solve the one that converges in ~37 iterations:
#
#   inner:  solve the 1597 system with s_1 PINNED at a trial value (sw19, no sw20)
#   outer:  s_1 <- min(Gamma_n, s_flex(w_1, tau_n_1)) evaluated AT that solution
#
# The map is DECREASING in s_1 (lower s_1 -> lower n_0 -> lower r_k* -> higher w_1 ->
# higher s_flex), so h(s) = s - g(s) has dh/ds = 1 - g'(s) > 1: strictly increasing,
# unique root.  h(Gamma_n) = 0.9023012 - 0.8914921 = +0.0108 > 0, so the root lies
# BELOW Gamma_n, bracketed by [g(Gamma_n), Gamma_n] = [0.8914921, 0.9023012].
#
# To keep the number of full solves small, STAGE A locates s_1* with the CHEAP
# period-0 GBC map (seconds per evaluation, r_b and y' held at a reference run), then
# STAGE B does full solves updating s_1 from the SOLVED s_flex until the outer
# condition holds.  A run is written only if the outer condition holds to OUTER_TOL
# AND the gate passes AND the branch is genuinely s_flex.
# ============================================================
using LinearAlgebra, ProgressMeter, Roots, NLsolve, StructArrays, StatsBase, Tullio,
      Polyester, PrettyTables, CairoMakie, LaTeXStrings, DelimitedFiles, XLSX, CSV,
      DataFrames, ColorTypes, Printf, JLD2, MAT
const ROOT = normpath(joinpath(@__DIR__, ".."))
for f in ("types","preferences","technology_and_fiscal","calibration","auxiliary_functions",
          "household_problem_and_aggregation","stationary_equilibria","transition",
          "union_block","transition_sw","counterfactuals_sw","statistics","time_series")
    include(joinpath(ROOT, "src", f * ".jl"))
end
const _TOL = 1e-12; const _ZERO_FTOL = 1e-5
const _MAX_ITERS = 20_000; const _ZERO_MAX_ITERS = 450

KAPPA_W=0.05; PHI_PI=1.5; NKPC_SCALE=240.0
TK_RULE2=-0.47520201; T1_RULE2=0.32446192; TU1_RULE2=0.31282997
SAVE_TAG="cf_cbo_debt_rule2"; OUTER_TOL=1e-7
# [SW21] ARGS[1] = s_1 to pin (skips Stage A); ARGS[2] = MAX_OUTER.
# The outer map has |g'| ~ 0.97, so the plain s_1 <- g(s_1) update converges at only
# ~3% per step -- useless at 7.6 min per full solve.  h(s) = s - g(s) is strictly
# increasing, so drive it with a SECANT step computed from measured full-solution
# points instead, passing each new trial in as ARGS[1].
S1_ARG   = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : nothing
MAX_OUTER = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4
_sw_runs_allowed[] = MAX_OUTER + 1

Base.@kwdef struct Calib
    θ::Float64=1.0; ν::Float64=1.0; ies::Float64=1.0; crra::Float64=4.0
    β::Float64=0.9880431080390226
    a_min::Float64=-0.02; a_max::Float64=120.0; n_a::Int=1000; z_bar::Int=12
    Pee::Float64=0.975; Puu::Float64=0.4
    Λ_t::Float64=0.0; Λ_tu::Float64=0.005; Λ_a::Float64=5.06*4.0
    Λ_g::Float64=0.05; Λ_f::Float64=0.04; Λ_b::Float64=1.05*4.0
    Λ_bstar::Float64=0.31*4.0; Λ_rk::Float64=0.33
    r_target::Float64=(1.0+0.01)^(1/4)-1.0; π_bar::Float64=(1.0+0.02)^(1/4)-1.0
    mc::Float64=5.0*4.0; τ::Float64=0.181; nu_target::Float64=0.7
    A::Float64=0.035; μ::Float64=1.0
    d_::Float64=1.0; Γ_n::Float64=1.0; Γ_c::Float64=1.0; covid::Float64=0.0
    χ::Float64=0.44
end
c = Calib()
_, _, z_scale = calibration_unemploy_mckay_(ROOT, c.Pee, c.Puu, 0.0)

_inp = load_object(joinpath(ROOT,"output","jld2_sw","steady_states.jld2"))
e_init, e_final_base = _inp[1], _inp[2]
h = e_init.h; t = e_init.t; g = e_init.g
δ = t.δ; π_bar = e_init.π_
i_ss_base = (1 + e_final_base.r)*(1 + π_bar) - 1
e_final = load_object(joinpath(ROOT,"output","jld2_sw","steady_states_cf_cbo_debt_rule2.jld2"))
i_ss = (1 + e_final.r)*(1 + π_bar) - 1
y_ss = y(e_init)
@printf("e_final r_f = %.4f%% ann   i_taylor_ss = %.4f%% ann\n",
        100*((1+e_final.r)^4-1), 100*((1+i_ss)^4-1))

b_target = e_init.b / y(e_init)
b_path = let
    T=700; H=100; ρB=0.08; bl=Array{Float64,1}(undef,T+H)
    bl[1]=b_target; bl[2]=b_target; bl[T:end].=b_target
    for i in 3:T-1; bl[i]=bl[1]*ρB^(i-1)+(1-ρB^(i-1))*b_target; end
    bl
end
ts = make_time_series_sw(c, e_init, e_final, b_path, 0.287*4.0, g)
T_u  = copy(ts.T_path);   T_u[1]  = T1_RULE2
Tu_u = copy(ts.Tu_path);  Tu_u[1] = TU1_RULE2
τk_u = copy(ts.τ_k_path); τk_u[1] = TK_RULE2
qb_u, qbs_u = cf_paths_cbo_debt_sw(ts, e_final; start_index=19, end_index=143,
    target=1.84*4.0, foreign_to_gdp=0.287*4.0, foreign_mode=:as_run)
pars_h = (θ=c.θ, ν=c.ν, ies=c.ies, crra=c.crra, β=c.β, a_min=c.a_min, a_max=c.a_max,
          z_u=0.0, n_a=c.n_a, init_state=e_init.h.Pss, z_scale=z_scale)
h_path, _, _ = h_path_transition_mckay_init(ROOT, ts.Pee_path, ts.Puu_path; pars_h...)
z_state = zeros(length(h_path), length(e_init.h.z_grid))
z_state[1,:] = transpose(e_init.h.Pss)*h_path[1].P
for i=2:lastindex(h_path); z_state[i,:] = transpose(z_state[i-1,:])*h_path[i].P; end
h0 = h_path[2]; z_dist0 = z_state[1,:]; n_ss = e_init.ws.n_pol[1,:]
r_min = -δ + eps(); r_max = 1/h_path[1].u.β
Γ = ts.Γ_n

_o = load_object(joinpath(ROOT,"output","jld2_sw","rule2_stage1.jld2"))
epsm = epsm_baseline_sw(copy(_o[3].r), copy(_o[3].π_), i_ss_base, π_bar, PHI_PI)
_o = nothing; GC.gc()
println("loading rule2_epsm (baseline comparison) ..."); flush(stdout)
_o = load_object(joinpath(ROOT,"output","jld2_sw","rule2_epsm.jld2"))
pB=_o[3]; π_base=copy(_o[4]); r_base=copy(pB.r); τn_base=copy(pB.τ_n)
yB=[pB.y[1],pB.y[2]]; CB=[pB.C[1],pB.C[2]]; nB=[pB.n[1],pB.n[2]]; kB=[pB.k[1],pB.k[2]]
_o=nothing; pB=nothing; GC.gc()
println("loading cf_cbo_debt (warm start + reference r_b, y') ..."); flush(stdout)
_o = load_object(joinpath(ROOT,"output","jld2_sw","cf_cbo_debt.jld2"))
pR=_o[3]; rb_ref=pR.r[1]; yp_ref=pR.y[2]
r0, s0 = warm_start_sw(pR, e_init, t; method=1, sw19=true)
_o=nothing; pR=nothing; GC.gc()

# ---- STAGE A: cheap period-0 map g(s_1) ----------------------------------
function g_cheap(s1)
    n0 = labor_supply_union_(h0; s=s1, n_ss=n_ss, z_dist=z_dist0)
    r_k = _rk_star_sw19(t, e_init; τ_k=τk_u[1], n0=n0, d=ts.d_path[1],
                        covid=ts.covid_path[1], r_min, r_max)
    pr = ProgressUnknown()
    f = (τn) -> _τ_n_excess_t0_port_sw_!(h0, t, g; τ_n=τn, τ_k=τk_u[1], T=T_u[1],
        Tu=Tu_u[1], G=ts.G_path[1], qbprime=qb_u[2], qb=qb_u[1], yprime=yp_ref,
        d=ts.d_path[1], s=s1, n_ss, covid=ts.covid_path[1], z_dist=z_dist0,
        r_b=rb_ref, r_k=r_k, verbose=false, p=pr).excess_g
    τn = _solve_tau_n_sw(f, 0.0, 8.0)
    out = _τ_n_excess_t0_port_sw_!(h0, t, g; τ_n=τn, τ_k=τk_u[1], T=T_u[1], Tu=Tu_u[1],
        G=ts.G_path[1], qbprime=qb_u[2], qb=qb_u[1], yprime=yp_ref, d=ts.d_path[1],
        s=s1, n_ss, covid=ts.covid_path[1], z_dist=z_dist0, r_b=rb_ref, r_k=r_k,
        verbose=false, p=pr)
    return min(Γ, s_flex_sw(g, c.ν; w=out.w, τ_n=τn, w_init=e_init.w, τ_n_init=e_init.τ_n))
end
println("\n", "="^100); println("STAGE A -- locate s_1* with the cheap period-0 map"); println("="^100)
@printf("Gamma_n = %.10f   g(Gamma_n) = %.10f   h(Gamma_n) = %+.10f\n",
        Γ, g_cheap(Γ), Γ - g_cheap(Γ))
lo = g_cheap(Γ)
@printf("g(lo=%.10f) = %.10f   h(lo) = %+.10f\n", lo, g_cheap(lo), lo - g_cheap(lo))
s1_star = !isnothing(S1_ARG) ? S1_ARG : try
    find_zero((s) -> s - g_cheap(s), (lo*0.98, Γ), Bisection(); xatol=1e-12)
catch e
    println("bisection failed (", e, "); falling back to fixed-point iteration")
    sx = Γ; for _ in 1:80; sx = 0.5*sx + 0.5*g_cheap(sx); end; sx
end
@printf("\nSTAGE A s_1* = %.12f   (Gamma_n = %.12f, margin %+.4f%%)\n",
        s1_star, Γ, 100*(s1_star/Γ - 1))

# ---- STAGE B: full solves, updating s_1 from the SOLVED s_flex ------------
n_per=length(τk_u); n_r=_n_r_sw(n_per,1,true); n_s=_n_s_sw(n_per,false)
DY_ss = union_aggregates_(e_init.h, g, e_init.ws; s=1.0, w=e_init.w, τ_n=e_init.τ_n, n_ss=n_ss)
wY_ss = e_init.w*DY_ss.Y
s1 = s1_star; local path_c, sol_c, π_c, info_c, sflex_sol
for it in 1:MAX_OUTER
    println("\n", "="^100)
    @printf("STAGE B outer iteration %d -- s_1 pinned at %.12f\n", it, s1)
    println("="^100); flush(stdout)
    global path_c, sol_c, π_c, info_c = solve_transition_sw_c(e_init, e_final;
        τ_k_path=τk_u, T_path=T_u, Tu_path=Tu_u, G_path=ts.G_path,
        qb_path=qb_u, qbstar_path=qbs_u, h_path=h_path, π_path=ts.π_path,
        d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=Γ, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
        method=1, kappa_w=KAPPA_W, epsm_path=epsm, phi_pi=PHI_PI, i_taylor_ss=i_ss,
        nkpc_scale=NKPC_SCALE, s_mode=:nkpc, r_path=r0, s_path=s0, sw19=true, s1_pin=s1)
    @printf("inner: iterations %d  converged %s  wall %.2f min\n",
            sol_c.iterations, sol_c.f_converged, info_c.wall/60)
    sol_c.f_converged || (println("*** inner solve failed; stopping ***"); break)
    global sflex_sol = s_flex_sw(g, c.ν; w=path_c.w[1], τ_n=path_c.τ_n[1],
                                 w_init=e_init.w, τ_n_init=e_init.τ_n)
    tgt = min(Γ, sflex_sol)
    @printf("OUTER: s_1 pinned %.12f   s_flex at solution %.12f   target %.12f   gap %+.3e\n",
            s1, sflex_sol, tgt, s1 - tgt)
    if abs(s1 - tgt) < OUTER_TOL
        println("OUTER CONDITION HOLDS to ", OUTER_TOL, " -- s_1 = min(Gamma_n, s_flex) at the solution.")
        break
    end
    global r0, s0 = warm_start_sw(path_c, e_init, t; method=1, sw19=true)
    global s1 = tgt
end

# ---- report ---------------------------------------------------------------
println("\n", "="^100); println("RESULT -- CBO off rule 2, outer iteration on s_1"); println("="^100)
tgt = min(Γ, sflex_sol)
outer_ok = abs(s1 - tgt) < OUTER_TOL
@printf("s_1 = %.12f   Gamma_n = %.12f   margin %+.4f%%\n", s1, Γ, 100*(s1/Γ-1))
@printf("s_flex at the solution = %.12f   outer gap %+.3e   holds? %s\n",
        sflex_sol, s1-tgt, outer_ok ? "YES" : "NO")
@printf("active branch = %s\n", s1 < Γ - 1e-10 ? "s_flex (lockdown SLACK)" : "Gamma_n (lockdown BINDS)")
@printf("[SW19] residuals[1] = %.3e   r_k_hat = %.12f\n", e_init.k-path_c.k[1], path_c.r_k[1]-δ)
gate = gate_lab_tax_sw(e_init, path_c, g); gate_ok = print_gate_sw(gate)
println("\n---- 2020Q2 OUTPUT vs the 8.4% Gamma_n calibration target ----")
@printf("y[1] CBO = %.10f   fall %.4f%%   |  baseline y[1] = %.10f   fall %.4f%%\n",
        path_c.y[1], 100*(1-path_c.y[1]/y_ss), yB[1], 100*(1-yB[1]/y_ss))
@printf("miss vs the 8.4%% target: %+.4f pp\n", 100*(1-path_c.y[1]/y_ss) - 8.4)
println("\n---- y, C, N, K at j = 1 and j = 2 vs the rule-2 baseline ----")
@printf("%-6s %4s %16s %16s %14s %11s\n","series","j","rule-2 base","CBO","diff","% diff")
for (nm,bb,cc) in (("y",yB,[path_c.y[1],path_c.y[2]]), ("C",CB,[path_c.C[1],path_c.C[2]]),
                   ("N",nB,[path_c.n[1],path_c.n[2]]), ("K",kB,[path_c.k[1],path_c.k[2]]))
    for j in 1:2
        @printf("%-6s %4d %16.10f %16.10f %+14.3e %+10.4f%%\n", nm, j, bb[j], cc[j],
                cc[j]-bb[j], 100*(cc[j]/bb[j]-1))
    end
end
@printf("\ntau_n j=1 = %.8f (1-tau_n = %+.8f)   j=2 = %.8f\n",
        path_c.τ_n[1], 1-path_c.τ_n[1], path_c.τ_n[2])
ic = intercept_channel_sw(i_ss, i_ss_base, PHI_PI, e_final.r)
@printf("\n%6s %12s %12s %14s %14s %14s\n","j","pi base","pi cbo","solved d(pi)","intercept","residual")
for j in (20,40,100,143)
    sv = 400*(π_c[j]-π_base[j]); @printf("%6d %12.4f %12.4f %+14.4f %+14.4f %+14.4f\n",
        j, 100*((1+π_base[j])^4-1), 100*((1+π_c[j])^4-1), sv, 400*ic.d_pi, sv-400*ic.d_pi)
end
ce = ce_summary_sw(e_init, path_c)
@printf("\nCE: mean %.4f  poor %.4f  middle %.4f  rich %.4f  min %.4f\n",
        ce.mean, ce.poor, ce.middle, ce.rich, ce.min)
println("deciles: ", join([@sprintf("%.4f", v) for v in ce.deciles], ", "))
ok = sol_c.f_converged && gate_ok && outer_ok && abs(e_init.k-path_c.k[1]) < 1e-12
@printf("\nconverged %s  gate %s  outer %s  |res[1]| %.3e\n",
        sol_c.f_converged, gate_ok, outer_ok, abs(e_init.k-path_c.k[1]))
if ok
    save_object(joinpath(ROOT,"output","jld2_sw",SAVE_TAG*".jld2"),(e_init,e_final,path_c,π_c))
    println("wrote output/jld2_sw/", SAVE_TAG, ".jld2")
else
    println("NOT WRITTEN.")
end
println("\n=== done ===")
