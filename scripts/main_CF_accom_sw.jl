# ============================================================
# [SW23] INTERCEPT-CHANNEL DECOMPOSITION: the two terminal-SS-moving counterfactuals
#        re-solved with the intercept channel switched OFF.
#
#   julia --project -t auto scripts/main_CF_accom_sw.jl cbo_debt
#   julia --project -t auto scripts/main_CF_accom_sw.jl lower_foreign
#
# THE SPECIFICATION, and why it is not the one first written down.
#
# i_taylor_ss and epsm enter solve_transition_sw_c in EXACTLY ONE place
# (transition_sw.jl:1092) and only as their SUM:
#
#   pi_j = pi_bar + ((1+r_{j+1})(1+pi_{j+1}) - 1 - i_taylor_ss - epsm_j)/phi
#
# Neither is passed to _residuals_h_r_sw!.  So "rebuild i_taylor_ss from r_final_cf,
# then shift epsm by -(i_ss_cf - i_ss_base) at EVERY date" restores the constant to
# i_ss_base + epsm_base[j] -- bit-identical (+0.000e+00, measured) to the earlier
# attempt that held i_taylor_ss at the baseline and diverged.  With the tail
# r -> r_final_cf and epsm -> 0 the recursion's own fixed point is then
#
#   pi* = pi_bar + d(i_ss)/(phi - 1 - r_f)
#       = pi_bar + 0.6885 pp ann  (CBO) ,  pi_bar + 0.1057 pp ann  (lower foreign)
#
# against the hard pin pi[800] = pi_bar.  Same boundary conflict, same magnitude.
# Restoring the constant is neutral only if the TAIL r is also the baseline's, and in
# these two counterfactuals it is not -- that is what they change.
#
# WHAT IS DONE HERE.  i_taylor_ss stays rebuilt from each counterfactual's own
# r_final, so Fisher holds at the terminal SS.  epsm is shifted by -d(i_ss) over the
# transition and returned to baseline for j > J_STOP = 600:
#
#   epsm_accom[j] = epsm_base[j] - d(i_ss)   for j <= 600
#                 = epsm_base[j]             for j >  600
#
# The backward recursion contracts at (1+r_f)/phi = 0.6690 per period, so
# 0.669^100 = 3.5e-18: dates j <= 500 are numerically untouched by the switch point,
# while the terminal boundary stays exactly consistent (pi = pi_bar on (600, 800]).
# With beta = 0.98804, beta^560 = 1.2e-03, so at the ANNOUNCEMENT DATE -- which is
# where the CE welfare object lives -- this is welfare-equivalent to a permanent
# accommodation.  Economically: a pre-announced return to the original rule, 150 years
# out.  d(i_ss) is applied in QUARTERLY units.
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

WHICH = isempty(ARGS) ? "cbo_debt" : ARGS[1]
WHICH in ("cbo_debt","lower_foreign") || error("ARGS[1] must be cbo_debt or lower_foreign")
J_STOP     = 600
KAPPA_W    = 0.05
PHI_PI     = 1.5
NKPC_SCALE = 240.0
TK_RULE2   = -0.47520201
T1_RULE2   = 0.32446192
TU1_RULE2  = 0.31282997
# forced hours ([SW22]) ONLY on CBO, matching the un-accommodated run it is compared to
FORCED   = (WHICH == "cbo_debt")
SW19     = true
UNACC    = WHICH == "cbo_debt" ? "cf_cbo_debt_rule2_forced" : "cf_lower_foreign_rule2"
SS_TAG   = WHICH == "cbo_debt" ? "steady_states_cf_cbo_debt_rule2" :
                                 "steady_states_cf_lower_foreign_rule2"
SAVE_TAG = UNACC * "_accom"
REPORT_J = [2, 6, 10, 20, 40, 100, 143, 400, 800]

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
Jf(f) = joinpath(ROOT,"output","jld2_sw",f*".jld2")

println("=== [SW23] intercept channel OFF: ", WHICH, " ===")
_inp = load_object(Jf("steady_states")); e_init, e_final_base = _inp[1], _inp[2]
h = e_init.h; t = e_init.t; g = e_init.g; δ = t.δ; π_bar = e_init.π_
i_ss_base = (1 + e_final_base.r)*(1 + π_bar) - 1
e_final   = load_object(Jf(SS_TAG))
i_ss      = (1 + e_final.r)*(1 + π_bar) - 1
d_iss     = i_ss - i_ss_base
y_ss      = y(e_init)
@printf("baseline r_f %.4f%% ann -> i_ss %.4f%% ann\n", 100*((1+e_final_base.r)^4-1), 100*((1+i_ss_base)^4-1))
@printf("THIS RUN r_f %.4f%% ann -> i_ss %.4f%% ann\n", 100*((1+e_final.r)^4-1), 100*((1+i_ss)^4-1))
@printf("d(i_ss) = %.12f quarterly = %+.4f pp (simple x4) = %+.4f pp (compounded ann)\n",
        d_iss, 400*d_iss, 100*((1+i_ss)^4 - (1+i_ss_base)^4))
@printf("forced hours [SW22] = %s   [SW19] = %s   J_STOP = %d\n", FORCED, SW19, J_STOP)

# ---- fiscal / foreign paths: identical to the un-accommodated run ---------------
b_target = e_init.b / y(e_init)
b_path = let
    T=700; H=100; ρB=0.08; bl=Array{Float64,1}(undef,T+H)
    bl[1]=b_target; bl[2]=b_target; bl[T:end].=b_target
    for i in 3:T-1; bl[i]=bl[1]*ρB^(i-1)+(1-ρB^(i-1))*b_target; end
    bl
end
Λ_bstar_f = WHICH == "lower_foreign" ? 0.18*4.0 : 0.287*4.0
ts = make_time_series_sw(c, e_init, e_final, b_path, Λ_bstar_f, g)
T_u  = copy(ts.T_path);   T_u[1]  = T1_RULE2
Tu_u = copy(ts.Tu_path);  Tu_u[1] = TU1_RULE2
τk_u = copy(ts.τ_k_path); τk_u[1] = TK_RULE2
@assert all(τk_u[j] == e_init.τ_k for j in 2:length(τk_u)) "tau_k not flat for j>=2"
if WHICH == "cbo_debt"
    qb_u, qbs_u = cf_paths_cbo_debt_sw(ts, e_final; start_index=19, end_index=143,
        target=1.84*4.0, foreign_to_gdp=0.287*4.0, foreign_mode=:as_run)
    @printf("qb: j=19 %.8f  j=143 %.8f   qbstar: j=20 %.8f  j=200 %.8f\n",
            qb_u[19], qb_u[143], qbs_u[20], qbs_u[200])
else
    qb_u  = ts.qb_path
    qbs_u = cf_qbstar_lower_foreign_sw(ts; drop=0.10, jlast=20)
    @printf("qbstar: j=2 %.8f  j=20 %.8f  j>=21 %.8f\n", qbs_u[2], qbs_u[20], qbs_u[21])
end
pars_h = (θ=c.θ, ν=c.ν, ies=c.ies, crra=c.crra, β=c.β, a_min=c.a_min, a_max=c.a_max,
          z_u=0.0, n_a=c.n_a, init_state=e_init.h.Pss, z_scale=z_scale)
h_path, _, _ = h_path_transition_mckay_init(ROOT, ts.Pee_path, ts.Puu_path; pars_h...)

# ---- epsm: baseline, then the accommodation ------------------------------------
_o = load_object(Jf("rule2_stage1"))
epsm_base = epsm_baseline_sw(copy(_o[3].r), copy(_o[3].π_), i_ss_base, π_bar, PHI_PI)
_o = nothing; GC.gc()
n_per = length(epsm_base)
epsm = copy(epsm_base)
for j in 1:min(J_STOP, n_per); epsm[j] -= d_iss; end

println("\n---- TERMINAL CONSISTENCY, checked BEFORE solving ----")
rf = e_final.r
fp(C) = π_bar + ((1+rf)*(1+π_bar) - 1 - C)/(PHI_PI - 1 - rf)
@printf("recursion constant (i_ss + epsm_j):  baseline %.14f\n", i_ss_base + epsm_base[400])
@printf("   un-accommodated  j=400 %.14f  (%+.3e vs baseline)\n",
        i_ss + epsm_base[400], i_ss + epsm_base[400] - i_ss_base - epsm_base[400])
@printf("   THIS RUN         j=400 %.14f  (%+.3e vs baseline)  <- channel off here\n",
        i_ss + epsm[400], i_ss + epsm[400] - i_ss_base - epsm_base[400])
@printf("   THIS RUN         j=700 %.14f  (%+.3e vs baseline)  <- channel back on\n",
        i_ss + epsm[700], i_ss + epsm[700] - i_ss_base - epsm_base[700])
@printf("terminal fixed point with the TAIL constant (j>%d, epsm=epsm_base):\n", J_STOP)
@printf("   pi* = %.12f = %.6f%% ann   dev from pi_bar %+.6f pp ann\n",
        fp(i_ss + epsm_base[n_per-1]), 100*((1+fp(i_ss + epsm_base[n_per-1]))^4-1),
        400*(fp(i_ss + epsm_base[n_per-1]) - π_bar))
@printf("   solver pin pi[%d] = pi_bar = %.6f%% ann  ->  CONSISTENT: conflict %+.3e\n",
        n_per, 100*((1+π_bar)^4-1), fp(i_ss + epsm_base[n_per-1]) - π_bar)
@printf("had the shift run to the terminal date instead, pi* = %.6f%% ann (dev %+.4f pp) -- REJECTED\n",
        100*((1+fp(i_ss - d_iss + epsm_base[n_per-1]))^4-1),
        400*(fp(i_ss - d_iss + epsm_base[n_per-1]) - π_bar))
cf_ = (1+rf)/PHI_PI
@printf("backward contraction (1+r_f)/phi = %.6f ; ^100 = %.3e -> j<=%d untouched by the switch\n",
        cf_, cf_^100, J_STOP-100)
@printf("beta^%d = %.3e -> announcement-date welfare weight past the switch\n", J_STOP-40, c.β^(J_STOP-40))
# tolerance 1e-7 quarterly = 4e-5 pp annualised.  The floor is the baseline epsm tail
# itself (epsm_base[799] = -6.5e-10, so 1.3e-09 of conflict), which the un-accommodated
# run carries too.  This still rejects the run-to-terminal-date variant (0.0017) by 1e4.
if abs(fp(i_ss + epsm_base[n_per-1]) - π_bar) > 1e-7
    error("terminal condition INCONSISTENT -- refusing to spend a solve")
end
println("terminal condition consistent. proceeding.")

# ---- warm start: the un-accommodated run, whose real side is nearly right -------
println("\nloading ", UNACC, " for the warm start and the comparison ..."); flush(stdout)
_o = load_object(Jf(UNACC))
pU = _o[3]
π_un = copy(_o[4]); r_un = copy(pU.r); τn_un = copy(pU.τ_n)
ce_un = ce_summary_sw(e_init, pU)
r0, s0 = warm_start_sw(pU, e_init, t; method=1, sw19=SW19)
_o = nothing; pU = nothing; GC.gc()

println("\n", "="^100)
println("SOLVE -- intercept channel OFF, ", WHICH)
println("="^100); flush(stdout)
path_c, sol_c, π_c, info_c = solve_transition_sw_c(e_init, e_final;
    τ_k_path=τk_u, T_path=T_u, Tu_path=Tu_u, G_path=ts.G_path,
    qb_path=qb_u, qbstar_path=qbs_u, h_path=h_path, π_path=ts.π_path,
    d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
    method=1, kappa_w=KAPPA_W, epsm_path=epsm, phi_pi=PHI_PI, i_taylor_ss=i_ss,
    nkpc_scale=NKPC_SCALE, s_mode=:nkpc, r_path=r0, s_path=s0,
    sw19=SW19, forced_hours=FORCED)

println("\n", "="^100)
println("RESULT -- ", WHICH, ", intercept channel OFF")
println("="^100)
@printf("iterations %d   converged %s   wall %.2f min   final inf-norm %.6e\n",
        sol_c.iterations, sol_c.f_converged, info_c.wall/60, sol_c.residual_norm)
guard_report_sw()
@printf("[SW19] %s   k[1] - e_init.k = %.3e   r_k_hat = %.12f\n",
        SW19 ? "ON" : "OFF", path_c.k[1]-e_init.k, path_c.r_k[1]-t.δ)
sol_c.f_converged || (println("\n*** NOT CONVERGED -- nothing written ***"); exit(0))

gate = gate_lab_tax_sw(e_init, path_c, g); gate_ok = print_gate_sw(gate)

println("\n---- PERIOD-0 BRANCH ----")
nt  = e_init.ws.n_pol[1,:]
emp = [z for z in eachindex(e_init.h.z_grid) if e_init.h.z_grid[z] != 0]
s1b = sum(path_c[1].n_pol[1,z]/nt[z] for z in emp)/length(emp)
sfx = s_flex_sw(g, c.ν; w=path_c.w[1], τ_n=path_c.τ_n[1], w_init=e_init.w, τ_n_init=e_init.τ_n)
@printf("s_1 booked by the household block = %.14f   Gamma_n = %.14f   diff %.3e\n",
        s1b, ts.Γ_n, s1b - ts.Γ_n)
@printf("s_flex_1 = %.10f   margin vs Gamma_n %+.4f%%   -> %s\n", sfx, 100*(sfx/ts.Γ_n-1),
        sfx < ts.Γ_n ? "ceiling would NOT bind" : "ceiling BINDS naturally")
branch_ok = FORCED ? true : (sfx >= ts.Γ_n)
FORCED && println("forced hours on: the blocks agree by construction, branch check not required")
!FORCED && println(branch_ok ? "BRANCH OK: Gamma_n is the active branch, as required" :
                               "*** BRANCH FAILED: s_flex_1 < Gamma_n and hours are not forced -> INVALID ***")

# ---- the baseline, for the three-way comparison --------------------------------
println("\nloading rule2_epsm for the baseline column ..."); flush(stdout)
_o = load_object(Jf("rule2_epsm"))
pB = _o[3]; π_B = copy(_o[4]); r_B = copy(pB.r); τn_B = copy(pB.τ_n); yB1 = pB.y[1]
ce_B = ce_summary_sw(e_init, pB)
_o = nothing; pB = nothing; GC.gc()

println("\n", "="^116)
println("INFLATION (annualised %):  baseline | un-accommodated | channel OFF")
println("="^116)
@printf("%6s %12s %12s %12s   %14s %14s %14s\n","j","baseline","un-accom","chan OFF",
        "unacc - base","OFF - base","OFF - unacc")
for j in REPORT_J
    j > length(π_c) && continue
    @printf("%6d %12.4f %12.4f %12.4f   %+14.4f %+14.4f %+14.4f\n", j,
            100*((1+π_B[j])^4-1), 100*((1+π_un[j])^4-1), 100*((1+π_c[j])^4-1),
            400*(π_un[j]-π_B[j]), 400*(π_c[j]-π_B[j]), 400*(π_c[j]-π_un[j]))
end
@printf("\nterminal pi: baseline %.6f%%   un-accom %.6f%%   channel OFF %.6f%% ann\n",
        100*((1+π_B[end])^4-1), 100*((1+π_un[end])^4-1), 100*((1+π_c[end])^4-1))
@printf("closed-form intercept effect -d(i_ss)/(phi-1-r_f) = %+.4f pp ann\n", -400*d_iss/(PHI_PI-1-rf))
@printf("mean (OFF - unacc) over j=2..%d = %+.4f pp ann   over j=2..60 = %+.4f pp ann\n",
        J_STOP-100, 400*sum(π_c[2:(J_STOP-100)].-π_un[2:(J_STOP-100)])/(J_STOP-101),
        400*sum(π_c[2:60].-π_un[2:60])/59)

println("\n", "="^116)
println("REAL RATE (annualised %):  baseline | un-accommodated | channel OFF")
println("="^116)
@printf("%6s %12s %12s %12s   %14s %14s %14s\n","j","baseline","un-accom","chan OFF",
        "unacc - base","OFF - base","OFF - unacc")
for j in REPORT_J
    j > length(path_c.r) && continue
    @printf("%6d %12.4f %12.4f %12.4f   %+14.4f %+14.4f %+14.4f\n", j,
            100*((1+r_B[j])^4-1), 100*((1+r_un[j])^4-1), 100*((1+path_c.r[j])^4-1),
            400*(r_un[j]-r_B[j]), 400*(path_c.r[j]-r_B[j]), 400*(path_c.r[j]-r_un[j]))
end

println("\n", "="^116)
println("WELFARE by decile (CE %, paper bin convention)")
println("   last column = THE INTERCEPT CHANNEL'S CONTRIBUTION")
println("="^116)
ce_c = ce_summary_sw(e_init, path_c)
@printf("%6s %12s %12s %12s   %14s %14s %16s\n","bin","baseline","un-accom","chan OFF",
        "unacc - base","OFF - base","OFF - unacc")
for k in 1:10
    @printf("%6d %12.4f %12.4f %12.4f   %+14.4f %+14.4f %+16.4f\n", k,
            ce_B.deciles[k], ce_un.deciles[k], ce_c.deciles[k],
            ce_un.deciles[k]-ce_B.deciles[k], ce_c.deciles[k]-ce_B.deciles[k],
            ce_c.deciles[k]-ce_un.deciles[k])
end
println("-"^116)
for (nm,kk) in (("mean",:mean),("poor",:poor),("middle",:middle),("rich",:rich),("min",:min))
    b_,u_,o_ = getfield(ce_B,kk), getfield(ce_un,kk), getfield(ce_c,kk)
    @printf("%6s %12.4f %12.4f %12.4f   %+14.4f %+14.4f %+16.4f\n", nm, b_,u_,o_, u_-b_, o_-b_, o_-u_)
end

println("\n---- supporting: tau_n and the period-0 aggregates ----")
@printf("%-14s %14.8f %14.8f %14.8f\n","tau_n j=1", τn_B[1], τn_un[1], path_c.τ_n[1])
@printf("%-14s %14.8f %14.8f %14.8f\n","tau_n j=2", τn_B[2], τn_un[2], path_c.τ_n[2])
@printf("%-14s %14.4f %14s %14.4f\n","y j=1 fall %", 100*(1-yB1/y_ss), "-", 100*(1-path_c.y[1]/y_ss))

ok = sol_c.f_converged && gate_ok && branch_ok && (!SW19 || abs(e_init.k-path_c.k[1]) < 1e-12)
@printf("\nconverged %s  gate %s  branch %s  |res[1]| %.3e  ->  %s\n",
        sol_c.f_converged, gate_ok, branch_ok, abs(e_init.k-path_c.k[1]),
        ok ? "WRITING" : "NOT WRITING")
if ok
    isfile(Jf(SAVE_TAG)) && error("$(SAVE_TAG).jld2 already exists -- refusing to overwrite")
    save_object(Jf(SAVE_TAG), (e_init, e_final, path_c, π_c, epsm))
    println("wrote output/jld2_sw/", SAVE_TAG, ".jld2")
end
println("\n=== ", WHICH, " (channel off) done ===")
