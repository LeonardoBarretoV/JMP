# ============================================================
# [SW20] INERTNESS CHECK / CBO SOLVE
# ============================================================
#   julia --project -t auto scripts/check_sw20_sw.jl <run>
#   <run> in: baseline | no_capital_subsidy | only_checks | lower_foreign | cbo_debt
#
# The four runs where Gamma_n binds must REPRODUCE with [SW20] on: the lifted row
# reads s_1 - Gamma_n there, is satisfied exactly at s_1 = Gamma_n, and the
# coordinate decouples.  Reports the period-0 branch MARGIN before and after, since
# lower_foreign sits at only +4.55% and is the one that could flip.
# cbo_debt is the run the lift is FOR: it sat at -1.20%, i.e. on the s_flex branch.
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

WHICH = isempty(ARGS) ? "baseline" : ARGS[1]
KAPPA_W = 0.05; PHI_PI = 1.5; NKPC_SCALE = 240.0
TK_RULE2 = -0.47520201; T1_RULE2 = 0.32446192; TU1_RULE2 = 0.31282997
# per-run: (saved solution to compare against, tau_k[1], T[1], Tu[1], Lambda_bstar_f, Lambda_b_f, own e_final?)
CFG = Dict(
 "baseline"           => (src="rule2_epsm",                  τk=TK_RULE2,     T1=T1_RULE2,   Tu1=TU1_RULE2,  Λbs=0.287*4.0, Λb=nothing,   ss=nothing),
 "no_capital_subsidy" => (src="cf_no_capital_subsidy_rule2", τk=nothing,      T1=T1_RULE2,   Tu1=TU1_RULE2,  Λbs=0.287*4.0, Λb=nothing,   ss=nothing),
 "only_checks"        => (src="cf_only_checks_rule2",        τk=nothing,      T1=0.97110527, Tu1=0.005,      Λbs=0.287*4.0, Λb=nothing,   ss=nothing),
 "lower_foreign"      => (src="cf_lower_foreign_rule2",      τk=TK_RULE2,     T1=T1_RULE2,   Tu1=TU1_RULE2,  Λbs=0.18*4.0,  Λb=nothing,   ss="steady_states_cf_lower_foreign_rule2"),
 "cbo_debt"           => (src="cf_cbo_debt",                 τk=TK_RULE2,     T1=T1_RULE2,   Tu1=TU1_RULE2,  Λbs=0.287*4.0, Λb=1.84*4.0,  ss="steady_states_cf_cbo_debt_rule2"))
haskey(CFG, WHICH) || error("unknown run: $WHICH")
cf = CFG[WHICH]
SAVE_TAG = WHICH == "cbo_debt" ? "cf_cbo_debt_rule2" : ""

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

println("=== [SW20] run: ", WHICH, " ===")
_inp = load_object(joinpath(ROOT,"output","jld2_sw","steady_states.jld2"))
e_init, e_final_base = _inp[1], _inp[2]
h = e_init.h; t = e_init.t; g = e_init.g
δ = t.δ; π_bar = e_init.π_
i_ss_base = (1 + e_final_base.r)*(1 + π_bar) - 1
e_final = isnothing(cf.ss) ? e_final_base :
          load_object(joinpath(ROOT,"output","jld2_sw",cf.ss*".jld2"))
i_ss = (1 + e_final.r)*(1 + π_bar) - 1
@printf("e_final r_f = %.4f%% ann   i_taylor_ss = %.4f%% ann\n",
        100*((1+e_final.r)^4-1), 100*((1+i_ss)^4-1))

b_target = e_init.b / y(e_init)
b_path = let
    T=700; H=100; ρB=0.08; bl = Array{Float64,1}(undef, T+H)
    bl[1]=b_target; bl[2]=b_target; bl[T:end] .= b_target
    for i in 3:T-1; bl[i] = bl[1]*ρB^(i-1) + (1-ρB^(i-1))*b_target; end
    bl
end
ts = make_time_series_sw(c, e_init, e_final, b_path, cf.Λbs, g)
T_u  = copy(ts.T_path);   T_u[1]  = cf.T1
Tu_u = copy(ts.Tu_path);  Tu_u[1] = cf.Tu1
τk_u = copy(ts.τ_k_path); τk_u[1] = isnothing(cf.τk) ? e_init.τ_k : cf.τk
if WHICH == "lower_foreign"
    qb_u = ts.qb_path; qbs_u = cf_qbstar_lower_foreign_sw(ts; drop = 0.10, jlast = 20)
elseif WHICH == "cbo_debt"
    qb_u, qbs_u = cf_paths_cbo_debt_sw(ts, e_final; start_index = 19, end_index = 143,
        target = 1.84*4.0, foreign_to_gdp = 0.287*4.0, foreign_mode = :as_run)
else
    qb_u = ts.qb_path; qbs_u = ts.bstar_path_r
end

pars_h = (θ=c.θ, ν=c.ν, ies=c.ies, crra=c.crra, β=c.β, a_min=c.a_min, a_max=c.a_max,
          z_u=0.0, n_a=c.n_a, init_state=e_init.h.Pss, z_scale=z_scale)
h_path, _, _ = h_path_transition_mckay_init(ROOT, ts.Pee_path, ts.Puu_path; pars_h...)
_o = load_object(joinpath(ROOT,"output","jld2_sw","rule2_stage1.jld2"))
epsm = epsm_baseline_sw(copy(_o[3].r), copy(_o[3].π_), i_ss_base, π_bar, PHI_PI)
_o = nothing; GC.gc()

println("loading ", cf.src, " for the warm start and the BEFORE comparison ..."); flush(stdout)
_o = load_object(joinpath(ROOT,"output","jld2_sw",cf.src*".jld2"))
pOLD = _o[3]
πOLD = length(_o) >= 4 ? copy(_o[4]) : copy(pOLD.π_)
rOLD = copy(pOLD.r); τnOLD = copy(pOLD.τ_n); wOLD = copy(pOLD.w); yOLD = copy(pOLD.y)
sOLD1 = let nt = e_init.ws.n_pol[1,:], emp = [z for z in eachindex(e_init.h.z_grid) if e_init.h.z_grid[z] != 0]
    sum(pOLD[1].n_pol[1,z]/nt[z] for z in emp)/length(emp)
end
r0, s0 = warm_start_sw(pOLD, e_init, t; method = 1, sw19 = true, sw20 = true)
_o = nothing; GC.gc()

sflex_before = s_flex_sw(g, c.ν; w = wOLD[1], τ_n = τnOLD[1], w_init = e_init.w, τ_n_init = e_init.τ_n)
@printf("\nBEFORE [SW20]: Gamma_n = %.10f   s_flex_1 = %.10f   margin %+.4f%%   branch %s\n",
        ts.Γ_n, sflex_before, 100*(sflex_before/ts.Γ_n - 1),
        sflex_before >= ts.Γ_n ? "Gamma_n" : "s_flex")
@printf("               s_1 implied by the saved hours policy = %.10f\n", sOLD1)
@printf("warm start s_1 = %.10f  (length %d, expected %d)\n",
        s0[1], length(s0), _n_s_sw(length(τk_u), true))

println("\n", "="^100); println("SOLVE with [SW19] + [SW20]"); println("="^100)
path_c, sol_c, π_c, info_c = solve_transition_sw_c(e_init, e_final;
    τ_k_path=τk_u, T_path=T_u, Tu_path=Tu_u, G_path=ts.G_path,
    qb_path=qb_u, qbstar_path=qbs_u, h_path=h_path, π_path=ts.π_path,
    d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
    method=1, kappa_w=KAPPA_W, epsm_path=epsm, phi_pi=PHI_PI, i_taylor_ss=i_ss,
    nkpc_scale=NKPC_SCALE, s_mode=:nkpc, r_path=r0, s_path=s0, sw19=true, sw20=true)

println("\n", "="^100); println("RESULT -- ", WHICH); println("="^100)
@printf("iterations %d   converged %s   wall %.2f min\n",
        sol_c.iterations, sol_c.f_converged, info_c.wall/60)
guard_report_sw()
n_per = length(τk_u); n_r = _n_r_sw(n_per, 1, true); n_s = _n_s_sw(n_per, true)
s1_new = sol_c.zero[n_r + 1]
sflex_after = s_flex_sw(g, c.ν; w = path_c.w[1], τ_n = path_c.τ_n[1],
                        w_init = e_init.w, τ_n_init = e_init.τ_n)
println("\n---- PERIOD-0 BRANCH, before vs after ----")
@printf("Gamma_n            = %.10f\n", ts.Γ_n)
@printf("s_flex_1  BEFORE   = %.10f   margin %+.4f%%\n", sflex_before, 100*(sflex_before/ts.Γ_n-1))
@printf("s_flex_1  AFTER    = %.10f   margin %+.4f%%   (moved %+.4e)\n",
        sflex_after, 100*(sflex_after/ts.Γ_n-1), sflex_after - sflex_before)
@printf("s_1 solved         = %.10f   = min(Gamma_n, s_flex)? %s\n",
        s1_new, abs(s1_new - min(ts.Γ_n, sflex_after)) < 1e-9 ? "YES" : "NO")
@printf("active branch      = %s\n", s1_new < ts.Γ_n - 1e-10 ? "s_flex  (lockdown SLACK)" : "Gamma_n (lockdown BINDS)")
@printf("[SW19] residuals[1]= %.3e   r_k_hat = %.12f\n", e_init.k - path_c.k[1], path_c.r_k[1]-δ)

println("\n---- 2020Q2 OUTPUT: does the Gamma_n calibration still hold? ----")
y_ss = y(e_init)
@printf("y[1] = %.10f   y_ss = %.10f   fall = %.4f%%   (calibration target 8.4%%)\n",
        path_c.y[1], y_ss, 100*(1 - path_c.y[1]/y_ss))
@printf("y[1] BEFORE = %.10f   fall = %.4f%%\n", yOLD[1], 100*(1 - yOLD[1]/y_ss))

gate = gate_lab_tax_sw(e_init, path_c, g); gate_ok = print_gate_sw(gate)
println("\n---- reproduction vs ", cf.src, " ----")
@printf("%-10s %16s %16s\n", "series", "max|diff| j>=2", "at j=1")
for (nm,a,b) in (("r",path_c.r,rOLD),("tau_n",path_c.τ_n,τnOLD),("w",path_c.w,wOLD),
                 ("y",path_c.y,yOLD),("pi",π_c,πOLD))
    @printf("%-10s %16.6e %16.6e\n", nm, maximum(abs.(a[2:end].-b[2:end])), abs(a[1]-b[1]))
end
ce = ce_summary_sw(e_init, path_c)
@printf("\nCE: mean %.4f  poor %.4f  middle %.4f  rich %.4f  min %.4f\n",
        ce.mean, ce.poor, ce.middle, ce.rich, ce.min)
@printf("tau_n j=1 = %.8f (1-tau_n = %+.8f)   j=2 = %.8f\n",
        path_c.τ_n[1], 1-path_c.τ_n[1], path_c.τ_n[2])

if SAVE_TAG != ""
    branch_ok = true   # [SW20] either branch is now legitimate
    ok = sol_c.f_converged && gate_ok && abs(e_init.k - path_c.k[1]) < 1e-12
    @printf("\nconverged %s  gate %s  |res[1]| %.3e\n", sol_c.f_converged, gate_ok,
            abs(e_init.k - path_c.k[1]))
    if ok
        save_object(joinpath(ROOT,"output","jld2_sw",SAVE_TAG*".jld2"),(e_init,e_final,path_c,π_c))
        println("wrote output/jld2_sw/", SAVE_TAG, ".jld2")
    else
        println("NOT WRITTEN.")
    end
end
println("\n=== ", WHICH, " done ===")
