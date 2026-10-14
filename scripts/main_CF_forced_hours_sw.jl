# ============================================================
# [SW22] FORCED PERIOD-0 HOURS -- baseline inertness check, then CBO
# ============================================================
#   julia --project -t auto scripts/main_CF_forced_hours_sw.jl baseline
#   julia --project -t auto scripts/main_CF_forced_hours_sw.jl cbo_debt
#
# Hours at period 0 are set to Gamma_n * ntilde_z EXACTLY in BOTH blocks. Gamma_n is
# NOT changed. The lockdown becomes a MANDATED level of hours rather than a ceiling,
# for the CBO column only.
#
# BASELINE run = the inertness check: the ceiling already binds there (s_flex_1 =
# 1.0357, margin +14.78%), so min(Gamma_n, s_flex) == Gamma_n and forcing must change
# NOTHING. If it reproduces `rule2_epsm` exactly, the modification is inert wherever
# the ceiling binds, which makes CBO the only column where the two specifications
# differ.
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
WHICH in ("baseline","cbo_debt") || error("ARGS[1] must be baseline or cbo_debt")
KAPPA_W=0.05; PHI_PI=1.5; NKPC_SCALE=240.0
TK_RULE2=-0.47520201; T1_RULE2=0.32446192; TU1_RULE2=0.31282997
SAVE_TAG = WHICH == "cbo_debt" ? "cf_cbo_debt_rule2_forced" : ""
SW19 = (WHICH == "cbo_debt")
PAPER = (mean=2.51, poor=4.00, middle=2.62, rich=0.50)
PAPER_BASE = (mean=1.82, poor=3.03, middle=2.00, rich=-0.20)
CEB = (mean=1.910, poor=3.350, middle=2.043, rich=-0.325)

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

println("=== [SW22] forced period-0 hours: ", WHICH, " ===")
_inp = load_object(joinpath(ROOT,"output","jld2_sw","steady_states.jld2"))
e_init, e_final_base = _inp[1], _inp[2]
h = e_init.h; t = e_init.t; g = e_init.g
δ = t.δ; π_bar = e_init.π_
i_ss_base = (1 + e_final_base.r)*(1 + π_bar) - 1
e_final = WHICH == "cbo_debt" ?
    load_object(joinpath(ROOT,"output","jld2_sw","steady_states_cf_cbo_debt_rule2.jld2")) :
    e_final_base
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
if WHICH == "cbo_debt"
    qb_u, qbs_u = cf_paths_cbo_debt_sw(ts, e_final; start_index=19, end_index=143,
        target=1.84*4.0, foreign_to_gdp=0.287*4.0, foreign_mode=:as_run)
else
    qb_u = ts.qb_path; qbs_u = ts.bstar_path_r
end
pars_h = (θ=c.θ, ν=c.ν, ies=c.ies, crra=c.crra, β=c.β, a_min=c.a_min, a_max=c.a_max,
          z_u=0.0, n_a=c.n_a, init_state=e_init.h.Pss, z_scale=z_scale)
h_path, _, _ = h_path_transition_mckay_init(ROOT, ts.Pee_path, ts.Puu_path; pars_h...)
Γ = ts.Γ_n

_o = load_object(joinpath(ROOT,"output","jld2_sw","rule2_stage1.jld2"))
epsm = epsm_baseline_sw(copy(_o[3].r), copy(_o[3].π_), i_ss_base, π_bar, PHI_PI)
_o = nothing; GC.gc()
println("loading rule2_epsm (baseline comparison + warm start) ..."); flush(stdout)
_o = load_object(joinpath(ROOT,"output","jld2_sw","rule2_epsm.jld2"))
pB=_o[3]; π_base=copy(_o[4]); r_base=copy(pB.r); τn_base=copy(pB.τ_n)
wB=copy(pB.w); yBv=copy(pB.y)
yB=[pB.y[1],pB.y[2]]; CB=[pB.C[1],pB.C[2]]; nB=[pB.n[1],pB.n[2]]; kB=[pB.k[1],pB.k[2]]
r0, s0 = warm_start_sw(pB, e_init, t; method=1, sw19=SW19)
_o=nothing; pB=nothing; GC.gc()
if WHICH == "cbo_debt"
    println("warm-starting from cf_cbo_debt (terminal r path already right) ..."); flush(stdout)
    _o = load_object(joinpath(ROOT,"output","jld2_sw","cf_cbo_debt.jld2"))
    global r0, s0 = warm_start_sw(_o[3], e_init, t; method=1, sw19=SW19)
    _o=nothing; GC.gc()
end

println("\n", "="^100); println("SOLVE -- [SW19] + [SW22] forced hours"); println("="^100)
path_c, sol_c, π_c, info_c = solve_transition_sw_c(e_init, e_final;
    τ_k_path=τk_u, T_path=T_u, Tu_path=Tu_u, G_path=ts.G_path,
    qb_path=qb_u, qbstar_path=qbs_u, h_path=h_path, π_path=ts.π_path,
    d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=Γ, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
    method=1, kappa_w=KAPPA_W, epsm_path=epsm, phi_pi=PHI_PI, i_taylor_ss=i_ss,
    nkpc_scale=NKPC_SCALE, s_mode=:nkpc, r_path=r0, s_path=s0,
    sw19=SW19, forced_hours=true)

println("\n", "="^100); println("RESULT -- ", WHICH, ", forced period-0 hours"); println("="^100)
@printf("iterations %d   converged %s   wall %.2f min   final inf-norm %.6e\n",
        sol_c.iterations, sol_c.f_converged, info_c.wall/60, sol_c.residual_norm)
guard_report_sw()
@printf("[SW19] %s   k[1] - e_init.k = %.3e   r_k_hat = %.12f\n", SW19 ? "ON" : "OFF", path_c.k[1]-e_init.k, path_c.r_k[1]-t.δ)
sol_c.f_converged || (println("\n*** NOT CONVERGED -- nothing written ***"); exit(0))

# --- item (b) gate, including t=1 -----------------------------------------
gate = gate_lab_tax_sw(e_init, path_c, g); gate_ok = print_gate_sw(gate)
# hours actually booked at j=1: s_1 recovered from n_pol must equal Gamma_n exactly
nt = e_init.ws.n_pol[1,:]
emp = [z for z in eachindex(e_init.h.z_grid) if e_init.h.z_grid[z] != 0]
s1_booked = sum(path_c[1].n_pol[1,z]/nt[z] for z in emp)/length(emp)
@printf("\ns_1 booked by the HOUSEHOLD block = %.14f   Gamma_n = %.14f   diff %.3e\n",
        s1_booked, Γ, s1_booked - Γ)
sflex1 = s_flex_sw(g, c.ν; w=path_c.w[1], τ_n=path_c.τ_n[1], w_init=e_init.w, τ_n_init=e_init.τ_n)
@printf("s_flex_1 at the solution = %.10f   margin vs Gamma_n %+.4f%%  (%s)\n",
        sflex1, 100*(sflex1/Γ-1),
        sflex1 < Γ ? "ceiling would NOT bind -> the mandate is doing real work" :
                     "ceiling would bind anyway -> forcing is inert here")

# --- the wedge -------------------------------------------------------------
wd = forced_hours_wedge_sw(e_init, g, path_c, Γ; j=1); print_wedge_sw(wd)

# --- tau_n, output ---------------------------------------------------------
println("\n---- tau_n and 1 - tau_n ----")
@printf("%-22s %14s %14s %14s\n", "", "tau_n", "1 - tau_n", "vs SS (%)")
@printf("%-22s %14.8f %14.8f %14s\n", "initial SS", e_init.τ_n, 1-e_init.τ_n, "-")
for j in (1,2)
    @printf("%-22s %14.8f %14.8f %+14.2f\n", "rule-2 baseline j=$j", τn_base[j], 1-τn_base[j],
            100*(τn_base[j]/e_init.τ_n-1))
    @printf("%-22s %14.8f %14.8f %+14.2f\n", "this run        j=$j", path_c.τ_n[j], 1-path_c.τ_n[j],
            100*(path_c.τ_n[j]/e_init.τ_n-1))
end
println("\n---- 2020Q2 OUTPUT (should be EXACTLY the 8.4002% baseline value, s_1 = Gamma_n) ----")
@printf("y[1] = %.12f   fall %.4f%%   |  rule-2 baseline y[1] = %.12f   fall %.4f%%\n",
        path_c.y[1], 100*(1-path_c.y[1]/y_ss), yB[1], 100*(1-yB[1]/y_ss))
@printf("difference in y[1] = %.3e\n", path_c.y[1]-yB[1])

println("\n---- y, C, N, K at j = 1 and j = 2 vs the rule-2 baseline ----")
@printf("%-6s %4s %16s %16s %14s %11s\n","series","j","rule-2 base","this run","diff","% diff")
for (nm,bb,cc) in (("y",yB,[path_c.y[1],path_c.y[2]]), ("C",CB,[path_c.C[1],path_c.C[2]]),
                   ("N",nB,[path_c.n[1],path_c.n[2]]), ("K",kB,[path_c.k[1],path_c.k[2]]))
    for j in 1:2
        @printf("%-6s %4d %16.10f %16.10f %+14.3e %+10.4f%%\n", nm, j, bb[j], cc[j],
                cc[j]-bb[j], 100*(cc[j]/bb[j]-1))
    end
end

println("\n---- pi and r vs the rule-2 baseline (annualised %) ----")
@printf("%5s %12s %12s %12s   %12s %12s %12s\n","j","pi base","pi run","diff pp","r base","r run","diff pp")
for j in vcat(collect(1:6),[12,20,40,100,143,400])
    @printf("%5d %12.4f %12.4f %+12.4f   %12.4f %12.4f %+12.4f\n", j,
            100*((1+π_base[j])^4-1), 100*((1+π_c[j])^4-1), 400*(π_c[j]-π_base[j]),
            100*((1+r_base[j])^4-1), 100*((1+path_c.r[j])^4-1), 400*(path_c.r[j]-r_base[j]))
end
@printf("\nmax|pi diff| j=2..60 = %.4f ann pp   terminal pi = %.6f%% ann\n",
        400*maximum(abs.(π_c[2:60].-π_base[2:60])), 100*((1+π_c[end])^4-1))
if WHICH == "cbo_debt"
    ic = intercept_channel_sw(i_ss, i_ss_base, PHI_PI, e_final.r)
    @printf("\n%6s %14s %14s %14s\n","j","solved d(pi)","intercept","residual")
    for j in (20,40,100,143)
        sv = 400*(π_c[j]-π_base[j])
        @printf("%6d %+14.4f %+14.4f %+14.4f\n", j, sv, 400*ic.d_pi, sv-400*ic.d_pi)
    end
end

# --- welfare ---------------------------------------------------------------
ce = ce_summary_sw(e_init, path_c)
println("\n---- WELFARE (paper bin convention) ----")
@printf("%-8s %12s %12s %14s   %12s %14s\n",
        "group","rule-2 base","this run","DIFFERENCE","paper level","paper diff")
for (nm,bv,cv,pl,pb) in (("mean",CEB.mean,ce.mean,PAPER.mean,PAPER_BASE.mean),
                         ("poor",CEB.poor,ce.poor,PAPER.poor,PAPER_BASE.poor),
                         ("middle",CEB.middle,ce.middle,PAPER.middle,PAPER_BASE.middle),
                         ("rich",CEB.rich,ce.rich,PAPER.rich,PAPER_BASE.rich))
    @printf("%-8s %12.3f %12.3f %+14.3f   %12.2f %+14.2f\n", nm, bv, cv, cv-bv, pl, pl-pb)
end
@printf("%-8s %12s %12.3f\n","min","-5.383",ce.min)
println("\nfull decile profile (level, and difference from the rule-2 baseline):")
curB = [3.3498,2.9171,2.7590,2.5988,2.2468,1.8393,1.4451,0.9105,0.3609,-0.3249]
@printf("%6s %12s %12s %12s\n","bin","this run","rule-2 base","diff")
for k in 1:10
    @printf("%6d %12.4f %12.4f %+12.4f\n", k, ce.deciles[k], curB[k], ce.deciles[k]-curB[k])
end

if WHICH == "baseline"
    println("\n---- INERTNESS CHECK vs rule2_epsm ----")
    @printf("%-10s %16s %16s\n","series","max|diff| j>=2","at j=1")
    for (nm,a,b) in (("r",path_c.r,r_base),("tau_n",path_c.τ_n,τn_base),
                     ("w",path_c.w,wB),("y",path_c.y,yBv),("pi",π_c,π_base))
        @printf("%-10s %16.6e %16.6e\n", nm, maximum(abs.(a[2:end].-b[2:end])), abs(a[1]-b[1]))
    end
    println("if these are ~0, forcing is INERT wherever the ceiling already binds,")
    println("so CBO is the ONLY column where the two specifications differ.")
end

ok = sol_c.f_converged && gate_ok && (!SW19 || abs(e_init.k-path_c.k[1]) < 1e-12)
@printf("\nconverged %s  gate %s  |res[1]| %.3e\n", sol_c.f_converged, gate_ok,
        abs(e_init.k-path_c.k[1]))
if SAVE_TAG != "" && ok
    save_object(joinpath(ROOT,"output","jld2_sw",SAVE_TAG*".jld2"),(e_init,e_final,path_c,π_c))
    println("wrote output/jld2_sw/", SAVE_TAG, ".jld2")
elseif SAVE_TAG != ""
    println("NOT WRITTEN.")
end
println("\n=== ", WHICH, " done ===")
