# ============================================================
# DEBT-FINANCED FISCAL STIMULUS, HETEROGENEITY, AND WELFARE
# [SW17] RULE-2 BASELINE  --  "Other" allocated at j = 1, epsm RE-DERIVED
# ============================================================
#     julia --project -t auto scripts/main_JMP_rule2_sw.jl
#
# WHAT CHANGES vs reduced_epsm
#   Table E2's "Other" (state/local grants, education, public services, federal
#   administration, global assistance) has no counterpart in the model, so the fed
#   debt path finances spending the model never does and tau_n absorbs the residue.
#   Rule 2 allocates it at j = 1 ONLY, by Table E2's internal logic:
#       state/local 13.1 + education 5.5 + public services 1.8 -> labour
#           (split checks:UI in their existing 21.4:19.4 ratio)
#       federal admin 1.8 + global assistance 0.3               -> firms
#   giving shares checks .4756 / UI .4311 / firms .0933 of X = $400bn (1.8714% of
#   annual 2019 GDP), the size that closes the j = 1 budget with tau_n at its
#   steady-state value.
#
#   EVERYTHING ELSE IS UNCHANGED: qb_path, G_path, bstar_path_r, the lockdown block,
#   and tau_k for every j >= 2.  Only T_path[1], Tu_path[1] and tau_k[1] move.
#
# WHY epsm IS RE-DERIVED, NOT REUSED
#   The saved epsm was backed out of a solve with the OLD fiscal paths.  Feeding it
#   unchanged would hold fixed a monetary residual computed for a different
#   economy.  So this script repeats the original two-step:
#     STAGE 1  solve with pi FED at the data path  (solve_transition_h_r_sw, no
#              monetary block) -> gives the new equilibrium r
#     STAGE 2  back epsm out of that solution, feed it to the reduced solver
#              (solve_transition_sw_c, pi eliminated) and VERIFY pi returns to data
#   That is exactly how baseline_stage3_kw05 -> reduced_epsm was built.
#
# [SW8] IS LEFT ENABLED (the default), because that is how the current baseline was
#   produced and this must be like-for-like.  But the script REFUSES TO SAVE if the
#   [SW8] grid fallback ever fired, since the returned tau_n would then not solve
#   the period budget and the "solution" would be suspect.
#
# WRITES output/jld2_sw/rule2_epsm.jld2 and rule2_stage1.jld2 ONLY if converged AND
# the item (b) gate passes AND the period-0 branch is Gamma_n AND no [SW8] fallback.
# reduced_epsm.jld2 is NEVER touched.
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
E2         = (checks = 21.4, ui = 19.4, firms = 30.8, other = 22.9)
SAVE_TAG   = "rule2_epsm"
STAGE1_TAG = "rule2_stage1"
GDP_2019   = 21380.0
# [SW13] two solves in this process (Stage 1 then Stage 2), same e_init/t/g
_sw_runs_allowed[] = 2

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

println("loading cached steady states ..."); flush(stdout)
_inp = load_object(joinpath(ROOT, "output", "jld2_sw", "steady_states.jld2"))
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

# ---- load the current baseline: warm start + old epsm for the comparison ----
println("loading reduced_epsm for the warm start and the old epsm ..."); flush(stdout)
_o = load_object(joinpath(ROOT, "output", "jld2_sw", "reduced_epsm.jld2"))
r0, s0 = warm_start_sw(_o[3], e_init, t; method = 1)
π_base_old = copy(_o[4]); τn_old = copy(_o[3].τ_n); r_old = copy(_o[3].r)
ce_old_path_v1 = nothing
_o = nothing; GC.gc()
_o = load_object(joinpath(ROOT, "output", "jld2_sw", "baseline_stage3_kw05.jld2"))
epsm_old = epsm_baseline_sw(copy(_o[3].r), copy(_o[3].π_), i_ss, π_bar, PHI_PI)
_o = nothing; GC.gc()

# ============================================================
# RULE 2: build the fiscal paths
# ============================================================
lab_sh = (13.1 + 5.5 + 1.8)/22.5
frm_sh = (1.8 + 0.3)/22.5
PC = lab_sh*E2.checks/(E2.checks + E2.ui)
PU = lab_sh*E2.ui/(E2.checks + E2.ui)
PF = frm_sh
@printf("\nrule 2 shares: checks %.6f   UI %.6f   firms %.6f   (sum %.6f)\n",
        PC, PU, PF, PC+PU+PF)

# period-0 GBC, used only to size X and tau_k (r_k held at the current baseline)
_ob = load_object(joinpath(ROOT, "output", "jld2_sw", "reduced_epsm.jld2"))
pb = _ob[3]
h0 = h_path[2]; z_dist0 = z_state[1,:]; n_ss = e_init.ws.n_pol[1,:]
r_k0 = pb.r_k[1]; r_b0 = pb.r[1]; y1_0 = pb.y[1]
qb0 = ts.qb_path[1]; qbp0 = ts.qb_path[2]; yp0 = pb.y[2]
T0 = ts.T_path[1]; Tu0 = ts.Tu_path[1]; G0 = ts.G_path[1]
d0 = ts.d_path[1]; cov0 = ts.covid_path[1]
_ob = nothing; GC.gc()

function period0(; τ_k, dTlev = 0.0, dTulev = 0.0)
    pr = ProgressUnknown()
    Tr = T0 + dTlev/y1_0; Tur = Tu0 + dTulev/y1_0
    f = (τ_n) -> _τ_n_excess_t0_port_sw_!(h0, t, g; τ_n, τ_k, T = Tr, Tu = Tur, G = G0,
        qbprime = qbp0, qb = qb0, yprime = yp0, d = d0, s = ts.Γ_n, n_ss, covid = cov0,
        z_dist = z_dist0, r_b = r_b0, r_k = r_k0, verbose = false, p = pr).excess_g
    τ_n = _solve_tau_n_sw(f, 0.0, 8.0)
    out = _τ_n_excess_t0_port_sw_!(h0, t, g; τ_n, τ_k, T = Tr, Tu = Tur, G = G0,
        qbprime = qbp0, qb = qb0, yprime = yp0, d = d0, s = ts.Γ_n, n_ss, covid = cov0,
        z_dist = z_dist0, r_b = r_b0, r_k = r_k0, verbose = false, p = pr)
    return (; τ_n = out.τ_n, fsub = -(τ_k*r_k0*out.k - caprev_ss), Tr, Tur)
end
base_sub = period0(; τ_k = -0.47).fsub
τ_k_for_sub(target; dTlev, dTulev) =
    find_zero((τ_k) -> period0(; τ_k, dTlev, dTulev).fsub - target, (-0.95, 0.10),
              Bisection(); xatol = 1e-12)
function rule2(X)
    dTl = PC*X; dTul = PU*X
    τ_k = τ_k_for_sub(base_sub + PF*X; dTlev = dTl, dTulev = dTul)
    a = period0(; τ_k, dTlev = dTl, dTulev = dTul)
    return a.τ_n - e_init.τ_n, τ_k, a
end
X2 = find_zero((X) -> rule2(X)[1], (0.0, 0.05), Bisection(); xatol = 1e-14)
_, τk2, a2 = rule2(X2)
@printf("X = %.10f model units = %.4f%% ann GDP = \$%.1f bn\n", X2, pct(X2), dol(X2))
@printf("tau_k[1] : %+.8f -> %+.8f\n", ts.τ_k_path[1], τk2)
@printf("T_path[1]: %.8f -> %.8f   (level +\$%.1f bn)\n", T0, a2.Tr, dol(PC*X2))
@printf("Tu_path[1]: %.8f -> %.8f  (level +\$%.1f bn)\n", Tu0, a2.Tur, dol(PU*X2))
@printf("firms slice: +\$%.1f bn of subsidy\n", dol(PF*X2))
@printf("period-0 GBC preview: tau_n = %.8f (SS %.8f)\n", a2.τ_n, e_init.τ_n)

T_r2  = copy(ts.T_path);   T_r2[1]  = a2.Tr
Tu_r2 = copy(ts.Tu_path);  Tu_r2[1] = a2.Tur
τk_r2 = copy(ts.τ_k_path); τk_r2[1] = τk2
@assert all(τk_r2[j] == ts.τ_k_path[j] for j in 2:length(τk_r2)) "tau_k changed beyond j=1"
@assert all(T_r2[j]  == ts.T_path[j]  for j in 2:length(T_r2))  "T changed beyond j=1"
@assert all(Tu_r2[j] == ts.Tu_path[j] for j in 2:length(Tu_r2)) "Tu changed beyond j=1"
println("assertions pass: only j = 1 differs in T, Tu and tau_k.")

n_per = length(τk_r2); n_r = _n_r_sw(n_per, 1); n_s = n_per - 2
DY_ss = union_aggregates_(e_init.h, g, e_init.ws; s=1.0, w=e_init.w, τ_n=e_init.τ_n, n_ss=n_ss)
wY_ss = e_init.w*DY_ss.Y

# ------------------------------------------------------------
# [SW17] SOLUTION-SPECIFIC GUARD CHECK.
# The [SW8] counters are CUMULATIVE over every residual evaluation in the process,
# so a non-zero count does not by itself condemn the solution -- the fallback firing
# at a wild intermediate Anderson iterate is by design harmless.  What matters is
# whether it fires AT THE FINAL POINT, because a fallback returns the grid argmin of
# |excess_g| rather than a root, i.e. a tau_n that does NOT balance that period's
# budget -- and excess_g is NOT part of the residual vector, so the reported
# inf-norm cannot detect it.  These helpers reset the counters and evaluate the
# residual ONCE at a given point, so the count that comes back belongs to that point
# alone.  The save gate uses THESE, not the cumulative totals.
# ------------------------------------------------------------
function stage1_guard_at(x)
    _SW8_WIDEN[] = 0; _SW8_FALLBACK[] = 0; _SW9_FLOOR[] = 0
    path_tmp = _initialize_path(e_init, e_final; qb_path = ts.qb_path)
    F = similar(x)
    _residuals_h_r_sw!(F, x; path = path_tmp, e_init = e_init, e_final = e_final,
        τ_k_path = τk_r2, T_path = T_r2, Tu_path = Tu_r2, G_path = ts.G_path,
        qb_path = ts.qb_path, qbstar_path = ts.bstar_path_r, h_path = h_path,
        π_path = ts.π_path, d_path = ts.d_path, covid_path = ts.covid_path,
        Γ_n = ts.Γ_n, Γ_c = ts.Γ_c, χ = c.χ, t = t, g = g, method = 1,
        kappa_w = KAPPA_W, wY_ss = wY_ss, nkpc_scale = NKPC_SCALE, s_mode = :nkpc)
    return (fb = _SW8_FALLBACK[], wd = _SW8_WIDEN[], inf = maximum(abs, F))
end

function stage2_guard_at(x, epsm_use)
    _SW8_WIDEN[] = 0; _SW8_FALLBACK[] = 0
    path_tmp = _initialize_path(e_init, e_final; qb_path = ts.qb_path)
    πv = Vector{Float64}(undef, n_per); F = similar(x)
    _residuals_sw_c!(F, x; π_vec = πv, epsm_path = epsm_use, i_taylor_ss = i_ss,
        phi_pi = PHI_PI, path = path_tmp, e_init = e_init, e_final = e_final,
        τ_k_path = τk_r2, T_path = T_r2, Tu_path = Tu_r2, G_path = ts.G_path,
        qb_path = ts.qb_path, qbstar_path = ts.bstar_path_r, h_path = h_path,
        π_path = ts.π_path, d_path = ts.d_path, covid_path = ts.covid_path,
        Γ_n = ts.Γ_n, Γ_c = ts.Γ_c, χ = c.χ, t = t, g = g, method = 1,
        kappa_w = KAPPA_W, wY_ss = wY_ss, nkpc_scale = NKPC_SCALE, s_mode = :nkpc)
    return (fb = _SW8_FALLBACK[], wd = _SW8_WIDEN[], inf = maximum(abs, F))
end

function warm_resid(x, s_mode, epsm_use)
    path = _initialize_path(e_init, e_final; qb_path = ts.qb_path)
    πv = Vector{Float64}(undef, n_per); F = similar(x)
    _residuals_sw_c!(F, x; π_vec=πv, epsm_path=epsm_use, i_taylor_ss=i_ss, phi_pi=PHI_PI,
        path=path, e_init=e_init, e_final=e_final, τ_k_path=τk_r2, T_path=T_r2,
        Tu_path=Tu_r2, G_path=ts.G_path, qb_path=ts.qb_path, qbstar_path=ts.bstar_path_r,
        h_path=h_path, π_path=ts.π_path, d_path=ts.d_path, covid_path=ts.covid_path,
        Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g, method=1, kappa_w=KAPPA_W,
        wY_ss=wY_ss, nkpc_scale=NKPC_SCALE, s_mode=s_mode)
    return (asset=maximum(abs, view(F,1:n_r)), nkpc=maximum(abs, view(F,(n_r+1):(n_r+n_s))),
            inf=maximum(abs, F))
end

# ============================================================
# STAGE 1 -- pi FED at the data path, no monetary block
# ============================================================
println("\n", "="^100)
println("STAGE 1 -- pi FED at the data path (solve_transition_h_r_sw)")
println("="^100)
w1 = warm_resid(vcat(r0, s0), :nkpc, epsm_old)   # epsm irrelevant here, pi is fed
@printf("warm start from reduced_epsm: asset %.6e   nkpc %.6e   inf-norm %.6e\n",
        w1.asset, w1.nkpc, w1.inf); flush(stdout)
t1 = time()
path1, sol1 = solve_transition_h_r_sw(e_init, e_final;
    τ_k_path=τk_r2, T_path=T_r2, Tu_path=Tu_r2, G_path=ts.G_path,
    qb_path=ts.qb_path, qbstar_path=ts.bstar_path_r, h_path=h_path, π_path=ts.π_path,
    d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
    method=1, kappa_w=KAPPA_W, nkpc_scale=NKPC_SCALE, s_mode=:nkpc, r_path=r0, s_path=s0)
wall1 = time() - t1
@printf("\nSTAGE 1: iterations = %d   converged = %s   wall = %.1f min\n",
        sol1.iterations, sol1.f_converged, wall1/60)
println("cumulative guard counters over the whole Stage 1 solve:")
cum1 = guard_report_sw()
if !sol1.f_converged
    println("\n*** STAGE 1 DID NOT CONVERGE -- nothing written, epsm cannot be re-derived ***")
    exit(0)
end
s1chk = stage1_guard_at(sol1.zero)
@printf("\n[SW17] AT STAGE 1's SOLUTION: [SW8] fallback fired %d times, widened %d, inf-norm %.6e\n",
        s1chk.fb, s1chk.wd, s1chk.inf)
println(s1chk.fb == 0 ?
    "  -> the solution point is CLEAN: tau_n solves every period budget exactly, so the\n" *
    "     cumulative count above came from intermediate Anderson iterates, which is the\n" *
    "     harmless case, and epsm derived from this path is trustworthy." :
    "  -> *** CONTAMINATED: tau_n does NOT balance the budget in some period AT THE\n" *
    "     SOLUTION. excess_g is not in the residual vector, so the inf-norm cannot see\n" *
    "     it. epsm derived from this path is NOT trustworthy. ***")

# back epsm out of the Stage 1 solution
epsm_new = epsm_baseline_sw(copy(path1.r), copy(path1.π_), i_ss, π_bar, PHI_PI)
@printf("\nnew epsm: max|epsm| = %.8f   at 2022Q2 (j=9) = %+.4f pp ann\n",
        maximum(abs, epsm_new), 400*epsm_new[9])
@printf("old epsm: max|epsm| = %.8f   at 2022Q2 (j=9) = %+.4f pp ann\n",
        maximum(abs, epsm_old), 400*epsm_old[9])

# ============================================================
# STAGE 2 -- feed the new epsm, pi eliminated; pi must return to data
# ============================================================
r1, s1 = warm_start_sw(path1, e_init, t; method = 1)
println("\n", "="^100)
println("STAGE 2 -- reduced solver, new epsm fed (pi must return to the data path)")
println("="^100)
w2 = warm_resid(vcat(r1, s1), :nkpc, epsm_new)
@printf("warm start from STAGE 1: asset %.6e   nkpc %.6e   inf-norm %.6e\n",
        w2.asset, w2.nkpc, w2.inf); flush(stdout)
path2, sol2, π2, info2 = solve_transition_sw_c(e_init, e_final;
    τ_k_path=τk_r2, T_path=T_r2, Tu_path=Tu_r2, G_path=ts.G_path,
    qb_path=ts.qb_path, qbstar_path=ts.bstar_path_r, h_path=h_path, π_path=ts.π_path,
    d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=ts.Γ_n, Γ_c=ts.Γ_c, χ=c.χ, t=t, g=g,
    method=1, kappa_w=KAPPA_W, epsm_path=epsm_new, phi_pi=PHI_PI, i_taylor_ss=i_ss,
    nkpc_scale=NKPC_SCALE, s_mode=:nkpc, r_path=r1, s_path=s1)
@printf("\nSTAGE 2: iterations = %d   converged = %s   wall = %.1f min\n",
        sol2.iterations, sol2.f_converged, info2.wall/60)
println("cumulative guard counters (Stage 1 + Stage 2):")
sw8w, sw8f, sw9f = guard_report_sw()
s2chk = stage2_guard_at(sol2.zero, epsm_new)
@printf("\n[SW17] AT STAGE 2's SOLUTION: [SW8] fallback fired %d times, widened %d, inf-norm %.6e\n",
        s2chk.fb, s2chk.wd, s2chk.inf)
println(s2chk.fb == 0 ? "  -> the reported baseline point is CLEAN." :
                        "  -> *** CONTAMINATED, see the Stage 1 note. ***")

# ============================================================
# REPORT
# ============================================================
println("\n", "="^100)
println("1. CONVERGENCE")
println("="^100)
@printf("STAGE 1  warm start: asset %.6e  nkpc %.6e  inf %.6e\n", w1.asset, w1.nkpc, w1.inf)
@printf("STAGE 1  iterations %d   converged %s   wall %.1f min\n",
        sol1.iterations, sol1.f_converged, wall1/60)
@printf("STAGE 2  warm start: asset %.6e  nkpc %.6e  inf %.6e\n", w2.asset, w2.nkpc, w2.inf)
@printf("STAGE 2  iterations %d   converged %s   wall %.1f min\n",
        sol2.iterations, sol2.f_converged, info2.wall/60)

println("\n---- regression (b): does pi return to the data path? ----")
@printf("max |pi_solved - pi_data| over j=1..60 = %.6e  (= %.6f ann pp)\n",
        maximum(abs.(π2[1:60] .- ts.π_path[1:60])),
        400*maximum(abs.(π2[1:60] .- ts.π_path[1:60])))
print_pi_table_sw(ts.π_path, π2, π_bar; label = "rule 2")

println("\n", "="^100)
println("2. ITEM (b) GATE and the PERIOD-0 BRANCH")
println("="^100)
gate = gate_lab_tax_sw(e_init, path2, g)
gate_ok = print_gate_sw(gate)
s_flex_1 = s_flex_sw(g, c.ν; w = path2.w[1], τ_n = path2.τ_n[1],
                     w_init = e_init.w, τ_n_init = e_init.τ_n)
branch_ok = ts.Γ_n <= s_flex_1
@printf("\nGamma_n = %.10f   s_flex(w_1, tau_n_1) = %.10f\n", ts.Γ_n, s_flex_1)
println("min(Gamma_n, s_flex) selects: ", branch_ok ? "Gamma_n  (CORRECT)" :
        "s_flex  (*** RUN IS INVALID ***)")

println("\n", "="^100)
println("3. tau_n AND 1 - tau_n  (Figure D1(d) plots 1 - tau_n)")
println("="^100)
@printf("%-14s %14s %14s %14s\n", "", "tau_n", "1 - tau_n", "vs SS (%)")
@printf("%-14s %14.8f %14.8f %14s\n", "initial SS", e_init.τ_n, 1-e_init.τ_n, "-")
for j in (1, 2)
    @printf("%-14s %14.8f %14.8f %+14.2f\n", "rule 2  j=$j", path2.τ_n[j], 1-path2.τ_n[j],
            100*(path2.τ_n[j]/e_init.τ_n - 1))
    @printf("%-14s %14.8f %14.8f %+14.2f\n", "  old    j=$j", τn_old[j], 1-τn_old[j],
            100*(τn_old[j]/e_init.τ_n - 1))
end
@printf("\nj=1: 1 - tau_n = %.8f  ->  %s\n", 1-path2.τ_n[1],
        path2.τ_n[1] > 1 ? "*** STILL NEGATIVE (retention above 1) ***" : "POSITIVE (retention below 1)")

println("\n", "="^100)
println("4. epsm: NEW vs OLD")
println("="^100)
@printf("%5s %16s %16s %16s %14s\n", "j", "old (ann pp)", "new (ann pp)", "diff (ann pp)", "rel change")
for j in (2, 6, 9, 10, 12, 20)
    @printf("%5d %16.4f %16.4f %16.4f %13.2f%%\n", j, 400*epsm_old[j], 400*epsm_new[j],
            400*(epsm_new[j]-epsm_old[j]),
            100*(epsm_new[j]/epsm_old[j] - 1))
end
@printf("\nmax |epsm_new - epsm_old| over j=1..60 = %.4f ann pp\n",
        400*maximum(abs.(epsm_new[1:60] .- epsm_old[1:60])))

println("\n", "="^100)
println("5. WELFARE (CE, paper bin convention)")
println("="^100)
ce = ce_summary_sw(e_init, path2)
@printf("%-10s %12s %12s %12s\n", "", "rule 2", "current", "paper")
@printf("%-10s %12.3f %12.3f %12.3f\n", "mean",   ce.mean,   1.797, 1.82)
@printf("%-10s %12.3f %12s %12s\n",     "min",    ce.min,    "-2.436", "-")
@printf("%-10s %12.3f %12.3f %12.3f\n", "poor",   ce.poor,   3.004, 3.03)
@printf("%-10s %12.3f %12.3f %12.3f\n", "middle", ce.middle, 1.971, 2.00)
@printf("%-10s %12.3f %12.3f %12.3f\n", "rich",   ce.rich,  -0.260, -0.20)
println("\nfull decile profile (rule 2), against the current baseline:")
cur = [3.0038, 2.6353, 2.5576, 2.4667, 2.1531, 1.7885, 1.4297, 0.9336, 0.4260, -0.2598]
@printf("%6s %12s %12s %12s\n", "bin", "rule 2", "current", "diff")
for k in 1:10
    @printf("%6d %12.4f %12.4f %+12.4f\n", k, ce.deciles[k], cur[k], ce.deciles[k]-cur[k])
end

# ============================================================
# SAVE
# ============================================================
# Gate on the SOLUTION-SPECIFIC fallback counts, not the cumulative totals: the
# fallback firing at an intermediate Anderson iterate is harmless by design, while
# firing at the solution would mean tau_n does not balance some period's budget.
ok = sol2.f_converged && gate_ok && branch_ok && s1chk.fb == 0 && s2chk.fb == 0
println("\n", "="^100)
@printf("converged %s   gate %s   branch %s\n", sol2.f_converged, gate_ok, branch_ok)
@printf("[SW8] fallback  cumulative %d   AT Stage 1 solution %d   AT Stage 2 solution %d\n",
        sw8f, s1chk.fb, s2chk.fb)
if ok
    save_object(joinpath(ROOT,"output","jld2_sw",STAGE1_TAG*".jld2"), (e_init, e_final, path1))
    save_object(joinpath(ROOT,"output","jld2_sw",SAVE_TAG*".jld2"), (e_init, e_final, path2, π2))
    println("wrote output/jld2_sw/", STAGE1_TAG, ".jld2  and  ", SAVE_TAG, ".jld2")
    println("reduced_epsm.jld2 NOT touched.")
else
    println("NOT WRITTEN.")
end
println("\n=== rule-2 re-solve done. Next: scripts/main_CSI_rule2_sw.jl ===")
