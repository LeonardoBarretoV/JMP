# ============================================================
# [SW16] TWO CHEAP CHECKS BEFORE COMMITTING TO A RULE-2 RE-SOLVE
# ============================================================
#     julia --project -t auto scripts/check_rule2_sw.jl
#
# CHECK 1 -- does rule 2 actually unblock Stage 5?
#   The no-capital-subsidy run blew up in the ASSET block at 5.6e20, driven by the
#   SIZE of the tau_k shock.  Rule 2 makes that shock slightly LARGER (7.876% of
#   annual GDP against 7.702%), so it may fix retention-above-1 without unblocking
#   anything.  Measure it: the :flex warm-start residual by block, under rule 2's
#   fiscal paths, against the 0.835 measured at the current baseline.
#
#   Four evaluations, all at the SAME point (the current baseline solution), so the
#   only thing varying is the fiscal path:
#     A  current paths, baseline tau_k      -> isolates the :flex mode effect alone
#     B  current paths, counterfactual tau_k -> the known 0.835
#     C  rule-2 paths, rule-2 baseline tau_k -> how far rule 2 moves the baseline
#     D  rule-2 paths, counterfactual tau_k  -> the rule-2 counterfactual
#   The tau_k shock effect is B relative to A under the status quo, and D relative
#   to C under rule 2.  That is the comparison that answers the question.
#
#   CAVEAT stated up front: the true rule-2 warm start would be the rule-2 BASELINE
#   SOLUTION, which does not exist until the re-solve.  Evaluating at the current
#   baseline point conflates "distance the tau_k shock moves things" with "distance
#   rule 2 moves the baseline".  C is reported precisely so those can be separated.
#
# CHECK 2 -- does retention of 0.638 at j=2 break anything MECHANICALLY?
#   Reports the implied aggregate average labour tax rate path (comparable with
#   Figure D1(d)), the minimum GHH composite x = c - theta*n^(1+1/nu)/(1+1/nu) per
#   period (the object [SW9] floors at 1e-9), and where tau_n sits in its root-find
#   bracket.  Nothing is re-solved; nothing is written.
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

GDP_2019 = 21380.0
E2 = (checks = 21.4, ui = 19.4, firms = 30.8, other = 22.9)
KAPPA_W = 0.05; PHI_PI = 1.5; NKPC_SCALE = 240.0

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
e_init, e_final, p, π_base = o[1], o[2], o[3], o[4]
h = e_init.h; t = e_init.t; g = e_init.g
δ = t.δ; π_bar = e_init.π_
y_ss = y(e_init); Y_ANN = 4*y_ss
pct(x) = 100*x/Y_ANN
dol(x) = GDP_2019*x/Y_ANN
labtax_ss = e_init.lab_tax
caprev_ss = e_init.τ_k*(e_init.r + δ)*e_init.k
T_ss = e_init.lump_sum_T; Tu_ss = e_init.lump_sum_Tu

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
h0 = h_path[2]; z_dist0 = z_state[1,:]
n_ss = e_init.ws.n_pol[1,:]; s_0 = ts.Γ_n
r_k0 = p.r_k[1]; r_b0 = p.r[1]
qb0 = ts.qb_path[1]; qbp0 = ts.qb_path[2]; yp0 = p.y[2]
T0 = ts.T_path[1]; Tu0 = ts.Tu_path[1]; G0 = ts.G_path[1]
d0 = ts.d_path[1]; cov0 = ts.covid_path[1]

# ---- reconstruct rule 2's (X, tau_k) exactly ---------------------------------
function period0(; τ_k = -0.47, dTlev = 0.0, dTulev = 0.0)
    pr = ProgressUnknown()
    Tr = T0 + dTlev/p.y[1]; Tur = Tu0 + dTulev/p.y[1]
    f = (τ_n) -> _τ_n_excess_t0_port_sw_!(h0, t, g; τ_n, τ_k, T = Tr, Tu = Tur, G = G0,
        qbprime = qbp0, qb = qb0, yprime = yp0, d = d0, s = s_0, n_ss, covid = cov0,
        z_dist = z_dist0, r_b = r_b0, r_k = r_k0, verbose = false, p = pr).excess_g
    τ_n = _solve_tau_n_sw(f, 0.0, 8.0)
    out = _τ_n_excess_t0_port_sw_!(h0, t, g; τ_n, τ_k, T = Tr, Tu = Tur, G = G0,
        qbprime = qbp0, qb = qb0, yprime = yp0, d = d0, s = s_0, n_ss, covid = cov0,
        z_dist = z_dist0, r_b = r_b0, r_k = r_k0, verbose = false, p = pr)
    fsub = -(τ_k*r_k0*out.k - caprev_ss)
    return (; τ_n = out.τ_n, k = out.k, y = out.y, fsub, Tr, Tur)
end
base = period0()
@assert abs(base.τ_n - p.τ_n[1]) < 1e-8 "period-0 GBC not reproduced"

lab_sh = (13.1 + 5.5 + 1.8)/22.5; frm_sh = (1.8 + 0.3)/22.5
PC = lab_sh*E2.checks/(E2.checks + E2.ui)
PU = lab_sh*E2.ui/(E2.checks + E2.ui)
PF = frm_sh
function τ_k_for_subsidy(target; dTlev, dTulev)
    f = (τ_k) -> period0(; τ_k, dTlev, dTulev).fsub - target
    return find_zero(f, (-0.95, 0.10), Bisection(); xatol = 1e-12)
end
function rule2_resid(X)
    dTl = PC*X; dTul = PU*X
    τ_k = τ_k_for_subsidy(base.fsub + PF*X; dTlev = dTl, dTulev = dTul)
    a = period0(; τ_k, dTlev = dTl, dTulev = dTul)
    return a.τ_n - e_init.τ_n, τ_k, a, dTl, dTul
end
X2 = find_zero((X) -> rule2_resid(X)[1], (0.0, 0.05), Bisection(); xatol = 1e-14)
_, τk2, a2, dTl2, dTul2 = rule2_resid(X2)
@printf("\nRULE 2 reconstructed:  X = %.10f model units = %.4f%% ann GDP = \$%.0f bn\n",
        X2, pct(X2), dol(X2))
@printf("  tau_k = %.8f    tau_n = %.8f (SS %.8f)    subsidy = %.4f%% ann GDP\n",
        τk2, a2.τ_n, e_init.τ_n, pct(a2.fsub))
@printf("  T_path[1]: %.8f -> %.8f      Tu_path[1]: %.8f -> %.8f\n",
        T0, a2.Tr, Tu0, a2.Tur)

# rule-2 fiscal paths: only j = 1 changes
T_r2  = copy(ts.T_path);  T_r2[1]  = a2.Tr
Tu_r2 = copy(ts.Tu_path); Tu_r2[1] = a2.Tur
τk_r2 = copy(ts.τ_k_path); τk_r2[1] = τk2
τk_cf = copy(ts.τ_k_path); τk_cf[1] = e_init.τ_k    # no capital subsidies

# ============================================================
# CHECK 1
# ============================================================
DY_ss = union_aggregates_(e_init.h, g, e_init.ws; s=1.0, w=e_init.w, τ_n=e_init.τ_n,
                          n_ss=e_init.ws.n_pol[1,:])
wY_ss = e_init.w*DY_ss.Y
n_per = length(ts.qb_path); n_r = _n_r_sw(n_per, 1); n_s = n_per - 2
_o = load_object(joinpath(ROOT,"output","jld2_sw","baseline_stage3_kw05.jld2"))
i_ss = (1 + e_final.r)*(1 + π_bar) - 1
epsm = epsm_baseline_sw(copy(_o[3].r), copy(_o[3].π_), i_ss, π_bar, PHI_PI)
_o = nothing; GC.gc()
r0, s0 = warm_start_sw(p, e_init, t; method = 1)
x0 = vcat(r0, s0)

function resid_blocks(τ_k_path, T_path, Tu_path, s_mode)
    path = _initialize_path(e_init, e_final; qb_path = ts.qb_path)
    πv = Vector{Float64}(undef, n_per); F = similar(x0)
    _residuals_sw_c!(F, x0; π_vec=πv, epsm_path=epsm, i_taylor_ss=i_ss, phi_pi=PHI_PI,
        path=path, e_init=e_init, e_final=e_final, τ_k_path=τ_k_path,
        T_path=T_path, Tu_path=Tu_path, G_path=ts.G_path, qb_path=ts.qb_path,
        qbstar_path=ts.bstar_path_r, h_path=h_path, π_path=ts.π_path,
        d_path=ts.d_path, covid_path=ts.covid_path, Γ_n=ts.Γ_n, Γ_c=ts.Γ_c,
        χ=c.χ, t=t, g=g, method=1, kappa_w=KAPPA_W, wY_ss=wY_ss,
        nkpc_scale=NKPC_SCALE, s_mode=s_mode)
    return (asset = maximum(abs, view(F, 1:n_r)),
            sblk  = maximum(abs, view(F, (n_r+1):(n_r+n_s))),
            inf   = maximum(abs, F))
end

println("\n", "="^100)
println("CHECK 1 -- does rule 2 unblock Stage 5?  :flex warm-start residual by block")
println("="^100)
A = resid_blocks(ts.τ_k_path, ts.T_path, ts.Tu_path, :flex)
B = resid_blocks(τk_cf,       ts.T_path, ts.Tu_path, :flex)
C = resid_blocks(τk_r2,       T_r2,      Tu_r2,      :flex)
D = resid_blocks(τk_cf,       T_r2,      Tu_r2,      :flex)
@printf("%-46s %14s %14s %14s\n", "evaluation (all at the current baseline x)", "asset", "s-block", "inf-norm")
@printf("%-46s %14.6e %14.6e %14.6e\n", "A  current paths, baseline tau_k",     A.asset, A.sblk, A.inf)
@printf("%-46s %14.6e %14.6e %14.6e\n", "B  current paths, NO-SUBSIDY tau_k",   B.asset, B.sblk, B.inf)
@printf("%-46s %14.6e %14.6e %14.6e\n", "C  rule-2 paths, rule-2 tau_k",        C.asset, C.sblk, C.inf)
@printf("%-46s %14.6e %14.6e %14.6e\n", "D  rule-2 paths, NO-SUBSIDY tau_k",    D.asset, D.sblk, D.inf)
@printf("\nthe tau_k shock effect on the ASSET block:\n")
@printf("  status quo : A %.6e -> B %.6e   (x%.2f)\n", A.asset, B.asset, B.asset/A.asset)
@printf("  rule 2     : C %.6e -> D %.6e   (x%.2f)\n", C.asset, D.asset, D.asset/C.asset)
@printf("\nD vs B, the like-for-like counterfactual residual: %.6e vs %.6e  -> %+.2f%%\n",
        D.asset, B.asset, 100*(D.asset/B.asset - 1))
println(D.asset < 0.9*B.asset ? "RULE 2 MATERIALLY IMPROVES the warm start." :
        D.asset > 1.1*B.asset ? "RULE 2 MAKES IT WORSE." :
        "RULE 2 LEAVES IT ESSENTIALLY UNCHANGED -- it does NOT unblock Stage 5.")

# ============================================================
# CHECK 2
# ============================================================
println("\n", "="^100)
println("CHECK 2 -- does retention of 0.638 at j=2 break anything MECHANICALLY?")
println("="^100)
θ_ = h_path[1].v.θ; ν_ = h_path[1].v.ν
atr_ss = labtax_ss/(e_init.w*e_init.n)
@printf("initial SS: tau_n = %.6f, aggregate average labour tax rate = %.4f%%\n",
        e_init.τ_n, 100*atr_ss)
@printf("Figure D1(d) reports the labour tax dipping to about 10%% in 2020Q2.\n\n")
@printf("%5s %10s %14s %16s %14s %12s\n",
        "j", "tau_n", "avg tax rate", "min GHH x", "bracket pos", "net inc > 0")
for j in 1:20
    atr = p.lab_tax[j]/(p.w[j]*p.n[j])
    xmin = Inf
    for i in axes(p[j].c_pol, 1), s in axes(p[j].c_pol, 2)
        nn = p[j].n_pol[i,s]
        xx = p[j].c_pol[i,s] - θ_*nn^(1 + 1/ν_)/(1 + 1/ν_)
        (e_init.ws.pdf[i,s] > 0 || true) && (xx < xmin) && (xmin = xx)
    end
    hi = j == 1 ? 5.0 : 2.0
    @printf("%5d %10.6f %13.4f%% %16.6e %13.1f%% %12s\n",
            j, p.τ_n[j], 100*atr, xmin, 100*p.τ_n[j]/hi, xmin > 0 ? "yes" : "NO")
end
@printf("\n[SW9] GHH floor (1e-9) fired %d times in this evaluation set.\n", _SW9_FLOOR[])
println("\nMECHANICAL VERDICT:")
println("  HSV net income is tau_n * y^(1-tau), positive for ANY tau_n > 0, so 0.638 is")
println("  interior and no expression is undefined. The root-find bracket for interior")
println("  periods is (0, 2) and 0.638 sits at 32% of it. The only object that can")
println("  actually break is the GHH composite x (log utility needs x > 0), which [SW9]")
println("  floors -- see the column above for whether it came close.")
println("\nDONE (nothing written)")
