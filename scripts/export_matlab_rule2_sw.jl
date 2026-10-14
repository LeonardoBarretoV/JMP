# ============================================================
# [SW26] EXPORT THE RULE-2 BASELINE TO MATLAB
#        READ-ONLY on the .jld2 files.  No solve.
#
#   julia --project -t auto scripts/export_matlab_rule2_sw.jl
#
# Writes into output/matlab_rule2/  (a NEW folder -- output/matlab_sw/ holds the
# PRE-rule-2 baseline and is left untouched, so the old figures stay reproducible):
#
#   model_time_series_ratio.mat   same 21-column layout as main_JMP_sw.jl
#   model_new_series_sw.mat       s_index, piw, wedge, kappa_w
#                                 + pi_endog  (inflation with the shock switched off)
#                                 + epsm      (the monetary shock itself)
#   welfare_time_paths.mat        w_gains_path, same as the sw track
#
# THE TWO NEW SERIES.
#   epsm      is the Taylor-rule residual backed out from the rule-2 Stage-1 solve,
#             i.e. exactly the shock sequence the baseline was solved with.
#   pi_endog  re-runs the inflation recursion on the SOLVED baseline r path with
#             epsm set to zero throughout.  It answers "what inflation would this
#             rule deliver, given this real allocation, with no monetary shocks?".
#             It is a decomposition of the inflation equation holding r fixed, NOT
#             a re-solved general-equilibrium counterfactual -- with epsm = 0 the
#             real side would itself move.  Labelled accordingly in the figure.
#             j = 1 is fed from the data in BOTH series (pi_1 is not produced by
#             the recursion), so the two lines coincide at j = 1 by construction
#             and only separate from j = 2 on.
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
Jf(f) = joinpath(ROOT, "output", "jld2_sw", f * ".jld2")

TAG        = get(ENV, "SW_TAG", "rule2_epsm")
OUTDIR     = joinpath(ROOT, "output", get(ENV, "SW_OUTDIR", "matlab_rule2"))
KAPPA_W    = 0.05
PHI_PI     = 1.5
β_         = 0.9880431080390226

println("=== [SW26] exporting ", TAG, " -> ", OUTDIR, " ===")
mkpath(OUTDIR)

_i = load_object(Jf("steady_states")); e_init0, e_fb = _i[1], _i[2]
g = e_init0.g; t = e_init0.t
π_bar  = e_init0.π_
i_base = (1 + e_fb.r) * (1 + π_bar) - 1
r_min  = -t.δ + eps(); r_max = 1 / β_

_o = load_object(Jf("rule2_stage1"))
epsm = epsm_baseline_sw(copy(_o[3].r), copy(_o[3].π_), i_base, π_bar, PHI_PI)
_o = nothing; GC.gc()

_b = load_object(Jf(TAG))
e_init, e_final, path = _b[1], _b[2], _b[3]
π_solved = length(_b) >= 4 ? copy(_b[4]) : copy(path.π_)
n_per = length(path.r)
i_ss  = (1 + e_final.r) * (1 + π_bar) - 1
@printf("n_per = %d   i_taylor_ss = %.4f%% ann   max|epsm| = %.8f\n",
        n_per, 100*((1+i_ss)^4 - 1), maximum(abs, epsm))

# ---- consistency: path.pi_ must equal the solved pi vector ------------------
d_pi = maximum(abs, path.π_ .- π_solved)
@printf("check: max|path.pi_ - pi_solved| = %.3e  %s\n", d_pi,
        d_pi < 1e-12 ? "(path carries the SOLVED inflation, good)" : "*** MISMATCH ***")

# ---- s_index: recovered from the stored hours policy ------------------------
nt  = e_init.ws.n_pol[1, :]
emp = [z for z in eachindex(e_init.h.z_grid) if e_init.h.z_grid[z] != 0]
s_index = [sum(path[j].n_pol[1, z] / nt[z] for z in emp) / length(emp) for j in 1:n_per]
spread  = maximum(j -> (v = [path[j].n_pol[1,z]/nt[z] for z in emp];
                        maximum(v) - minimum(v)), 1:n_per)
@printf("s_index: s_1 = %.14f   max z-spread across all j = %.3e (0 => exact rationing)\n",
        s_index[1], spread)

# ---- piw: nominal wage inflation -------------------------------------------
piw = [j == 1 ? NaN : piw_at_sw(path, j, n_per, π_bar) for j in 1:n_per]

# ---- union labor wedge ------------------------------------------------------
hh  = e_init.h
θ_  = hh.v.θ; ν_ = hh.v.ν
wedge = Vector{Float64}(undef, n_per)
for j in 1:n_per
    DY = union_aggregates_(hh, g, path[j]; s = s_index[j], w = path.w[j],
                           τ_n = path.τ_n[j], n_ss = nt)
    wedge[j] = s_index[j] * (θ_ * s_index[j]^(1/ν_) * DY.D - path.w[j] * DY.Y)
end
@printf("wedge: j=1 %+.8e   j=2 %+.8e   max|.| over j<=40 %.3e\n",
        wedge[1], wedge[2], maximum(abs, wedge[1:40]))

# ---- pi with the monetary shock switched off --------------------------------
rblk = collect(path.r[2:end])
π_endog = Vector{Float64}(undef, n_per)
_pi_from_r_recursion_sw!(π_endog, rblk; π_path = π_solved, π_bar = π_bar,
    i_taylor_ss = i_ss, phi_pi = PHI_PI, epsm_path = zeros(n_per), n_per = n_per,
    r_min = r_min, r_max = r_max, sw19 = true)

# sanity: feeding the ACTUAL epsm back must reproduce the solved path exactly
π_chk = Vector{Float64}(undef, n_per)
_pi_from_r_recursion_sw!(π_chk, rblk; π_path = π_solved, π_bar = π_bar,
    i_taylor_ss = i_ss, phi_pi = PHI_PI, epsm_path = epsm, n_per = n_per,
    r_min = r_min, r_max = r_max, sw19 = true)
@printf("check: recursion with the true epsm reproduces the solved pi to %.3e\n",
        maximum(abs, π_chk .- π_solved))

println("\n---- inflation decomposition, annualised %% ----")
@printf("%6s %10s %14s %14s %14s\n", "j", "date", "pi solved", "pi, no shock", "epsm (pp ann)")
qlab(j) = (k = 2020*4 + 1 + (j-1); string(k ÷ 4, "Q", k % 4 + 1))
for j in vcat(collect(1:12), [16, 20, 24, 40])
    @printf("%6d %10s %14.4f %14.4f %+14.4f\n", j, qlab(j),
            100*((1+π_solved[j])^4 - 1), 100*((1+π_endog[j])^4 - 1), 400*epsm[j])
end
pk  = argmax(abs.(π_solved[2:end] .- π_bar)) + 1
@printf("\npeak solved pi   : %.4f%% ann at %s\n", 100*((1+π_solved[pk])^4-1), qlab(pk))
@printf("pi with no shock : %.4f%% ann at the same date\n", 100*((1+π_endog[pk])^4-1))
@printf("share of the peak deviation from target explained without shocks: %.2f%%\n",
        100*(π_endog[pk] - π_bar)/(π_solved[pk] - π_bar))

# ============================================================
# write
# ============================================================
mat1 = joinpath(OUTDIR, "model_time_series_ratio.mat")
matopen(mat1, "w") do f
    write(f, "time_series", [path.y path.π_ path.un path.q path.b path.bd path.bstar path.l path.lump_sum_T path.lump_sum_Tu path.r path.r_k path.k path.n path.w path.C path.lab_tax path.τ_k path.τ_n path.X path.G])
    write(f, "ss_init", [y(e_init) e_init.π_ e_init.h.Pss[1] e_init.q e_init.b e_init.bd e_init.bstar (e_init.b - ((1-g.ϑ)/(1+e_init.π_))*e_init.b) e_init.lump_sum_T e_init.lump_sum_Tu e_init.r e_init.r e_init.k e_init.n e_init.w (y(e_init) - t.δ*e_init.k - e_init.G) e_init.lab_tax e_init.τ_k e_init.τ_n (t.δ*e_init.k) e_init.G])
    write(f, "ss_final", [y(e_final) e_final.π_ e_final.h.Pss[1] e_final.q e_final.b e_final.bd e_final.bstar (e_final.b - ((1-g.ϑ)/(1+e_final.π_))*e_final.b) e_final.lump_sum_T e_final.lump_sum_Tu e_final.r e_final.r e_final.k e_final.n e_final.w (y(e_final) - t.δ*e_final.k - e_final.G) e_final.lab_tax e_final.τ_k e_final.τ_n (t.δ*e_final.k) e_final.G])
end
println("\nwrote ", mat1)

mat2 = joinpath(OUTDIR, "model_new_series_sw.mat")
matopen(mat2, "w") do f
    write(f, "s_index",  collect(s_index))
    write(f, "piw",      piw)
    write(f, "wedge",    wedge)
    write(f, "kappa_w",  KAPPA_W)
    write(f, "phi_pi",   PHI_PI)
    write(f, "pi_solved", collect(π_solved))
    write(f, "pi_endog", collect(π_endog))   # [SW26] shock switched off
    write(f, "epsm",     collect(epsm))      # [SW26] the shock itself
    write(f, "pi_bar",   π_bar)
    write(f, "i_taylor_ss", i_ss)
    write(f, "tag",      TAG)
end
println("wrote ", mat2)

mat3 = joinpath(OUTDIR, "welfare_time_paths.mat")
matopen(mat3, "w") do f
    write(f, "welfare_mat", w_gains_path(e_init, path))
end
println("wrote ", mat3)
println("\n=== done.  output/matlab_sw/ was NOT touched ===")
