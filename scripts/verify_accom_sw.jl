using LinearAlgebra, ProgressMeter, Roots, NLsolve, StructArrays, StatsBase, Tullio,
      Polyester, PrettyTables, CairoMakie, LaTeXStrings, DelimitedFiles, XLSX, CSV,
      DataFrames, ColorTypes, Printf, JLD2, MAT
const ROOT = normpath(joinpath(@__DIR__, ".."))
for f in ("types","preferences","technology_and_fiscal","calibration","auxiliary_functions",
          "household_problem_and_aggregation","stationary_equilibria","transition",
          "union_block","transition_sw","counterfactuals_sw","statistics","time_series")
    include(joinpath(ROOT, "src", f * ".jl"))
end
J(f) = joinpath(ROOT,"output","jld2_sw",f*".jld2")
PHI = 1.5
_i = load_object(J("steady_states")); e_init = _i[1]; e_fb = _i[2]
pib = e_init.π_
i_base = (1+e_fb.r)*(1+pib) - 1
e_cbo = load_object(J("steady_states_cf_cbo_debt_rule2"))
e_lf  = load_object(J("steady_states_cf_lower_foreign_rule2"))
_o = load_object(J("rule2_stage1"))
epsm = epsm_baseline_sw(copy(_o[3].r), copy(_o[3].π_),
                        i_base, pib, PHI)
_o = nothing; GC.gc()
n = length(epsm)
@printf("pi_bar = %.12f (%.4f%% ann)\n", pib, 100*((1+pib)^4-1))
@printf("epsm: max|.| = %.10f   tail: epsm[700] = %.3e  epsm[799] = %.3e  epsm[800] = %.3e\n",
        maximum(abs,epsm), epsm[700], epsm[799], epsm[n])
println()
for (nm, ef) in (("CBO debt", e_cbo), ("lower foreign", e_lf))
    i_cf = (1+ef.r)*(1+pib) - 1
    d    = i_cf - i_base
    rf   = ef.r
    @printf("=== %s ===\n", nm)
    @printf("  r_final = %.10f (%.4f%% ann)   i_ss_cf = %.10f (%.4f%% ann)\n",
            rf, 100*((1+rf)^4-1), i_cf, 100*((1+i_cf)^4-1))
    @printf("  d(i_ss) = %.12f = %+.4f pp ann\n", d, 400*d)
    # the recursion uses ONLY the sum (i_taylor_ss + epsm_j).  Compare the three specs.
    c_base  = i_base + epsm[400]                     # baseline
    c_unacc = i_cf   + epsm[400]                     # un-accommodated (i_ss rebuilt)
    c_accom = i_cf   + (epsm[400] - d)               # the requested accommodation
    c_old   = i_base + epsm[400]                     # the earlier divergent attempt
    @printf("  recursion constant (i_ss + epsm_j) at j=400:\n")
    @printf("     baseline                 %.14f\n", c_base)
    @printf("     un-accommodated          %.14f   (diff from baseline %+.3e)\n", c_unacc, c_unacc-c_base)
    @printf("     REQUESTED accommodation  %.14f   (diff from baseline %+.3e)\n", c_accom, c_accom-c_base)
    @printf("     earlier divergent try    %.14f   (diff from REQUESTED %+.3e)\n", c_old, c_old-c_accom)
    # terminal fixed point of pi_j = pi_bar + ((1+rf)(1+pi_{j+1}) - 1 - C)/phi with tail C
    fp(C) = pib + ( (1+rf)*(1+pib) - 1 - C ) / (PHI - 1 - rf)
    @printf("  terminal fixed point of the recursion (tail r = r_final_cf, tail epsm = 0):\n")
    @printf("     un-accommodated  pi* = %.12f = %.4f%% ann   (dev from pi_bar %+.4f pp ann)\n",
            fp(i_cf), 100*((1+fp(i_cf))^4-1), 400*(fp(i_cf)-pib))
    @printf("     REQUESTED        pi* = %.12f = %.4f%% ann   (dev from pi_bar %+.4f pp ann)\n",
            fp(i_cf-d), 100*((1+fp(i_cf-d))^4-1), 400*(fp(i_cf-d)-pib))
    @printf("     but the solver PINS pi[800] = pi_bar = %.4f%% ann  -> boundary conflict %+.4f pp ann\n",
            100*((1+pib)^4-1), 400*(fp(i_cf-d)-pib))
    @printf("  closed form d(i_ss)/(phi-1-r_f) = %+.4f pp ann\n", 400*d/(PHI-1-rf))
    # backward decay per period
    @printf("  backward contraction factor (1+r_f)/phi = %.6f ; ^100 = %.3e ; ^200 = %.3e\n",
            (1+rf)/PHI, ((1+rf)/PHI)^100, ((1+rf)/PHI)^200)
    println()
end
