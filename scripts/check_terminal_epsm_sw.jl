# ============================================================
# [SW24] Does a terminal-date-only epsm correction close the system?
#        READ-ONLY.  No solve.  Nothing written.
#
# The proposed specification:
#     i_taylor_ss = i_ss_base           (baseline, unchanged)
#     epsm[j]     = epsm_base[j]        j = 1..799, untouched
#     epsm[800]   = (1+r_final_cf)*(1+pi_bar) - 1 - i_ss_base
#
# Two things are checked:
#   (1) whether epsm[800] is READ by the recursion at all;
#   (2) what the recursion actually does at the tail -- growth or decay.
#
# The tail numbers here are EXACT for any run: in the tail r -> r_final_cf in every
# specification, so pi_j for large j does not depend on which run's r path is used.
# The j <= 400 numbers are indicative (a full solve would move r), but the tail is
# what the question turns on.
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
Jf(f) = joinpath(ROOT,"output","jld2_sw",f*".jld2")
PHI = 1.5
β_  = 0.9880431080390226

_i = load_object(Jf("steady_states")); e_init = _i[1]; e_fb = _i[2]
t = e_init.t; π_bar = e_init.π_
i_base = (1+e_fb.r)*(1+π_bar) - 1
_o = load_object(Jf("rule2_stage1"))
epsm_base = epsm_baseline_sw(copy(_o[3].r), copy(_o[3].π_), i_base, π_bar, PHI)
_o = nothing; GC.gc()
n = length(epsm_base)
r_min = -t.δ + eps(); r_max = 1/β_

println("="^104)
println("(0) WHICH epsm INDICES DOES THE RECURSION READ?")
println("="^104)
println("the loop is  `for j in (n_per-1):-1:2`  reading epsm_path[j]")
@printf("  n_per = %d  ->  epsm indices READ = 2 .. %d\n", n, n-1)
@printf("  epsm[1]   is NOT read (pi_1 is fed from the data)\n")
@printf("  epsm[%d] is NOT read  <-- the proposed correction sits here\n", n)

REPORT = [100,200,300,400,500,600,700,750,780,790,795,798,799,800]
for (nm, sstag, untag) in (("CBO debt","steady_states_cf_cbo_debt_rule2","cf_cbo_debt_rule2_forced"),
                           ("lower foreign","steady_states_cf_lower_foreign_rule2","cf_lower_foreign_rule2"))
    ef   = load_object(Jf(sstag))
    rf   = ef.r
    i_cf = (1+rf)*(1+π_bar) - 1
    d    = i_cf - i_base
    println("\n", "="^104)
    println("=== ", nm, " ===")
    println("="^104)
    @printf("r_final_cf = %.10f (%.4f%% ann)   i_ss_cf = %.10f (%.4f%% ann)\n",
            rf, 100*((1+rf)^4-1), i_cf, 100*((1+i_cf)^4-1))
    @printf("i_ss_base  = %.10f (%.4f%% ann)   d(i_ss) = %.12f (%+.4f pp ann simple)\n",
            i_base, 100*((1+i_base)^4-1), d, 400*d)

    # ---- (1) the terminal block, exactly as proposed ------------------------
    epsm_T = i_cf - i_base                     # = (1+r_f)(1+pi_bar)-1 - i_ss_base
    rule   = i_base + PHI*(π_bar - π_bar) + epsm_T
    fisher = (1+rf)*(1+π_bar) - 1
    println("\n---- (1) TERMINAL BLOCK at j=800: pi=pi_bar, r=r_final_cf ----")
    @printf("  epsm[800] set to %.14f (%+.4f pp ann)\n", epsm_T, 400*epsm_T)
    @printf("  rule   i = i_ss_base + phi*(pi_bar-pi_bar) + epsm[800] = %.16f\n", rule)
    @printf("  Fisher i = (1+r_final_cf)*(1+pi_bar) - 1              = %.16f\n", fisher)
    @printf("  difference = %.3e   -> CLOSES to machine precision: %s\n",
            rule-fisher, abs(rule-fisher) < 1e-16 ? "YES" : "no")
    println("  BUT epsm[800] is never read by the recursion (see (0)), so this")
    println("  identity has NO effect on the solved path.")

    # ---- (2) what the recursion actually does -------------------------------
    # analytical: x_j = pi_j - pi_bar,  x_800 = 0,  x_j = d/phi + ((1+r_f)/phi) x_{j+1}
    g = (1+rf)/PHI
    xfix = (d/PHI)/(1-g)
    println("\n---- (2) THE BACKWARD RECURSION, analytically ----")
    @printf("  x_j = pi_j - pi_bar   satisfies   x_j = d/phi + %.6f * x_{j+1},  x_800 = 0\n", g)
    @printf("  d/phi              = %.10f = %+.4f pp ann   <- your predicted x_799\n", d/PHI, 400*d/PHI)
    @printf("  FIXED POINT x*     = (d/phi)/(1-%.6f) = %.10f = %+.4f pp ann\n", g, xfix, 400*xfix)
    @printf("  d(i_ss)/(phi-1-r_f) = %+.4f pp ann   (same number, closed form)\n", 400*d/(PHI-1-rf))
    println("  the factor 0.669 is the rate of convergence TO x*, not decay toward 0.")
    println("  running BACKWARD from x_800 = 0, |x| GROWS monotonically to x*.")

    # numeric: run the actual recursion with the un-accommodated r path
    _u = load_object(Jf(untag)); pu = _u[3]
    rblk = collect(pu.r[2:end])                      # r_block[j] = r_{j+1}, sw19 layout
    πdat = copy(_u[4])
    _u = nothing; pu = nothing; GC.gc()
    epsm_u = copy(epsm_base); epsm_u[n] = epsm_T     # the proposed spec
    πv = Vector{Float64}(undef, n)
    _pi_from_r_recursion_sw!(πv, rblk; π_path = πdat, π_bar = π_bar,
        i_taylor_ss = i_base, phi_pi = PHI, epsm_path = epsm_u, n_per = n,
        r_min = r_min, r_max = r_max, sw19 = true)
    println("\n---- pi under the PROPOSED spec (i_ss = i_ss_base, epsm untouched on 2..799) ----")
    @printf("%6s %14s %16s %14s %16s\n","j","pi (ann %)","dev from target","analytic x_j","beta^j weight")
    for j in REPORT
        xa = xfix*(1 - g^(n-j))                       # closed form from x_800 = 0
        @printf("%6d %14.6f %+16.4f %+14.4f %16.3e\n", j,
                100*((1+πv[j])^4-1), 400*(πv[j]-π_bar), 400*xa, β_^j)
    end
    @printf("\nterminal pi = %.6f%% ann ; pi at j=400 sits %+.4f pp from target\n",
            100*((1+πv[n])^4-1), 400*(πv[400]-π_bar))
    @printf("periods for |x| to reach 99%% of x*: %d  (0.669^k < 0.01)\n",
            ceil(Int, log(0.01)/log(g)))

    # ---- (3) the impossibility, and the one relabelling that works ----------
    println("\n---- (3) CONSTANT COMPARISON (the recursion only ever sees i_ss + epsm_j) ----")
    @printf("%-42s %20s %14s\n","specification","constant at j=400","vs baseline")
    cb = i_base + epsm_base[400]
    for (lab, c_) in (("baseline", cb),
                      ("un-accommodated (i_ss_cf, epsm_base)", i_cf + epsm_base[400]),
                      ("PROPOSED (i_ss_base, epsm_base)", i_base + epsm_base[400]),
                      ("j<=600 accommodated", i_cf + epsm_base[400] - d),
                      ("RELABEL (i_ss_base, epsm_base + d)", i_base + epsm_base[400] + d))
        @printf("%-42s %20.14f %+14.3e\n", lab, c_, c_-cb)
    end
    @printf("\nPROPOSED minus j<=600-accommodated constant = %+.3e  (same system)\n",
            (i_base + epsm_base[400]) - (i_cf + epsm_base[400] - d))
    @printf("RELABEL minus un-accommodated constant      = %+.3e  (same system)\n",
            (i_base + epsm_base[400] + d) - (i_cf + epsm_base[400]))
    println("\nTERMINAL STEADY STATE REQUIREMENT.  At j=800 with pi=pi_bar and r=r_f the rule")
    println("gives i = i_ss + epsm_inf and Fisher gives (1+r_f)(1+pi_bar)-1, so")
    println("     i_ss + epsm_inf = (1+r_f)(1+pi_bar) - 1")
    println("whose RHS differs across columns whenever r_f does.  i_ss and epsm_inf are the")
    println("only two objects, so AT LEAST ONE must differ across columns.  Holding both")
    println("fixed while terminal pi = pi_bar is not available: one equation, no free slot.")
    @printf("     RHS: baseline %.14f   this column %.14f   gap %+.3e\n",
            i_base, i_cf, d)
end
println("\n", "="^104)
println("nothing solved, nothing written")

# ============================================================
# (4) THE FALLBACK: spread the correction over the last K periods instead.
#     epsm[j] = epsm_base[j] + d  for j in (800-K+1)..800, epsm_base elsewhere.
# ============================================================
println("\n", "="^104)
println("(4) FALLBACK: spreading the correction over the last K periods")
println("="^104)
for (nm, sstag, untag) in (("CBO debt","steady_states_cf_cbo_debt_rule2","cf_cbo_debt_rule2_forced"),
                           ("lower foreign","steady_states_cf_lower_foreign_rule2","cf_lower_foreign_rule2"))
    ef = load_object(Jf(sstag)); rf = ef.r
    d  = (1+rf)*(1+π_bar) - 1 - i_base
    _u = load_object(Jf(untag)); pu = _u[3]
    rblk = collect(pu.r[2:end]); πdat = copy(_u[4]); _u = nothing; pu = nothing; GC.gc()
    println("\n--- ", nm, " ---")
    @printf("%6s %16s %16s %16s %16s\n","K","pi[400] dev pp","pi[100] dev pp","pi[770] dev pp","pi[790] dev pp")
    for K in (1, 20, 30, 50, 100, 200, 400, 799)
        em = copy(epsm_base)
        for j in max(1, n-K+1):n; em[j] += d; end
        πv = Vector{Float64}(undef, n)
        _pi_from_r_recursion_sw!(πv, rblk; π_path = πdat, π_bar = π_bar, i_taylor_ss = i_base,
            phi_pi = PHI, epsm_path = em, n_per = n, r_min = r_min, r_max = r_max, sw19 = true)
        @printf("%6d %+16.4f %+16.4f %+16.4f %+16.4f\n", K,
                400*(πv[400]-π_bar), 400*(πv[100]-π_bar), 400*(πv[770]-π_bar), 400*(πv[790]-π_bar))
    end
    println("K = 799 restores the constant EVERYWHERE the recursion reads -> the un-accommodated run.")
end
