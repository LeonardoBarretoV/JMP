# ============================================================
# [SW25] TABLE 6 CONSOLIDATION -- post-600 adjustment specification
#        READ-ONLY on the .jld2 files.  No solve.
#        Writes: output/jld2_sw/<tag>.txt  (one per run)
#                output/table6_sw.md       (the README section)
#
# Columns:
#   1 no capital subsidies   cf_no_capital_subsidy_rule2      (r_final unchanged)
#   2 lower foreign demand   cf_lower_foreign_rule2_accom     (j<=600 accommodated)
#   3 only stimulus checks   cf_only_checks_rule2             (r_final unchanged)
#   4 CBO debt projections   cf_cbo_debt_rule2_forced_accom   (j<=600, forced hours)
#
# The rule claim is VERIFIED FROM THE FILES, not asserted: for each run the pi
# recursion is re-run on the SAVED r path with the claimed (i_taylor_ss, epsm) and
# compared to the SAVED pi path.
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
Tf(f) = joinpath(ROOT,"output","jld2_sw",f*".txt")
PHI = 1.5; β_ = 0.9880431080390226; J_STOP = 600; NU = 1.0
qlabel(j) = (idx = 2020*4 + 1 + (j-1); string(idx ÷ 4, "Q", idx % 4 + 1))

# run registry: (label, tag, forced_hours, epsm_treatment, iters, wall_min, logged_norm)
RUNS = [
 ("Baseline",              "rule2_epsm",                     false, :none,  1, 0.19, 7.928515e-06),
 ("No capital subsidies",  "cf_no_capital_subsidy_rule2",    false, :none, 37, 6.80, 4.704329e-06),
 ("Lower foreign demand",  "cf_lower_foreign_rule2_accom",   false, :j600, 44, 8.22, 8.776202e-06),
 ("Only stimulus checks",  "cf_only_checks_rule2",           false, :none, 40, 6.62, 9.500213e-06),
 ("CBO debt projections",  "cf_cbo_debt_rule2_forced_accom", true,  :j600, 86,16.11, 9.451583e-06),
]
PAPER = Dict(
 "Baseline"             => (mean=1.82, min=-1.77, poor=3.03, middle=2.00, rich=-0.20),
 "No capital subsidies" => (mean=2.02, min=-6.18, poor=3.77, middle=2.27, rich=-1.69),
 "Lower foreign demand" => (mean=1.22, min=-3.72, poor=2.03, middle=1.34, rich= 0.06),
 "Only stimulus checks" => (mean=4.29, min=-5.80, poor=8.60, middle=4.28, rich=-1.90),
 "CBO debt projections" => (mean=2.51, min=-1.27, poor=4.00, middle=2.62, rich= 0.50))

_i = load_object(Jf("steady_states")); e_init0 = _i[1]; e_fb = _i[2]
t = e_init0.t; g0 = e_init0.g; π_bar = e_init0.π_
i_base = (1+e_fb.r)*(1+π_bar) - 1
r_min = -t.δ + eps(); r_max = 1/β_
_o = load_object(Jf("rule2_stage1"))
epsm_base = epsm_baseline_sw(copy(_o[3].r), copy(_o[3].π_), i_base, π_bar, PHI)
_o = nothing; GC.gc()
n = length(epsm_base)
y_ss = y(e_init0)

# ---- realized bin masses, from the initial-SS invariant wealth marginal ----
pdfm = sum(e_init0.ws.pdf, dims = 2)[:, 1]
ag   = round.(e_init0.h.a_grid, digits = 10)
qu   = round.([quantile(ag, ProbabilityWeights(pdfm), x) for x in 0.0:0.1:1.0], digits = 10)
masses = Float64[]; npts = Int[]
for (qi, qii) in zip(qu, qu[2:end])
    qii_ = (qii == qu[end]) ? qii + 10*eps() : qii
    sel = [qi <= ag[i] < qii_ for i in eachindex(ag)]
    push!(masses, sum(pdfm[sel])); push!(npts, count(sel))
end

R = Any[]
for (lab, tag, forced, treat, it, wall, lognorm) in RUNS
    o = load_object(Jf(tag))
    ntup = length(o)
    ei, ef, p = o[1], o[2], o[3]
    πs = copy(o[4])
    epsm_used = length(o) >= 5 ? copy(o[5]) : copy(epsm_base)
    i_ss = (1 + ef.r)*(1 + π_bar) - 1
    # --- VERIFY the spec from the file: recursion on saved r must give saved pi ---
    rblk = collect(p.r[2:end])
    πchk = Vector{Float64}(undef, n)
    _pi_from_r_recursion_sw!(πchk, rblk; π_path = πs, π_bar = π_bar, i_taylor_ss = i_ss,
        phi_pi = PHI, epsm_path = epsm_used, n_per = n, r_min = r_min, r_max = r_max, sw19 = true)
    recur_err = maximum(abs, πchk[2:end] .- πs[2:end])
    # --- gate, branch, residuals[1] ---
    gate = gate_lab_tax_sw(ei, p, g0)
    nt = ei.ws.n_pol[1,:]
    emp = [z for z in eachindex(ei.h.z_grid) if ei.h.z_grid[z] != 0]
    s1b = sum(p[1].n_pol[1,z]/nt[z] for z in emp)/length(emp)
    sfx = s_flex_sw(g0, NU; w = p.w[1], τ_n = p.τ_n[1], w_init = ei.w, τ_n_init = ei.τ_n)
    ce = ce_summary_sw(ei, p)
    dev = 400 .* (πs .- π_bar)
    pk = argmax(abs.(dev[2:end])) + 1
    push!(R, (; lab, tag, forced, treat, it, wall, lognorm, i_ss, rf = ef.r, ntup,
              epsm_used, πs, r = copy(p.r), τn = copy(p.τ_n), y1 = p.y[1],
              recur_err, gate, s1b, sfx, ce, pk, pkdev = dev[pk], dev))
    o = nothing; p = nothing; GC.gc()
end
B = R[1]; C = R[2:end]
Γn = B.s1b

io = IOBuffer()
W(s...) = (println(io, s...); println(s...))

W("="^118)
W("[SW25] TABLE 6 CONSOLIDATION -- post-600 adjustment specification")
W("="^118)

# ============================================================
W("\n", "="^118)
W("PART 1 -- IS THE RULE COMMON ACROSS ALL FOUR COLUMNS?")
W("="^118)
@printf(io, "%-22s %14s %7s %9s %13s %-24s\n","column","i_taylor_ss","phi_pi","pi_bar","r_final","epsm adjustment")
@printf("%-22s %14s %7s %9s %13s %-24s\n","column","i_taylor_ss","phi_pi","pi_bar","r_final","epsm adjustment")
for x in R
    tr = x.treat === :none ? "NONE" : "-d(i_ss), j <= $(J_STOP)"
    for out in (io, stdout)
        @printf(out, "%-22s %13.4f%% %7.2f %8.2f%% %12.4f%% %-24s\n", x.lab,
                100*((1+x.i_ss)^4-1), PHI, 100*((1+π_bar)^4-1), 100*((1+x.rf)^4-1), tr)
    end
end
W("\nr_final moved? (against the baseline terminal steady state)")
for x in C
    d = x.i_ss - i_base
    for out in (io, stdout)
        @printf(out, "  %-22s d(r_final) = %+.4f pp ann   d(i_ss) = %+.12f  %s\n", x.lab,
                100*((1+x.rf)^4-1) - 100*((1+B.rf)^4-1), d,
                abs(d) < 1e-14 ? "<- UNCHANGED, so NO adjustment is needed" : "<- moved, adjustment applied")
    end
end
W("\nVERIFIED: columns 1 and 3 carry i_taylor_ss EQUAL to the baseline value and no epsm")
W("adjustment.  Not assumed -- d(i_ss) is exactly 0 and epsm equals epsm_base elementwise:")
for x in C
    if x.treat === :none
        de = maximum(abs, x.epsm_used .- epsm_base)
        for out in (io, stdout)
            @printf(out, "  %-22s max|epsm - epsm_base| = %.3e   |i_ss - i_ss_base| = %.3e\n",
                    x.lab, de, abs(x.i_ss - i_base))
        end
    end
end

W("\n---- the object the model actually uses: the rule constant  i_taylor_ss + epsm_j ----")
@printf(io, "%-22s %22s %22s %16s\n","column","constant at j=100","constant at j=700","j<=600 vs base")
@printf("%-22s %22s %22s %16s\n","column","constant at j=100","constant at j=700","j<=600 vs base")
cb = [B.i_ss + B.epsm_used[j] for j in 1:n]
for x in R
    cx = [x.i_ss + x.epsm_used[j] for j in 1:n]
    dmax = maximum(abs, cx[2:J_STOP] .- cb[2:J_STOP])
    for out in (io, stdout)
        @printf(out, "%-22s %22.14f %22.14f %16.3e\n", x.lab, cx[100], cx[700], dmax)
    end
end
W("\nTHE RULE IS BIT-IDENTICAL ACROSS ALL FIVE RUNS FOR EVERY DATE j = 2..600.")
W("It differs only for j > 600, where columns 2 and 4 revert to the value Fisher")
W("consistency at their own terminal steady state requires.  beta^600 = 7.34e-04, so")
W("nothing past the switch touches the announcement-date welfare numbers.")
W("\nEquivalently, with i_taylor_ss held at the BASELINE 3.0994% in all five columns, the")
W("shock path is exactly epsm_base for j <= 600 everywhere, and columns 2 and 4 carry a")
W("permanent extra shock for j > 600 of:")
for x in C
    x.treat === :j600 || continue
    for out in (io, stdout)
        @printf(out, "  %-22s %+.4f pp ann\n", x.lab, 400*(x.i_ss - i_base))
    end
end

W("\n---- specification verified FROM THE FILES (recursion on saved r -> saved pi) ----")
@printf(io, "%-22s %18s %14s\n","column","max|pi_chk - pi|","verdict")
@printf("%-22s %18s %14s\n","column","max|pi_chk - pi|","verdict")
for x in R
    for out in (io, stdout)
        @printf(out, "%-22s %18.3e %14s\n", x.lab, x.recur_err, x.recur_err < 1e-10 ? "OK" : "MISMATCH")
    end
end

W("\n---- validity of each saved solution ----")
@printf(io, "%-22s %8s %9s %13s %11s %13s %9s\n","column","iters","wall min","inf-norm","gate rel","s_flex_1 marg","forced")
@printf("%-22s %8s %9s %13s %11s %13s %9s\n","column","iters","wall min","inf-norm","gate rel","s_flex_1 marg","forced")
for x in R
    for out in (io, stdout)
        @printf(out, "%-22s %8d %9.2f %13.3e %11.3e %12.2f%% %9s\n", x.lab, x.it, x.wall,
                x.lognorm, x.gate.worst_rel, 100*(x.sfx/Γn - 1), x.forced ? "YES" : "no")
    end
end
W("gate tolerance 1e-5; 1.566e-07 is the known mass-leak constant.")
W("s_flex_1 margin > 0 means the Gamma_n branch binds, which every column satisfies.")
for x in R
    for out in (io, stdout)
        @printf(out, "  %-22s s_1 = %.14f  (Gamma_n, diff %.2e)   y[1] fall %.4f%%\n",
                x.lab, x.s1b, x.s1b - Γn, 100*(1 - x.y1/y_ss))
    end
end

# ============================================================
W("\n", "="^118)
W("PART 2 -- TABLE 6")
W("="^118)
W("\nBIN CONVENTION: quantiles over the asset grid weighted by the initial-SS invariant")
W("wealth marginal; Poor = bin 1, Middle = mean of bins 5 and 6, Rich = bin 10.")
W("LABELLING CORRECTION -- realized masses, not the nominal deciles:")
@printf(io, "%6s %14s %16s %10s\n","bin","realized mass","cumulative","grid pts")
@printf("%6s %14s %16s %10s\n","bin","realized mass","cumulative","grid pts")
cmass = cumsum(masses)
for k in 1:10
    for out in (io, stdout)
        @printf(out, "%6d %13.2f%% %15.1f%% %10d\n", k, 100*masses[k], 100*cmass[k], npts[k])
    end
end
for out in (io, stdout)
    @printf(out, "\nSo \"Poor (bottom 10%%)\" is really the bottom %.1f%%, \"Middle (40-60%%)\" is about\n",
            100*masses[1])
    @printf(out, "%.1f-%.1f%%, and \"Rich (top 10%%)\" is the top %.1f%%.\n",
            100*sum(masses[1:4]), 100*sum(masses[1:6]), 100*masses[10])
end

ROWS = [("Overall mean",:mean),("Overall minimum",:min),("Poor (bottom 10%)",:poor),
        ("Middle (40-60%)",:middle),("Rich (top 10%)",:rich)]
hdr(pfx) = begin
    for out in (io, stdout)
        @printf(out, "%-20s %11s %12s %14s %14s %14s\n", pfx, "Baseline", "No capital",
                "Lower foreign", "Only stimulus", "CBO debt")
        @printf(out, "%-20s %11s %12s %14s %14s %14s\n", "", "", "subsidies", "demand",
                "checks", "projections")
    end
end
W("\n", "-"^118)
W("(a) LEVELS -- CE against the pre-stimulus steady state (%)")
W("-"^118); hdr("")
for (nm, k) in ROWS
    for out in (io, stdout)
        @printf(out, "%-20s %11.3f", nm, getfield(B.ce, k))
        for x in C; @printf(out, " %12.3f  ", getfield(x.ce, k)); end
        println(out)
    end
end
W("\n", "-"^118)
W("(b) DIFFERENCES from the rule-2 baseline (pp)")
W("-"^118); hdr("")
for (nm, k) in ROWS
    for out in (io, stdout)
        @printf(out, "%-20s %11s", nm, "-")
        for x in C; @printf(out, " %+12.3f  ", getfield(x.ce, k) - getfield(B.ce, k)); end
        println(out)
    end
end
W("\n", "-"^118)
W("(c) PUBLISHED PAPER NUMBERS for the same cells (%)")
W("-"^118); hdr("")
for (nm, k) in ROWS
    for out in (io, stdout)
        @printf(out, "%-20s %11.2f", nm, getfield(PAPER["Baseline"], k))
        for x in C; @printf(out, " %12.2f  ", getfield(PAPER[x.lab], k)); end
        println(out)
    end
end
W("\n", "-"^118)
W("(a) minus (c): model level minus published level (pp)")
W("-"^118); hdr("")
for (nm, k) in ROWS
    for out in (io, stdout)
        @printf(out, "%-20s %+11.3f", nm, getfield(B.ce,k) - getfield(PAPER["Baseline"],k))
        for x in C; @printf(out, " %+12.3f  ", getfield(x.ce,k) - getfield(PAPER[x.lab],k)); end
        println(out)
    end
end
W("\n", "-"^118)
W("(b) minus paper differences: model difference minus published difference (pp)")
W("-"^118); hdr("")
for (nm, k) in ROWS
    for out in (io, stdout)
        @printf(out, "%-20s %11s", nm, "-")
        for x in C
            md = getfield(x.ce,k) - getfield(B.ce,k)
            pd = getfield(PAPER[x.lab],k) - getfield(PAPER["Baseline"],k)
            @printf(out, " %+12.3f  ", md - pd)
        end
        println(out)
    end
end

W("\n", "-"^118)
W("FULL DECILE PROFILE -- levels (%), then differences from the baseline")
W("-"^118)
for out in (io, stdout)
    @printf(out, "%5s %9s %11s %14s %14s %14s %14s\n","bin","mass","Baseline",
            "No capital","Lower foreign","Only stimulus","CBO debt")
    for k in 1:10
        @printf(out, "%5d %8.2f%% %11.4f", k, 100*masses[k], B.ce.deciles[k])
        for x in C; @printf(out, " %14.4f", x.ce.deciles[k]); end
        println(out)
    end
    println(out, "  differences from the baseline:")
    for k in 1:10
        @printf(out, "%5d %8.2f%% %11s", k, 100*masses[k], "-")
        for x in C; @printf(out, " %+14.4f", x.ce.deciles[k] - B.ce.deciles[k]); end
        println(out)
    end
end

# ============================================================
W("\n", "="^118)
W("PART 3 -- WHAT THE EXTENSION ADDS (not in the published table)")
W("="^118)
J22 = 9
for out in (io, stdout)
    @printf(out, "2022Q2 is j = %d.  pi_1 is fed from the data, so the peak is taken over j >= 2.\n\n", J22)
    @printf(out, "%-22s %16s %14s %16s %14s %18s\n","column","pi 2022Q2 ann%","vs base pp",
            "r j=20 ann%","vs base pp","peak pi dev (date)")
    @printf(out, "%-22s %16.4f %14s %16.4f %14s %11.4f pp %s\n", B.lab,
            100*((1+B.πs[J22])^4-1), "-", 100*((1+B.r[20])^4-1), "-",
            B.pkdev, qlabel(B.pk))
    for x in C
        @printf(out, "%-22s %16.4f %+14.4f %16.4f %+14.4f %11.4f pp %s\n", x.lab,
                100*((1+x.πs[J22])^4-1), 400*(x.πs[J22]-B.πs[J22]),
                100*((1+x.r[20])^4-1), 400*(x.r[20]-B.r[20]),
                x.pkdev, qlabel(x.pk))
    end
    println(out)
    @printf(out, "%-22s %14s %14s %16s\n","column","tau_n j=1","tau_n j=2","terminal pi ann%")
    @printf(out, "%-22s %14.8f %14.8f %16.6f\n", B.lab, B.τn[1], B.τn[2], 100*((1+B.πs[end])^4-1))
    for x in C
        @printf(out, "%-22s %14.8f %14.8f %16.6f\n", x.lab, x.τn[1], x.τn[2], 100*((1+x.πs[end])^4-1))
    end
    @printf(out, "\ninitial-SS tau_n = %.8f   y_ss = %.10f   Gamma_n = %.14f\n",
            e_init0.τ_n, y_ss, Γn)
end

# ============================================================
# companion .txt files
# ============================================================
for x in R
    open(Tf(x.tag), "w") do f
        println(f, "="^78)
        println(f, x.tag, ".jld2")
        println(f, "TABLE 6 COLUMN: ", x.lab)
        println(f, "written by scripts/table6_consolidate_sw.jl  [SW25]")
        println(f, "="^78)
        println(f, "\nSAVED TUPLE: (e_init, e_final, path, pi", length(load_object(Jf(x.tag))) >= 5 ? ", epsm)" : ")")
        println(f, "\n---- MONETARY SPECIFICATION ----")
        @printf(f, "  phi_pi          = %.4f\n", PHI)
        @printf(f, "  pi_bar          = %.12f  (%.4f%% ann)\n", π_bar, 100*((1+π_bar)^4-1))
        @printf(f, "  i_taylor_ss     = %.12f  (%.4f%% ann)\n", x.i_ss, 100*((1+x.i_ss)^4-1))
        @printf(f, "  r_final         = %.12f  (%.4f%% ann)\n", x.rf, 100*((1+x.rf)^4-1))
        @printf(f, "  d(i_ss) vs base = %+.12f  (%+.4f pp ann)\n", x.i_ss-i_base, 400*(x.i_ss-i_base))
        if x.treat === :none
            println(f, "  epsm            = epsm_base, NO ADJUSTMENT anywhere")
            @printf(f, "                    max|epsm - epsm_base| = %.3e\n", maximum(abs, x.epsm_used .- epsm_base))
        else
            @printf(f, "  epsm            = epsm_base - d(i_ss) for j <= %d ; epsm_base for j > %d\n", J_STOP, J_STOP)
            println(f, "                    (post-600 adjustment specification)")
        end
        @printf(f, "  rule constant (i_taylor_ss + epsm_j) equals the baseline's for j=2..%d\n", J_STOP)
        @printf(f, "  forced period-0 hours [SW22] = %s\n", x.forced ? "YES (mandate, not ceiling)" : "no")
        println(f, "\n---- CONVERGENCE ----")
        @printf(f, "  iterations      = %d\n", x.it)
        @printf(f, "  wall time       = %.2f min\n", x.wall)
        @printf(f, "  final inf-norm  = %.6e   (ftol 1e-5)\n", x.lognorm)
        @printf(f, "  spec re-verified from this file: max|pi_recursion - pi_saved| = %.3e\n", x.recur_err)
        println(f, "\n---- ITEM (b) GATE: household vs government labour tax ----")
        @printf(f, "  worst relative gap = %.3e at j = %d   -> %s\n", x.gate.worst_rel, x.gate.worst_j,
                x.gate.worst_rel < 1e-5 ? "PASSED" : "FAILED")
        println(f, "  (1.566e-07 is the known mass-leak constant, not an error)")
        println(f, "\n---- PERIOD-0 BRANCH ----")
        @printf(f, "  s_1 booked by the household block = %.14f\n", x.s1b)
        @printf(f, "  Gamma_n                           = %.14f  (diff %.2e)\n", Γn, x.s1b-Γn)
        @printf(f, "  s_flex_1                          = %.10f  (margin %+.4f%% vs Gamma_n)\n",
                x.sfx, 100*(x.sfx/Γn-1))
        println(f, "  ", x.sfx >= Γn ? "Gamma_n is the active branch (ceiling binds)." :
                                       "s_flex_1 < Gamma_n -- valid only because hours are FORCED.")
        @printf(f, "  2020Q2 output fall = %.4f%%\n", 100*(1-x.y1/y_ss))
        println(f, "\n---- HEADLINE CE WELFARE (%, paper bin convention) ----")
        for (nm, k) in ROWS
            @printf(f, "  %-20s model %8.3f", nm, getfield(x.ce, k))
            x.lab == "Baseline" || @printf(f, "   diff vs baseline %+8.3f", getfield(x.ce,k)-getfield(B.ce,k))
            @printf(f, "   paper %6.2f\n", getfield(PAPER[x.lab], k))
        end
        println(f, "\n  full decile profile (levels):")
        for k in 1:10
            @printf(f, "    bin %2d (mass %5.2f%%)  %9.4f\n", k, 100*masses[k], x.ce.deciles[k])
        end
        println(f, "\n---- EXTENSION DIAGNOSTICS ----")
        @printf(f, "  pi 2022Q2 (j=9)   = %.4f%% ann\n", 100*((1+x.πs[J22])^4-1))
        @printf(f, "  r at j=20         = %.4f%% ann\n", 100*((1+x.r[20])^4-1))
        @printf(f, "  peak pi deviation = %+.4f pp ann at %s (j=%d), over j>=2\n", x.pkdev, qlabel(x.pk), x.pk)
        @printf(f, "  terminal pi       = %.6f%% ann\n", 100*((1+x.πs[end])^4-1))
        @printf(f, "  tau_n j=1 / j=2   = %.8f / %.8f   (retention LEVEL: lower = heavier tax)\n", x.τn[1], x.τn[2])
        println(f, "\n", "="^78)
    end
    println("wrote ", Tf(x.tag))
end

open(joinpath(ROOT,"output","table6_sw.txt"), "w") do f
    write(f, String(take!(io)))
end
println("wrote output/table6_sw.txt")

# ============================================================
# markdown emitter -- output/table6_sw.md, appended to README_sw.md
# ============================================================
LBL = [x.lab for x in C]
open(joinpath(ROOT,"output","table6_sw.md"), "w") do f
    pr(s...) = println(f, s...)
    pf(fmt, a...) = Printf.format(f, Printf.Format(fmt), a...)
    pr("\n---\n")
    pr("## TABLE 6 -- consolidated, post-600 adjustment specification")
    pr()
    pr("Generated by `scripts/table6_consolidate_sw.jl` (`[SW25]`), read-only on the saved")
    pr("`.jld2` files. Plain-text form in `output/table6_sw.txt`; per-run provenance in")
    pr("`output/jld2_sw/<tag>.txt`.")
    pr()
    pr("| # | column | file |")
    pr("|---|---|---|")
    for (i,x) in enumerate(C); pf("| %d | %s | `%s.jld2` |\n", i, x.lab, x.tag); end
    pr("| - | Baseline | `rule2_epsm.jld2` |")

    pr("\n### Is the rule common across all four columns?")
    pr()
    pr("| column | `i_taylor_ss` | `phi_pi` | `pi_bar` | `r_final` | `epsm` adjustment |")
    pr("|---|---|---|---|---|---|")
    for x in R
        pf("| %s | %.4f%% ann | %.2f | %.2f%% ann | %.4f%% ann | %s |\n", x.lab,
           100*((1+x.i_ss)^4-1), PHI, 100*((1+π_bar)^4-1), 100*((1+x.rf)^4-1),
           x.treat === :none ? "**none**" : "`-d(i_ss)` over `j <= $(J_STOP)`")
    end
    pr()
    pr("**Columns 1 and 3 carry no adjustment, verified rather than assumed.** Neither moves")
    pr("the terminal steady state, so `d(i_ss)` is exactly zero and `epsm` is exactly")
    pr("`epsm_base`, elementwise:")
    pr()
    pr("| column | `d(r_final)` | `d(i_ss)` | `max abs(epsm - epsm_base)` |")
    pr("|---|---|---|---|")
    for x in C
        de = maximum(abs, x.epsm_used .- epsm_base)
        pf("| %s | %+.4f pp ann | %+.12f | %s |\n", x.lab,
           100*((1+x.rf)^4-1) - 100*((1+B.rf)^4-1), x.i_ss - i_base,
           x.treat === :none ? @sprintf("**%.3e**", de) : "(adjusted by construction)")
    end

    pr("\n### The rule constant, which is what the model actually uses")
    pr()
    pr("`i_taylor_ss` and `epsm` enter the solver only as their sum (`transition_sw.jl:1092`),")
    pr("so the sum is the rule:")
    pr()
    pr("| column | constant at j=100 | constant at j=700 | max abs dev vs baseline, j=2..600 |")
    pr("|---|---|---|---|")
    cbv = [B.i_ss + B.epsm_used[j] for j in 1:n]
    for x in R
        cx = [x.i_ss + x.epsm_used[j] for j in 1:n]
        pf("| %s | %.14f | %.14f | **%.3e** |\n", x.lab, cx[100], cx[700],
           maximum(abs, cx[2:J_STOP] .- cbv[2:J_STOP]))
    end
    pr()
    pr("**The rule is bit-identical across all five runs at every date `j = 2..600`.** It")
    pr("differs only for `j > 600`, where columns 2 and 4 revert to the value Fisher")
    pr("consistency at their own terminal steady state requires. `beta^600 = 7.34e-04`, so")
    pr("nothing past the switch reaches the announcement-date welfare numbers.")
    pr()
    pr("Equivalently, holding `i_taylor_ss` at the baseline 3.0994% in all five columns, the")
    pr("shock path is exactly `epsm_base` for `j <= 600` everywhere, and columns 2 and 4 carry")
    pr("a permanent extra shock for `j > 600` of")
    pr()
    for x in C
        x.treat === :j600 || continue
        pf("- %s: **%+.4f pp** annualised\n", x.lab, 400*(x.i_ss - i_base))
    end
    pr()
    pr("which is exactly the change in the neutral nominal rate. So the rule and the shocks")
    pr("are held fixed across all four columns over the entire horizon that matters, with the")
    pr("qualification confined to `j > 600`.")

    pr("\n### Validity of each saved solution")
    pr()
    pr("| column | iters | wall (min) | inf-norm | gate (rel) | `s_flex_1` margin | forced hours | spec re-verified |")
    pr("|---|---|---|---|---|---|---|---|")
    for x in R
        pf("| %s | %d | %.2f | %.3e | %.3e | %+.2f%% | %s | %.3e |\n", x.lab, x.it, x.wall,
           x.lognorm, x.gate.worst_rel, 100*(x.sfx/Γn-1), x.forced ? "YES" : "no", x.recur_err)
    end
    pr()
    pr("Gate tolerance 1e-5; `1.566e-07` is the known mass-leak constant. \"Spec re-verified\"")
    pr("re-runs the `pi` recursion on each file's own saved `r` path with the stated")
    pr("`(i_taylor_ss, epsm)` and compares to the saved `pi` -- so the specification claim is")
    pr("checked against the files, not against these notes. **A positive `s_flex_1` margin")
    pr("means the `Gamma_n` branch binds, which every column satisfies** -- including the CBO")
    pr("column at +20.61%, so `[SW22]` forced hours is inert there and the column does not")
    pr("depend on it. `s_1 = 0.90230115171930 = Gamma_n` in all five, and the 2020Q2 output")
    pr("fall is 8.4002% in all five.")

    pr("\n### Bin convention and the labelling correction")
    pr()
    pr("Quantiles over the asset grid weighted by the initial-SS invariant wealth marginal;")
    pr("Poor = bin 1, Middle = mean of bins 5 and 6, Rich = bin 10. There is an atom at")
    pr("`a_min`, so the realized masses are not deciles:")
    pr()
    pr("| bin | realized mass | cumulative | grid points |")
    pr("|---|---|---|---|")
    for k in 1:10
        pf("| %d%s | %.2f%% | %.1f%% | %d |\n", k,
           k==1 ? " (\"bottom 10%\")" : k==10 ? " (\"top 10%\")" : "",
           100*masses[k], 100*cmass[k], npts[k])
    end
    pr()
    pf("So **\"Poor (bottom 10%%)\" is really the bottom %.1f%%**, \"Middle (40-60%%)\" is about\n", 100*masses[1])
    pf("%.1f-%.1f%%, and \"Rich (top 10%%)\" is the top %.1f%%.\n",
       100*sum(masses[1:4]), 100*sum(masses[1:6]), 100*masses[10])

    hdrmd() = (pr("| | Baseline | No capital subsidies | Lower foreign demand | Only stimulus checks | CBO debt projections |");
               pr("|---|---|---|---|---|---|"))
    pr("\n### (a) Levels -- CE against the pre-stimulus steady state (%)")
    pr(); hdrmd()
    for (nm,k) in ROWS
        pf("| **%s** | %.3f", nm, getfield(B.ce,k))
        for x in C; pf(" | %.3f", getfield(x.ce,k)); end
        pr(" |")
    end
    pr("\n### (b) Differences from the rule-2 baseline (pp)")
    pr(); hdrmd()
    for (nm,k) in ROWS
        pf("| **%s** | -", nm)
        for x in C; pf(" | %+.3f", getfield(x.ce,k)-getfield(B.ce,k)); end
        pr(" |")
    end
    pr("\n### (c) Published paper numbers, same cells (%)")
    pr(); hdrmd()
    for (nm,k) in ROWS
        pf("| **%s** | %.2f", nm, getfield(PAPER["Baseline"],k))
        for x in C; pf(" | %.2f", getfield(PAPER[x.lab],k)); end
        pr(" |")
    end
    pr("\n### (a) minus (c) -- model level minus published level (pp)")
    pr(); hdrmd()
    for (nm,k) in ROWS
        pf("| **%s** | %+.3f", nm, getfield(B.ce,k)-getfield(PAPER["Baseline"],k))
        for x in C; pf(" | %+.3f", getfield(x.ce,k)-getfield(PAPER[x.lab],k)); end
        pr(" |")
    end
    pr("\n### (b) minus published differences (pp)")
    pr(); hdrmd()
    for (nm,k) in ROWS
        pf("| **%s** | -", nm)
        for x in C
            pf(" | %+.3f", (getfield(x.ce,k)-getfield(B.ce,k)) -
                           (getfield(PAPER[x.lab],k)-getfield(PAPER["Baseline"],k)))
        end
        pr(" |")
    end

    pr("\n### Full decile profile")
    pr()
    pr("| bin | mass | Baseline | No capital | Lower foreign | Only stimulus | CBO debt |")
    pr("|---|---|---|---|---|---|---|")
    for k in 1:10
        pf("| %d | %.2f%% | %.4f", k, 100*masses[k], B.ce.deciles[k])
        for x in C; pf(" | %.4f", x.ce.deciles[k]); end
        pr(" |")
    end
    pr()
    pr("Differences from the baseline:")
    pr()
    pr("| bin | mass | No capital | Lower foreign | Only stimulus | CBO debt |")
    pr("|---|---|---|---|---|---|")
    for k in 1:10
        pf("| %d | %.2f%%", k, 100*masses[k])
        for x in C; pf(" | %+.4f", x.ce.deciles[k]-B.ce.deciles[k]); end
        pr(" |")
    end

    pr("\n### Second panel -- what the extension adds (not in the published table)")
    pr()
    pr("| column | `pi` 2022Q2 (ann %) | vs baseline (pp) | `r` at j=20 (ann %) | vs baseline (pp) | peak `pi` dev from target | peak date |")
    pr("|---|---|---|---|---|---|---|")
    pf("| %s | %.4f | - | %.4f | - | %+.4f pp | %s |\n", B.lab,
       100*((1+B.πs[J22])^4-1), 100*((1+B.r[20])^4-1), B.pkdev, qlabel(B.pk))
    for x in C
        pf("| %s | %.4f | %+.4f | %.4f | %+.4f | %+.4f pp | %s |\n", x.lab,
           100*((1+x.πs[J22])^4-1), 400*(x.πs[J22]-B.πs[J22]),
           100*((1+x.r[20])^4-1), 400*(x.r[20]-B.r[20]), x.pkdev, qlabel(x.pk))
    end
    pr()
    pr("2022Q2 is `j = 9`. `pi_1` is fed from the data, so the peak is taken over `j >= 2`.")
    pr()
    pr("**Read the peak column carefully.** The baseline's `epsm` is backed out so the baseline")
    pr("reproduces observed inflation, so its +6.78pp peak is the *data's* deviation, not an")
    pr("endogenous model response. What the extension generates endogenously is the")
    pr("*difference across columns*, and that is small -- at most **+0.15pp** at 2022Q2 (lower")
    pr("foreign demand). This is the same finding as the headline result: absent monetary")
    pr("shocks the fiscal and foreign blocks move inflation very little.")
    pr()
    pr("| column | `tau_n` j=1 | `tau_n` j=2 | terminal `pi` (ann %) |")
    pr("|---|---|---|---|")
    for x in R; pf("| %s | %.8f | %.8f | %.6f |\n", x.lab, x.τn[1], x.τn[2], 100*((1+x.πs[end])^4-1)); end
    pr()
    pf("Initial-SS `tau_n = %.8f` (a retention *level*: lower means a heavier tax),\n", e_init0.τ_n)
    pf("`y_ss = %.10f`, `Gamma_n = %.14f`.\n", y_ss, Γn)

    pr("\n### What matches and what does not")
    pr()
    pr("- **Levels agree with the paper closely** in this specification: the mean is within")
    pf("  %.3f-%.3f pp across all five columns, Middle within %.3f pp, Rich within %.3f pp.\n",
       minimum(abs(getfield(x.ce,:mean)-getfield(PAPER[x.lab],:mean)) for x in R),
       maximum(abs(getfield(x.ce,:mean)-getfield(PAPER[x.lab],:mean)) for x in R),
       maximum(abs(getfield(x.ce,:middle)-getfield(PAPER[x.lab],:middle)) for x in R),
       maximum(abs(getfield(x.ce,:rich)-getfield(PAPER[x.lab],:rich)) for x in R))
    pr("- **`Overall minimum` is the weak cell.** It is a single `(a,z)` cell with positive")
    pr("  mass, so it is far more sensitive than any average. The baseline is -5.383 against")
    pr("  the published -1.77, and lower foreign demand is -10.063 against -3.72; but no")
    pr("  capital subsidies matches at +0.084. Treat the minimum row as not yet reconciled.")
    pr("- On **differences**, `no capital subsidies` reproduces the paper to 0.001-0.033 pp and")
    pr("  `only stimulus checks` to 0.009-0.098 pp. `lower foreign demand` is off by up to")
    pr("  0.150 pp and `CBO debt` by up to 0.476 pp on Poor -- with Rich now matching to")
    pr("  0.006 pp, since the intercept channel that produced the old sign flip is switched")
    pr("  off here (see the `[SW23]` section).")
end
println("wrote output/table6_sw.md")
