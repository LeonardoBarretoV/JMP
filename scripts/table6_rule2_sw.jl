# ============================================================
# TABLE 6, all four columns on ONE fiscal baseline (rule 2), read-only.
#   baseline                 rule2_epsm.jld2
#   no capital subsidies     cf_no_capital_subsidy_rule2.jld2    [SW19]
#   only stimulus checks     cf_only_checks_rule2.jld2           option (ii)
#   lower foreign demand     cf_lower_foreign_rule2.jld2
#   CBO debt                 cf_cbo_debt_rule2_forced.jld2       [SW19]+[SW22] mandate
# Paper bin convention throughout (10/20/10 over initial-SS asset deciles).
# No solving. Nothing written.
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
J(f) = joinpath(ROOT,"output","jld2_sw",f*".jld2")

COLS = [("no capital subsidies","cf_no_capital_subsidy_rule2",(mean=2.02,poor=3.77,middle=2.27,rich=-1.69)),
        ("only stimulus checks","cf_only_checks_rule2",       nothing),
        ("lower foreign demand","cf_lower_foreign_rule2",     (mean=1.22,poor=2.03,middle=1.34,rich=0.06)),
        ("CBO debt (mandate)",  "cf_cbo_debt_rule2_forced",   (mean=2.51,poor=4.00,middle=2.62,rich=0.50))]
PB = (mean=1.82, poor=3.03, middle=2.00, rich=-0.20)

_b = load_object(J("rule2_epsm")); e_init = _b[1]
ceB = ce_summary_sw(e_init, _b[3]); rB = copy(_b[3].r); τB = copy(_b[3].τ_n)
yB1 = _b[3].y[1]; _b = nothing; GC.gc()
y_ss = y(e_init)

res = Any[]
for (nm, tag, pap) in COLS
    isfile(J(tag)) || (println("MISSING: ", tag); continue)
    o = load_object(J(tag))
    push!(res, (nm=nm, tag=tag, pap=pap, ce=ce_summary_sw(o[1], o[3]),
                τ1=o[3].τ_n[1], τ2=o[3].τ_n[2], y1=o[3].y[1], r1=o[3].r[1],
                π5=o[4][5], π20=o[4][20]))
    o = nothing; GC.gc()
end

println("\n", "="^118)
println("TABLE 6 -- CE welfare gains (%), all columns on the rule-2 fiscal baseline")
println("="^118)
@printf("%-24s %9s %9s %9s %9s %9s\n","","mean","poor","middle","rich","min")
@printf("%-24s %9.3f %9.3f %9.3f %9.3f %9.3f   <- BASELINE (rule2_epsm)\n",
        "baseline  model", ceB.mean, ceB.poor, ceB.middle, ceB.rich, ceB.min)
@printf("%-24s %9.2f %9.2f %9.2f %9.2f %9s   <- paper baseline\n",
        "baseline  paper", PB.mean, PB.poor, PB.middle, PB.rich, "-")
for x in res
    println("-"^118)
    @printf("%-24s %9.3f %9.3f %9.3f %9.3f %9.3f\n","  "*x.nm*"  model",
            x.ce.mean, x.ce.poor, x.ce.middle, x.ce.rich, x.ce.min)
    if isnothing(x.pap)
        @printf("%-24s %9s %9s %9s %9s  (no like-for-like paper column)\n",
                "  "*x.nm*"  paper","-","-","-","-")
    else
        @printf("%-24s %9.2f %9.2f %9.2f %9.2f\n","  "*x.nm*"  paper",
                x.pap.mean, x.pap.poor, x.pap.middle, x.pap.rich)
    end
end

println("\n", "="^118)
println("THE OBJECT OF INTEREST: the DIFFERENCE from the baseline (counterfactual minus baseline)")
println("="^118)
@printf("%-24s %19s %19s %19s %19s\n","column","mean","poor","middle","rich")
@printf("%-24s %8s %10s %8s %10s %8s %10s %8s %10s\n","",
        "model","paper","model","paper","model","paper","model","paper")
for x in res
    ms = [getfield(x.ce,k)-getfield(ceB,k) for k in (:mean,:poor,:middle,:rich)]
    if isnothing(x.pap)
        @printf("%-24s %+8.3f %10s %+8.3f %10s %+8.3f %10s %+8.3f %10s\n", x.nm,
                ms[1],"-",ms[2],"-",ms[3],"-",ms[4],"-")
    else
        ps = [getfield(x.pap,k)-getfield(PB,k) for k in (:mean,:poor,:middle,:rich)]
        @printf("%-24s %+8.3f %+10.2f %+8.3f %+10.2f %+8.3f %+10.2f %+8.3f %+10.2f\n", x.nm,
                ms[1],ps[1],ms[2],ps[2],ms[3],ps[3],ms[4],ps[4])
    end
end
println("\ngap between the model difference and the paper difference (pp):")
@printf("%-24s %10s %10s %10s %10s\n","column","mean","poor","middle","rich")
for x in res
    isnothing(x.pap) && continue
    @printf("%-24s", x.nm)
    for k in (:mean,:poor,:middle,:rich)
        @printf(" %+10.3f", (getfield(x.ce,k)-getfield(ceB,k)) - (getfield(x.pap,k)-getfield(PB,k)))
    end
    println()
end

println("\n", "="^118)
println("FULL DECILE PROFILE (level, %)")
println("="^118)
@printf("%5s %12s", "bin", "baseline"); for x in res; @printf(" %20s", first(x.nm,20)); end; println()
for k in 1:10
    @printf("%5d %12.4f", k, ceB.deciles[k])
    for x in res; @printf(" %20.4f", x.ce.deciles[k]); end; println()
end
println("\ndifference from the baseline:")
@printf("%5s %12s","bin",""); for x in res; @printf(" %20s", first(x.nm,20)); end; println()
for k in 1:10
    @printf("%5d %12s", k, "")
    for x in res; @printf(" %+20.4f", x.ce.deciles[k]-ceB.deciles[k]); end; println()
end

println("\n", "="^118)
println("SUPPORTING PATHS")
println("="^118)
@printf("%-24s %12s %12s %12s %14s %12s %12s\n","column","tau_n j=1","tau_n j=2",
        "r j=1 (ann%)","y j=1 fall %","pi j=5 ann%","pi j=20 ann%")
@printf("%-24s %12.8f %12.8f %12.4f %14.4f %12.4f %12.4f\n","baseline",
        τB[1], τB[2], 100*((1+rB[1])^4-1), 100*(1-yB1/y_ss), NaN, NaN)
for x in res
    @printf("%-24s %12.8f %12.8f %12.4f %14.4f %12.4f %12.4f\n", x.nm,
            x.τ1, x.τ2, 100*((1+x.r1)^4-1), 100*(1-x.y1/y_ss),
            100*((1+x.π5)^4-1), 100*((1+x.π20)^4-1))
end
@printf("\ninitial-SS tau_n = %.8f   y_ss = %.10f\n", e_init.τ_n, y_ss)
