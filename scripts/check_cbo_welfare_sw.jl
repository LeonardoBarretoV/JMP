# [SW22] read-only: is the CBO welfare pattern intrinsic, or a rule-2/forcing artifact?
# Compares the CBO-minus-baseline welfare DIFFERENCE on two pairs:
#   pre-rule-2:  cf_cbo_debt        - reduced_epsm     (ceiling bound naturally)
#   rule-2:      cf_cbo_debt_rule2_forced - rule2_epsm (hours mandated)
# No solving. Nothing written.
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
function ce_of(tag)
    o = load_object(J(tag)); ce = ce_summary_sw(o[1], o[3]); o = nothing; GC.gc(); ce
end
PAPER      = (mean=2.51, poor=4.00, middle=2.62, rich=0.50)
PAPER_BASE = (mean=1.82, poor=3.03, middle=2.00, rich=-0.20)
println("=== CBO welfare: pre-rule-2 (ceiling) vs rule-2 (mandate) ===")
b0 = ce_of("reduced_epsm");             c0 = ce_of("cf_cbo_debt")
b2 = ce_of("rule2_epsm");               c2 = ce_of("cf_cbo_debt_rule2_forced")
@printf("\n%-8s | %8s %8s %9s | %8s %8s %9s | %8s %8s %9s\n","group",
        "base0","cbo0","DIFF0","base2","cbo2","DIFF2","paperB","paper","PDIFF")
for k in (:mean,:poor,:middle,:rich)
    d0 = getfield(c0,k)-getfield(b0,k); d2 = getfield(c2,k)-getfield(b2,k)
    pd = getfield(PAPER,k)-getfield(PAPER_BASE,k)
    @printf("%-8s | %8.3f %8.3f %+9.3f | %8.3f %8.3f %+9.3f | %8.2f %8.2f %+9.2f\n",
            String(k), getfield(b0,k), getfield(c0,k), d0,
            getfield(b2,k), getfield(c2,k), d2,
            getfield(PAPER_BASE,k), getfield(PAPER,k), pd)
end
@printf("%-8s | %8.3f %8.3f %+9.3f | %8.3f %8.3f %+9.3f |\n","min",
        b0.min,c0.min,c0.min-b0.min, b2.min,c2.min,c2.min-b2.min)
println("\ndecile DIFFERENCES (counterfactual minus its own baseline):")
@printf("%6s %12s %12s %12s\n","bin","pre-rule-2","rule-2 forced","gap")
for k in 1:10
    d0 = c0.deciles[k]-b0.deciles[k]; d2 = c2.deciles[k]-b2.deciles[k]
    @printf("%6d %+12.4f %+12.4f %+12.4f\n", k, d0, d2, d2-d0)
end
println("\nif DIFF0 and DIFF2 agree, the pattern is INTRINSIC to the CBO experiment;")
println("if they disagree, rule 2 and/or the mandate is responsible.")
