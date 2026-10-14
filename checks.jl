# ============================================================
# checks.jl  --  post-run diagnostics for the [Z-NORM] change
#
# HOW TO USE: run scripts/main_JMP.jl first (so e_init, e_final,
# path_r, h_path, c, t, g are in the workspace), then select a
# block below and hit Ctrl+Enter.  Each block is self-contained,
# so you can run them in any order.
#
# WHAT YOU ARE TESTING: the Z-normalization scales the whole
# economy by ~1.573 in LEVELS.  So raw levels should all move by
# that same factor, while every RATIO and every path normalized
# by its own steady state should be essentially unchanged from
# the un-normalized run.  Levels moving together = fine.
# Ratios moving = something real changed.
# ============================================================

using Printf

## ===========================================================
## BLOCK 0 -- is the workspace loaded?
## ===========================================================
for v in [:e_init, :e_final, :path_r, :h_path]
    @printf("%-10s %s\n", v, isdefined(Main, v) ? "OK" : "MISSING -- run scripts/main_JMP.jl first")
end


## ===========================================================
## BLOCK 1 -- steady states and the level factor
## ===========================================================
println("\n--- STEADY STATES ---")
@printf("Z_ss = %.12f   (should be exactly 1.0)\n", sum(e_init.h.Pss .* e_init.h.z_grid))
for (nm, es) in (("init", e_init), ("final", e_final))
    @printf("%-5s r=%7.4f%%a  tau_n=%.4f  w=%.5f  N=%.4f  Y=%.5f  b/y=%.3f  a/y=%.3f  k/y=%.3f  valid=%s\n",
        nm, 400*es.r, es.τ_n, es.w, es.n, y(es), es.b/y(es), es.a/y(es), es.k/y(es), is_pol_valid(es))
end
# Reference SS from the un-normalized grid (solved 2026-08-24, same calibration):
#   r = 0.00248356, tau_n = 0.816555, w = 0.095449, N = 0.316312, K = 0.779003,
#   Y = 0.045062, b/y = 4.2000, a/y = 20.2473, k/y = 17.2873, valid = true
@printf("\nlevel factor vs un-normalized run:  Y=%.4f  N=%.4f  K=%.4f  A=%.4f   (all should be ~1.573)\n",
        y(e_init)/0.045062, e_init.n/0.316312, e_init.k/0.779003, e_init.a/0.912387)
@printf("ratios (should be ~unchanged):      b/y %.4f vs 4.2000   a/y %.4f vs 20.2473   k/y %.4f vs 17.2873\n",
        e_init.b/y(e_init), e_init.a/y(e_init), e_init.k/y(e_init))
@printf("r:  %.8f vs 0.00248356   (%.2f%% change)\n", e_init.r, 100*(e_init.r/0.00248356 - 1))
@printf("mass at a_max = %.3e   (was 0.0 before the change -- known accepted cost)\n", sum(e_init.ws.pdf[end,:]))


## ===========================================================
## BLOCK 1B -- is the TOP of the asset grid binding?
## a_max is still in the OLD (un-normalized) units while incomes
## grew ~1.57x, so this is the diagnostic most exposed to the
## Z-NORM change. "% at a_max" is the share of households sitting
## exactly on the last grid point.
## ===========================================================
amax = last(e_init.h.a_grid)
amin = first(e_init.h.a_grid)
ntop = 10
println("\n--- TOP OF THE ASSET GRID ---")
@printf("a_grid: %d pts,  a_min = %.5f,  a_max = %.4f,  2nd-highest = %.4f\n\n",
        length(e_init.h.a_grid), amin, amax, e_init.h.a_grid[end-1])
println("                     % at a_max    % in top 10     max a_pol   a_pol/a_max    % at a_min")
for (nm, ws) in (("initial SS", e_init.ws), ("final SS", e_final.ws))
    @printf("%-18s %11.6f %14.6f %12.3f %13.4f %13.3f\n", nm,
        100*sum(ws.pdf[end,:]), 100*sum(ws.pdf[end-ntop+1:end,:]),
        maximum(ws.a_pol), maximum(ws.a_pol)/amax, 100*sum(ws.pdf[1,:]))
end

mass_top  = [sum(path_r[i].pdf[end,:])              for i in eachindex(path_r)]
mass_topN = [sum(path_r[i].pdf[end-ntop+1:end,:])   for i in eachindex(path_r)]
apol_max  = [maximum(path_r[i].a_pol)               for i in eachindex(path_r)]
mass_bot  = [sum(path_r[i].pdf[1,:])                for i in eachindex(path_r)]
it, ia = argmax(mass_top), argmax(apol_max)
@printf("%-18s %11.6f %14.6f %12.3f %13.4f %13.3f\n", "transition (worst)",
    100*mass_top[it], 100*maximum(mass_topN), apol_max[ia], apol_max[ia]/amax, 100*maximum(mass_bot))

@printf("\n   worst mass at a_max : %.4e (%.6f%% of hh) at t=%d\n", mass_top[it], 100*mass_top[it], it)
@printf("   worst max a_pol     : %.4f at t=%d   (a_max = %.4f)\n", apol_max[ia], ia, amax)
@printf("   periods with ANY mass at a_max            : %d of %d\n", count(>(0), mass_top), length(mass_top))
@printf("   periods where a_pol reaches/exceeds a_max : %d of %d\n", count(>=(amax), apol_max), length(apol_max))

binding = (maximum(mass_top) > 0) || (maximum(apol_max) >= amax) ||
          (sum(e_init.ws.pdf[end,:]) > 0) || (maximum(e_init.ws.a_pol) >= amax)
if binding
    @printf("\n   => a_max IS BINDING. Households want to save past the top of the grid,\n")
    @printf("      so the top of the wealth distribution is censored. Fix = raise a_max\n")
    @printf("      (scaling it by ~1.573 gives %.1f) or rescale the grid with the units-only variant.\n", 1.573*amax)
else
    println("\n   => a_max is NOT binding: no mass on the last grid point and a_pol stays inside the grid.")
end

## ===========================================================
## BLOCK 2 -- did the transition actually converge?
## RUN THIS BEFORE TRUSTING ANYTHING BELOW.
## ===========================================================
println("\n--- CONVERGENCE ---")
res = [path_r.a[j] - path_r.k[j] - path_r.q[j-1]*path_r.bd[j] for j in 2:length(path_r.a)]
@printf("asset-market residual   max|res| = %.3e  at t=%d   (nlsolve ftol = %.0e)\n",
        maximum(abs, res), argmax(abs.(res))+1, _ZERO_FTOL)
@printf("                        L2 norm  = %.3e\n", sqrt(sum(abs2, res)))
@printf("period-0 capital market resid    = %.3e\n", e_init.k - path_r.k[1])
@printf("period-0 bond market resid       = %.3e\n", e_init.b - path_r.b[1])
@printf("terminal match to e_final:  r %.3e   k %.3e   q %.3e\n",
        path_r.r[end]-e_final.r, path_r.k[end]-e_final.k, path_r.q[end]-e_final.q)
@printf("any NaN/Inf in r,k,n,y,q,C?  %s\n",
        any(!isfinite, vcat(path_r.r, path_r.k, path_r.n, path_r.y, path_r.q, path_r.C)))


## ===========================================================
## BLOCK 3 -- Z_t and the normalized transition paths
## This is the table to eyeball against your old figures.
## ===========================================================
z_state = zeros(length(h_path), length(e_init.h.z_grid))
z_state[1,:] = e_init.h.Pss' * h_path[1].P
for i = 2:lastindex(h_path)
    z_state[i,:] = z_state[i-1,:]' * h_path[i].P
end
Z_t = [sum(z_state[i,:] .* h_path[i].z_grid) for i in eachindex(h_path)]

ys, ns, ks, cs = y(e_init), e_init.n, e_init.k, e_init.c
println("\n--- TRANSITION (normalized by initial SS) ---")
println("  t    un%     Z_t     Y/Yss   N/Nss   K/Kss   C/Css    r(a%)      q      pi(a%)   b/y")
for i in vcat(1:8, 12, 20, 40, 80)
    @printf("%3d  %6.2f  %.4f  %.4f  %.4f  %.4f  %.4f  %7.3f  %.4f  %6.2f  %6.3f\n",
        i, 100*path_r.un[i], Z_t[i], path_r.y[i]/ys, path_r.n[i]/ns, path_r.k[i]/ks,
        path_r.C[i]/cs, 400*path_r.r[i], path_r.q[i], 400*path_r.π_[i],
        path_r.b[i]/path_r.y[i])
end
# Reference Z_t under this change:  t=1 -> 0.944835,  t=12 -> 0.978226,  t=40 -> 0.995530
# NOTE: path_r.C[1] and path_r.G[1] are built in LEVELS at transition.jl:523,541 while
# i>=2 uses ratios*y, so the t=1 consumption entry can look out of line. Pre-existing.
@printf("\nZ_t check:  t=1 %.6f (exp 0.944835)   t=12 %.6f (exp 0.978226)   t=40 %.6f (exp 0.995530)\n",
        Z_t[1], Z_t[12], Z_t[40])
@printf("hours per efficiency unit N_t/Z_t:  t=1 %.4f   t=40 %.4f   ss %.4f\n",
        path_r.n[1]/Z_t[1], path_r.n[40]/Z_t[40], e_init.n)


## ===========================================================
## BLOCK 4 -- headline numbers to compare against the paper
## ===========================================================
ys = y(e_init)
println("\n--- HEADLINES (first 40 quarters) ---")
@printf("trough Y    = %+.2f%%   at t=%d\n", 100*(minimum(path_r.y[1:40])/ys-1), argmin(path_r.y[1:40]))
@printf("trough C    = %+.2f%%   at t=%d\n", 100*(minimum(path_r.C[2:40])/e_init.c-1), argmin(path_r.C[2:40])+1)
@printf("peak un     =  %.2f%%   at t=%d\n", 100*maximum(path_r.un[1:40]), argmax(path_r.un[1:40]))
@printf("peak b/y    =  %.3f    at t=%d\n", maximum(path_r.b[1:40]./path_r.y[1:40]), argmax(path_r.b[1:40]./path_r.y[1:40]))
@printf("r  range    =  %.3f to %.3f  (ann %%)\n", 400*minimum(path_r.r[1:40]), 400*maximum(path_r.r[1:40]))
@printf("q  range    =  %.4f to %.4f\n", minimum(path_r.q[1:40]), maximum(path_r.q[1:40]))
@printf("pi range    =  %.2f to %.2f  (ann %%, exogenous -- must match old run exactly)\n",
        400*minimum(path_r.π_[1:40]), 400*maximum(path_r.π_[1:40]))


## ===========================================================
## BLOCK 5 -- welfare / distribution tables
## Most sensitive to the change, since the asset grid is the one
## object still denominated in the OLD units.
## ===========================================================
println(generate_tables(e_init, e_final, path_r, io=String)...)


## ===========================================================
## BLOCK 6 -- save these series so the NEXT change can be diffed
## against this run. Writes output/jld2/checks_znorm.jld2
## ===========================================================
let
    snap = Dict(
        "Z_ss"  => sum(e_init.h.Pss .* e_init.h.z_grid),
        "ss"    => [e_init.r, e_init.τ_n, e_init.w, e_init.n, e_init.k, y(e_init),
                    e_init.a, e_init.b, e_init.c],
        "y"=>collect(path_r.y), "n"=>collect(path_r.n), "k"=>collect(path_r.k),
        "C"=>collect(path_r.C), "r"=>collect(path_r.r), "q"=>collect(path_r.q),
        "b"=>collect(path_r.b), "un"=>collect(path_r.un),
    )
    d = joinpath(ROOT, "output", "jld2"); mkpath(d)
    save_object(joinpath(d, "checks_znorm.jld2"), snap)
    println("\nsaved -> ", joinpath(d, "checks_znorm.jld2"))
end
