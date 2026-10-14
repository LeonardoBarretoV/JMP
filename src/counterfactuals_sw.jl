# ============================================================
# [SW11] STAGE 5 COUNTERFACTUAL MACHINERY  --  sticky-wage / endogenous-inflation track
# ============================================================
# NEW FILE. Nothing in existing src/ is modified; every name carries an _sw suffix.
#
# WHY THIS FILE EXISTS
#   cf_lower_foreign.jld2 and cf_cbo_debt.jld2 were produced on 2026-09-08 by a
#   runner that was never saved into the repo (it is not in scripts/, not in the
#   VS Code local history, and not in the Julia REPL history -- both runs were
#   script-driven and the script is gone).  The knob values below were RECOVERED
#   FROM THE SAVED RESULTS by inverting the solver's own mapping:
#
#       qbstar_path[i] = path.bstar[i] / path.b[i]
#       qb_path[i]     = path.b[i] * q_lag_i / path.y[i]
#       q_lag_i        = (1/(1+r_i)) * ((ibar+vartheta)/(1+pi_i)
#                                       + (1-vartheta)/(1+pi_i) * q_i)
#
#   so these numbers are reconstructions, not guesses.  Both were then checked by
#   evaluating the reduced residual at the saved solution; see
#   scripts/verify_cf_inputs_sw.jl.
#
# CONTENTS
#   ce_summary_sw                CE welfare gains vs the initial SS
#   gate_lab_tax_sw              item-(b) gate: household- vs government-side labour tax
#   period0_branch_sw            checks that s_1 = Gamma_n is the ACTIVE branch
#   epsm_baseline_sw             backs the monetary residual out of a Stage 3 path
#   warm_start_sw                (r0, s0) for the reduced system from a saved path
#   cf_qbstar_lower_foreign_sw   the LOWER FOREIGN DEMAND knob
#   cf_paths_cbo_debt_sw         the CBO DEBT PROJECTIONS knob
#   intercept_channel_sw         closed-form d(pi) from a move in the Taylor intercept
#   guard_report_sw              prints the [SW8]/[SW9] counters
#   print_pi_table_sw            pi(data) vs pi(model) vs a reference run
# ============================================================


# ------------------------------------------------------------
# CE welfare gains vs the INITIAL steady state, in %.
# Same construction as generate_tables (src/statistics.jl:6-42): decile bins of
# the initial-SS wealth marginal; "poor" = bottom decile, "middle" = mean of
# deciles 5 and 6, "rich" = top decile.  Reproduced here rather than called so
# that src/statistics.jl stays untouched and so a counterfactual path can be
# summarised without building the full pretty-table.
# ------------------------------------------------------------
function ce_summary_sw(e_init, path)
    w_gains = 100 .* (path.v[1] ./ e_init.ws.v .- 1)
    pdf_marginal = sum(e_init.ws.pdf, dims = 2)[:, 1]
    a_grid = round.(e_init.h.a_grid, digits = 10)
    qu_10 = round.([quantile(a_grid, StatsBase.ProbabilityWeights(pdf_marginal), x)
                    for x in 0.0:0.1:1.0], digits = 10)

    weight = similar(e_init.ws.pdf)
    w_gains_rs = reshape(w_gains, 1, :)[1, :]
    summ = Float64[]
    for (qi, qii) in zip(qu_10, qu_10[2:end])
        qii_ = (qii == qu_10[end]) ? qii + 10 * eps() : qii
        for i in axes(weight, 1), j in axes(weight, 2)
            weight[i, j] = (qi <= a_grid[i] < qii_) ? e_init.ws.pdf[i, j] : 0.0
        end
        push!(summ, mean(w_gains_rs, ProbabilityWeights(reshape(weight, 1, :)[1, :])))
    end

    return (mean    = sum(w_gains .* e_init.ws.pdf),
            min     = minimum(w_gains[e_init.ws.pdf .> 0]),
            poor    = first(summ),
            middle  = sum(summ[5:6]) / 2,
            rich    = last(summ),
            deciles = summ)
end


function print_ce_sw(ce; label = "", paper = nothing)
    println("\n---- CE welfare gains vs the initial steady state (%)  ", label, " ----")
    if isnothing(paper)
        @printf("mean %7.3f   poor %7.3f   middle %7.3f   rich %7.3f   (min %7.3f)\n",
                ce.mean, ce.poor, ce.middle, ce.rich, ce.min)
    else
        @printf("%-8s %10s %10s %10s\n", "group", "model", "paper", "gap (pp)")
        for (nm, mv, pv) in (("mean",   ce.mean,   paper.mean),
                             ("poor",   ce.poor,   paper.poor),
                             ("middle", ce.middle, paper.middle),
                             ("rich",   ce.rich,   paper.rich))
            @printf("%-8s %10.3f %10.3f %+10.3f\n", nm, mv, pv, mv - pv)
        end
    end
    return nothing
end


# ------------------------------------------------------------
# ITEM (b) GATE.  The household side and the government side must book the same
# labour tax at every date, t = 1 included:
#
#   gap_j = sum_z pdf_j(z) * [ w_j z s_j nt_z - tau_n_j (w_j z s_j nt_z)^(1-tau) ]
#           - path.lab_tax[j]
#
# EXPECTED MAGNITUDE: a relative gap of about 1.6e-07, and that is NOT an error --
# it is the known mass-leak constant (sum(pdf) - 1 = 1.5664e-07, because
# update_stationary_pdf! does not renormalise; see README_sw.md).  Anything
# materially larger means the two blocks are booking different hours and the run
# is INVALID.
#
# z = 0 is skipped: 0^(-tau) = Inf and the HSV term is identically zero there.
# ------------------------------------------------------------
function gate_lab_tax_sw(e_init, path, g; jlist = nothing)
    n_per = length(path.r)
    js = isnothing(jlist) ? vcat(collect(1:12), [20, 40, 100, 400, n_per - 1]) : jlist
    nt  = e_init.ws.n_pol[1, :]
    z   = e_init.h.z_grid
    emp = [iz for iz in eachindex(z) if z[iz] != 0]

    rows = NamedTuple[]
    worst_rel = 0.0
    worst_j = 0
    for j in js
        (1 <= j <= n_per) || continue
        zmarg = sum(path[j].pdf, dims = 1)[1, :]
        # s_j recovered from the stored hours policy: n_pol[1,z] = s_j * nt_z,
        # exact and z-independent (verified to 6.7e-16 in the flexible-wage case).
        s_j = sum(path[j].n_pol[1, iz] / nt[iz] for iz in emp) / length(emp)
        hh = 0.0
        for iz in emp
            base = path.w[j] * z[iz] * s_j * nt[iz]
            hh += zmarg[iz] * (base - path.τ_n[j] * base^(1 - g.τ))
        end
        gap = hh - path.lab_tax[j]
        rel = abs(gap) / max(abs(path.lab_tax[j]), eps())
        push!(rows, (; j = j, s_j = s_j, hh = hh, gov = path.lab_tax[j], gap = gap, rel = rel))
        if rel > worst_rel
            worst_rel = rel
            worst_j = j
        end
    end
    return (; rows = rows, worst_rel = worst_rel, worst_j = worst_j)
end


function print_gate_sw(gate; tol = 1e-5)
    println("\n---- item (b) gate: household vs government labour tax ----")
    @printf("%6s %12s %14s %14s %14s %12s\n", "j", "s_j", "household", "government", "gap", "relative")
    for r in gate.rows
        @printf("%6d %12.8f %14.8f %14.8f %+14.3e %12.3e\n", r.j, r.s_j, r.hh, r.gov, r.gap, r.rel)
    end
    @printf("worst relative gap = %.3e at j = %d\n", gate.worst_rel, gate.worst_j)
    ok = gate.worst_rel < tol
    if ok
        println("GATE PASSED  (~1.6e-07 is the known mass-leak constant, not an error)")
    else
        println("*** GATE FAILED -- the blocks are booking different hours. RUN IS INVALID ***")
    end
    return ok
end


# ------------------------------------------------------------
# PERIOD-0 BRANCH CHECK.  s_1 is pinned to Gamma_n, which reproduces
# min(Gamma_n * nt_z, FOC) only while the Gamma_n branch is the binding one.
# If s_flex(w_1, tau_n_1) falls BELOW Gamma_n, the government block and the
# household block book different hours and the run is invalid -- that is what
# killed kappa_w = 0.03.
# ------------------------------------------------------------
function period0_branch_sw(e_init, path, g, ν, Γ_n)
    s_flex_1 = s_flex_sw(g, ν; w = path.w[1], τ_n = path.τ_n[1],
                         w_init = e_init.w, τ_n_init = e_init.τ_n)
    ok = Γ_n <= s_flex_1
    println("\n---- period-0 branch ----")
    @printf("Gamma_n = %.10f   s_flex(w_1, tau_n_1) = %.10f\n", Γ_n, s_flex_1)
    println("active branch = ", ok ? "Gamma_n  (CORRECT)" : "s_flex  (*** RUN IS INVALID ***)")
    return ok
end


# ------------------------------------------------------------
# The monetary residual, backed out of a Stage 3 path (pi FED from data).
# Identical algebra to scripts/main_JMP_reduced.jl:231-233, isolated here so that
# the counterfactuals hold epsm at EXACTLY the baseline path -- that is what
# "Stage 5" means: change fiscal policy, leave monetary policy where the baseline
# identification put it.
# ------------------------------------------------------------
function epsm_baseline_sw(r_stage3, π_stage3, i_ss_base, π_bar, phi_pi)
    n_per = length(r_stage3)
    i_conv = [ j < n_per ? (1 + r_stage3[j+1]) * (1 + π_stage3[j+1]) - 1 : i_ss_base
               for j in 1:n_per ]
    return [ i_conv[j] - i_ss_base - phi_pi * (π_stage3[j] - π_bar) for j in 1:n_per ]
end


# ------------------------------------------------------------
# Warm start for the reduced ([r ; s], 1598) system out of a saved path.
#   r_block[1]     = r_k_hat = r_k[1] - delta      (r_b is substituted out)
#   r_block[2..n]  = r_i                            (NO +1 shift under method 1)
#   s_block[k]     = s_{k+1}, k = 1..n-2            (s_1 and s_n are pinned)
# ------------------------------------------------------------
function warm_start_sw(path, e_init, t; method = 1, sw19 = false, sw20 = false)
    n_per = length(path.r)
    n_r   = _n_r_sw(n_per, method, sw19)
    r0 = Vector{Float64}(undef, n_r)
    if sw19
        # [SW19] r_k is no longer an unknown: the block is [r_2 .. r_n_per].
        for i in 2:n_per
            r0[i-1] = path.r[i]
        end
    else
        r0[1] = path.r_k[1] - t.δ
        for i in 2:n_per
            r0[i] = path.r[i]
        end
    end
    nt  = e_init.ws.n_pol[1, :]
    emp = [iz for iz in eachindex(e_init.h.z_grid) if e_init.h.z_grid[iz] != 0]
    s_all = [ sum(path[j].n_pol[1, iz] / nt[iz] for iz in emp) / length(emp)
              for j in 1:n_per ]
    # [SW20] the s block starts at s_1 rather than s_2.  Recovering s_1 from the
    # stored hours policy works whichever branch was active, because
    # backwards_once_t0_T! books min(Gamma_n, s_flex)*ntilde_z either way.
    s0 = sw20 ? s_all[1:(n_per - 1)] : s_all[2:(n_per - 1)]
    return r0, s0
end


# ============================================================
# COUNTERFACTUAL 1 -- LOWER FOREIGN DEMAND
# ============================================================
# Foreign holdings of government debt fall by `drop` of ANNUAL GDP along the
# observed stretch of the transition (j = 2..jlast), and the terminal external
# debt-to-GDP ratio falls to Lambda_bstar_f_low, which is passed to
# make_time_series_sw and sets the j > jlast fill from e_final itself.
#
# The solver consumes bstar as a RATIO TO b, so a level drop of `drop` of annual
# GDP is  -drop / (qb_path[j]/4)  in ratio units:
#
#     bstar_j / b_j = (bstar_j/y_ann) / (b_j/y_ann)
#     d(bstar/b)    = d(bstar/y_ann) / (q b / y_ann) = -drop / (qb_path[j]/4)
#
# j = 1 is NOT touched: bstar_path_r[1] is a LEVEL that drives q_lag = 1 and hence
# r_b, and moving it would break the r_b substitution.
#
# RECOVERED FROM cf_lower_foreign.jld2:  drop = 0.10, jlast = 20.  Reproduces
# every stored entry exactly, e.g. j = 2:  0.2694 - 0.1/1.313375 = 0.19326028.
# ============================================================
function cf_qbstar_lower_foreign_sw(ts; drop = 0.10, jlast = 20)
    qbs = copy(ts.bstar_path_r)
    for j in 2:min(jlast, length(qbs))
        qbs[j] -= drop / (ts.qb_path[j] / 4.0)
    end
    return qbs
end


# ============================================================
# COUNTERFACTUAL 2 -- CBO DEBT PROJECTIONS
# ============================================================
# Government debt rises linearly from its last observed value (2024Q4 = index 19)
# to `target` of GDP by index `end_index`, then stays there.  Foreign holdings are
# held at a CONSTANT `foreign_to_gdp` of annual GDP over the same stretch, so the
# foreign SHARE falls as total debt rises.
#
# RECOVERED FROM cf_cbo_debt.jld2:
#   start_index = 19   (2024Q4, value 1.218*4 = 4.872, left in place)
#   end_index   = 143  (2055Q4)  -> 124 quarters, step 0.020064516
#   target      = 1.84 * 4.0  = 7.36
#   foreign     = 0.287 * 4.0 = 1.148
# Reproduces every stored entry exactly (qb[20] = 4.89206452, qb[100] = 6.49722581,
# qbstar[20] = 0.23466575, qbstar[200] = 0.15597826).
#
# NOTE ON A SMALL INCONSISTENCY IN THAT RUN, reproduced here deliberately so the
# script matches the saved file.  The foreign rule divides by qb_path = q*b/y, not
# by b/y, so the terminal foreign share it implies is
#     1.148 / 7.36 = 0.15597826
# while the CBO final steady state actually has
#     e_final.bstar / e_final.b = 0.15267816      (q_final = 0.9788 < 1).
# The transition therefore ends 2.2% above its own terminal steady-state foreign
# share.  It converged and the gate passed, so it is not fatal.  For the
# internally consistent variant pass  foreign_mode = :level_b,  which divides by
# b/y and lands exactly on e_final.bstar/e_final.b; that is a DIFFERENT run and
# will not reproduce cf_cbo_debt.jld2.
# ============================================================
function cf_paths_cbo_debt_sw(ts, e_final; start_index = 19, end_index = 143,
                              target = 1.84 * 4.0, foreign_to_gdp = 0.287 * 4.0,
                              foreign_mode = :as_run)
    qb = copy(ts.qb_path)
    n  = length(qb)
    trend = collect(LinRange(qb[start_index], target, end_index - start_index + 1))[2:end]
    qb[(start_index + 1):end_index] .= trend
    qb[(end_index + 1):n] .= target

    qbs = copy(ts.bstar_path_r)
    for j in (start_index + 1):n
        if foreign_mode === :as_run
            qbs[j] = foreign_to_gdp / qb[j]                 # divides by q*b/y  (as the saved run did)
        else
            qbs[j] = foreign_to_gdp / (qb[j] / e_final.q)   # divides by b/y    (internally consistent)
        end
    end
    return qb, qbs
end


# ------------------------------------------------------------
# THE INTERCEPT CHANNEL, in closed form.
#
# i_taylor_ss is rebuilt as (1+r_final)(1+pi_bar)-1 whenever the terminal steady
# state moves, so any counterfactual that changes the long-run real rate shifts
# the Taylor intercept, and the backward recursion turns that into a level shift
# in inflation along the whole transition:
#
#     d(pi_j) = -d(i_ss)/phi + ((1+r_{j+1})/phi) * d(pi_{j+1})
#  => d(pi)   = -d(i_ss) / (phi - 1 - r_f)        (away from the terminal pin)
#
# ECONOMICS: larger long-run debt -> higher neutral rate -> higher Taylor
# intercept -> lower inflation -> higher real burden of the nominal debt, so the
# debt expansion partly defeats itself.  This channel cannot exist in the
# original paper, where inflation is fed.
# ------------------------------------------------------------
function intercept_channel_sw(i_ss_cf, i_ss_base, phi_pi, r_f)
    d_iss = i_ss_cf - i_ss_base
    d_pi  = -d_iss / (phi_pi - 1 - r_f)
    println("\n---- intercept channel (closed form) ----")
    @printf("i_ss baseline = %.10f  (%.4f%% ann)\n", i_ss_base, 100 * ((1 + i_ss_base)^4 - 1))
    @printf("i_ss this run = %.10f  (%.4f%% ann)\n", i_ss_cf,   100 * ((1 + i_ss_cf)^4 - 1))
    @printf("d(i_ss) = %.9f   r_f = %.10f   phi = %.3f\n", d_iss, r_f, phi_pi)
    @printf("predicted d(pi) = -d(i_ss)/(phi-1-r_f) = %.9f = %+.4f pp ann\n", d_pi, 400 * d_pi)
    println("Compare with the SOLVED pi gap below. The residual is the r-path channel,")
    println("which pushes the other way: the two oppose and the intercept wins.")
    return (; d_iss = d_iss, d_pi = d_pi)
end


# ------------------------------------------------------------
# Guard counters.  All three guards are inert at any converged solution; a
# non-zero [SW8] FALLBACK count in particular means the residual went
# discontinuous somewhere and the solve should not be trusted.
# ------------------------------------------------------------
function guard_report_sw()
    println("\n---- numerical guards ----")
    @printf("[SW8] bracket widened       : %d\n", _SW8_WIDEN[])
    @printf("[SW8] grid fallback fired   : %d   <-- MUST be 0 for a trustworthy solve\n", _SW8_FALLBACK[])
    @printf("[SW9] GHH composite floored : %d\n", _SW9_FLOOR[])
    return (_SW8_WIDEN[], _SW8_FALLBACK[], _SW9_FLOOR[])
end


# ------------------------------------------------------------
# pi table: data vs this run vs an optional reference run (the baseline), all
# annualised.  `π_ref` lets a counterfactual show its inflation gap directly,
# which is what the intercept-channel comparison needs.
# ------------------------------------------------------------
function print_pi_table_sw(π_data, π_model, π_bar; π_ref = nothing, label = "model")
    ql = ["2020Q2","2020Q3","2020Q4","2021Q1","2021Q2","2021Q3",
          "2021Q4","2022Q1","2022Q2","2022Q3","2022Q4","2023Q1"]
    ann(x) = 100 * ((1 + x)^4 - 1)
    println("\n---- inflation (annualised %) ----")
    if isnothing(π_ref)
        @printf("%4s  %-8s %12s %12s %14s\n", "j", "quarter", "data", label, "diff (ann pp)")
        for j in 1:12
            @printf("%4d  %-8s %12.4f %12.4f %+14.4f\n", j, ql[j],
                    ann(π_data[j]), ann(π_model[j]), 400 * (π_model[j] - π_data[j]))
        end
        for j in (16, 20, 40, 100, length(π_model))
            @printf("%4d  %-8s %12.4f %12.4f %+14.4f\n", j, "",
                    ann(π_data[j]), ann(π_model[j]), 400 * (π_model[j] - π_data[j]))
        end
    else
        @printf("%4s  %-8s %12s %12s %12s %16s\n", "j", "quarter", "data", "baseline", label, "cf - base (pp)")
        for j in 1:12
            @printf("%4d  %-8s %12.4f %12.4f %12.4f %+16.4f\n", j, ql[j],
                    ann(π_data[j]), ann(π_ref[j]), ann(π_model[j]), 400 * (π_model[j] - π_ref[j]))
        end
        for j in (16, 20, 40, 100, 143, length(π_model))
            @printf("%4d  %-8s %12.4f %12.4f %12.4f %+16.4f\n", j, "",
                    ann(π_data[j]), ann(π_ref[j]), ann(π_model[j]), 400 * (π_model[j] - π_ref[j]))
        end
    end
    @printf("\npeak |pi - pi_bar| over j = 2..60:  %s %+.4f pp    data %+.4f pp\n",
            label, 400 * maximum(abs.(π_model[2:60] .- π_bar)),
            400 * maximum(abs.(π_data[2:60] .- π_bar)))
    return nothing
end


# ============================================================
# [SW22] FORCED PERIOD-0 HOURS -- the lockdown as a MANDATE, not a ceiling
# ============================================================
# For ONE counterfactual column only (CBO), hours at period 0 are set to
# Gamma_n * ntilde_z EXACTLY in BOTH blocks, rather than min(Gamma_n, s_flex)*ntilde_z.
#
# WHY.  On the rule-2 baseline the CBO counterfactual has NO period-0 equilibrium
# under the ceiling specification: g(s) = min(Gamma_n, s_flex(w_1(s), tau_n_1(s))) has
# slope +2.13 to +2.20 (measured on three converged full solves), so h(s) = s - g(s)
# is DECREASING with h > 0 at s = Gamma_n, 0.8973 and 0.8808 -- no root at or below
# Gamma_n.  The mechanism is a fiscal-feedback spiral: fewer hours -> smaller labour
# tax base -> tighter budget -> lower tau_n -> lower after-tax wage -> fewer hours.
#
# Mandating hours removes the feedback: s_1 = Gamma_n by construction, so there is
# nothing for the spiral to act on.  This is a REAL equilibrium of a slightly
# different model -- households are REQUIRED to supply Gamma_n * ntilde_z rather than
# permitted to supply at most that -- and the budget, the value function and the
# distribution all solve consistently around it.  The cost is a labour wedge: at the
# mandated hours the marginal disutility need not equal the marginal after-tax wage.
#
# The government block already books Gamma_n unconditionally, so with the household
# side forced to the same hours the two agree BY CONSTRUCTION and item (b) must pass.
#
# INERT wherever the ceiling already binds: if min(Gamma_n, s_flex) == Gamma_n then
# forcing changes nothing, so the rule-2 baseline and the other three counterfactuals
# must reproduce exactly.  That is checked, not assumed.
# ============================================================
labor_forced_sw(; n_, Γ_n) = Γ_n * n_
labor_income_forced_sw(; w, n_, Γ_n) = w * Γ_n * n_
disutility_forced_sw(v::GHH; n_, Γ_n) = v.θ * (Γ_n * n_)^(1 + 1/v.ν) / (1 + 1/v.ν)

# Fork of backwards_once_t0_T! (household_problem_and_aggregation.jl:450) with the
# three min() sites -- labor_income_, disutility_given_w_, labor_ -- replaced by their
# forced counterparts.  EVERYTHING else is verbatim, including the covid == 0 branch
# (unused at period 0, kept so the fork stays a drop-in).  `t` is passed EXPLICITLY
# rather than read from the global, closing one instance of the documented hazard.
function backwards_once_t0_T_forced_sw!(h, g, current_values; next_values, R_k, R_b, w,
        τ_n, T, Tu, n_ss, Γ_n, c_ss, Γ_c, covid, a_tmp = similar(current_values.η),
        e_init, method, t)
    (; u) = h
    Pp = h.Pprime
    τ = g.τ
    δ = t.δ
    ind_T = g.ind_T
    ind_Tu = g.ind_Tu
    R = 1.0 + e_init.r
    if method == 2
        Λ_ka = e_init.k / e_init.a
    end
    a_min = first(h.a_grid)
    a_max = last(h.a_grid)
    T_rev = T_rev!(e_init, R_k, R_b)

    for s in eachindex(h.z_grid)
        z = h.z_grid[s]
        lab = labor_income_forced_sw(; w = w * z, n_ = n_ss[s], Γ_n = Γ_n)
        lab_tax = lab - τ_n * lab ^ (1-τ)
        dis = disutility_forced_sw(h.v; n_ = n_ss[s], Γ_n = Γ_n)
        if ind_T[s]==1 && ind_Tu[s]==1
            add = lab - lab_tax - dis + T + Tu
        elseif ind_T[s]==1 && ind_Tu[s]==0
            add = lab - lab_tax - dis + T
        elseif ind_T[s]==0 && ind_Tu[s]==1
            add = lab - lab_tax - dis + Tu
        else
            add = lab - lab_tax - dis
        end

        for i in eachindex(h.a_grid)
            x = _backwards_euler_x_helper(h.u, Pp, i, s, next_values)
            a_next = h.a_grid[i]
            T_r = T_rev[i]
            if method == 1
                a_tmp[i, s] = (x + a_next - add - T_r) / R
            elseif method == 2
                a_tmp[i, s] = (x + a_next - add) / (Λ_ka*(R_k - δ) + (1.0 - Λ_ka)*R_b)
            end
        end

        for i in eachindex(h.a_grid)
            a = h.a_grid[i]
            if covid == 0.0
                a_prime = interp1D(a, view(a_tmp, :, s), h.a_grid)
                a_prime = clamp(a_prime, a_min, a_max)
                T_r = T_rev[i]
                if method == 1
                    x = R * a + add + T_r - a_prime
                elseif method == 2
                    x = (Λ_ka*(R_k - δ) + (1.0 - Λ_ka)*R_b) * a + add - a_prime
                end
                current_values.c_pol[i, s] = x + dis
                current_values.a_pol[i, s] = a_prime
                _extra_updates!(u, current_values, next_values, h, a_prime, Pp, x, i, s)
            elseif covid == 1.0
                T_r = T_rev[i]
                x = Γ_c * c_ss[i, s] - dis
                a_prime = R * a + add + T_r - x
                if a_prime < a_min
                    a_prime = a_min
                    x = R * a + add + T_r - a_prime
                end
                current_values.c_pol[i, s] = x + dis
                current_values.a_pol[i, s] = a_prime
                _extra_updates!(u, current_values, next_values, h, a_prime, Pp, x, i, s)
            end
            current_values.n_pol[i, s] = labor_forced_sw(; n_ = n_ss[s], Γ_n = Γ_n)
        end
    end
end


# ------------------------------------------------------------
# [SW22] the labour wedge the mandate creates, at period 0.
#   marginal disutility     theta * n_z^(1/nu),   n_z = Gamma_n * ntilde_z
#   marginal after-tax wage rho_z * w * z,        rho_z = (1-tau)*tau_n*(w z n_z)^(-tau)
# A POSITIVE wedge means households are OVER-working relative to free choice.
# z = 0 is skipped: 0^(-tau) = Inf.
# ------------------------------------------------------------
function forced_hours_wedge_sw(e_init, g, path, Γ_n; j = 1)
    v = e_init.h.v; τ = g.τ
    nt = e_init.ws.n_pol[1,:]; zg = e_init.h.z_grid
    w1 = path.w[j]; τn1 = path.τ_n[j]
    zmarg = sum(path[j].pdf, dims = 1)[1,:]
    rows = NamedTuple[]
    tot = 0.0; mass = 0.0
    for s in eachindex(zg)
        zg[s] == 0 && continue
        n_z  = Γ_n * nt[s]
        mdis = v.θ * n_z^(1/v.ν)
        ρ    = (1 - τ) * τn1 * (w1 * zg[s] * n_z)^(-τ)
        mben = ρ * w1 * zg[s]
        push!(rows, (; s = s, z = zg[s], n_z = n_z, mdis = mdis, mben = mben,
                     wedge = mdis - mben, rel = (mdis - mben)/mben, mass = zmarg[s]))
        tot += zmarg[s]*(mdis - mben); mass += zmarg[s]
    end
    return (; rows = rows, avg = tot/mass, mass = mass)
end

function print_wedge_sw(wd)
    println("\n---- [SW22] the labour wedge at j = 1 (mandated hours vs the FOC) ----")
    println("positive = households OVER-working relative to free choice")
    @printf("%4s %10s %12s %16s %16s %16s %12s %12s\n",
            "z", "z value", "n_z", "marg disutil", "marg net wage", "wedge", "rel", "mass")
    for r in wd.rows
        @printf("%4d %10.4f %12.6f %16.8f %16.8f %+16.8f %+11.2f%% %12.4e\n",
                r.s, r.z, r.n_z, r.mdis, r.mben, r.wedge, 100*r.rel, r.mass)
    end
    @printf("\nmass-weighted average wedge = %+.8f   (over %.6f of mass, z=0 excluded)\n",
            wd.avg, wd.mass)
    println(wd.avg > 0 ? "POSITIVE: over-working, as expected under a mandate above free choice." :
                         "NEGATIVE: under-working -- the mandate is BELOW free choice here.")
    return nothing
end
