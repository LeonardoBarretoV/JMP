# ============================================================
# [SW] Forked transition solver  --  STICKY-WAGE TRACK ONLY
# ============================================================
# Included by scripts/main_JMP_sw.jl ONLY. Never by scripts/main_JMP.jl.
#
# STAGE 2: plumbing only. The unknown vector is restructured to [r ; s] and the
# labor block is routed through union_block.jl, but the ECONOMICS IS UNCHANGED:
# the interior s residual pins s[j] to s_flex,j, the flexible-wage FOC ratio, so
# hours are exactly what the original solver produces. This stage is a pure
# regression test. The NKPC is NOT on.
#
# NAMESPACE RULE
#   Every forked name carries a _sw suffix. Checked against the 127 top-level
#   names in existing src/ files; no collisions. The originals must keep
#   resolving to their own file:line after this file loads.
#
# ORIGINALS FORKED HERE
#   solve_transition_h_r     transition.jl:384            -> solve_transition_h_r_sw
#   _residuals_h_r!          transition.jl:414            -> _residuals_h_r_sw!
#   _τ_n_excess_port_!       stationary_equilibria.jl:145 -> _τ_n_excess_port_sw_!
#   _τ_n_excess_t0_port_!    stationary_equilibria.jl:159 -> _τ_n_excess_t0_port_sw_!
#
# The two GBC helpers live in stationary_equilibria.jl, not transition.jl, and
# both call labor_supply_ (:148, :161) and taxes_HSV_ (:151, :164). Forking only
# the solver would leave this track routing through tau_n-dependent hours.
#
# NOT used as templates (all dead; the first pair calls three helpers that do
# not exist anywhere in src/): solve_transition_h, _residuals_h!,
# solve_transition_h_, solve_transition_h_l, solve_transition,
# solve_transition_slow_but_robust, transition_partial.
# ============================================================


# ------------------------------------------------------------
# s_flex: the flexible-wage FOC hours ratio relative to ntilde_z.
#
#   ntilde_z = A(tau_n_init, w_init) * z^((1-tau)nu/(1+tau*nu))
#   n_z(w, tau_n) = A(tau_n, w) * z^((1-tau)nu/(1+tau*nu))
#   => n_z(w, tau_n) / ntilde_z = A(tau_n,w)/A(tau_n_init,w_init), z-INDEPENDENT
#
# The z^(1-tau) factor cancels exactly, which is why proportional rationing spans
# the whole family of flexible-wage hours profiles. Measured z-independence of
# n_j(z)/ntilde_z along the original solution: 6.7e-16.
#
# Normalizing constants are read from e_init; nothing is hard-coded.
# ------------------------------------------------------------
s_flex_sw(g, ν; w, τ_n, w_init, τ_n_init) =
    ((τ_n / τ_n_init) * (w / w_init)^(1 - g.τ))^(ν / (1 + g.τ * ν))

# Nominal wage inflation. Terminal condition: piw at the LAST transition period
# is pi_bar, NOT the formula -- otherwise the j+1 term at the final interior
# index (j = 799) would read past the intended terminal block.
piw_at_sw(path, j, n_per, π_bar) = j == n_per ? π_bar : (1 + path.π_[j]) * path.w[j] / path.w[j-1] - 1


# ============================================================
# [SW] GLOBAL-STATE GUARD  and  ONE-RUN-PER-PROCESS GUARD
# ============================================================
# Two functions in existing src/ read GLOBALS rather than arguments:
#   time_series.jl:25                      make_time_series -> global g  (g.pi_bar)
#   household_problem_and_aggregation.jl:455  backwards_once_t0_T! -> global t  (t.delta)
# If two runs coexist in one Julia process, these resolve to whichever was
# assigned LAST, and make_time_series would silently build pi_path from the
# wrong pi_bar. That is the silent failure mode; this guard makes it loud.
_sw_run_count = Ref(0)
# [SW13] How many solves this process may perform.  Stays 1 by default, which is the
# original one-run-per-process rule.  A TWO-STAGE warm start (solve in s_mode = :flex,
# then use that as the warm start for :nkpc) needs 2, and it is safe there ONLY because
# both stages share the SAME e_init/t/g -- the hazard the guard exists for is two runs
# with DIFFERENT globals coexisting, which is not what a staged warm start does.
# Raise it deliberately in the script; never as a way past an unexpected guard trip.
_sw_runs_allowed = Ref(1)

function _assert_globals_sw(e_init, t, g)
    isdefined(Main, :t) || error("[SW] global `t` undefined; backwards_once_t0_T! (household_problem_and_aggregation.jl:455) reads it")
    isdefined(Main, :g) || error("[SW] global `g` undefined; make_time_series (time_series.jl:25) reads it")
    gt = getfield(Main, :t); gg = getfield(Main, :g)
    for (nm, a, b, c) in (("delta", gt.δ, t.δ, e_init.t.δ), ("A", gt.A, t.A, e_init.t.A), ("alpha", gt.α, t.α, e_init.t.α))
        (a == b == c) || error("[SW] technology mismatch on `", nm, "`: global=", a, " arg=", b, " e_init=", c)
    end
    for (nm, a, b, c) in (("pi_bar", gg.π_bar, g.π_bar, e_init.g.π_bar), ("tau", gg.τ, g.τ, e_init.g.τ), ("vartheta", gg.ϑ, g.ϑ, e_init.g.ϑ))
        (a == b == c) || error("[SW] tax-system mismatch on `", nm, "`: global=", a, " arg=", b, " e_init=", c)
    end
    (g.π_bar == e_init.π_) || error("[SW] g.pi_bar (", g.π_bar, ") != e_init.pi_ (", e_init.π_, ")")
    _sw_run_count[] += 1
    _sw_run_count[] <= _sw_runs_allowed[] || error("[SW] solve_transition_h_r_sw called ", _sw_run_count[],
        " times in this process, limit ", _sw_runs_allowed[],
        ". ONE RUN PER PROCESS: start a fresh julia for each counterfactual.")
    println("[SW] global-state guard OK (pi_bar=", g.π_bar, ", delta=", t.δ, "); run 1 of this process")
    return nothing
end

# ------------------------------------------------------------
# Forked GBC helper, interior periods.  Original: stationary_equilibria.jl:145
# Changes vs original: labor_supply_ -> labor_supply_union_(s),
#                      taxes_HSV_    -> w*n - τ_n*hsv_base_union_(s).
# Everything else verbatim. Γ_n dropped from the signature: hours are s*ntilde.
# ------------------------------------------------------------
function _τ_n_excess_port_sw_!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, yprime, d, s, n_ss, covid, z_dist, r, verbose, p)
    w = w_from_r_(t; r, τ_k, d, covid)
    n = labor_supply_union_(h; s, n_ss, z_dist)
    k = k_from_r_n_(t; r, n = n, τ_k, d, covid)
    y = t.A * (d * k)^t.α * n^(1 - t.α)
    Ξ = hsv_base_union_(h, g; w, s, n_ss, z_dist)
    lab_tax = w * n - τ_n * Ξ
    lump_sum_T = T * y
    lump_sum_Tu = Tu * y
    excess_g = lab_tax + τ_k * (r + t.δ) * k + qbprime * yprime - (1 + r) * qb * y - lump_sum_T - lump_sum_Tu - G * y
    verbose && next!(p, showvalues = [(:r, r), (:error, excess_g)])
    return (; excess_g, w, y, n, k, τ_n, τ_k, lab_tax, lump_sum_T, lump_sum_Tu, G)
end


# ------------------------------------------------------------
# Forked GBC helper, period 0.  Original: stationary_equilibria.jl:159
# Period 0 is otherwise UNTOUCHED: r_k/r_b split, backwards_once_t0_T!,
# residuals[1] and residuals[2] are all unchanged. s[1] is pinned to Γ_n, which
# reproduces min(Γ_n*ntilde_z, FOC) exactly because the min binds for all 14
# employed states there (verified: implied s_1 has min = max = 0.902301).
# ------------------------------------------------------------
function _τ_n_excess_t0_port_sw_!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, yprime, d, s, n_ss, covid, z_dist, r_b, r_k, verbose, p)
    w = w_from_r_k_(t; r_k, τ_k, d, covid)
    n = labor_supply_union_(h; s, n_ss, z_dist)
    k = k_from_r_k_n_(t; r_k, n = n, τ_k, d, covid)
    y = t.A * (d * k)^t.α * n^(1 - t.α)
    Ξ = hsv_base_union_(h, g; w, s, n_ss, z_dist)
    lab_tax = w * n - τ_n * Ξ
    lump_sum_T = T * y
    lump_sum_Tu = Tu * y
    excess_g = lab_tax + τ_k * r_k * k + qbprime * yprime - (1 + r_b) * qb - lump_sum_T - lump_sum_Tu - G
    verbose && next!(p, showvalues = [(:r_b, r_b), (:error, excess_g)])
    return (; excess_g, w, y, n, k, τ_n, τ_k, lab_tax, lump_sum_T, lump_sum_Tu, G)
end


# ------------------------------------------------------------
# Forked residual.  Original: transition.jl:414
#
# UNKNOWN VECTOR  x = [r_block ; s_block]                       length 1599
#
#   r_block  length 801 = length(tau_k_path) + 1   (UNCHANGED from the original)
#       [1]    r_k_hat  period-0 capital return net of delta
#       [2]    r_b      period-0 bond return
#       [i+1]  r_j      for j = 2..800  ->  indices 3..801
#
#   s_block  length 798 = length(tau_k_path) - 2
#       s is aligned 1:1 with `path`, so the full s vector has 800 entries.
#       s[1]   pinned to Gamma_n                          (first transition period)
#       s[800] pinned to s_flex(e_final.w, e_final.tau_n) (last transition period)
#       s_block[k] = s[k+1]  for k = 1..798
#
#   OFF-BY-ONE: r carries a +1 shift because r_path_rest holds TWO entries for
#   the first transition period (the MIT-shock r_k/r_b split), so r_path_rest[i+1]
#   is period i from i=2 on. Hours have no such split: one s per period, so
#   s[i] is period i with NO shift. That asymmetry is the trap.
#
# RESIDUALS  length 1599
#   [1]        k_agg - path.k[1]                     period-0 capital market
#   [2]        b_agg - path.b[1]                     period-0 bond market
#   [i]        a - k - q*bd, i = 3..801              asset clearing for period i-1
#   [801 + k]  s[k+1] - s_flex(w_{k+1}, tau_n_{k+1}) k = 1..798, NO shift
# ------------------------------------------------------------
const _SW9_XMIN = 1e-9                    # [SW9] floor on the GHH composite x
const _SW9_FLOOR = Threads.Atomic{Int}(0)  # times it fired (atomic: the loop is @batch-parallel)


# [SW8] Robust tau_n root-find.  `maximum(find_zeros(f, lo, hi))` THROWS on an empty
# result, which killed every counterfactual run: an intermediate Anderson iterate
# wanders somewhere the period budget has no root in the bracket and the whole
# solve dies instead of backing off.  Here we widen once, and if there is still no
# sign change we return the grid point that minimises |excess| -- a finite, well
# defined fallback that leaves a large asset residual, which is exactly the signal
# the solver needs to retreat.  At any true solution a root exists, so this changes
# nothing about the converged answer.
const _SW8_WIDEN = Ref(0)      # times the widened bracket was needed
const _SW8_FALLBACK = Ref(0)   # times NO root existed and the grid fallback fired
# [SW13] [SW8] KILL SWITCH.  Set _SW8_ENABLED[] = false to restore the ORIGINAL
# behaviour: `maximum(find_zeros(...))` on an empty result, which throws.  The
# widen-then-grid-argmin fallback makes the residual DISCONTINUOUS in the unknowns,
# and Anderson acceleration cannot handle that -- it is the prime suspect for the
# 5.6e20 / step-norm-3.1e37 blow-up seen while chasing the tau_k counterfactuals.
# Reverting trades a graceful degradation for a loud failure, which is what you
# want while diagnosing.  Default stays true so nothing else on the track changes.
const _SW8_ENABLED = Ref(true)
function _solve_tau_n_sw(f, lo, hi; widen = 4.0, ngrid = 96)
    r = find_zeros(f, lo, hi)
    isempty(r) || return maximum(r)
    # [SW13] reverted path: same ArgumentError the pre-[SW8] code raised.
    _SW8_ENABLED[] || return maximum(r)
    _SW8_WIDEN[] += 1
    hi2 = lo + widen * (hi - lo)
    r = find_zeros(f, lo, hi2)
    isempty(r) || return maximum(r)
    _SW8_FALLBACK[] += 1
    xs = range(lo + 1e-8, hi2; length = ngrid)
    best_x = first(xs); best_v = Inf
    for x in xs
        v = abs(f(x))
        if isfinite(v) && v < best_v; best_v = v; best_x = x; end
    end
    return best_x
end


# [SW4] length of the r block.  method == 1 substitutes r_b out via q_lag == 1,
# so the block holds [r_k_hat ; r_2 .. r_800] = n_per entries and the +1 shift is
# gone.  method == 2 keeps the original [r_k_hat ; r_b ; r_2 .. r_800] = n_per + 1.
_n_r_sw(n_per, method) = method == 1 ? n_per : n_per + 1

# ============================================================
# [SW19] EXACT ELIMINATION OF r_k AND residuals[1]
# ============================================================
# residuals[1] = k_agg - path.k[1], and under method == 1 k_agg = e_init.k exactly
# (see the k_agg assignment in the forward loop).  Period-0 capital is
# predetermined, so this is one equation in r_k -- and under the [SW] fork it is
# INVERTIBLE IN CLOSED FORM, because period-0 hours no longer depend on r_k:
#
#   k[1] = [ (alpha*A*dfac) / (mu*(1+tau_k)*r_k) ]^(1/(1-alpha)) * n
#   n    = labor_supply_union_(s = Gamma_n, n_ss, z_dist)   <- no w, no r_k, no tau_n
#   set k[1] = e_init.k:
#   r_k* = (alpha*A*dfac) / [ mu*(1+tau_k)*(e_init.k/n)^(1-alpha) ],  dfac = d^alpha if covid
#
# In the ORIGINAL model period-0 hours came from the household labour FOC and so
# depended on w, hence on r_k: genuinely simultaneous, which is why r_k had to be an
# unknown.  The fork replaced that FOC with proportional rationing at s*ntilde and
# pinned s_1 = Gamma_n, making n exogenous at period 0.  So this is the same move as
# retiring residuals[2] via q_lag == 1, and it is available for the same reason: a
# modelling change made an equation trivial and nobody revisited it.
#
# k is strictly decreasing in r_k (k ~ r_k^(-1/(1-alpha))), so residuals[1] is
# strictly increasing in r_k and the root is unique.
#
# Unknowns 1598 -> 1597: the r block becomes [r_2 .. r_n_per] and residuals[1] is
# removed rather than left near-zero.  OPT-IN via `sw19 = true`; every existing
# solver path is bit-identical when it is false.
# ============================================================
_n_r_sw(n_per, method, sw19::Bool) = sw19 ? _n_r_sw(n_per, method) - 1 : _n_r_sw(n_per, method)

# ============================================================
# [SW20] LIFT s_1 INTO THE OUTER UNKNOWN VECTOR
# ============================================================
# The lockdown is a CEILING on hours, not a requirement:
#     s_1 = min(Gamma_n, s_flex(w_1, tau_n_1))
# backwards_once_t0_T! ALREADY delivers min(Gamma_n*ntilde, FOC) and both branches
# are proportional to ntilde_z, so the HOUSEHOLD side has always been right.  The
# bug was the GOVERNMENT block booking Gamma_n unconditionally, so the two book
# different hours whenever s_flex_1 < Gamma_n (the item (b) gate then fails at j=1).
#
# This is NOT invertible: tau_n_1 comes from a root-find of the period-0 GBC.  It is
# a genuine scalar fixed point.  But given s_1 the whole chain is ONE FORWARD PASS --
# n_0 -> r_k* (closed form, [SW19]) -> w_1 (closed form) -> tau_n_1 (the scalar
# root-find every period already does) -> s_flex_1 -> residual -- so lifting s_1 into
# x lets Anderson see it instead of hiding a damped inner loop inside the residual.
# The earlier MINFIX attempt used a 60-iteration damped inner loop and that
# non-smoothness is what made Anderson diverge even where the branch never flipped.
#
# NOTE [SW19] does NOT make n_0 exogenous per se -- that holds only while s_1 is
# PINNED.  With s_1 lifted, n_0 = labor_supply_union_(s_1,...) depends on it, hence so
# do r_k* and w_1.  That second channel is STABILISING: lower s_1 -> lower n_0 ->
# lower r_k* -> HIGHER w_1 -> higher s_flex, so the map is decreasing in s_1.
#
# Layout: s_block[k] = s_vec[k] for k = 1..n_per-1, paired with residuals[n_r + k];
# k = 1 carries the min condition, k = 2..n_per-1 the NKPC.  1597 -> 1598.
# Where Gamma_n binds the residual is s_1 - Gamma_n, its Jacobian row is 1, the
# coordinate decouples and stays at Gamma_n -- so it is INERT by construction.
# ============================================================
_n_s_sw(n_per, sw20::Bool) = sw20 ? n_per - 1 : n_per - 2

function _rk_star_sw19(t, e_init; τ_k, n0, d, covid, r_min, r_max)
    dfac = covid == 1 ? d^t.α : one(d)
    r_k = (t.α * t.A * dfac) / (t.μ * (1 + τ_k) * (e_init.k / n0)^(1 - t.α))
    # If this fell outside the [SW10] clamp the clamp would silently override the
    # substitution and residuals[1] would NOT be zero -- so fail loudly instead.
    (r_min <= r_k - t.δ <= r_max) || error(
        "[SW19] r_k_hat* = ", r_k - t.δ, " lies outside the [SW10] clamp [",
        r_min, ", ", r_max, "]; the elimination would be silently overridden")
    return r_k
end


function _residuals_h_r_sw!(residuals, x; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, χ, t, g, method, kappa_w, wY_ss, nkpc_scale = 1.0, s_mode = :nkpc, π_over = nothing, sw19 = false, sw20 = false, s1_pin = nothing, forced_hours = false)   # [SW4] π_over: nothing => read π_path (Stage 3, bit-identical); else read the solved π vector
    π_bar = e_init.π_
    ī = (1 + e_init.r)*(1 + π_bar) - 1
    A = t.A
    α = t.α
    δ = t.δ
    β = h_path[1].u.β
    r_min = - δ + eps()
    r_max = 1 / β
    n_ss = e_init.ws.n_pol[1,:]
    c_ss = e_init.ws.c_pol

    n_per = length(τ_k_path)
    n_r   = _n_r_sw(n_per, method, sw19)   # [SW4] method==1 drops r_b; [SW19] also drops r_k
    n_s   = _n_s_sw(n_per, sw20)   # [SW20] s_1 lifted into x
    r_path_rest = view(x, 1:n_r)
    s_block     = view(x, (n_r + 1):(n_r + n_s))

    ν_ = h_path[1].v.ν
    θ_ = h_path[1].v.θ
    hh_agg = e_init.h    # θ, ν, z_grid, a_grid, u are identical across h_path; only P varies
    w_init = e_init.w
    τ_n_init = e_init.τ_n
    s_last = s_flex_sw(g, ν_; w = e_final.w, τ_n = e_final.τ_n, w_init = w_init, τ_n_init = τ_n_init)

    s_vec = Vector{eltype(x)}(undef, n_per)
    s_vec[n_per] = s_last
    if sw20
        @views s_vec[1:(n_per - 1)] .= s_block   # [SW20] s_1 is x, not Gamma_n
    else
        # [SW21] s1_pin lets the OUTER iteration place s_1 on the s_flex branch.
        # backwards_once_t0_T! still books min(Gamma_n, s_flex)*ntilde on the household
        # side, so a run is only VALID when s1_pin == min(Gamma_n, s_flex(w_1, tau_n_1))
        # at the solution -- that equality is the outer condition, checked ex post.
        s_vec[1] = isnothing(s1_pin) ? Γ_n : s1_pin
        @views s_vec[2:(n_per - 1)] .= s_block
    end

    z_state = zeros(length(h_path), length(e_init.h.z_grid))
    z_state[1,:] = transpose(e_init.h.Pss) * h_path[1].P
    for i = 2:lastindex(h_path)
        z_state[i,:] = transpose(z_state[i - 1,:]) * h_path[i].P
    end

    for i in reverse(eachindex(path, qb_path, τ_k_path))
        if i > firstindex(path)
            z_dist = z_state[i, :]
            h = i == lastindex(path) ? h_path[end] : h_path[i + 1]
            qb, qbstar, τ_k, T, Tu, G, π_, d, covid = qb_path[i], qbstar_path[i], τ_k_path[i], T_path[i], Tu_path[i], G_path[i], (π_over === nothing ? π_path[i] : π_over[i]), d_path[i], covid_path[i]
            qbprime = i == lastindex(path) ? e_final.q * e_final.b / output_(t; e_final.k, e_final.n, d, covid) : qb_path[i + 1]
            qbstar_prime = i == lastindex(path) ?  e_final.bstar / e_final.b : qbstar_path[i + 1]
            r = r_path_rest[sw19 ? i-1 : (method == 1 ? i : i+1)]   # [SW4] no +1 shift once r_b is gone; [SW19] r block starts at r_2
            r = clamp(r, r_min, r_max)
            R = 1 + r

            yprime = i == lastindex(path) ? output_(t; e_final.k, e_final.n, d, covid) : path.y[i + 1]
            qprime = i == lastindex(path) ? e_final.q : path.q[i + 1]
            rprime = i == lastindex(path) ? e_final.r : path.r[i + 1]
            π_prime = i == lastindex(path) ? e_final.π_ : path.π_[i + 1]
            q = (1/(1 + rprime))*((ī + g.ϑ)/(1 + π_prime) + (1 - g.ϑ)/(1 + π_prime) * qprime)
            bprime_level = (qbprime / q) * yprime
            q_lag = (1/(1 + r))*((ī + g.ϑ)/(1 + π_) + (1 - g.ϑ)/(1 + π_) * q)

            s_i = s_vec[i]
            p = ProgressUnknown()
            f = (τ_n) -> _τ_n_excess_port_sw_!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, yprime, d, s = s_i, n_ss, covid, z_dist, r, verbose = false, p).excess_g
            τ_n = _solve_tau_n_sw(f, 0.0, 2.0)   # [SW8]
            (; w, y, n, k, lab_tax, lump_sum_T, lump_sum_Tu) = _τ_n_excess_port_sw_!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, yprime, d, s = s_i, n_ss, covid, z_dist, r, verbose = false, p)
            X = i == lastindex(path) ? δ * k : path.k[i + 1] - (1.0 - δ) * k

            b_level = (qb / q_lag)* y
            path.G[i] = G * y
            path.l[i] = bprime_level - (1 - g.ϑ)/(1 + π_) * b_level
            path.h[i] = h
            path.q[i] = q
            path.π_[i] = π_
            path.τ_k[i] = τ_k
            path.τ_n[i] = τ_n
            path.r[i] = r
            path.w[i] = w
            path.k[i] = k
            path.b[i] = b_level
            path.bstar[i] = qbstar * path.b[i]
            path.Ψ[i] = path.bstar[i] / q_lag^(-χ)
            path.bd[i] = path.b[i] - path.bstar[i]
            path.n[i] = n
            path.y[i] = y
            path.X[i] = X
            path.C[i] = path.y[i] - path.X[i] - path.G[i] - (1 + r) * qbstar * qb * y + qbstar_prime * qbprime * yprime
            path.lab_tax[i] = lab_tax
            path.lump_sum_T[i] = lump_sum_T
            path.lump_sum_Tu[i] = lump_sum_Tu
            next_values = i == lastindex(path) ? e_final.ws : path[i+1]
            ind_T_ = e_init.g.ind_T[:] .== 1.0
            ind_Tu_ = e_init.g.ind_Tu[:] .== 1.0
            T_ = T * y / sum(z_dist[ind_T_])
            Tu_ = Tu * y / sum(z_dist[ind_Tu_])
            path.T[i] = T_
            path.Tu[i] =Tu_
            backwards_once_sw_!(h, g, path[i]; next_values = next_values, R, w, τ_n, T = T_, Tu = Tu_, s = s_i, n_ss, c_ss, Γ_c, covid, path[i].a_tmp)
            asset_policy_given_η_sw_!(path[i].a_pol; h, g, path[i].η, R, w, τ_n, T = T_, Tu = Tu_, s = s_i, n_ss)
            interpolate_asset_policy!(path[i].lower_index, path[i].lower_weight; h.a_grid, path[i].a_pol)
        end
    end

    forward_pdf!(path[1].pdf; h = h_path[1], e_init.ws.pdf, e_init.ws.lower_index, e_init.ws.lower_weight)
    h = h_path[2]
    z_dist = z_state[1,:]
    # [SW10] clamp r_k_hat exactly as the interior clamps r (line 237).  Without this,
    # an overshooting Anderson iterate can drive the gross capital return negative and
    # w_from_r_k_ raises a negative base to a fractional power.  Same bounds, same intent.
    if sw19
        # [SW19] r_k is no longer an unknown: period-0 hours are independent of it, so
        # the period-0 capital market inverts in closed form.  n0 is recomputed here
        # rather than read from the GBC helper because it must be known BEFORE r_k.
        n0_sw19 = labor_supply_union_(h; s = s_vec[1], n_ss = n_ss, z_dist = z_dist)
        # READ tau_k / d / covid STRAIGHT FROM THE PATHS, not from the local names.
        # `tau_k`, `d` and `covid` are function-locals that the interior loop also
        # assigns, and they are not (re)set to their period-0 values until the
        # `qb, qbstar, tau_k, ... = qb_path[1], ...` line BELOW this block.  Reading
        # them here picked up the loop's last iteration (period 2), which made r_k*
        # a factor (1-0.4752)/(1+0.1379) too small and k[1] 3.17x too large.  Indexing
        # the paths directly makes this independent of statement order.
        r_k = _rk_star_sw19(t, e_init; τ_k = τ_k_path[1], n0 = n0_sw19,
                            d = d_path[1], covid = covid_path[1], r_min, r_max)
        r_k_hat = r_k - t.δ
    else
        r_k_hat = clamp(r_path_rest[1], r_min, r_max)
        r_k = r_k_hat + t.δ
    end
    r_b = method == 1 ? 0.0 : r_path_rest[2]   # [SW4] method==1: overwritten below from q_lag == 1
    R_k = 1 + r_k
    R_b = 1 + r_b
    qb, qbstar, τ_k, T, Tu, G, π_, d, covid = qb_path[1], qbstar_path[1], τ_k_path[1], T_path[1], Tu_path[1], G_path[1], (π_over === nothing ? π_path[1] : π_over[1]), d_path[1], covid_path[1]
    qbprime = qb_path[2]
    qbstar_prime = qbstar_path[2]
    qprime = path.q[2]
    yprime = path.y[2]

    rprime = path.r[2]
    π_prime = path.π_[2]
    q = (1/(1 + rprime))*((ī + g.ϑ)/(1 + π_prime) + (1 - g.ϑ)/(1 + π_prime) * qprime)
    bprime_level = (qbprime / q) * yprime
    # ---- [SW4] SUBSTITUTION: q_lag == 1 pins r_b in CLOSED FORM ----------------
    # q_lag is the PRE-announcement bond price, an initial-steady-state fact, so it
    # is normalised to 1.  r_b enters only as the factor 1/(1+r_b) multiplying a term
    # free of r_b, so the condition inverts exactly:
    #        1 + r_b = [ (ibar + vartheta) + (1 - vartheta) * q_1 ] / (1 + pi_1)
    # q_1 is already computed above and r_b is not consumed until the tau_n solve
    # below, so this is well ordered.  Under method == 1, b_agg = e_init.b is a
    # constant and e_init.b == qb_path[1] BITWISE (verified), so residuals[2] =
    # e_init.b - qb/q_lag vanishes identically; the row is removed from the system
    # rather than left degenerate.  Under method == 2, b_agg comes from path[1].pdf,
    # residuals[2] is genuine market clearing, and r_b stays an unknown.
    local q_lag
    if method == 1
        r_b   = ((ī + g.ϑ)/(1 + π_) + (1 - g.ϑ)/(1 + π_) * q) - 1.0
        R_b   = 1 + r_b
        q_lag = 1.0
    else
        q_lag = (1/(1 + r_b))*((ī + g.ϑ)/(1 + π_) + (1 - g.ϑ)/(1 + π_) * q)
    end
    println("q_lag = ", q_lag); flush(stdout)   # [SW] flush: stdout is block-buffered when redirected to a file
    s_0 = s_vec[1]    # [SW] period-0 hours = Gamma_n (pre-fix). The Gamma_n branch is active
    # at kappa_w = 0.05, so both blocks use the same hours and item (b) reconciles.
    p = ProgressUnknown()
    f = (τ_n) -> _τ_n_excess_t0_port_sw_!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, yprime, d, s = s_0, n_ss, covid, z_dist, r_b, r_k, verbose = false, p).excess_g
    τ_n = _solve_tau_n_sw(f, 0.0, 5.0)   # [SW8]
    (; w, y, n, k, lab_tax, lump_sum_T, lump_sum_Tu) = _τ_n_excess_t0_port_sw_!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, yprime, d, s = s_0, n_ss, covid, z_dist, r_b, r_k, verbose = false, p)
    # [SW20] the lifted s_1 residual.  `w` here IS w_1 = w_from_r_k_(r_k, tau_k, d,
    # covid), so s_flex_1 is a closed-form function of the quantities just computed --
    # one forward pass, no inner loop.  Where Gamma_n binds this reads s_1 - Gamma_n
    # and is satisfied exactly at s_1 = Gamma_n, so the row is inert by construction.
    if sw20
        s_flex_1 = s_flex_sw(g, ν_; w = w, τ_n = τ_n, w_init = w_init, τ_n_init = τ_n_init)
        residuals[n_r + 1] = s_vec[1] - min(Γ_n, s_flex_1)
    end
    X = path.k[2] - (1.0 - δ) * k

    b_level = qb / q_lag
    path.G[1] = G
    path.l[1] = bprime_level - (1 - g.ϑ)/(1 + π_) * b_level
    path.h[1] = h
    path.q[1] = q
    path.π_[1] = π_
    path.τ_k[1] = τ_k
    path.τ_n[1] = τ_n
    path.r[1] = r_b
    path.r_k[1] = r_k
    path.w[1] = w
    path.k[1] = k
    path.b[1] = b_level
    path.bstar[1] = qbstar / q_lag
    path.Ψ[1] = path.bstar[1] / q_lag^(-χ)
    path.bd[1] = path.b[1] - path.bstar[1]
    path.n[1] = n
    path.y[1] = y
    path.X[1] = X
    path.C[1] = path.y[1] - path.X[1] - path.G[1] - (1 + r_b) * qbstar + qbstar_prime * qbprime * yprime
    path.lab_tax[1] = lab_tax
    path.lump_sum_T[1] = lump_sum_T
    path.lump_sum_Tu[1] = lump_sum_Tu
    next_values = path[2]
    ind_T_ = e_init.g.ind_T[:] .== 1.0
    ind_Tu_ = e_init.g.ind_Tu[:] .== 1.0
    T_ = T * y / sum(z_dist[ind_T_])
    Tu_ = Tu * y / sum(z_dist[ind_Tu_])
    path.T[1] = T_
    path.Tu[1] = Tu_
    # [SW22] forced_hours: the household block books Gamma_n*ntilde_z EXACTLY rather than
    # min(Gamma_n, s_flex)*ntilde_z, so the lockdown is a MANDATE not a ceiling.  The
    # government block already books Gamma_n, so the two then agree by construction.
    if forced_hours
        backwards_once_t0_T_forced_sw!(h, g, path[1]; next_values = next_values, R_k, R_b, w, τ_n, T = T_, Tu = Tu_, n_ss, Γ_n, c_ss, Γ_c, covid, path[1].a_tmp, e_init, method, t)
    else
        backwards_once_t0_T!(h, g, path[1]; next_values = next_values, R_k, R_b, w, τ_n, T = T_, Tu = Tu_, n_ss, Γ_n, c_ss, Γ_c, covid, path[1].a_tmp, e_init, method)
    end
    interpolate_asset_policy!(path[1].lower_index, path[1].lower_weight; h.a_grid, path[1].a_pol)

    if method ==1
        k_agg = e_init.k
        b_agg = e_init.b
    elseif method == 2
        Λ_ka = e_init.k / e_init.a
        k_agg = sum((Λ_ka.*h.a_grid)  .* path[1].pdf)
        b_agg = sum(((1.0 - Λ_ka).*h.a_grid) .* path[1].pdf)
    end
    for i in 1:length(path.r)+1
        if i==1
            path.un[1] = sum(path[1].pdf, dims = 1)[1]
            a = asset_supply(h_path[i].a_grid, path[i].pdf)
            path.a[i] = a
            # [SW19] identically 0 when sw19: r_k was chosen to make k[1] == k_agg == e_init.k.
            # The row is REMOVED from the system rather than left at machine zero.
            sw19 || (residuals[1] =  k_agg - path.k[1])
        elseif i==2
            method == 1 || (residuals[2] = b_agg - path.b[1])   # [SW4] identically 0 when method==1; row removed
        else
            a = asset_supply(h_path[i-1].a_grid, path[i-1].pdf)
            path.a[i-1] = a
            residuals[sw19 ? i-2 : (method == 1 ? i-1 : i)] =  a - path.k[i-1] - path.q[i - 2] * path.bd[i-1]   # [SW4] index shift; [SW19] one further

            # ---- [SW] wage Phillips curve residual for period j = i-1 ----
            # In the FORWARD loop because D and Y need path[j].pdf, which only
            # exists after forward iteration j. Period j sits at i = j+1, next
            # to the asset residual for the same j.
            j_nk = i - 1
            if 2 <= j_nk <= n_per - 1
                piw_j   = piw_at_sw(path, j_nk,     n_per, π_bar)
                piw_jp1 = piw_at_sw(path, j_nk + 1, n_per, π_bar)
                DYj = union_aggregates_(hh_agg, g, path[j_nk]; s = s_vec[j_nk], w = path.w[j_nk], τ_n = path.τ_n[j_nk], n_ss = n_ss)
                wedge = s_vec[j_nk] * (θ_ * s_vec[j_nk]^(1/ν_) * DYj.D - path.w[j_nk] * DYj.Y)
                Rj = (piw_j - π_bar)*(1 + piw_j) - (kappa_w / wY_ss) * wedge - β * (piw_jp1 - π_bar)*(1 + piw_jp1)
                if s_mode === :flex
                    # verification (a): pin interior s to s_flex, forked household block.
                    # Must reproduce Stage 2 to solver tolerance -- isolates the fork from the NKPC.
                    residuals[n_r + (sw20 ? j_nk : j_nk - 1)] = s_vec[j_nk] - s_flex_sw(g, ν_; w = path.w[j_nk], τ_n = path.τ_n[j_nk], w_init = w_init, τ_n_init = τ_n_init)
                else
                    residuals[n_r + (sw20 ? j_nk : j_nk - 1)] = -nkpc_scale * Rj    # leading minus: match asset-block Jacobian sign
                end
            end
        end
        (i == length(path.r)+1) && break

        if i > 1
            forward_pdf!(path[i].pdf; h = h_path[i], path[i-1].pdf, path[i-1].lower_index, path[i-1].lower_weight)
            path.un[i] = sum(path[i].pdf, dims = 1)[1]
        end
    end

    return nothing
end


# ------------------------------------------------------------
# Forked driver.  Original: transition.jl:384
# Returns (path, sol) so the caller can report iterations / residual norm.
# ------------------------------------------------------------
function solve_transition_h_r_sw(
    e_init, e_final; τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, χ, t, g, method, kappa_w, nkpc_scale = 1.0, s_mode = :nkpc, r_path = nothing, s_path = nothing, verbose = true,
    nlsolve_kwargs = nothing
)
    nlsolve_baseline_kwargs = (
        ftol = _ZERO_FTOL,
        show_trace = verbose,
        method = :anderson,
        m = 10,
        iterations = _ZERO_MAX_ITERS,
        beta = -0.01
    )

    nlsolve_kwargs_merged = isnothing(nlsolve_kwargs) ?
        nlsolve_baseline_kwargs :
        merge(nlsolve_baseline_kwargs, nlsolve_kwargs)

    _assert_globals_sw(e_init, t, g)
    path = _initialize_path(e_init, e_final; qb_path)

    n_per = length(τ_k_path)
    n_r   = _n_r_sw(n_per, method)   # [SW4]
    n_s   = n_per - 2

    r0 = isnothing(r_path) ? fill(e_final.r, n_r) : collect(r_path)
    ν_ = h_path[1].v.ν
    s_last = s_flex_sw(g, ν_; w = e_final.w, τ_n = e_final.τ_n, w_init = e_init.w, τ_n_init = e_init.τ_n)
    DY_ss = union_aggregates_(e_init.h, g, e_init.ws; s = 1.0, w = e_init.w, τ_n = e_init.τ_n, n_ss = e_init.ws.n_pol[1,:])
    wY_ss = e_init.w * DY_ss.Y    # computed in-solver, not hard-coded
    println("[SW] w_ss*Y_ss = ", wY_ss, "   kappa_w = ", kappa_w, "   nkpc_scale = ", nkpc_scale)
    s0 = isnothing(s_path) ? fill(s_last, n_s) : collect(s_path)

    length(r0) == n_r || error("r_path must have length $(n_r), got $(length(r0))")
    length(s0) == n_s || error("s_path must have length $(n_s), got $(length(s0))")

    x0 = vcat(r0, s0)
    f! = (F, x) -> _residuals_h_r_sw!(F, x; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, χ, t, g, method, kappa_w, wY_ss, nkpc_scale, s_mode)
    sol = nlsolve(f!, x0; nlsolve_kwargs_merged...)
    Ffin = similar(sol.zero)
    _residuals_h_r_sw!(Ffin, sol.zero; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, χ, t, g, method, kappa_w, wY_ss, nkpc_scale, s_mode)
    Ra = view(Ffin, 1:n_r)
    Rn = view(Ffin, (n_r+1):length(Ffin))
    println("[SW] max|R_asset| = ", maximum(abs, Ra), "   max|R_s/nkpc| = ", maximum(abs, Rn))
    println("[SW] med|R_asset| = ", median(abs.(Ra)), "   med|R_s/nkpc| = ", median(abs.(Rn)))
    println("[SW] block ratio med|R_asset|/med|R_s/nkpc| = ", median(abs.(Ra))/median(abs.(Rn)))
    return path, sol
end


# ============================================================
# [SW] Forked household block -- hours are RATIONED at s_j * ntilde_z
# ============================================================
# Originals: backwards_once_!        household_problem_and_aggregation.jl:317
#            asset_policy_given_η_!  household_problem_and_aggregation.jl:368
#
# WHY THIS FORK EXISTS
#   Once the NKPC is on, s_j != s_flex,j by construction. The originals call
#   labor_income_ / disutility_given_w_, which on the non-covid branch return the
#   flexible-wage FOC. Households would then face FOC hours while the government
#   block booked s_j*ntilde_z hours: government revenue and household tax payments
#   would not reconcile, and union_aggregates_ would build lambda from a c_pol
#   generated under different hours than the disutility subtracted from it.
#
# WHY NOT REUSE THE min(Γ_n*n_, FOC) PATH
#   Passing Γ_n := s into labor_ gives min(s*ntilde_z, FOC). When s > s_flex --
#   precisely the high-inflation quarters of interest -- the min clamps back to
#   the FOC and the rationing silently stops binding. It would fail exactly where
#   the economics matters, with no error.
#
# PERIOD 0 IS NOT TOUCHED. backwards_once_t0_T! stays as-is: there s[1] = Γ_n and
# the min genuinely binds for all 14 employed states (verified, implied s_1 has
# min = max = 0.902301).
#
# NOTE the loop variable over productivity is `sx`, not `s`: `s` is now the
# scalar hours index.
# ============================================================

function backwards_once_sw_!(h, g, current_values; next_values, R, w, τ_n, T, Tu, s, n_ss, c_ss, Γ_c, covid, a_tmp = similar(current_values.η))
    η_now = current_values.η
    (; u) = h
    P′ = h.Pprime
    ξ = get_inverse_ies(u)
    R_ = R^(-1/ξ)
    τ = g.τ
    ind_T = g.ind_T
    ind_Tu = g.ind_Tu
    θ = h.v.θ
    ν = h.v.ν

    a_min = first(h.a_grid)
    a_max = last(h.a_grid)

    @batch for sx in eachindex(h.z_grid)
        z = h.z_grid[sx]
        n_iz = s * n_ss[sx]
        lab  = w * z * n_iz
        lab_tax = lab - τ_n * lab ^ (1-τ)
        dis  = θ * n_iz^(1 + 1/ν) / (1 + 1/ν)
        if ind_T[sx]==1 && ind_Tu[sx]==1
            add = lab - lab_tax - dis + T + Tu
        elseif ind_T[sx]==1 && ind_Tu[sx]==0
            add = lab - lab_tax - dis + T
        elseif ind_T[sx]==0 && ind_Tu[sx]==1
            add = lab - lab_tax - dis + Tu
        else
            add = lab - lab_tax - dis
        end

        for i in eachindex(h.a_grid)
            x = _backwards_euler_x_helper(h.u, P′, i, sx, next_values)
            a_next = h.a_grid[i]
            a_tmp[i, sx]  = (x + a_next - add) / R
        end

        for i in eachindex(h.a_grid)
            a = h.a_grid[i]
            a_prime = interp1D(a, view(a_tmp, :, sx), h.a_grid)
            a_prime = clamp(a_prime, a_min, a_max)
            if covid == 0.0
                x = R * a + add - a_prime
            elseif covid == 1.0
                x = Γ_c * c_ss[i, sx] - dis
            end
            # [SW9] floor the GHH composite.  x <= 0 makes log(x) throw and kills the whole
            # solve at an intermediate Anderson iterate.  Flooring returns a large-but-finite
            # utility instead, so the solver can retreat.  At any true solution x > 0 for every
            # household, so this must fire ZERO times at a converged point -- counted and reported.
            if x <= _SW9_XMIN
                Threads.atomic_add!(_SW9_FLOOR, 1)
                x = _SW9_XMIN
            end
            current_values.c_pol[i, sx] = x + dis
            current_values.n_pol[i, sx] = n_iz
            η_now[i, sx] = R_ * x
            _extra_updates!(u, current_values, next_values, h, a_prime, P′, x, i, sx)
        end
    end
end


function asset_policy_given_η_sw_!(a_pol; h, g, η, R, w, τ_n, T, Tu, s, n_ss)
    ξ = get_inverse_ies(h.u)
    τ = g.τ
    ind_T = g.ind_T
    ind_Tu = g.ind_Tu
    θ = h.v.θ
    ν = h.v.ν
    @batch for i in eachindex(h.a_grid)
        for sx in eachindex(h.z_grid)
            n_iz = s * n_ss[sx]
            lab  = w * h.z_grid[sx] * n_iz
            dis  = θ * n_iz^(1 + 1/ν) / (1 + 1/ν)
            base = R * h.a_grid[i] + τ_n * lab^(1-τ) - η[i, sx] * R^(1/ξ) - dis
            if ind_T[sx]==1 && ind_Tu[sx]==1
                a_pol[i, sx] = base + T + Tu
            elseif ind_T[sx]==1 && ind_Tu[sx]==0
                a_pol[i, sx] = base + T
            elseif ind_T[sx]==0 && ind_Tu[sx]==1
                a_pol[i, sx] = base + Tu
            else
                a_pol[i, sx] = base
            end
        end
    end
end


# ============================================================
# [SW] make_time_series_sw -- fork of make_time_series (time_series.jl:4)
# ============================================================
# TWO changes vs the original, both first-transition-period debt targets:
#   b_path[2]        4.9940 -> 5.1949 -> 5.2535   (two Newton steps)
#   bstar_path_r[2]  0.2660 -> 0.2771 -> 0.2694  (step1 overshot; qbstar is Bstar-over-B)
# Newton step new = old * data/model, from the measured 2020Q2 gaps:
#   debt-to-GDP      data 132.815  model 127.678  (-5.137 pp)
#   ext debt-to-GDP  data  35.375  model  33.962  (-1.412 pp)
# NOTE the figure plots q(t)*b(t+1)/(4*y(t)), so the point at date t is driven by
# qb_path[t+1]. b_path[1] and bstar_path_r[1] are LEVELS and drive q_lag = 1 and
# hence r_b; they do NOT move any plotted debt point. Left untouched.
#
# Also takes `g` as an explicit ARGUMENT. The original reads a global g at
# time_series.jl:25 (g.π_bar), one of the two global leaks documented in
# README_sw.md. This fork closes it for the _sw track.
# ============================================================
function make_time_series_sw(c, e_init, e_final, b_path_in, Λ_bstar_f, g)
    Tsim = length(b_path_in)

    Pee_path = fill(c.Pee, Tsim)
    Pee_path[1] = 0.89
    Pee_path[2] = 0.96
    Pee_path[3] = 0.962
    Pee_path[3] = 0.964
    Pee_path[4] = 0.962
    Pee_path[5] = 0.964
    Pee_path[6] = 0.971

    Puu_path = fill(c.Puu, Tsim)
    Puu_path[1] = 0.55

    π_path = fill(g.π_bar, Tsim)
    π_path[1] = (1.0065)^0.25-1; π_path[2] = (1.0137)^0.25-1; π_path[3] = (1.0136)^0.25-1;
    π_path[4] = (1.0262)^0.25-1; π_path[5] = (1.0539)^0.25-1; π_path[6] = (1.0539)^0.25-1; π_path[7] = (1.0704)^0.25-1;
    π_path[8] = (1.0854)^0.25-1; π_path[9] = (1.0906)^0.25-1; π_path[10] = (1.0820)^0.25-1; π_path[11] = (1.0645)^0.25-1;
    π_path[12] = (1.0498)^0.25-1; π_path[13] = (1.0297)^0.25-1; π_path[14] = (1.0370)^0.25-1; π_path[15] = (1.0335)^0.25-1;
    π_path[16] = (1.0348)^0.25-1; π_path[17] = (1.0297)^0.25-1; π_path[18] = (1.0244)^0.25-1; π_path[19] = (1.0289)^0.25-1;

    b_path = fill(e_final.q * e_final.b, Tsim) ./ y(e_final)
    b_path[1] = e_init.b ; b_path[2] = 5.2535; b_path[3] = 1.26 * 4.0      # [SW] b_path[2] raised
    b_path[4] = 1.24 * 4.0; b_path[5] = 1.22 * 4.0; b_path[6] = 1.19 * 4.0 ; b_path[7] = 1.20 * 4.0
    b_path[8] = 1.21 * 4.0; b_path[9] = 1.18 * 4.0; b_path[10] = 1.18 * 4.0; b_path[11] = 1.18 * 4.0
    b_path[12] = 1.16 * 4.0; b_path[13] = 1.177 * 4.0; b_path[14] = 1.186 * 4.0; b_path[15] = 1.20 * 4.0
    b_path[16] = 1.208 * 4.0; b_path[17] = 1.20 * 4.0; b_path[18] = 1.207 * 4.0; b_path[19] = 1.218 * 4.0
    b_path[20] = 1.208 * 4.0

    covid_path = fill(c.covid, Tsim)
    covid_path[1] = 1.0
    d_path = fill(c.d_, Tsim)
    Γ_n = 0.9023011517193031
    Γ_c = 0.9229604357071445

    T_path = fill(e_init.T, Tsim)
    T_path[1] = 1.4 * 4 * 0.051

    Tu_path = fill(e_init.lump_sum_Tu, Tsim) ./ y(e_init)
    Tu_path[1] = 1.4 * 4 * 0.046 + 4 * 0.005

    τ_k_path = fill(e_init.τ_k, Tsim)
    τ_k_path[1] = -0.47

    G_path = fill(e_init.G, Tsim) ./ y(e_init)
    G_path[1]  = e_init.G * 1.1; G_path[2]  = 0.0533;   G_path[3]  = 0.0526;
    G_path[4]  = 0.0548;    G_path[5]  = 0.0523;   G_path[6]  = 0.0503;
    G_path[7]  = 0.0494;    G_path[8]  = 0.0481;   G_path[9]  = 0.0471;
    G_path[10] = 0.0470;    G_path[11] = 0.0474;   G_path[12] = 0.0478;
    G_path[13] = 0.0473;    G_path[14] = 0.0478;   G_path[15] = 0.0476;
    G_path[16] = 0.0478;    G_path[17] = 0.0478;   G_path[18] = 0.0484;
    G_path[19] = 0.0487;

    bstar_path = fill(Λ_bstar_f, Tsim)
    bstar_path[1] = e_init.bstar ; bstar_path[2] = 0.326 * 4.0; bstar_path[3] = 0.320 * 4.0
    bstar_path[4] = 0.310 * 4.0; bstar_path[5] = 0.322 * 4.0; bstar_path[6] = 0.316 * 4.0 ; bstar_path[7] = 0.312 * 4.0
    bstar_path[8] = 0.30 * 4.0; bstar_path[9] = 0.29 * 4.0; bstar_path[10] = 0.276 * 4.0; bstar_path[11] = 0.274 * 4.0
    bstar_path[12] = 0.275 * 4.0; bstar_path[13] = 0.275 * 4.0; bstar_path[14] = 0.268 * 4.0; bstar_path[15] = 0.28 * 4.0
    bstar_path[16] = 0.283 * 4.0; bstar_path[17] = 0.284 * 4.0; bstar_path[18] = 0.296 * 4.0; bstar_path[19] = 0.287 * 4.0

    bstar_path_r = fill(e_final.bstar, Tsim) ./ e_final.b
    bstar_path_r[1]  = e_init.bstar;
    bstar_path_r[2]  = 0.2694;   bstar_path_r[3]  = 0.262;    bstar_path_r[4]  = 0.255;    bstar_path_r[5]  = 0.250;   # [SW] [2] raised
    bstar_path_r[6]  = 0.264;    bstar_path_r[7]  = 0.266;    bstar_path_r[8]  = 0.261;    bstar_path_r[9]  = 0.250;
    bstar_path_r[10] = 0.243;    bstar_path_r[11] = 0.234;    bstar_path_r[12] = 0.233;    bstar_path_r[13] = 0.237;
    bstar_path_r[14] = 0.234;    bstar_path_r[15] = 0.226;    bstar_path_r[16] = 0.234;    bstar_path_r[17] = 0.234;
    bstar_path_r[18] = 0.237;    bstar_path_r[19] = 0.246;    bstar_path_r[20] = 0.236;

    return (
        Pee_path = Pee_path, Puu_path = Puu_path, π_path = π_path, qb_path = b_path,
        covid_path = covid_path, d_path = d_path, Γ_n = Γ_n, Γ_c = Γ_c,
        T_path = T_path, Tu_path = Tu_path, τ_k_path = τ_k_path, G_path = G_path,
        bstar_path = bstar_path, bstar_path_r = bstar_path_r,
    )
end


# ============================================================
# [SW4] STAGE 4 STEP 2 -- 3x3 joint solve on (r, s, pi)
# ============================================================
# UNKNOWNS   x = [ r_block(801) ; s_block(798) ; pi_block(798) ]   length 2397
#   r_block   [1] r_k_hat, [2] r_b, [i+1] r_i for i=2..800     (+1 SHIFT)
#   s_block   s_block[k] = s[k+1], k=1..798                    (NO shift)
#   pi_block  pi_block[k] = pi[k+1], k=1..798                  (NO shift)
#             pi[1] FED from data (NOT an unknown) <- the j=2 guard
#             pi[800] pinned at pi_bar
#
# RESIDUALS                                                    length 2397
#   [1..801]        asset block   (unchanged, from _residuals_h_r_sw!)
#   [802..1599]     NKPC at period k+1  (unchanged)
#   [1600..2397]    monetary at period j=k+1, k=1..798:
#                     (1+r_{j+1}) - (1+i_j)/(1+pi_{j+1})
#                     i_j = i_ss + phi_pi*(pi_j - pi_bar) + epsm_j
#
# The first monetary residual (k=1, j=2) ties pi_2, pi_3, r_3.  The earliest r
# it touches is r_3 -- two periods clear of the period-0 r_k/r_b split, so no
# instance ever involves r_1 (= r_b, pinned by q_lag = 1 rather than a policy
# rate).  The dlogw boundary at j=2 enters only through the NKPC via
# piw_2 = (1+pi_2)*w_2/w_1 - 1, which reads the stored w_1 directly and never
# the -alpha/(1-alpha) r-based identity.
# ------------------------------------------------------------
function _residuals_sw3!(residuals, x; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, χ, t, g, method, kappa_w, wY_ss, epsm_path, i_taylor_ss, phi_pi, nkpc_scale = 1.0, mon_scale = 1.0, s_mode = :nkpc)
    n_per = length(τ_k_path)
    n_r   = _n_r_sw(n_per, method)   # [SW4]
    n_s   = n_per - 2
    n_p   = n_per - 2
    π_bar = e_init.π_

    # ---- assemble the full pi vector: pi_1 fed, pi_800 pinned, interior solved
    π_vec = Vector{eltype(x)}(undef, n_per)
    π_vec[1]     = π_path[1]        # FED -- not an unknown
    π_vec[n_per] = π_bar            # pinned
    @inbounds for k in 1:n_p
        π_vec[k + 1] = x[n_r + n_s + k]
    end

    # ---- blocks 1 and 2: asset + NKPC, with pi read from the solved vector
    _residuals_h_r_sw!(residuals, x; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, χ, t, g, method, kappa_w, wY_ss, nkpc_scale, s_mode, π_over = π_vec)

    # ---- block 3: monetary / Fisher, periods j = 2..799
    @inbounds for k in 1:n_p
        j   = k + 1
        i_j = i_taylor_ss + phi_pi * (π_vec[j] - π_bar) + epsm_path[j]
        residuals[n_r + n_s + k] = mon_scale * ((1 + path.r[j + 1]) - (1 + i_j)/(1 + π_vec[j + 1]))
    end
    _SW4_CALLS[] += 1
    if _SW4_REPORT_EVERY[] > 0 && _SW4_CALLS[] % _SW4_REPORT_EVERY[] == 0
        _report_blocks_sw3("call " * string(_SW4_CALLS[]), residuals, n_r, n_s, n_p)
    end
    return nothing
end


const _SW4_CALLS = Ref(0)          # [SW4] residual-call counter, for mid-solve block reports
const _SW4_REPORT_EVERY = Ref(25)  # [SW4] print the 3-block split every N calls; 0 disables


function solve_transition_sw3(
    e_init, e_final; τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, χ, t, g, method, kappa_w, epsm_path, phi_pi,
    i_taylor_ss = nothing, nkpc_scale = 1.0, mon_scale = 1.0, s_mode = :nkpc,
    r_path = nothing, s_path = nothing, pi_path0 = nothing, verbose = true, nlsolve_kwargs = nothing
)
    nlsolve_baseline_kwargs = (ftol = _ZERO_FTOL, show_trace = verbose, method = :anderson, m = 10, iterations = _ZERO_MAX_ITERS, beta = -0.01)
    nlsolve_kwargs_merged = isnothing(nlsolve_kwargs) ? nlsolve_baseline_kwargs : merge(nlsolve_baseline_kwargs, nlsolve_kwargs)

    _assert_globals_sw(e_init, t, g)
    _SW4_CALLS[] = 0
    path = _initialize_path(e_init, e_final; qb_path)

    n_per = length(τ_k_path); n_r = _n_r_sw(n_per, method); n_s = n_per - 2; n_p = n_per - 2   # [SW4]
    π_bar = e_init.π_
    i_ss  = isnothing(i_taylor_ss) ? (1 + e_final.r)*(1 + π_bar) - 1 : i_taylor_ss

    ν_ = h_path[1].v.ν
    s_last = s_flex_sw(g, ν_; w = e_final.w, τ_n = e_final.τ_n, w_init = e_init.w, τ_n_init = e_init.τ_n)
    DY_ss  = union_aggregates_(e_init.h, g, e_init.ws; s = 1.0, w = e_init.w, τ_n = e_init.τ_n, n_ss = e_init.ws.n_pol[1,:])
    wY_ss  = e_init.w * DY_ss.Y

    r0  = isnothing(r_path)   ? fill(e_final.r, n_r) : collect(r_path)
    s0  = isnothing(s_path)   ? fill(s_last, n_s)    : collect(s_path)
    p0  = isnothing(pi_path0) ? collect(π_path[2:(n_per-1)]) : collect(pi_path0)

    length(r0) == n_r || error("r_path must have length $(n_r), got $(length(r0))")
    length(s0) == n_s || error("s_path must have length $(n_s), got $(length(s0))")
    length(p0) == n_p || error("pi_path0 must have length $(n_p), got $(length(p0))")
    length(epsm_path) == n_per || error("epsm_path must have length $(n_per), got $(length(epsm_path))")

    println("[SW4] w_ss*Y_ss = ", wY_ss, "   kappa_w = ", kappa_w, "   phi_pi = ", phi_pi)
    println("[SW4] i_taylor_ss = ", i_ss, "  (", 100*((1+i_ss)^4-1), "% ann)   pi_bar = ", π_bar, "  (", 100*((1+π_bar)^4-1), "% ann)")
    println("[SW4] nkpc_scale = ", nkpc_scale, "   mon_scale = ", mon_scale, "   max|epsm| = ", maximum(abs, epsm_path))
    println("[SW4] unknowns: ", n_r, " r + ", n_s, " s + ", n_p, " pi = ", n_r+n_s+n_p)

    x0 = vcat(r0, s0, p0)
    kw = (; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, χ, t, g, method, kappa_w, wY_ss, epsm_path, i_taylor_ss = i_ss, phi_pi, nkpc_scale, mon_scale, s_mode)

    # ---- residual blocks AT THE WARM START, before solving
    F0 = similar(x0)
    _residuals_sw3!(F0, x0; kw...)
    _report_blocks_sw3("WARM START", F0, n_r, n_s, n_p)

    f! = (F, x) -> _residuals_sw3!(F, x; kw...)
    sol = nlsolve(f!, x0; nlsolve_kwargs_merged...)

    Ffin = similar(sol.zero)
    _residuals_sw3!(Ffin, sol.zero; kw...)
    _report_blocks_sw3("CONVERGED", Ffin, n_r, n_s, n_p)
    return path, sol
end


function _report_blocks_sw3(tag, F, n_r, n_s, n_p)
    Ra = view(F, 1:n_r)
    Rn = view(F, (n_r+1):(n_r+n_s))
    Rm = view(F, (n_r+n_s+1):(n_r+n_s+n_p))
    println("[SW4] ---- residual blocks @ ", tag, " ----")
    for (nm, R) in (("asset  ", Ra), ("nkpc   ", Rn), ("monetary", Rm))
        println("[SW4]   ", nm, "  max|R| = ", maximum(abs, R), "   med|R| = ", median(abs.(R)), "   argmax = ", argmax(abs.(R)))
    end
    println("[SW4]   inf-norm all = ", maximum(abs, F))
    flush(stdout)
    return nothing
end


# ============================================================
# [SW5] NESTED SOLVER -- outer on pi (monetary residual), inner 2x2 on (r, s)
# ============================================================
# WHY NESTED.  Row scaling equalises diagonals but multiplies cross-block
# couplings too, so with (nkpc=240, mon=69) the pi column's largest entry moved
# onto the NKPC row and the unknown-to-equation assignment inverted; regression
# (b) diverged.  Nesting sidesteps the coupling instead of balancing it: the
# outer problem never sees the asset or NKPC rows at all.
#
#   OUTER   unknown pi_block(798) = pi_2 .. pi_799   (pi_1 fed, pi_800 = pi_bar)
#           residual  (1+r_{j+1}) - (1+i_j)/(1+pi_{j+1}),  j = 2..799
#           evaluated at the INNER solution's r path.  No scaling: one block.
#   INNER   the Stage 3 2x2 on (r, s) at nkpc_scale, with pi held fixed via
#           pi_over.  Tolerance TIGHTER than the outer (default 1e-7 vs 1e-5)
#           so the outer sees a smooth map rather than inner-solve noise.
#
# The inner solve is warm-started from the previous outer iterate's solution, so
# only the first one pays full price.  If any inner solve fails to converge the
# routine THROWS rather than returning a partial answer -- an unconverged inner
# solve is precisely the noise source the tolerance gap is there to avoid.
# ------------------------------------------------------------
struct InnerSolveFailure <: Exception
    outer_call::Int
    inner_iters::Int
    inner_norm::Float64
end
Base.showerror(io::IO, e::InnerSolveFailure) = print(io,
    "InnerSolveFailure: inner 2x2 did NOT converge at outer call ", e.outer_call,
    " (inner iterations = ", e.inner_iters, ", inner residual inf-norm = ", e.inner_norm,
    "). Stopping rather than returning a partial answer.")

function solve_transition_nested_sw(
    e_init, e_final; τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path,
    d_path, covid_path, Γ_n, Γ_c, χ, t, g, method, kappa_w, epsm_path, phi_pi,
    i_taylor_ss = nothing, nkpc_scale = 240.0, s_mode = :nkpc,
    r_path, s_path, pi_path0 = nothing,
    inner_ftol = 1e-7, outer_ftol = _ZERO_FTOL,
    inner_iters = 450, outer_iters = 200, verbose = true
)
    _assert_globals_sw(e_init, t, g)
    path = _initialize_path(e_init, e_final; qb_path)
    n_per = length(τ_k_path); n_r = _n_r_sw(n_per, method); n_s = n_per - 2; n_p = n_per - 2
    π_bar = e_init.π_
    i_ss  = isnothing(i_taylor_ss) ? (1 + e_final.r)*(1 + π_bar) - 1 : i_taylor_ss

    DY_ss = union_aggregates_(e_init.h, g, e_init.ws; s = 1.0, w = e_init.w, τ_n = e_init.τ_n,
                              n_ss = e_init.ws.n_pol[1,:])
    wY_ss = e_init.w * DY_ss.Y

    r_warm = collect(r_path); s_warm = collect(s_path)
    length(r_warm) == n_r || error("r_path must have length $(n_r), got $(length(r_warm))")
    length(s_warm) == n_s || error("s_path must have length $(n_s), got $(length(s_warm))")
    p0 = isnothing(pi_path0) ? collect(π_path[2:(n_per-1)]) : collect(pi_path0)
    length(p0) == n_p || error("pi_path0 must have length $(n_p), got $(length(p0))")

    println("[SW5] NESTED   inner ftol = ", inner_ftol, "   outer ftol = ", outer_ftol,
            "   nkpc_scale = ", nkpc_scale, "   phi_pi = ", phi_pi)
    println("[SW5] outer unknowns = ", n_p, "   inner unknowns = ", n_r + n_s,
            "   max|epsm| = ", maximum(abs, epsm_path)); flush(stdout)

    inner_iter_log = Int[]; inner_norm_log = Float64[]; outer_call = Ref(0)
    π_vec = Vector{Float64}(undef, n_per)
    F_in  = Vector{Float64}(undef, n_r + n_s)
    kw_in = (; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path,
             h_path, π_path, d_path, covid_path, Γ_n, Γ_c, χ, t, g, method, kappa_w, wY_ss,
             nkpc_scale, s_mode)

    function outer_residual!(Fout, pib)
        outer_call[] += 1
        π_vec[1] = π_path[1]; π_vec[n_per] = π_bar
        @inbounds for k in 1:n_p; π_vec[k+1] = pib[k]; end

        f_in! = (F, x) -> _residuals_h_r_sw!(F, x; kw_in..., π_over = π_vec)
        sol_in = nlsolve(f_in!, vcat(r_warm, s_warm); ftol = inner_ftol, method = :anderson,
                         m = 10, beta = -0.01, iterations = inner_iters, show_trace = false)
        # leave `path` at the inner SOLUTION: nlsolve's last probe need not be there
        _residuals_h_r_sw!(F_in, sol_in.zero; kw_in..., π_over = π_vec)
        innorm = maximum(abs, F_in)
        sol_in.f_converged || throw(InnerSolveFailure(outer_call[], sol_in.iterations, innorm))
        r_warm .= view(sol_in.zero, 1:n_r); s_warm .= view(sol_in.zero, (n_r+1):(n_r+n_s))
        push!(inner_iter_log, sol_in.iterations); push!(inner_norm_log, innorm)

        @inbounds for k in 1:n_p
            j = k + 1
            i_j = i_ss + phi_pi * (π_vec[j] - π_bar) + epsm_path[j]
            Fout[k] = (1 + path.r[j+1]) - (1 + i_j)/(1 + π_vec[j+1])
        end
        if verbose
            @printf("[SW5] outer %3d  |R_mon| = %.6e   inner: %3d iters, |R_in| = %.3e\n",
                    outer_call[], maximum(abs, Fout), sol_in.iterations, innorm); flush(stdout)
        end
        return nothing
    end

    t0 = time()
    sol_out = nlsolve(outer_residual!, p0; ftol = outer_ftol, method = :anderson, m = 10,
                      beta = -0.01, iterations = outer_iters, show_trace = false)
    wall = time() - t0

    println("\n[SW5] ---- nested summary ----")
    @printf("[SW5] outer iterations      = %d   converged = %s   |R_mon| = %.6e\n",
            sol_out.iterations, sol_out.f_converged, sol_out.residual_norm)
    @printf("[SW5] outer residual calls  = %d\n", outer_call[])
    if !isempty(inner_iter_log)
        @printf("[SW5] inner iters/outer     : first %d, median %.1f, max %d, total %d\n",
                inner_iter_log[1], median(inner_iter_log), maximum(inner_iter_log), sum(inner_iter_log))
        @printf("[SW5] inner final |R|       : median %.3e, max %.3e   (inner ftol %.1e)\n",
                median(inner_norm_log), maximum(inner_norm_log), inner_ftol)
    end
    @printf("[SW5] wall time             = %.1f s  (%.2f min)\n", wall, wall/60)
    return path, sol_out, (; inner_iter_log, inner_norm_log, outer_calls = outer_call[], wall,
                           r = copy(r_warm), s = copy(s_warm))
end


# ============================================================
# [SW6] pi(r) BACKWARD RECURSION  --  exact inverse of the monetary block
# ============================================================
# Setting the monetary residual to zero,
#     (1 + r_{j+1})(1 + pi_{j+1}) = 1 + i_ss + phi*(pi_j - pi_bar) + epsm_j
# solves for pi_j given pi_{j+1}:
#     pi_j = pi_bar + [ (1+r_{j+1})(1+pi_{j+1}) - 1 - i_ss - epsm_j ] / phi
#
# DIRECTION.  d pi_j / d pi_{j+1} = (1 + r_{j+1}) / phi ~ 0.67 < 1, so the
# recursion CONTRACTS backward; forward it would be phi/(1+r) ~ 1.49 and
# explode.  That is the Taylor principle: phi > 1 makes the path
# forward-explosive and backward-contracting, so the bounded solution is
# UNIQUE and this recursion constructs it rather than searching for it.  No
# multi-start determinacy check is needed.
#
# TERMINAL CONDITION.  pi_800 = pi_bar reproduces itself exactly: at r = r_f,
# epsm = 0, pi_{j+1} = pi_bar the bracket is (1+r_f)(1+pi_bar) - 1 - i_ss = 0
# because i_ss is DEFINED as (1+r_f)(1+pi_bar) - 1.  No drift from the tail.
#
# pi_1 is FED from data and is not produced here.  r_2 never enters (the
# earliest rate the monetary block touches is r_3), consistent with the j=2
# guard.  Factored out so the Gauss-Seidel variant can reuse it unchanged.
# ------------------------------------------------------------
function _pi_from_r_recursion_sw!(π_vec, r_block; π_path, π_bar, i_taylor_ss, phi_pi,
                                  epsm_path, n_per, r_min, r_max, sw19 = false)
    π_vec[n_per] = π_bar
    @inbounds for j in (n_per - 1):-1:2
        # [SW19] the r block starts at r_2, so r_{j+1} sits at index j, not j+1.
        rj1 = clamp(r_block[sw19 ? j : j + 1], r_min, r_max)   # same clamp the residual applies
        π_vec[j] = π_bar + ((1 + rj1) * (1 + π_vec[j + 1]) - 1 - i_taylor_ss - epsm_path[j]) / phi_pi
    end
    π_vec[1] = π_path[1]     # FED, not produced by the recursion
    return π_vec
end


# ============================================================
# [SW6] REDUCED SYSTEM  --  pi eliminated, unknowns back to [r ; s]
# ============================================================
# Exact elimination, the same move that retired residual[2]: the monetary block
# is satisfied IDENTICALLY at every evaluation, so it leaves the residual vector
# entirely and the system collapses from 2396 unknowns to the 1598 the Stage 3
# machinery already solves in 44-71 Anderson iterations.  The solution set is
# unchanged.
# ------------------------------------------------------------
function _residuals_sw_c!(residuals, x; π_vec, epsm_path, i_taylor_ss, phi_pi, path, e_init,
                          e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path,
                          h_path, π_path, d_path, covid_path, Γ_n, Γ_c, χ, t, g, method,
                          kappa_w, wY_ss, nkpc_scale = 1.0, s_mode = :nkpc, sw19 = false, sw20 = false,
                          s1_pin = nothing, forced_hours = false)
    n_per = length(τ_k_path)
    n_r   = _n_r_sw(n_per, method, sw19)
    β     = h_path[1].u.β
    r_min = -t.δ + eps()
    r_max = 1 / β
    _pi_from_r_recursion_sw!(π_vec, view(x, 1:n_r); π_path, π_bar = e_init.π_,
                             i_taylor_ss, phi_pi, epsm_path, n_per, r_min, r_max, sw19)
    _residuals_h_r_sw!(residuals, x; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path,
                       qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, χ,
                       t, g, method, kappa_w, wY_ss, nkpc_scale, s_mode, π_over = π_vec, sw19, sw20, s1_pin, forced_hours)
    return nothing
end


function solve_transition_sw_c(
    e_init, e_final; τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path,
    d_path, covid_path, Γ_n, Γ_c, χ, t, g, method, kappa_w, epsm_path, phi_pi,
    i_taylor_ss = nothing, nkpc_scale = 240.0, s_mode = :nkpc, r_path, s_path,
    verbose = true, nlsolve_kwargs = nothing, sw19 = false, sw20 = false, s1_pin = nothing,
    forced_hours = false
)
    _assert_globals_sw(e_init, t, g)
    path  = _initialize_path(e_init, e_final; qb_path)
    n_per = length(τ_k_path); n_r = _n_r_sw(n_per, method, sw19); n_s = _n_s_sw(n_per, sw20)
    π_bar = e_init.π_
    i_ss  = isnothing(i_taylor_ss) ? (1 + e_final.r)*(1 + π_bar) - 1 : i_taylor_ss
    DY_ss = union_aggregates_(e_init.h, g, e_init.ws; s = 1.0, w = e_init.w, τ_n = e_init.τ_n,
                              n_ss = e_init.ws.n_pol[1,:])
    wY_ss = e_init.w * DY_ss.Y

    r0 = collect(r_path); s0 = collect(s_path)
    length(r0) == n_r || error("r_path must have length $(n_r), got $(length(r0))")
    length(s0) == n_s || error("s_path must have length $(n_s), got $(length(s0))")
    length(epsm_path) == n_per || error("epsm_path must have length $(n_per)")

    println("[SW6] REDUCED  unknowns = ", n_r + n_s, "  (pi eliminated by the backward recursion)",
            sw19 ? "  [SW19] r_k and residuals[1] also eliminated" : "",
            sw20 ? "  [SW20] s_1 lifted into x" : "")
    println("[SW6] nkpc_scale = ", nkpc_scale, "   phi_pi = ", phi_pi, "   max|epsm| = ", maximum(abs, epsm_path))
    println("[SW6] i_taylor_ss = ", i_ss, "   pi_bar = ", π_bar); flush(stdout)

    π_vec = Vector{Float64}(undef, n_per)
    kw = (; π_vec, epsm_path, i_taylor_ss = i_ss, phi_pi, path, e_init, e_final, τ_k_path, T_path,
          Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c,
          χ, t, g, method, kappa_w, wY_ss, nkpc_scale, s_mode, sw19, sw20, s1_pin, forced_hours)

    base = (ftol = _ZERO_FTOL, show_trace = verbose, method = :anderson, m = 10,
            iterations = _ZERO_MAX_ITERS, beta = -0.01)
    kwm = isnothing(nlsolve_kwargs) ? base : merge(base, nlsolve_kwargs)

    x0 = vcat(r0, s0)
    F0 = similar(x0); _residuals_sw_c!(F0, x0; kw...)
    @printf("[SW6] warm start: asset max|R| = %.6e   nkpc max|R| = %.6e   inf-norm = %.6e\n",
            maximum(abs, view(F0, 1:n_r)), maximum(abs, view(F0, (n_r+1):(n_r+n_s))), maximum(abs, F0))
    flush(stdout)

    t0 = time()
    sol = nlsolve((F, x) -> _residuals_sw_c!(F, x; kw...), x0; kwm...)
    wall = time() - t0

    Ffin = similar(sol.zero); _residuals_sw_c!(Ffin, sol.zero; kw...)
    Ra = view(Ffin, 1:n_r); Rn = view(Ffin, (n_r+1):(n_r+n_s))
    println("\n[SW6] ---- reduced solve summary ----")
    @printf("[SW6] iterations = %d   converged = %s   inf-norm = %.6e\n",
            sol.iterations, sol.f_converged, maximum(abs, Ffin))
    @printf("[SW6] asset  max|R| = %.6e  med|R| = %.6e\n", maximum(abs, Ra), median(abs.(Ra)))
    @printf("[SW6] nkpc   max|R| = %.6e  med|R| = %.6e   (= %.6e unscaled)\n",
            maximum(abs, Rn), median(abs.(Rn)), maximum(abs, Rn)/nkpc_scale)
    @printf("[SW6] wall time  = %.1f s  (%.2f min)   <-- Stage 5 viability number\n", wall, wall/60)
    return path, sol, copy(π_vec), (; wall,)
end
