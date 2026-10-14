# ============================================================
# [SW] Union block  --  STICKY-WAGE TRACK ONLY
# ============================================================
# Included by scripts/main_JMP_sw.jl ONLY. Never by scripts/main_JMP.jl.
#
# STAGE 1: aggregates only. Nothing here is wired into any solver yet.
#
# NAMESPACE RULE
#   Julia has no file-level namespacing without modules. A function defined here
#   with the same name as one in an existing src/ file would OVERWRITE it in the
#   same session, even though that file is byte-identical on disk. Every name in
#   this file carries a _union suffix and was checked against the 127 top-level
#   names defined across existing src/ files. No collisions.
#
# CONVENTIONS THIS FILE RELIES ON (all verified against the repo)
#   - tau_n is the HSV RETENTION rate (1 - tax rate), not the tax rate.
#   - pdf / c_pol / v / n_pol are indexed (a, z): rows = assets, cols = productivity.
#   - ntilde_z = n_ss[z] = e_init.ws.n_pol[1,:]. Row 1 is the lowest asset grid
#     point; exact because GHH makes hours a-independent (measured spread 0.0).
#   - z_grid[1] = 0 is the unemployed state, with n_ss[1] = 0.
# ============================================================


# ------------------------------------------------------------
# Aggregate effective labor under proportional rationing.
#   N_t = s_t * sum_z( z_dist[z] * z_grid[z] * ntilde_z )
# NOTE: this is NOT s * sum(z_dist .* z_grid). At the initial steady state the
# first equals e_init.n = 0.49721019; the second equals 1.0. Factor of ~2 apart.
# ------------------------------------------------------------
function labor_supply_union_(h; s, n_ss, z_dist)
    n = zero(eltype(h.z_grid))
    for i in eachindex(z_dist)
        n += z_dist[i] * h.z_grid[i] * n_ss[i]
    end
    return s * n
end


# ------------------------------------------------------------
# HSV tax base under proportional rationing.
#   Xi = sum_{z>0} z_dist[z] * (w * z_grid[z] * s * ntilde_z)^(1 - tau)
# Net labor tax revenue is then  w*N - tau_n*Xi, which is AFFINE in tau_n once
# hours no longer depend on it.
# ------------------------------------------------------------
function hsv_base_union_(h, g; w, s, n_ss, z_dist)
    base = zero(eltype(h.z_grid))
    for i in eachindex(z_dist)
        h.z_grid[i] == 0 && continue
        base += z_dist[i] * (w * h.z_grid[i] * s * n_ss[i])^(1 - g.τ)
    end
    return base
end


# ------------------------------------------------------------
# lambda-weighted union aggregates D and Y (= Ycal), over EMPLOYED (z > 0) only.
#
#   ntilde = n_ss[z]           (steady-state hours profile)
#   n_iz   = s * ntilde        (actual hours this period)
#   x      = c_pol - theta*n_iz^(1+1/nu)/(1+1/nu)
#   lambda = (1-beta) * v^xi_inv * x^(-xi_inv)
#   rho    = (1-tau) * tau_n * (w * z * n_iz)^(-tau)
#   D     += pdf * lambda * ntilde^(1+1/nu)
#   Y     += pdf * lambda * rho * z * ntilde
#
# CRITICAL: D and Y are built from ntilde, NOT from n_iz. Only rho and x are
# evaluated at actual hours. The s-dependence enters the Phillips curve
# explicitly as s and s^(1/nu). Building D and Y from n_iz instead double-counts
# s; it is exactly right at s = 1, so a steady-state check CANNOT catch it.
#
# Skipping z = 0 is MANDATORY, not stylistic: rho at z = 0 is 0^(-tau) = Inf,
# and Inf * 0 = NaN, which propagates silently instead of erroring.
#
# `item` is duck-typed: it needs .c_pol, .v, .pdf indexed (a, z). Both a
# HouseholdWorkspace (e_init.ws) and a transition path row (path[i]) satisfy it.
# ------------------------------------------------------------
function union_aggregates_(h, g, item; s, w, τ_n, n_ss)
    θ  = h.v.θ
    ν  = h.v.ν
    τ  = g.τ
    β  = get_β(h.u)
    ξi = get_inverse_ies(h.u)
    z_grid = h.z_grid

    D = zero(eltype(h.a_grid))
    Y = zero(eltype(h.a_grid))

    for iz in eachindex(z_grid)
        z = z_grid[iz]
        z == 0 && continue                      # MANDATORY, see header
        ñ    = n_ss[iz]
        n_iz = s * ñ
        ρ    = (1 - τ) * τ_n * (w * z * n_iz)^(-τ)
        dis  = θ * n_iz^(1 + 1/ν) / (1 + 1/ν)
        ñ_p  = ñ^(1 + 1/ν)
        for ia in axes(item.pdf, 1)
            x = item.c_pol[ia, iz] - dis
            λ = (1 - β) * item.v[ia, iz]^ξi * x^(-ξi)
            μ = item.pdf[ia, iz]
            D += μ * λ * ñ_p
            Y += μ * λ * ρ * z * ñ
        end
    end
    return (D = D, Y = Y)
end
