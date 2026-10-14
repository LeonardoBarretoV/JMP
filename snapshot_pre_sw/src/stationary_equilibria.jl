# ============================================================
# Steady State Functions
# ============================================================

function stationary_equilibrium_init(h, t, g, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, Λ_bstar; r_range,
    root_max, 
    verbose = true,   
    hh_problem_kwargs = nothing,
    solver = A42(), 
    find_zero_kwargs = (; atol = _ZERO_FTOL)
)
    hh_problem_options_baseline = (; 
        verbose = false
    ) 

    stationary_kwargs = isnothing(hh_problem_kwargs) ?
        hh_problem_options_baseline :
        merge(hh_problem_options_baseline, hh_problem_kwargs)

    r0 = r_range[1]
    π_ = g.π_bar
    p = ProgressUnknown()
    f = (τ_n) -> _τ_n_excess_init!(h, t, g; τ_n, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, r = r0, verbose, p).excess_g
    τ_n = root_max == 1 ? maximum(find_zeros(f, 0.0, 5.0)) : minimum(find_zeros(f, 0.0, 5.0))
    (; excess_g, w, n, k, T, Tu) = _τ_n_excess_init!(h, t, g; τ_n, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, r = r0, verbose = false, p)
    verbose && finish!(p, showvalues = [(:τ_n, τ_n), (:error, excess_g)])
    pars = (R = 1 + r0, w = w, τ_n = τ_n, T = T, Tu = Tu)

    ws = HouseholdWorkspace(; h, g, pars...)
    p = ProgressUnknown()    
    f = (r) -> _k_b_excess_init!(ws, r, h, t, g, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, Λ_bstar, root_max, verbose, p, stationary_kwargs).excess
    r = find_zero(f, r_range, solver; find_zero_kwargs...)
    (; w, a, n, k, y, τ_n, lab_tax, lump_sum_T, lump_sum_Tu, T, Tu, b, bd, bstar, c, G, q) = _k_b_excess_init!(ws, r, h, t, g, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, Λ_bstar, root_max, false, p, stationary_kwargs)
    verbose && finish!(p, showvalues = [(:r, r), (:error, a - k - q*bd)])
    e = StationaryEquilibrium(; h, t, g, ws, G, w, n, k, a, r, τ_n, τ_k, lab_tax, lump_sum_T, lump_sum_Tu, T, Tu, b, bd, bstar, c, q, π_)
    return e
end

function _τ_n_excess_init!(h, t, g; τ_n, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, r, verbose, p)
    π_ = g.π_bar
    q = 1.0 # (1/(1 + r))*((ī + g.ϑ)/(1 + π_) + (1 - g.ϑ)/(1 + π_) * qprime)
    ī= (1 + r) *(1 + π_) - 1
    w = w_from_r(t; r, τ_k)
    n = labor_supply(h, g; w = w, τ_n)
    k = k_from_r_n(t; r, n = n, τ_k)
    y = t.A * k^t.α * n^(1- t.α)
    lab_tax = taxes_HSV(h, g; w, τ_n)
    Tu = Λ_tu * y / h.Pss[1]
    T = Λ_t * y # Could put here /sum(h.Pss[1:z_bar])
    b = Λ_b * y
    G = Λ_g * y
    lump_sum_T = agg_T(h, g; T)
    lump_sum_Tu = agg_Tu(h, g; Tu)
    excess_g = lab_tax + τ_k * (r + t.δ)*k + q*(b - (1 - g.ϑ)/(1 + π_) * b) - (ī + g.ϑ)/(1 + π_) * b - lump_sum_T - lump_sum_Tu - G  #excess_g = lab_tax + τ_k * (r + t.δ)*k + q*(bprime - (1 - g.ϑ)*b) - (r + g.ϑ)*b - lump_sum
    verbose && next!(p, showvalues = [(:r, r), (:error, excess_g)])
    return (; excess_g, w, n, k, y, τ_n, τ_k, lab_tax, lump_sum_T, lump_sum_Tu, T, Tu, b, G, q)
end

function _k_b_excess_init!(ws, r, h, t, g, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, Λ_bstar, root_max, verbose, p, stationary_kwargs)
    f = (τ_n) -> _τ_n_excess_init!(h, t, g; τ_n, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, r, verbose, p).excess_g
    τ_n = root_max == 1 ? maximum(find_zeros(f, 0.0, 5.0)) : minimum(find_zeros(f, 0.0, 5.0))
    (; w, n, k, y, lab_tax, lump_sum_T, lump_sum_Tu, T, Tu, b, G, q) = _τ_n_excess_init!(h, t, g; τ_n, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, r, verbose, p)
    stationary!(ws, g; R = 1 + r, w = w, τ_n = τ_n, T = T, Tu = Tu, stationary_kwargs...)
    a = asset_supply(h.a_grid, ws.pdf)    
    bstar = Λ_bstar * y
    c = y - t.δ * k - G - (1 + r) * q * bstar + q * bstar
    bd = b - bstar
    excess = a - k - q * bd
    verbose && next!(p, showvalues = [(:r, r), (:error, excess)])
    return (; excess, w, a, n, k, y, τ_n, lab_tax, lump_sum_T, lump_sum_Tu, T, Tu, b, bd, bstar, c, G, q)
end

function stationary_equilibrium_final(h, t, g, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, Λ_bstar, ī ; r_range,
    root_max, 
    verbose = true,   
    hh_problem_kwargs = nothing,
    solver = A42(), 
    find_zero_kwargs = (; atol = _ZERO_FTOL)
)
    hh_problem_options_baseline = (; 
        verbose = false
    ) 

    stationary_kwargs = isnothing(hh_problem_kwargs) ?
        hh_problem_options_baseline :
        merge(hh_problem_options_baseline, hh_problem_kwargs)
    
        r0 = r_range[1]
    π_ = g.π_bar
    p = ProgressUnknown()
    f = (τ_n) -> _τ_n_excess_final!(h, t, g; τ_n, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, r = r0, ī, verbose, p).excess_g
    τ_n = root_max == 1 ? maximum(find_zeros(f, 0.0, 5.0)) : minimum(find_zeros(f, 0.0, 5.0))
    (; excess_g, w, n, k, T, Tu) = _τ_n_excess_final!(h, t, g; τ_n, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, r = r0, ī, verbose = false, p)
    verbose && finish!(p, showvalues = [(:τ_n, τ_n), (:error, excess_g)])
    pars = (R = 1 + r0, w = w, τ_n = τ_n, T = T, Tu = Tu)

    ws = HouseholdWorkspace(; h, g, pars...)
    p = ProgressUnknown()    
    f = (r) -> _k_b_excess_final!(ws, r, h, t, g, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, Λ_bstar, ī, root_max, verbose, p, stationary_kwargs).excess
    r = find_zero(f, r_range, solver; find_zero_kwargs...)
    (; w, a, n, k, τ_n, lab_tax, lump_sum_T, lump_sum_Tu, T, Tu, b, bd, bstar, c, G, q) = _k_b_excess_final!(ws, r, h, t, g, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, Λ_bstar, ī, root_max, false, p, stationary_kwargs)
    verbose && finish!(p, showvalues = [(:r, r), (:error, a - k - q*b)])
    e = StationaryEquilibrium(; h, t, g, ws, G, w, n, k, a, r, τ_n, τ_k, lab_tax, lump_sum_T, lump_sum_Tu, T, Tu, b, bd, bstar, c, q, π_)
    return e
end

function _τ_n_excess_final!(h, t, g; τ_n, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, r, ī, verbose, p)
    π_ = g.π_bar
    q = (ī + g.ϑ) / ((1 + r)*(1 + π_) - (1 - g.ϑ)) # q in the SS after imposing q = qprime. Would be equal to 1 in the initial SS, but no in the final
    w = w_from_r(t; r, τ_k)
    n = labor_supply(h, g; w = w, τ_n)
    k = k_from_r_n(t; r, n = n, τ_k)

    y = t.A * k^t.α * n^(1- t.α)
    lab_tax = taxes_HSV(h, g; w, τ_n)
    Tu = Λ_tu * y / h.Pss[1]
    T = Λ_t * y # Could put here /sum(h.Pss[1:z_bar])
    b = Λ_b * y / q
    G = Λ_g * y
    lump_sum_T = agg_T(h, g; T)
    lump_sum_Tu = agg_Tu(h, g; Tu)
    excess_g = lab_tax + τ_k * (r + t.δ)*k + q*b - (1 + r) * q*b - lump_sum_T - lump_sum_Tu - G   #excess_g = lab_tax + τ_k * (r + t.δ)*k + q*(bprime - (1 - g.ϑ)*b) - (r + g.ϑ)*b - lump_sum
    verbose && next!(p, showvalues = [(:r, r), (:error, excess_g)])
    return (; excess_g, w, n, k, y, τ_n, τ_k, lab_tax, lump_sum_T, lump_sum_Tu, T, Tu, b, G, q)
end

function _k_b_excess_final!(ws, r, h, t, g, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, Λ_bstar, ī, root_max, verbose, p, stationary_kwargs)
    f = (τ_n) -> _τ_n_excess_final!(h, t, g; τ_n, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, r, ī, verbose, p).excess_g
    τ_n = root_max == 1 ? maximum(find_zeros(f, 0.0, 5.0)) : minimum(find_zeros(f, 0.0, 5.0))
    (; w, n, k, y, lab_tax, lump_sum_T, lump_sum_Tu, T, Tu, b, G, q) = _τ_n_excess_final!(h, t, g; τ_n, τ_k, Λ_t, Λ_tu, Λ_g, Λ_b, r, ī, verbose, p)
    stationary!(ws, g; R = 1 + r, w = w, τ_n = τ_n, T = T, Tu = Tu, stationary_kwargs...)
    a = asset_supply(h.a_grid, ws.pdf)
    bstar = Λ_bstar * y
    c = y - t.δ * k - G - (1 + r) * q * bstar + q * bstar
    bd = b - bstar
    excess = a - k - q * bd
    verbose && next!(p, showvalues = [(:r, r), (:error, excess)])
    return (; excess, w, a, n, k, y, τ_n, lab_tax, lump_sum_T, lump_sum_Tu, T, Tu, b, bd, bstar, c, G, q)
end

# ============================================================
# Government Budget Constraint and Portfolios
# ============================================================

function _τ_n_excess_port_!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, yprime, d, n_ss, Γ_n, covid, z_dist, r, verbose, p)
    #Here qbprime = q_t b_t+1 and qb = q_t-1 b_t. This is what I feed in the paths
    w = w_from_r_(t; r, τ_k, d, covid)
    n = labor_supply_(h, g; w = w, τ_n, n_ss, Γ_n, covid, z_dist)
    k = k_from_r_n_(t; r, n = n, τ_k, d, covid)
    y = t.A *(d * k)^t.α * n^(1- t.α)
    lab_tax = taxes_HSV_(h, g; w = w, τ_n, n_ss, Γ_n, covid, z_dist)
    lump_sum_T = T * y #agg_T_(g; T, z_dist)
    lump_sum_Tu = Tu * y #agg_Tu_(g; Tu, z_dist)
    excess_g = lab_tax + τ_k * (r + t.δ)*k + qbprime * yprime - (1 + r) * qb * y - lump_sum_T - lump_sum_Tu - G * y #Everything in levels 
    verbose && next!(p, showvalues = [(:r, r), (:error, excess_g)])
    return (; excess_g, w, y, n, k, τ_n, τ_k, lab_tax, lump_sum_T, lump_sum_Tu, G)
end

function _τ_n_excess_t0_port_!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, yprime, d, n_ss, Γ_n, covid, z_dist, r_b, r_k, verbose, p)
    w = w_from_r_k_(t; r_k, τ_k, d, covid)
    n = labor_supply_(h, g; w = w, τ_n, n_ss, Γ_n, covid, z_dist)
    k = k_from_r_k_n_(t; r_k, n = n, τ_k, d, covid)
    y = t.A *(d * k)^t.α * n^(1- t.α)
    lab_tax = taxes_HSV_(h, g; w = w, τ_n, n_ss, Γ_n, covid, z_dist)
    lump_sum_T = T * y #agg_T_(g; T, z_dist)
    lump_sum_Tu = Tu * y #agg_Tu_(g; Tu, z_dist)
    excess_g = lab_tax + τ_k * r_k * k + qbprime * yprime - (1 + r_b) * qb - lump_sum_T - lump_sum_Tu - G
    verbose && next!(p, showvalues = [(:r_b, r_b), (:error, excess_g)])
    return (; excess_g, w, y, n, k, τ_n, τ_k, lab_tax, lump_sum_T, lump_sum_Tu, G)
end

function _τ_n_excess_port_l!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, d, n_ss, Γ_n, covid, z_dist, r, verbose, p)
    w = w_from_r_(t; r, τ_k, d, covid)
    n = labor_supply_(h, g; w = w, τ_n, n_ss, Γ_n, covid, z_dist)
    k = k_from_r_n_(t; r, n = n, τ_k, d, covid)
    y = t.A *(d * k)^t.α * n^(1- t.α)
    lab_tax = taxes_HSV_(h, g; w = w, τ_n, n_ss, Γ_n, covid, z_dist)
    lump_sum_T = T * y #agg_T_(g; T, z_dist)
    lump_sum_Tu = Tu * y #agg_Tu_(g; Tu, z_dist)
    excess_g = lab_tax + τ_k * (r + t.δ)*k + qbprime - (1 + r) * qb - lump_sum_T - lump_sum_Tu - G * y #Everything in levels 
    verbose && next!(p, showvalues = [(:r, r), (:error, excess_g)])
    return (; excess_g, w, y, n, k, τ_n, τ_k, lab_tax, lump_sum_T, lump_sum_Tu, G)
end

function _τ_n_excess_t0_port_l!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, d, n_ss, Γ_n, covid, z_dist, r_b, r_k, verbose, p)
    w = w_from_r_k_(t; r_k, τ_k, d, covid)
    n = labor_supply_(h, g; w = w, τ_n, n_ss, Γ_n, covid, z_dist)
    k = k_from_r_k_n_(t; r_k, n = n, τ_k, d, covid)
    y = t.A *(d * k)^t.α * n^(1- t.α)
    lab_tax = taxes_HSV_(h, g; w = w, τ_n, n_ss, Γ_n, covid, z_dist)
    lump_sum_T = T * y #agg_T_(g; T, z_dist)
    lump_sum_Tu = Tu * y #agg_Tu_(g; Tu, z_dist)
    excess_g = lab_tax + τ_k * r_k * k + qbprime - (1 + r_b) * qb - lump_sum_T - lump_sum_Tu - G
    verbose && next!(p, showvalues = [(:r_b, r_b), (:error, excess_g)])
    return (; excess_g, w, y, n, k, τ_n, τ_k, lab_tax, lump_sum_T, lump_sum_Tu, G)
end
