# ============================================================
# Transition path functions
# ============================================================

function h_path_transition(Pee_path, Puu_path; θ, ν, ies, crra, β, nP, nIID , ar1, sigmaP, sigmaIID, a_min, a_max, z_u, n_a)
    out = Array{Any}(undef, length(Pee_path))
    Threads.@threads for i in eachindex(Pee_path)
        Pee = Pee_path[i]
        Puu = Puu_path[i]
        v = GHH(θ = θ, ν = ν)
        ar1 = ar1
        P, z_vals = calibration_unemploy(nP, nIID , ar1, sigmaP, sigmaIID, Pee, Puu, z_u)
        u = EZ(ies = ies, ra = crra, β = β)
        out[i] = Household(u = u, a_grid = grid(start = a_min, stop = a_max, length = n_a, scale = :log),
            v = v, P = P, z_grid = z_vals)
    end
    return out
end

function h_path_transition_mckay(ROOT, Pee_path, Puu_path; θ, ν, ies, crra, β, a_min, a_max, z_u, n_a)
    out = Array{Any}(undef, length(Pee_path))
    Threads.@threads for i in eachindex(Pee_path)
        Pee = Pee_path[i]
        Puu = Puu_path[i]
        v = GHH(θ = θ, ν = ν)
        P, z_vals = calibration_unemploy_mckay(ROOT, Pee, Puu, z_u)
        u = EZ(ies = ies, ra = crra, β = β)
        out[i] = Household(u = u, a_grid = grid(start = a_min, stop = a_max, length = n_a, scale = :log),
            v = v, P = P, z_grid = z_vals)
    end
    return out
end

function h_path_transition_mckay_(ROOT, Pee_path, Puu_path; θ, ν, ies, crra, β, a_min, a_max, z_u, n_a)
    out = Array{Any}(undef, length(Pee_path))
    Threads.@threads for i in eachindex(Pee_path)
        Pee = Pee_path[i]
        Puu = Puu_path[i]
        v = GHH(θ = θ, ν = ν)
        P, z_vals, _ = calibration_unemploy_mckay_(ROOT, Pee, Puu, z_u)   # [Z-NORM] third return value is the scale factor (was: P, z_vals = ...)
        u = EZ(ies = ies, ra = crra, β = β)
        out[i] = Household(u = u, a_grid = grid(start = a_min, stop = a_max, length = n_a, scale = :log),
            v = v, P = P, z_grid = z_vals)
    end
    return out
end

function h_path_transition_mckay_init(ROOT, Pee_path, Puu_path; θ, ν, ies, crra, β, a_min, a_max, z_u, n_a, init_state, z_scale = 1.0)   # [Z-NORM] added z_scale kwarg
    chain_long = zeros(length(init_state),length(Pee_path))
    chain = zeros(2,length(Pee_path))
    init_state_short = zeros(2,1); init_state_short[1] = init_state[1]; init_state_short[2] = 1.0 - init_state[1]; #init_state needs to add up to 1. Always happening in a given distribution
    for i in eachindex(Pee_path)
        Pee = Pee_path[i]
        Puu = Puu_path[i]
        Pprime_u = [Puu 1-Pee; 1-Puu Pee]
        chain[:,i] = i ==1 ? init_state_short' * Pprime_u' : chain[:,i - 1]' * Pprime_u'
    end

    out = Array{Any}(undef, length(Pee_path))
    init_state_updated = similar(init_state)
    for i in eachindex(Pee_path) #Threads.@threads
        Pee = Pee_path[i]
        Puu = Puu_path[i]
        v = GHH(θ = θ, ν = ν)
        init_state_new = i == 1 ? init_state : init_state_updated'
        P, z_vals = calibration_unemploy_mckay_init(ROOT, Pee, Puu, z_u; unemp_rate_target = chain[1,i], Big_init_state = init_state_new, z_scale = z_scale)   # [Z-NORM] pass the initial-SS scale factor through
        init_state_updated = init_state_new' * P
        chain_long[:,i] = init_state_updated
        u = EZ(ies = ies, ra = crra, β = β)
        out[i] = Household(u = u, a_grid = grid(start = a_min, stop = a_max, length = n_a, scale = :log),
            v = v, P = P, z_grid = z_vals)
    end
    return out, chain, chain_long
end

function _transition_item(h::Household)  
    (; v, η, a_pol, pdf, lower_index, lower_weight, a_tmp, c_pol, n_pol) = _generate_base_workspace_matrices(h)

    r = w = τ_n = τ_k = lab_tax = lump_sum_T = lump_sum_Tu = T = Tu = a = k = b = bd =  bstar = n = y = C = X = q = r_k = π_ = un = l = G = Ψ = zero(eltype(h.a_grid)) 
    return (; h, v, η, a_pol, pdf, lower_index, lower_weight, a_tmp, c_pol, n_pol, r, w, τ_n, τ_k, lab_tax, lump_sum_T, lump_sum_Tu, T, Tu, a, k, b, bd, bstar, n, y, C, X, q, r_k, π_, un, l, G, Ψ)
end 


function _initialize_path(e_init, e_final; qb_path)
    h = e_init.h

    lst = StructArray(_transition_item(h) for _ in qb_path)
    lst[1].pdf .= e_init.ws.pdf

    lst[end].v .= e_final.ws.v 
    lst[end].η .= e_final.ws.η
    lst[end].a_pol .= e_final.ws.a_pol
    lst[end].lower_index .= e_final.ws.lower_index
    lst[end].lower_weight .= e_final.ws.lower_weight
    lst[end].c_pol .= e_final.ws.c_pol
    lst[end].n_pol .= e_final.ws.n_pol

    lst.r[end] = e_final.r
    lst.q[end] = e_final.q
    lst.π_[end] = e_final.π_
    lst.w[end] = e_final.w
    lst.τ_n[end] = e_final.τ_n
    lst.τ_k[end] = e_final.τ_k
    lst.lab_tax[end] = e_final.lab_tax
    lst.lump_sum_T[end] = e_final.lump_sum_T
    lst.lump_sum_Tu[end] = e_final.lump_sum_Tu
    lst.T[end] = e_final.T
    lst.Tu[end] = e_final.Tu
    lst.b[end] = e_final.b
    lst.k[end] = e_final.k
    lst.a[end] = e_final.a
    lst.n[end] = e_final.n

    return lst 
end 


function solve_transition_slow_but_robust(
    e_init, e_final; k_path, b_path, h, t, g, r_path = nothing, verbose = true,
    nlsolve_kwargs = nothing
)
    nlsolve_baseline_kwargs = (
        ftol = _ZERO_FTOL,
        show_trace = verbose,
        method = :anderson,
        m = 10,
        iterations = _ZERO_MAX_ITERS,
        beta = -0.1 #-0.1  # adjustment parameter for the fixed point algorithm
    )
    
    nlsolve_kwargs_merged = isnothing(nlsolve_kwargs) ?             
        nlsolve_baseline_kwargs :
        merge(nlsolve_baseline_kwargs, nlsolve_kwargs)

    path = _initialize_path(e_init, e_final; b_path)

    if isnothing(r_path) 
        r_path_rest = fill(e_final.r, length(k_path) - 1)
    else 
        r_path_rest = r_path[2:end]
    end 
    r_path_rest[end] = e_final.r

    f! = (F, x) -> _residuals!(F, x; path, e_init, e_final, k_path, b_path, h, t, g)
    r_sol = nlsolve(f!, r_path_rest; nlsolve_kwargs_merged...)
    _residuals!(similar(r_sol.zero), r_sol.zero; path, e_init, e_final, k_path, b_path, h, t, g)
    
    return path 
end

function solve_transition(
    e_init, e_final; τ_k_path, T_path, Tu_path, b_path, h, t, g, r_path = nothing, verbose = true,
    nlsolve_kwargs = nothing
)
    nlsolve_baseline_kwargs = (
        ftol = _ZERO_FTOL,
        show_trace = verbose,
        method = :anderson,
        m = 10,
        iterations = _ZERO_MAX_ITERS,
        beta = -0.1 #-0.1  # adjustment parameter for the fixed point algorithm
    )
    
    nlsolve_kwargs_merged = isnothing(nlsolve_kwargs) ?             
        nlsolve_baseline_kwargs :
        merge(nlsolve_baseline_kwargs, nlsolve_kwargs)

    path = _initialize_path(e_init, e_final; b_path)

    if isnothing(r_path) 
        r_path_rest = fill(e_final.r, length(τ_k_path))
    else 
        r_path_rest = r_path
    end 
    r_path_rest[end] = e_final.r

    f! = (F, x) -> _residuals!(F, x; path, e_init, e_final, τ_k_path, T_path, Tu_path, b_path, h, t, g)
    r_sol = nlsolve(f!, r_path_rest; nlsolve_kwargs_merged...)
    _residuals!(similar(r_sol.zero), r_sol.zero; path, e_init, e_final, τ_k_path, T_path, Tu_path, b_path, h, t, g)
    
    return path 
end

# ============================================================
# Transition with revaluation of assets
# ============================================================

function solve_transition_h(
    e_init, e_final; τ_k_path, T_path, Tu_path, b_path, h_path, π_path, t, g, method, r_path = nothing, verbose = true,
    nlsolve_kwargs = nothing
)
    nlsolve_baseline_kwargs = (
        ftol = _ZERO_FTOL,
        show_trace = verbose,
        method = :anderson,
        m = 10,
        iterations = _ZERO_MAX_ITERS,
        beta = -0.01 # adjustment parameter for the fixed point algorithm
    )
    
    nlsolve_kwargs_merged = isnothing(nlsolve_kwargs) ?             
        nlsolve_baseline_kwargs :
        merge(nlsolve_baseline_kwargs, nlsolve_kwargs)

    path = _initialize_path(e_init, e_final; b_path)

    if isnothing(r_path) 
        r_path_rest = fill(e_final.r, length(τ_k_path)+1) # The length is increased by 1, as we need to solve for r_b and r_k in period 0
    else 
        r_path_rest = r_path
    end 
    f! = (F, x) -> _residuals_h!(F, x; path, e_init, e_final, τ_k_path, T_path, Tu_path, b_path, h_path, π_path, t, g, method)
    r_sol = nlsolve(f!, r_path_rest; nlsolve_kwargs_merged...)
    _residuals_h!(similar(r_sol.zero), r_sol.zero; path, e_init, e_final, τ_k_path, T_path, Tu_path, b_path, h_path, π_path, t, g, method)
    
    return path 
end

function _residuals_h!(residuals, r_path_rest; path, e_init, e_final, τ_k_path, T_path, Tu_path, b_path, h_path, π_path, t, g, method)
    #r̄ = e_init.r # Return on bonds in the initial SS that makes the bond price equal to 1 (normalization of the coupons)
    π_bar = e_init.π_
    ī = (1 + e_init.r)*(1 + π_bar) - 1
    A = t.A
    α = t.α
    δ = t.δ
    β = h_path[1].u.β
    r_min = - δ + eps()  # minimum interest rate
    r_max = 1 / β  # a maximum -- could be relaxed
    
    for i in reverse(eachindex(path, b_path, τ_k_path))
        if i > firstindex(path)
            #println("i= ", i)
            h = h_path[i]
            b, τ_k, T, Tu, π_ = b_path[i], τ_k_path[i], T_path[i], Tu_path[i], π_path[i]
            bprime = i == lastindex(b_path) ? e_final.b : b_path[i + 1] 
            r = r_path_rest[i+1] # Note r_path_rest has one more entry, as we need to solve for two returns in the first period. However [i+1] here refers to period t
            r = clamp(r, r_min, r_max)
            R = 1 + r

            qprime = i == lastindex(b_path) ? e_final.q : path.q[i + 1] 
            p = ProgressUnknown()
            f = (τ_n) -> _τ_n_excess_port!(h, t, g; τ_n, τ_k, T, Tu, π_, bprime, b, qprime, r, ī, verbose = false, p).excess_g
            τ_n = i == lastindex(b_path) ? e_final.τ_n : fzero(f, path.τ_n[i+1]) # ; Guess of τ_n in period t is τ_n in period t+1 
            (; w, k, n, lab_tax, lump_sum_T, lump_sum_Tu, q) = _τ_n_excess_port!(h, t, g; τ_n, τ_k, T, Tu, π_, bprime, b, qprime, r, ī, verbose = false, p)
            
            path.G[i] = 0.0
            path.l[i] = bprime - (1 - g.ϑ)/(1 + π_) * b
            path.h[i] = h
            path.q[i] = q
            path.π_[i] = π_
            path.τ_k[i] = τ_k
            path.τ_n[i] = τ_n
            path.r[i] = r 
            path.w[i] = w
            path.k[i] = k
            path.n[i] = n
            path.y[i] = A * path.k[i]^α * path.n[i]^(1- α)
            path.lab_tax[i] = lab_tax 
            path.lump_sum_T[i] = lump_sum_T 
            path.lump_sum_Tu[i] = lump_sum_Tu 
            next_values = i == lastindex(path) ? e_final.ws : path[i+1]
            backwards_once!(h, g, path[i]; next_values = next_values, R, w, τ_n, T, Tu, path[i].a_tmp)
            asset_policy_given_η!(path[i].a_pol; h, g, path[i].η, R, w, τ_n, T, Tu)
            interpolate_asset_policy!(path[i].lower_index, path[i].lower_weight; h.a_grid, path[i].a_pol)
        end
    end
    # Period 1
    h = h_path[1]
    r_k_hat = r_path_rest[1] # Note here that we are storing the return on capital in period 1 (after the MIT shock) in r_path_rest[2]
    r_k = r_k_hat + t.δ
    r_b = r_path_rest[2] # Note here that we are storing the return on bonds in period 1 (after the MIT shock) in r_path_rest[2]
    R_k = 1 + r_k
    R_b = 1 + r_b
    b, τ_k, T, Tu, π_ = b_path[1], τ_k_path[1], T_path[1], Tu_path[1], π_path[1]
    bprime = b_path[2]
    qprime = path.q[2]

    p = ProgressUnknown()
    f = (τ_n) -> _τ_n_excess_t0_port!(h, t, g; τ_n, τ_k, T, Tu, π_, bprime, b, qprime, r_b, r_k, ī, verbose = false, p).excess_g
    τ_n = fzero(f, path.τ_n[2]) #; path.τ_n[2] # Guess of τ_n in period t is τ_n in period t+1 
    (; w, n, k, lab_tax, lump_sum_T, lump_sum_Tu, q) = _τ_n_excess_t0_port!(h, t, g; τ_n, τ_k, T, Tu, π_, bprime, b, qprime, r_b, r_k, ī, verbose = false, p)
        
    path.G[1] = 0.0 
    path.l[1] = bprime - (1 - g.ϑ)/(1 + π_) * b
    path.h[1] = h
    path.q[1] = q
    path.π_[1] = π_
    path.τ_k[1] = τ_k
    path.τ_n[1] = τ_n
    path.r[1] = r_b 
    path.r_k[1] = r_k 
    path.w[1] = w
    path.k[1] = k
    path.b[1] = b
    path.n[1] = n
    path.y[1] = A * path.k[1]^α * path.n[1]^(1- α)
    path.lab_tax[1] = lab_tax 
    path.lump_sum_T[1] = lump_sum_T
    path.lump_sum_Tu[1] = lump_sum_Tu
    next_values = path[2]
    backwards_once_t0_port!(h, g, path[1]; next_values = next_values, R_k, R_b, w, τ_n, T, Tu, path[1].a_tmp, e_init, method)
    interpolate_asset_policy!(path[1].lower_index, path[1].lower_weight; h.a_grid, path[1].a_pol)

    # Forward iterations
    # Initial portfolio
    if method ==1
        Λ_ka = e_init.k / e_init.a
        k_agg = sum((Λ_ka.*h.a_grid)  .* path[1].pdf)
        b_agg = sum(((1.0 - Λ_ka).*h.a_grid) .* path[1].pdf)
    elseif method == 2
        (; Θ_k, Λ_ka, Λ_ba) = set_portfolio_t0(h; e_init, method = 2)
        k_agg = sum((Λ_ka.*h.a_grid)  .* path[1].pdf)
        b_agg = sum(((1.0 .- Λ_ka).*h.a_grid) .* path[1].pdf)
    elseif method == 3
        (; Θ_k, Λ_ka, Λ_ba) = set_portfolio_t0(h; e_init, method = 3)
        k_agg = sum((Λ_ka.*h.a_grid)  .* path[1].pdf)
        b_agg = sum(((1.0 .- Λ_ka).*h.a_grid) .* path[1].pdf)
    end
    #k_agg = aggregate_capital(h.a_grid, path[1].pdf) # Aggregate capital in period 0
    #b_agg = aggregate_bonds(h.a_grid, path[1].pdf) # Aggregate bonds in period 0
    #for i in eachindex(path)
    for i in 1:length(path.r)+1
        if i==1 # Find r_k to clear the capital market in period 0 after the MIT shock
            path.un[1] = sum(path[1].pdf, dims = 1)[1]
            a = asset_supply(h_path[i].a_grid, path[i].pdf)
            path.a[i] = a # Store path of assets in period 0 (predetermined by the initial SS)          
            residuals[1] =  k_agg - path.k[1] # Period 0, index 1
        elseif i==2 # Find r_b to clear the bonds market in period 0 after MIT shock
            residuals[2] =  b_agg - path.b[1] # Period 0, index 1
        else
            a = asset_supply(h_path[i-1].a_grid, path[i-1].pdf)
            path.a[i-1] = a
            path.b[i-1] = b_path[i-1]     
            residuals[i] =  a - path.k[i-1] - path.q[i - 2] * path.b[i-1] # Period 1, index 2
        end
        (i == length(path.r)+1) && break

        if i > 1
            forward_pdf!(path[i].pdf; h = h_path[i], path[i-1].pdf, path[i-1].lower_index, path[i-1].lower_weight)
            path.un[i] = sum(path[i].pdf, dims = 1)[1]
        end
    end
end 

# ============================================================
# Adding Pandemic Paths
# ============================================================

function solve_transition_h_(
    e_init, e_final; τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, t, g, method, r_path = nothing, verbose = true,
    nlsolve_kwargs = nothing
)
    nlsolve_baseline_kwargs = (
        ftol = _ZERO_FTOL,
        show_trace = verbose,
        method = :anderson,
        m = 10,
        iterations = _ZERO_MAX_ITERS,
        beta = -0.01 # adjustment parameter for the fixed point algorithm
    )
    
    nlsolve_kwargs_merged = isnothing(nlsolve_kwargs) ?             
        nlsolve_baseline_kwargs :
        merge(nlsolve_baseline_kwargs, nlsolve_kwargs)

    path = _initialize_path(e_init, e_final; qb_path)

    if isnothing(r_path) 
        r_path_rest = fill(e_final.r, length(τ_k_path)+1) # The length is increased by 1, as we need to solve for r_b and r_k in period 0
    else 
        r_path_rest = r_path
    end 
    f! = (F, x) -> _residuals_h_!(F, x; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, t, g, method)
    r_sol = nlsolve(f!, r_path_rest; nlsolve_kwargs_merged...)
    _residuals_h_!(similar(r_sol.zero), r_sol.zero; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, t, g, method)
    
    return path 
end

# ============================================================
# Transition with Open Economy: Adding Foreign Debt ratio
# ============================================================

function solve_transition_h_r(
    e_init, e_final; τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, χ, t, g, method, r_path = nothing, verbose = true,
    nlsolve_kwargs = nothing
)
    nlsolve_baseline_kwargs = (
        ftol = _ZERO_FTOL,
        show_trace = verbose,
        method = :anderson,
        m = 10,
        iterations = _ZERO_MAX_ITERS,
        beta = -0.01 #-0.002#-0.01 #-0.1  # adjustment parameter for the fixed point algorithm
    )
    
    nlsolve_kwargs_merged = isnothing(nlsolve_kwargs) ?             
        nlsolve_baseline_kwargs :
        merge(nlsolve_baseline_kwargs, nlsolve_kwargs)

    path = _initialize_path(e_init, e_final; qb_path)

    if isnothing(r_path) 
        r_path_rest = fill(e_final.r, length(τ_k_path)+1) # The length is increased by 1, as we need to solve for r_b and r_k in period 0
    else 
        r_path_rest = r_path
    end 
    f! = (F, x) -> _residuals_h_r!(F, x; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, χ, t, g, method)
    r_sol = nlsolve(f!, r_path_rest; nlsolve_kwargs_merged...)
    _residuals_h_r!(similar(r_sol.zero), r_sol.zero; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, χ, t, g, method)
    return path 
end

function _residuals_h_r!(residuals, r_path_rest; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, χ, t, g, method)
    #r̄ = e_init.r # Return on bonds in the initial SS that makes the bond price equal to 1 (normalization of the coupons)
    π_bar = e_init.π_
    ī = (1 + e_init.r)*(1 + π_bar) - 1
    A = t.A
    α = t.α
    δ = t.δ
    β = h_path[1].u.β
    r_min = - δ + eps()  # minimum interest rate
    r_max = 1 / β  # a maximum -- could be relaxed
    n_ss = e_init.ws.n_pol[1,:]
    c_ss = e_init.ws.c_pol

    #Transition of employment states
    #From period 0 to period 1 I do the transition using h_path[2]
    z_state = zeros(length(h_path), length(e_init.h.z_grid))
    z_state[1,:] = e_init.h.Pss'*h_path[1].P
    for i = 2:lastindex(h_path)
        z_state[i,:] = z_state[i - 1,:]' * h_path[i].P 
    end
    #Backward iteration
    for i in reverse(eachindex(path, qb_path, τ_k_path))
        if i > firstindex(path)
            #println("i= ", i)
            z_dist = z_state[i, :]
            h = i == lastindex(path) ? h_path[end] : h_path[i + 1] # In period 599, I use h_path[600] to compute transition from states in 599 (z_it) to 600 (z_it+1)
            qb, qbstar, τ_k, T, Tu, G, π_, d, covid = qb_path[i], qbstar_path[i], τ_k_path[i], T_path[i], Tu_path[i], G_path[i], π_path[i], d_path[i], covid_path[i]
            qbprime = i == lastindex(path) ? e_final.q * e_final.b / output_(t; e_final.k, e_final.n, d, covid) : qb_path[i + 1] 
            qbstar_prime = i == lastindex(path) ?  e_final.bstar / e_final.b : qbstar_path[i + 1] # qbstar_path[end] : qbstar_path[i + 1] 
            r = r_path_rest[i+1] # Note r_path_rest has one more entry, as we need to solve for two returns in the first period. However [i+1] here refers to period t
            r = clamp(r, r_min, r_max) # This is r_t. The return in t of the bonds that I bought in t-1
            R = 1 + r

            yprime = i == lastindex(path) ? output_(t; e_final.k, e_final.n, d, covid) : path.y[i + 1] # Here [i+1] refers to period t
            qprime = i == lastindex(path) ? e_final.q : path.q[i + 1] # Here [i+1] refers to period t
            rprime = i == lastindex(path) ? e_final.r : path.r[i + 1] # This refers to r_t+1
            π_prime = i == lastindex(path) ? e_final.π_ : path.π_[i + 1] # Here [i+1] refers to period t
            q = (1/(1 + rprime))*((ī + g.ϑ)/(1 + π_prime) + (1 - g.ϑ)/(1 + π_prime) * qprime) # Only used to compute q_t recover the face value of b_t+1
            bprime_level = (qbprime / q) * yprime
            q_lag = (1/(1 + r))*((ī + g.ϑ)/(1 + π_) + (1 - g.ϑ)/(1 + π_) * q) # This is q_t-1, which is a function of r_t. Needed to recover b_t

            p = ProgressUnknown()
            f = (τ_n) -> _τ_n_excess_port_!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, yprime, d, n_ss, Γ_n, covid, z_dist, r, verbose = false, p).excess_g
            τ_n = maximum(find_zeros(f, 0.0, 2.0))
            (; w, y, n, k, lab_tax, lump_sum_T, lump_sum_Tu) = _τ_n_excess_port_!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, yprime, d, n_ss, Γ_n, covid, z_dist, r, verbose = false, p)
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
            path.bstar[i] = qbstar * path.b[i]# Note with qstar as a ratio the prices cancels (qbstar / q_lag)* y. Here qbstar = B^star / B, so qbstar*path.b[i] = B^star. 
            path.Ψ[i] = path.bstar[i] / q_lag^(-χ)
            path.bd[i] = path.b[i] - path.bstar[i] 
            path.n[i] = n
            path.y[i] = y
            path.X[i] = X
            path.C[i] = path.y[i] - path.X[i] - path.G[i] - (1 + r) * qbstar * qb * y + qbstar_prime * qbprime * yprime # Recall I feed ratios of next period GDP (q_t-1 B_t / Y_t)
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
            backwards_once_!(h, g, path[i]; next_values = next_values, R, w, τ_n, T = T_, Tu = Tu_, n_ss, Γ_n, c_ss, Γ_c, covid, path[i].a_tmp)
            asset_policy_given_η_!(path[i].a_pol; h, g, path[i].η, R, w, τ_n, T = T_, Tu = Tu_, n_ss, Γ_n, covid)
            interpolate_asset_policy!(path[i].lower_index, path[i].lower_weight; h.a_grid, path[i].a_pol)
        end
    end
    # Period 1
    forward_pdf!(path[1].pdf; h = h_path[1], e_init.ws.pdf, e_init.ws.lower_index, e_init.ws.lower_weight)
    h = h_path[2] # Given distribution in period 0 (with high unemployment), transition matrix here (to period 1) is h_path[2].P
    z_dist = z_state[1,:]
    r_k_hat = r_path_rest[1] # Note here that we are storing the return on capital in period 1 (after the MIT shock) in r_path_rest[2]
    r_k = r_k_hat + t.δ
    r_b = r_path_rest[2] # Note here that we are storing the return on bonds in period 1 (after the MIT shock) in r_path_rest[2]
    R_k = 1 + r_k
    R_b = 1 + r_b
    qb, qbstar, τ_k, T, Tu, G, π_, d, covid = qb_path[1], qbstar_path[1], τ_k_path[1], T_path[1], Tu_path[1], G_path[1], π_path[1], d_path[1], covid_path[1] #Here qfb is in levels, but the rest of qb_path is a ratio
    qbprime = qb_path[2]
    qbstar_prime = qbstar_path[2]
    qprime = path.q[2]
    yprime = path.y[2]

    rprime = path.r[2] # This refers to r_t+1
    π_prime = path.π_[2] # Here [i+1] refers to period t
    q = (1/(1 + rprime))*((ī + g.ϑ)/(1 + π_prime) + (1 - g.ϑ)/(1 + π_prime) * qprime) # Only used to compute q_t recover the face value of b_t+1
    bprime_level = (qbprime / q) * yprime
    q_lag = (1/(1 + r_b))*((ī + g.ϑ)/(1 + π_) + (1 - g.ϑ)/(1 + π_) * q) # This is q_t-1, which is a function of r_t. Needed to recover b_t
    println("q_lag = ", q_lag)
    p = ProgressUnknown()
    f = (τ_n) -> _τ_n_excess_t0_port_!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, yprime, d, n_ss, Γ_n, covid, z_dist, r_b, r_k, verbose = false, p).excess_g
    τ_n = maximum(find_zeros(f, 0.0, 5.0))
    (; w, y, n, k, lab_tax, lump_sum_T, lump_sum_Tu) = _τ_n_excess_t0_port_!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, yprime, d, n_ss, Γ_n, covid, z_dist, r_b, r_k, verbose = false, p)
    X = path.k[2] - (1.0 - δ) * k

    b_level = qb / q_lag #Note here it is the level, and not a ratio
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
    backwards_once_t0_T!(h, g, path[1]; next_values = next_values, R_k, R_b, w, τ_n, T = T_, Tu = Tu_, n_ss, Γ_n, c_ss, Γ_c, covid, path[1].a_tmp, e_init, method)
    interpolate_asset_policy!(path[1].lower_index, path[1].lower_weight; h.a_grid, path[1].a_pol)

    # Forward iterations
    # Initial portfolio
    if method ==1
        k_agg = e_init.k
        b_agg = e_init.b # This is aggregate demand for bonds (e_init.b - b_star is the domestic demand)
    elseif method == 2
        Λ_ka = e_init.k / e_init.a
        k_agg = sum((Λ_ka.*h.a_grid)  .* path[1].pdf)
        b_agg = sum(((1.0 - Λ_ka).*h.a_grid) .* path[1].pdf)
    end
    #k_agg = aggregate_capital(h.a_grid, path[1].pdf) # Aggregate capital in period 0
    #b_agg = aggregate_bonds(h.a_grid, path[1].pdf) # Aggregate bonds in period 0
    for i in 1:length(path.r)+1
        if i==1 # Find r_k to clear the capital market in period 0 after the MIT shock
            path.un[1] = sum(path[1].pdf, dims = 1)[1]
            a = asset_supply(h_path[i].a_grid, path[i].pdf)
            path.a[i] = a # Store path of assets in period 0 (predetermined by the initial SS)          
            residuals[1] =  k_agg - path.k[1] # Period 0, index 1
        elseif i==2 # Find r_b to clear the bonds market in period 0 after MIT shock
            residuals[2] =  b_agg - path.b[1] # Period 0, index 1
        else
            a = asset_supply(h_path[i-1].a_grid, path[i-1].pdf)
            path.a[i-1] = a
            residuals[i] =  a - path.k[i-1] - path.q[i - 2] * path.bd[i-1] # Period 1, index 2
        end
        (i == length(path.r)+1) && break

        if i > 1
            forward_pdf!(path[i].pdf; h = h_path[i], path[i-1].pdf, path[i-1].lower_index, path[i-1].lower_weight)
            path.un[i] = sum(path[i].pdf, dims = 1)[1]
        end
    end
end 

# ============================================================
# Checking previous code is correct. Running the transition in levels
# ============================================================
function solve_transition_h_l(
    e_init, e_final; τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, t, g, method, r_path = nothing, verbose = true,
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

    path = _initialize_path(e_init, e_final; qb_path)

    if isnothing(r_path) 
        r_path_rest = fill(e_final.r, length(τ_k_path)+1) # The length is increased by 1, as we need to solve for r_b and r_k in period 0
    else 
        r_path_rest = r_path
    end 
    f! = (F, x) -> _residuals_h_l!(F, x; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, t, g, method)
    r_sol = nlsolve(f!, r_path_rest; nlsolve_kwargs_merged...)
    _residuals_h_l!(similar(r_sol.zero), r_sol.zero; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, t, g, method)
    
    return path 
end


function _residuals_h_l!(residuals, r_path_rest; path, e_init, e_final, τ_k_path, T_path, Tu_path, G_path, qb_path, qbstar_path, h_path, π_path, d_path, covid_path, Γ_n, Γ_c, t, g, method)
    #r̄ = e_init.r # Return on bonds in the initial SS that makes the bond price equal to 1 (normalization of the coupons)
    π_bar = e_init.π_
    ī = (1 + e_init.r)*(1 + π_bar) - 1
    A = t.A
    α = t.α
    δ = t.δ
    β = h_path[1].u.β
    r_min = - δ + eps()  # minimum interest rate
    r_max = 1 / β  # a maximum -- could be relaxed
    n_ss = e_init.ws.n_pol[1,:]
    c_ss = e_init.ws.c_pol

    z_state = zeros(length(h_path), length(e_init.h.z_grid))
    z_state[1,:] = e_init.h.Pss'*h_path[1].P
    for i = 2:lastindex(h_path)
        z_state[i,:] = z_state[i - 1,:]' * h_path[i].P 
    end
    #Backward iteration
    for i in reverse(eachindex(path, qb_path, τ_k_path))
        if i > firstindex(path)
            z_dist = z_state[i, :]
            h = i == lastindex(path) ? h_path[end] : h_path[i + 1] # In period 599, I use h_path[600] to compute transition from states in 599 (z_it) to 600 (z_it+1)
            qb, qbstar, τ_k, T, Tu, G, π_, d, covid = qb_path[i], qbstar_path[i], τ_k_path[i], T_path[i], Tu_path[i], G_path[i], π_path[i], d_path[i], covid_path[i]
            qbprime = i == lastindex(path) ? e_final.q * e_final.b : qb_path[i + 1] 
            qbstar_prime = i == lastindex(path) ?  e_final.bstar / e_final.b : qbstar_path[i + 1] # qbstar_path[end] : qbstar_path[i + 1] 
            r = r_path_rest[i+1] # Note r_path_rest has one more entry, as we need to solve for two returns in the first period. However [i+1] here refers to period t
            r = clamp(r, r_min, r_max) # This is r_t. The return in t of the bonds that I bought in t-1
            R = 1 + r

            qprime = i == lastindex(path) ? e_final.q : path.q[i + 1] # Here [i+1] refers to period t
            rprime = i == lastindex(path) ? e_final.r : path.r[i + 1] # This refers to r_t+1
            π_prime = i == lastindex(path) ? e_final.π_ : path.π_[i + 1] # Here [i+1] refers to period t
            q = (1/(1 + rprime))*((ī + g.ϑ)/(1 + π_prime) + (1 - g.ϑ)/(1 + π_prime) * qprime) # Only used to compute q_t recover the face value of b_t+1
            bprime_level = qbprime / q
            q_lag = (1/(1 + r))*((ī + g.ϑ)/(1 + π_) + (1 - g.ϑ)/(1 + π_) * q) # This is q_t-1, which is a function of r_t. Needed to recover b_t

            p = ProgressUnknown()
            f = (τ_n) -> _τ_n_excess_port_l!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, d, n_ss, Γ_n, covid, z_dist, r, verbose = false, p).excess_g
            τ_n = maximum(find_zeros(f, 0.0, 2.0))
            (; w, y, n, k, lab_tax, lump_sum_T, lump_sum_Tu) = _τ_n_excess_port_l!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, d, n_ss, Γ_n, covid, z_dist, r, verbose = false, p)
            X = i == lastindex(path) ? δ * k : path.k[i + 1] - (1.0 - δ) * k

            b_level = qb / q_lag
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
            path.bd[i] = path.b[i] - path.bstar[i] 
            path.n[i] = n
            path.y[i] = y 
            path.X[i] = X
            path.C[i] = path.y[i] - path.X[i] - path.G[i] - (1 + r) * qbstar * qb + qbstar_prime * qbprime # Recall I feed ratios of next period GDP (q_t-1 B_t / Y_t)
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
            backwards_once_!(h, g, path[i]; next_values = next_values, R, w, τ_n, T = T_, Tu = Tu_, n_ss, Γ_n, c_ss, Γ_c, covid, path[i].a_tmp)
            asset_policy_given_η_!(path[i].a_pol; h, g, path[i].η, R, w, τ_n, T = T_, Tu = Tu_, n_ss, Γ_n, covid)
            interpolate_asset_policy!(path[i].lower_index, path[i].lower_weight; h.a_grid, path[i].a_pol)
        end
    end
    # Period 1
    forward_pdf!(path[1].pdf; h = h_path[1], e_init.ws.pdf, e_init.ws.lower_index, e_init.ws.lower_weight)
    h = h_path[2] # Given distribution in period 0 (with high unemployment), transition matrix here (to period 1) is h_path[2].P
    z_dist = z_state[1,:]
    r_k_hat = r_path_rest[1] # Note here that we are storing the return on capital in period 1 (after the MIT shock) in r_path_rest[2]
    r_k = r_k_hat + t.δ
    r_b = r_path_rest[2] # Note here that we are storing the return on bonds in period 1 (after the MIT shock) in r_path_rest[2]
    R_k = 1 + r_k
    R_b = 1 + r_b
    qb, qbstar, τ_k, T, Tu, G, π_, d, covid = qb_path[1], qbstar_path[1], τ_k_path[1], T_path[1], Tu_path[1], G_path[1], π_path[1], d_path[1], covid_path[1] #Here qfb is in levels, but the rest of qb_path is a ratio
    qbprime = qb_path[2]
    qbstar_prime = qbstar_path[2]
    qprime = path.q[2]

    rprime = path.r[2] # This refers to r_t+1
    π_prime = path.π_[2] # Here [i+1] refers to period t
    q = (1/(1 + rprime))*((ī + g.ϑ)/(1 + π_prime) + (1 - g.ϑ)/(1 + π_prime) * qprime) # Only used to compute q_t recover the face value of b_t+1
    #bprime = qbprime / q
    bprime_level = qbprime / q
    q_lag = (1/(1 + r_b))*((ī + g.ϑ)/(1 + π_) + (1 - g.ϑ)/(1 + π_) * q) # This is q_t-1, which is a function of r_t. Needed to recover b_t
    println("q_lag = ", q_lag)
    p = ProgressUnknown()
    f = (τ_n) -> _τ_n_excess_t0_port_l!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, d, n_ss, Γ_n, covid, z_dist, r_b, r_k, verbose = false, p).excess_g
    τ_n = maximum(find_zeros(f, 0.0, 5.0))
    (; w, y, n, k, lab_tax, lump_sum_T, lump_sum_Tu) = _τ_n_excess_t0_port_l!(h, t, g; τ_n, τ_k, T, Tu, G, qbprime, qb, d, n_ss, Γ_n, covid, z_dist, r_b, r_k, verbose = false, p)
    X = path.k[2] - (1.0 - δ) * k

    b_level = qb / q_lag #Note here it is the level, and not a ratio

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
    path.bd[1] = path.b[1] - path.bstar[1] 
    path.n[1] = n
    path.y[1] = y
    path.X[1] = X
    path.C[1] = path.y[1] - path.X[1] - path.G[1] - (1 + r_b) * qbstar + qbstar_prime * qbprime
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
    backwards_once_t0_T!(h, g, path[1]; next_values = next_values, R_k, R_b, w, τ_n, T = T_, Tu = Tu_, n_ss, Γ_n, c_ss, Γ_c, covid, path[1].a_tmp, e_init, method)
    interpolate_asset_policy!(path[1].lower_index, path[1].lower_weight; h.a_grid, path[1].a_pol)

    # Forward iterations
    # Initial portfolio
    if method ==1
        k_agg = e_init.k
        b_agg = e_init.b # This is aggregate demand for bonds (e_init.b - b_star is the domestic demand)
    elseif method == 2
        Λ_ka = e_init.k / e_init.a
        k_agg = sum((Λ_ka.*h.a_grid)  .* path[1].pdf)
        b_agg = sum(((1.0 - Λ_ka).*h.a_grid) .* path[1].pdf)
    end
    for i in 1:length(path.r)+1
        if i==1 # Find r_k to clear the capital market in period 0 after the MIT shock
            path.un[1] = sum(path[1].pdf, dims = 1)[1]
            a = asset_supply(h_path[i].a_grid, path[i].pdf)
            path.a[i] = a # Store path of assets in period 0 (predetermined by the initial SS)          
            residuals[1] =  k_agg - path.k[1] # Period 0, index 1
        elseif i==2 # Find r_b to clear the bonds market in period 0 after MIT shock
            residuals[2] =  b_agg - path.b[1] # Period 0, index 1
        else
            a = asset_supply(h_path[i-1].a_grid, path[i-1].pdf)
            path.a[i-1] = a
            #path.b[i-1] = b_path[i-1]     
            residuals[i] =  a - path.k[i-1] - path.q[i - 2] * path.bd[i-1] # Period 1, index 2
        end
        (i == length(path.r)+1) && break

        if i > 1
            forward_pdf!(path[i].pdf; h = h_path[i], path[i-1].pdf, path[i-1].lower_index, path[i-1].lower_weight)
            path.un[i] = sum(path[i].pdf, dims = 1)[1]
        end
    end
end 

# ============================================================
# Decomposing the Transition: Partial Equilibrium
# ============================================================

function transition_partial(e_init, e_final; r_k, r_path, w_path, τ_n_path, T_path, Tu_path, h_path, covid_path, Γ_n, Γ_c, t, g, method)
    path = _initialize_path(e_init, e_final; qb_path = T_path)
    _residuals_partial!(; path, e_init, e_final, r_k, r_path, w_path, τ_n_path, T_path, Tu_path, h_path, covid_path, Γ_n, Γ_c, t, g, method)
    return path 
end

function _residuals_partial!(; path, e_init, e_final, r_k, r_path, w_path, τ_n_path, T_path, Tu_path, h_path, covid_path, Γ_n, Γ_c, t, g, method)
    π_bar = e_init.π_
    ī = (1 + e_init.r)*(1 + π_bar) - 1
    A = t.A
    α = t.α
    δ = t.δ
    β = h_path[1].u.β
    r_min = - δ + eps()  # minimum interest rate
    r_max = 1 / β  # a maximum -- could be relaxed
    n_ss = e_init.ws.n_pol[1,:]
    c_ss = e_init.ws.c_pol

    #Transition of employment states
    #From period 0 to period 1 I do the transition using h_path[2]
    z_state = zeros(length(h_path), length(e_init.h.z_grid))
    z_state[1,:] = e_init.h.Pss'*h_path[1].P
    for i = 2:lastindex(h_path)
        z_state[i,:] = z_state[i - 1,:]' * h_path[i].P 
    end
    #Backward iteration
    for i in reverse(eachindex(path))
        if i > firstindex(path)
            #println("i= ", i)
            z_dist = z_state[i, :]
            h = i == lastindex(path) ? h_path[end] : h_path[i + 1] # In period 599, I use h_path[600] to compute transition from states in 599 (z_it) to 600 (z_it+1)
            T, Tu, covid =  T_path[i], Tu_path[i], covid_path[i]
            τ_n = τ_n_path[i] #path.τ_n[i]
            r = r_path[i] #path.r[i] 
            w = w_path[i] #path.w[i]
            R = 1 + r
            path.h[i] = h
            next_values = i == lastindex(path) ? e_final.ws : path[i+1]
            backwards_once_!(h, g, path[i]; next_values = next_values, R, w, τ_n, T, Tu, n_ss, Γ_n, c_ss, Γ_c, covid, path[i].a_tmp)
            asset_policy_given_η_!(path[i].a_pol; h, g, path[i].η, R, w, τ_n, T, Tu, n_ss, Γ_n, covid)
            interpolate_asset_policy!(path[i].lower_index, path[i].lower_weight; h.a_grid, path[i].a_pol)
        end
    end
    # Period 1
    forward_pdf!(path[1].pdf; h = h_path[1], e_init.ws.pdf, e_init.ws.lower_index, e_init.ws.lower_weight)
    h = h_path[2] # Given distribution in period 0 (with high unemployment), transition matrix here (to period 1) is h_path[2].P
    z_dist = z_state[1,:]

    τ_n = τ_n_path[1] 
    r_b = r_path[1] 
    w = w_path[1] 
    R_k = 1 + r_k
    R_b = 1 + r_b
    T, Tu,covid = T_path[1], Tu_path[1], covid_path[1]
    next_values = path[2]
    backwards_once_t0_T!(h, g, path[1]; next_values = next_values, R_k, R_b, w, τ_n, T, Tu, n_ss, Γ_n, c_ss, Γ_c, covid, path[1].a_tmp, e_init, method)
    interpolate_asset_policy!(path[1].lower_index, path[1].lower_weight; h.a_grid, path[1].a_pol)
end