# ============================================================
# Solving the Household Problem
# ============================================================

# Household's euler equation 

function backwards_euler_x(u::CRRA, p_vec, η_vec)
    ξ = get_inverse_ies(u)
    β = get_β(u)
    Evx = zero(eltype(η_vec))
    for (p, η) in zip(p_vec, η_vec)
        Evx += p * (η)^(-ξ)
    end 
    x = (β * Evx)^(-1/ξ)
    return x
end


function backwards_euler_x(u::EZ, p_vec, v_vec, η_vec)
    ξ = get_inverse_ies(u)
    γ = get_ra(u)
    β = get_β(u)
    Ev = Evx = zero(eltype(v_vec))
    for (v, η, p) in zip(v_vec, η_vec, p_vec)
        Evx += p * (v)^((ξ - γ)) * (η)^(-ξ)
        Ev += p * u.risk(v)
    end 
    x = (β * inverse(u.risk, Ev^((γ - ξ))) * Evx)^(-1/ξ)
    return x
end


_backwards_euler_x_helper(u::CRRA, P′, i, s, next_values) = backwards_euler_x(u, view(P′, :, s), view(next_values.η, i, :))
_backwards_euler_x_helper(u::EZ, P′, i, s, next_values) = backwards_euler_x(u, view(P′, :, s), view(next_values.v, i, :), view(next_values.η, i, :))


function backwards_once!(h, g, current_values; next_values, R, w, τ_n, T, Tu, a_tmp = similar(current_values.η))   
    η_now = current_values.η
    (; u) = h
    P′ = h.Pprime
    ξ = get_inverse_ies(u)
    R_ = R^(-1/ξ)
    τ = g.τ
    ind_T = g.ind_T
    ind_Tu = g.ind_Tu

    a_min = first(h.a_grid)
    a_max = last(h.a_grid)

    @batch for s in eachindex(h.z_grid)
        z = h.z_grid[s]
        lab = labor_income(h.v, g; w = w * z, τ_n = τ_n)
        lab_tax = lab - τ_n * lab ^ (1-τ)  # HSV(2017)
        if ind_T[s]==1 && ind_Tu[s]==1
            add = lab - lab_tax - disutility_given_w(h.v, g; w = w * z, τ_n = τ_n) + T + Tu
        elseif ind_T[s]==1 && ind_Tu[s]==0
            add = lab - lab_tax - disutility_given_w(h.v, g; w = w * z, τ_n = τ_n) + T
        elseif ind_T[s]==0 && ind_Tu[s]==1
            add = lab - lab_tax - disutility_given_w(h.v, g; w = w * z, τ_n = τ_n) + Tu
        else
            add = lab - lab_tax - disutility_given_w(h.v, g; w = w * z, τ_n = τ_n)
        end

        for i in eachindex(h.a_grid)
            x = _backwards_euler_x_helper(h.u, P′, i, s, next_values)
            a_next = h.a_grid[i]
            a_tmp[i, s]  = (x + a_next - add) / R # a_2(a_3,z_2)
        end

        for i in eachindex(h.a_grid)
            a = h.a_grid[i] # R_2 a_2
            a_prime = interp1D(a, view(a_tmp, :, s), h.a_grid) #a_3(a_2,z_2)
            a_prime = clamp(a_prime, a_min, a_max)
            x = R * a + add - a_prime
            η_now[i, s] = R_ * x
            _extra_updates!(u, current_values, next_values, h, a_prime, P′, x, i, s)
        end
    end    
end

_extra_updates!(::CRRA, args...) = nothing

function _extra_updates!(u::EZ, current_values, next_values, h, a_prime, P′, x, i, s)
    β = u.β

    Ev = zero(a_prime)
    for s2 in axes(next_values.v, 2)
        v_prime = interp1D(a_prime, h.a_grid, view(next_values.v, :, s2))
        Ev += P′[s2, s] * u.risk(v_prime)
    end
    current_values.v[i, s] =
        inverse(u.temporal, (1 - β) * u.temporal(x) + β * u.temporal(inverse(u.risk, Ev)))
end 


function asset_policy_given_η!(a_pol; h, g, η, R, w, τ_n, T, Tu)
    ξ = get_inverse_ies(h.u)
    τ = g.τ
    ind_T = g.ind_T
    ind_Tu = g.ind_Tu
    @batch for i in eachindex(h.a_grid)
        for s in eachindex(h.z_grid)
            if ind_T[s]==1 && ind_Tu[s]==1
                a_pol[i, s] = R * h.a_grid[i] + τ_n * labor_income(h.v, g; w = w * h.z_grid[s], τ_n = τ_n)^(1-τ) -  η[i, s] * R^(1/ξ) - disutility_given_w(h.v, g; w = w * h.z_grid[s], τ_n = τ_n) + T + Tu
            elseif ind_T[s]==1 && ind_Tu[s]==0
                a_pol[i, s] = R * h.a_grid[i] + τ_n * labor_income(h.v, g; w = w * h.z_grid[s], τ_n = τ_n)^(1-τ) -  η[i, s] * R^(1/ξ) - disutility_given_w(h.v, g; w = w * h.z_grid[s], τ_n = τ_n) + T
            elseif ind_T[s]==0 && ind_Tu[s]==1
                a_pol[i, s] = R * h.a_grid[i] + τ_n * labor_income(h.v, g; w = w * h.z_grid[s], τ_n = τ_n)^(1-τ) -  η[i, s] * R^(1/ξ) - disutility_given_w(h.v, g; w = w * h.z_grid[s], τ_n = τ_n) + Tu
            else
                a_pol[i, s] = R * h.a_grid[i] + τ_n * labor_income(h.v, g; w = w * h.z_grid[s], τ_n = τ_n)^(1-τ) -  η[i, s] * R^(1/ξ) - disutility_given_w(h.v, g; w = w * h.z_grid[s], τ_n = τ_n)
            end
        end
    end
end
asset_policy_given_η!(ws::HouseholdWorkspace, g ; R, w, τ_n, T, Tu) = asset_policy_given_η!(ws.a_pol; ws.h, g, ws.η, R, w, τ_n, T, Tu) 


function interpolate_asset_policy!(lower_index, lower_weight; a_grid, a_pol)
    @batch for i in eachindex(a_pol, lower_index, lower_weight)
        a = a_pol[i]
        j = searchsortedfirst(a_grid, a)
        ind = clamp(j, firstindex(a_grid) + 1, lastindex(a_grid))
        l_weight =  (a_grid[ind] - a) / (a_grid[ind] - a_grid[ind - 1]) 
        lower_index[i] = ind - 1
        lower_weight[i] = clamp(l_weight, zero(l_weight), one(l_weight))
    end 
end
interpolate_asset_policy!(ws::HouseholdWorkspace) = interpolate_asset_policy!(ws.lower_index, ws.lower_weight; ws.h.a_grid, ws.a_pol)


function forward_pdf!(pdf_next; h, pdf, lower_index, lower_weight)
    P′ = h.Pprime
    fill!(pdf_next, zero(eltype(pdf_next)))
    for s in axes(pdf, 2)
        for i in axes(pdf, 1)
            ind = lower_index[i, s]
            weight = lower_weight[i, s]
            for s′ in axes(pdf, 2)
                mass = P′[s′, s] * pdf[i, s]
                p = mass * weight 
                pdf_next[ind, s′] += p
                pdf_next[ind + 1, s′] += mass - p
            end 
        end 
    end 
end 

function update_stationary_pdf!(ws; tol = _TOL, max_iters = 10_000, print_every = 100, verbose = false)
    (; h, pdf, pdf_tmp, lower_index, lower_weight) = ws 
    pdf_0, pdf_1 = pdf, pdf_tmp
    iter = 1 
    while true 
        forward_pdf!(pdf_1; h, pdf = pdf_0, lower_index, lower_weight)
        distance = chebyshev(pdf_0, pdf_1)
        verbose && (iter % print_every == 1) && println("iter $iter: $distance")
        if (distance < tol) 
            verbose && println("Converged. Iter $iter: $distance")
            break
        end 
        iter += 1 
        (iter > max_iters) && (@warn "PDF did not converge!"; break) 
        pdf_0, pdf_1 = pdf_1, pdf_0 
    end 
    pdf .= pdf_0
end 


asset_supply(a_grid, pdf) = @tullio s := a_grid[i] * pdf[i, s]

function aggregate_c(h::Household, g::TaxSys; R, w, τ_n, T, a_pol, pdf)
    ind_T = g.ind_T
    agg_c = zero(eltype(z_grid))
    @batch for i in eachindex(h.a_grid)
        for s in eachindex(h.z_grid)
            if ind_T[s]==1
                agg_c += (R * h.a_grid[i] + τ_n * labor_income(h.v, g; w = w * h.z_grid[s], τ_n = τ_n)^(1-g.τ) - a_pol[i, s] + T) * pdf[i, s]
            else
                agg_c += (R * h.a_grid[i] + τ_n * labor_income(h.v, g; w = w * h.z_grid[s], τ_n = τ_n)^(1-g.τ) - a_pol[i, s]) * pdf[i, s]
            end
        end
    end
    return agg_c
end

function labor_supply(h::Household, g::TaxSys; w, τ_n)
    (; v, z_grid, Pss) = h
    n = zero(eltype(z_grid))
    for (p, z) in zip(Pss, z_grid)
        n += p * z * labor(v, g; w = w * z, τ_n = τ_n)
    end
    return n
end

function cons_labor_policy!(c_pol, n_pol, lab_income, income, lab_tax, a_pol, h, g; R, w, τ_n, T, Tu)
    (; z_grid,  a_grid) = h
    ind_T = g.ind_T
    ind_Tu = g.ind_Tu
    @batch for s in axes(c_pol, 2)
        for i in axes(c_pol, 1)
            lab_income[i, s] = labor_income(h.v, g; w = w * z_grid[s], τ_n = τ_n)
            lab_tax[i, s] = τ_n * lab_income[i, s] ^ (1 - g.τ)
            n_pol[i, s] = labor(h.v, g; w = w * z_grid[s], τ_n = τ_n)
            if ind_T[s]==1 && ind_Tu[s]==1
                income[i, s] = R * a_grid[i] + τ_n * labor_income(h.v, g; w = w * z_grid[s], τ_n = τ_n) ^ (1 - g.τ) + T + Tu
            elseif ind_T[s]==1 && ind_Tu[s]==0
                income[i, s] = R * a_grid[i] + τ_n * labor_income(h.v, g; w = w * z_grid[s], τ_n = τ_n) ^ (1 - g.τ) + T
            elseif ind_T[s]==0 && ind_Tu[s]==1
                income[i, s] = R * a_grid[i] + τ_n * labor_income(h.v, g; w = w * z_grid[s], τ_n = τ_n) ^ (1 - g.τ) + Tu
            else
                income[i, s] = R * a_grid[i] + τ_n * labor_income(h.v, g; w = w * z_grid[s], τ_n = τ_n) ^ (1 - g.τ)
            end
            c_pol[i, s]  = income[i, s] - a_pol[i, s]
        end
    end
end
cons_labor_policy!(ws::HouseholdWorkspace, g; R, w, τ_n, T, Tu) = cons_labor_policy!(ws.c_pol, ws.n_pol, ws.lab_income, ws.income, ws.lab_tax, ws.a_pol, ws.h, g; R, w, τ_n, T, Tu) 

is_pol_valid(a_pol, h) = all(a_pol .< last(h.a_grid))
is_pol_valid(e::StationaryEquilibrium) = is_pol_valid(e.ws.a_pol, e.h)



function stationary(h, g; R, w, τ_n, T, Tu, kwargs...)
    ws = HouseholdWorkspace(;h, g, R, w, τ_n, T, Tu)
    stationary!(ws, g; R, w, τ_n, T, Tu, kwargs...)
    return ws
end 


stationary!(ws, g; R, w, τ_n, T, Tu, kwargs...) = stationary!(ws.h.u, ws, g; R, w, τ_n, T, Tu, kwargs...)


function stationary!(
    ::EZ, ws, g; 
    R, w, τ_n, T, Tu,
    max_iters = _MAX_ITERS, print_every = 10, 
    value_tol = _TOL,
    policy_tol = _TOL,
    pdf_tol = _TOL, 
    verbose = true
)

    (; h, v, η, v_tmp, η_tmp, a_tmp) = ws
    v_η_0 = (; v = v, η = η)
    v_η_1 = (; v = v_tmp, η = η_tmp)
    i = 1
    while true
        backwards_once!(h, g, v_η_1; next_values = v_η_0, R, w, τ_n, T, Tu, a_tmp)
        dis1 = chebyshev(v_η_0.v, v_η_1.v)
        if dis1 < value_tol 
            # checked that policy converged too 
            dis2 = chebyshev(v_η_0.η, v_η_1.η)
            if dis2 < policy_tol
                verbose && println("Converged. iter $i;  R = $R, τ_n = $τ_n, w = $w, |v|: $dis1 , |η| : $dis2")
                break
            end
        end 
        if i > max_iters
            dis2 = (v_η_0.η, v_η_1.η)
            @warn("Max iters reached. DID NOT CONVERGE! R = $R, τ_n = $τ_n, w = $w,  |v|: $dis1 , |η| : $dis2")
            break
        end 
        verbose && (i % print_every == 1) && println("iter $i: $dis1") 
        i += 1
        v_η_0, v_η_1 = v_η_1, v_η_0
    end
    ws.v .= v_η_0.v
    ws.η .= v_η_0.η

    asset_policy_given_η!(ws, g; R, w, τ_n, T, Tu)
    cons_labor_policy!(ws, g; R, w, τ_n, T, Tu)
    interpolate_asset_policy!(ws)
    update_stationary_pdf!(ws; tol = pdf_tol)
end 


function stationary!(
    ::CRRA, ws, g; 
    R, w, τ_n, T, Tu,
    max_iters = _MAX_ITERS, print_every = 100, 
    policy_tol = _TOL,
    pdf_tol = _TOL, 
    verbose = true
)

    (; h, η, η_tmp, a_tmp) = ws
    η_0 = η
    η_1 = η_tmp
    i = 1
    while true
        backwards_once!(h, g, (; η = η_1); next_values = (; η = η_0), R, w, τ_n, T, Tu, a_tmp)
        dis1 = chebyshev(η_0, η_1)
        if dis1 < policy_tol 
            verbose && println("Converged. iter $i;  R = $R, τ_n = $τ_n, w = $w, |η|: $dis1")
            break
        end 
        if i > max_iters
            @warn("Max iters reached. DID NOT CONVERGE! R = $R, τ_n = $τ_n, w = $w,  |η|: $dis1")
            break
        end 
        verbose && (i % print_every == 1) && println("iter $i: $dis1") 
        i += 1
        η_0, η_1  = η_1, η_0
    end
    # η stores the final answer
    ws.η .= η_0

    asset_policy_given_η!(ws, g; R, w, τ_n, T, Tu)
    interpolate_asset_policy!(ws)
    update_stationary_pdf!(ws; tol = pdf_tol)
end

# ============================================================
# Household problem during the Pandemic
# ============================================================

function backwards_once_!(h, g, current_values; next_values, R, w, τ_n, T, Tu, n_ss, Γ_n, c_ss, Γ_c, covid, a_tmp = similar(current_values.η))   
    η_now = current_values.η
    (; u) = h
    P′ = h.Pprime
    ξ = get_inverse_ies(u)
    R_ = R^(-1/ξ)
    τ = g.τ
    ind_T = g.ind_T
    ind_Tu = g.ind_Tu

    a_min = first(h.a_grid)
    a_max = last(h.a_grid)

    @batch for s in eachindex(h.z_grid) #
        z = h.z_grid[s]
        lab = labor_income_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)
        lab_tax = lab - τ_n * lab ^ (1-τ)  # HSV(2017)
        if ind_T[s]==1 && ind_Tu[s]==1
            add = lab - lab_tax - disutility_given_w_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid) + T + Tu
        elseif ind_T[s]==1 && ind_Tu[s]==0
            add = lab - lab_tax - disutility_given_w_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid) + T
        elseif ind_T[s]==0 && ind_Tu[s]==1
            add = lab - lab_tax - disutility_given_w_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid) + Tu
        else
            add = lab - lab_tax - disutility_given_w_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)
        end
        
        
        for i in eachindex(h.a_grid)
            x = _backwards_euler_x_helper(h.u, P′, i, s, next_values)
            a_next = h.a_grid[i]
            a_tmp[i, s]  = (x + a_next - add) / R # a_2(a_3,z_2)
        end

        for i in eachindex(h.a_grid)
            a = h.a_grid[i] # R_2 a_2
            a_prime = interp1D(a, view(a_tmp, :, s), h.a_grid) #a_3(a_2,z_2)
            a_prime = clamp(a_prime, a_min, a_max)
            if covid == 0.0
                x = R * a + add - a_prime
            elseif covid == 1.0
                x =  Γ_c * c_ss[i, s] - disutility_given_w_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)
            end
            current_values.c_pol[i, s] = x + disutility_given_w_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)
            current_values.n_pol[i, s] = labor_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)
            η_now[i, s] = R_ * x
            _extra_updates!(u, current_values, next_values, h, a_prime, P′, x, i, s)
        end
    end    
end

function asset_policy_given_η_!(a_pol; h, g, η, R, w, τ_n, T, Tu, n_ss, Γ_n, covid)
    ξ = get_inverse_ies(h.u)
    τ = g.τ
    ind_T = g.ind_T
    ind_Tu = g.ind_Tu
    @batch for i in eachindex(h.a_grid)
        for s in eachindex(h.z_grid)
            if ind_T[s]==1 && ind_Tu[s]==1
                a_pol[i, s] = R * h.a_grid[i] + τ_n * labor_income_(h.v, g; w = w * h.z_grid[s], τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)^(1-τ) -  η[i, s] * R^(1/ξ) - disutility_given_w_(h.v, g; w = w * h.z_grid[s], τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid) + T + Tu
            elseif ind_T[s]==1 && ind_Tu[s]==0
                a_pol[i, s] = R * h.a_grid[i] + τ_n * labor_income_(h.v, g; w = w * h.z_grid[s], τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)^(1-τ) -  η[i, s] * R^(1/ξ) - disutility_given_w_(h.v, g; w = w * h.z_grid[s], τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid) + T
            elseif ind_T[s]==0 && ind_Tu[s]==1
                a_pol[i, s] = R * h.a_grid[i] + τ_n * labor_income_(h.v, g; w = w * h.z_grid[s], τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)^(1-τ) -  η[i, s] * R^(1/ξ) - disutility_given_w_(h.v, g; w = w * h.z_grid[s], τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid) + Tu
            else
                a_pol[i, s] = R * h.a_grid[i] + τ_n * labor_income_(h.v, g; w = w * h.z_grid[s], τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)^(1-τ) -  η[i, s] * R^(1/ξ) - disutility_given_w_(h.v, g; w = w * h.z_grid[s], τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)
            end
        end
    end
end
asset_policy_given_η_!(ws::HouseholdWorkspace, g ; R, w, τ_n, T, Tu, n_ss, Γ_n, covid) = asset_policy_given_η_!(ws.a_pol; ws.h, g, ws.η, R, w, τ_n, T, Tu, n_ss, Γ_n, covid) 

function labor_supply_(h::Household, g::TaxSys; w, τ_n, n_ss, Γ_n, covid, z_dist)
    (; v, z_grid) = h
    n = zero(eltype(z_grid))
    for (p, z, n_g) in zip(z_dist, z_grid, n_ss)
        n += p * z * labor_(v, g; w = w * z, τ_n = τ_n, n_ = n_g, Γ_n = Γ_n, covid = covid)
    end
    return n
end

function cons_labor_policy_!(c_pol, n_pol, lab_income, lab_tax, a_pol, h, g; R, w, τ_n, T, Tu, n_ss, Γ_n, covid)
    (; z_grid,  a_grid) = h
    ind_T = g.ind_T
    ind_Tu = g.ind_Tu
    @batch for s in axes(c_pol, 2)
        for i in axes(c_pol, 1)
            lab_income[i, s] = labor_income_(h.v, g; w = w * z_grid[s], τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid) 
            lab_tax[i, s] = τ_n * lab_income[i, s] ^ (1 - g.τ)
            n_pol[i, s] = labor_(h.v, g; w = w * z_grid[s], τ_n = τ_n , n_ = n_ss[s], Γ_n = Γ_n, covid = covid)
            if ind_T[s]==1 && ind_Tu[s]==1
                c_pol[i, s] = R * a_grid[i] + τ_n * labor_income_(h.v, g; w = w * z_grid[s], τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid) ^ (1 - g.τ) - a_pol[i, s] + T + Tu
            elseif ind_T[s]==1 && ind_Tu[s]==0
                c_pol[i, s] = R * a_grid[i] + τ_n * labor_income_(h.v, g; w = w * z_grid[s], τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid) ^ (1 - g.τ) - a_pol[i, s] + T
            elseif ind_T[s]==0 && ind_Tu[s]==1
                c_pol[i, s] = R * a_grid[i] + τ_n * labor_income_(h.v, g; w = w * z_grid[s], τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid) ^ (1 - g.τ) - a_pol[i, s] + Tu
            else
                c_pol[i, s] = R * a_grid[i] + τ_n * labor_income_(h.v, g; w = w * z_grid[s], τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid) ^ (1 - g.τ) - a_pol[i, s]
            end
        end
    end
end
cons_labor_policy_!(ws::HouseholdWorkspace, g; R, w, τ_n, T, Tu, n_ss, Γ_n, covid) = cons_labor_policy_!(ws.c_pol, ws.n_pol, ws.lab_income, ws.lab_tax, ws.a_pol, ws.h, g; R, w, τ_n, T, Tu, n_ss, Γ_n, covid) 

function T_rev!(e_init, R_k, R_b)
    δ = e_init.t.δ
    r_ss = e_init.r
    a_ss = e_init.a
    a_grid = round.(e_init.h.a_grid, digits = 10)
    T_rev = similar(a_grid)

    pdf_marginal = sum(e_init.ws.pdf, dims=2)[:, 1]
    qu_10 = round.([quantile(a_grid, StatsBase.ProbabilityWeights(pdf_marginal), x) for x in 0.0:0.1:1.0], digits = 10)

    Λ_bi = [-0.015; -0.027; -0.072; -0.11; -0.12; -0.12; -0.11; -0.06; 0.06; NaN]; # b_i / a_ss for each decile 
    dist_hhs = [sum(pdf_marginal[a_grid.<=qu_10[qi]]) for qi in 2:lastindex(qu_10)] # Starts from 2 to remove the zero. This has size 10, roughly [0.1 0.2 0.3 ... 1]
    mass_ = [dist_hhs[1] ; dist_hhs |> diff]
    Λ_bi[end] = (e_init.bd - sum(mass_[1:end-1] .* Λ_bi[1:end-1] .* a_ss)) / (mass_[end] * a_ss) # Make sure that the last decile guarantees SS debt 
    
    #Make capital shares consistent with wealth of each decile
    a_decile = [sum((a_grid .<= qi) .* (a_grid .* pdf_marginal)) for qi in qu_10] |> diff # This has size 10. 
    a_levels = a_decile ./ mass_
    k_levels = a_levels .- (Λ_bi .* a_ss)
    Λ_ki = k_levels ./ a_ss

    for i in eachindex(a_grid)
        idx = min(searchsortedfirst(qu_10, a_grid[i]), length(qu_10)) # I added this because for high values of a_grid, they were greater than the max of qu_10. So searchsorted gives 12, and can be at most 11.
        #println(idx)
        T_rev[i] = idx == 1 ? (R_k - δ - (1 + r_ss))* Λ_ki[1] * a_ss + (R_b - (1 + r_ss)) * Λ_bi[1] * a_ss : (R_k - δ - (1 + r_ss))* Λ_ki[idx - 1] * a_ss + (R_b - (1 + r_ss)) * Λ_bi[idx -  1] * a_ss
    end
    return T_rev
end

function backwards_once_t0_T!(h, g, current_values; next_values, R_k, R_b, w, τ_n, T, Tu, n_ss, Γ_n, c_ss, Γ_c, covid, a_tmp = similar(current_values.η), e_init, method)   
    η_now = current_values.η
    (; u) = h
    P′ = h.Pprime
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

    for s in eachindex(h.z_grid) #@batch 
        z = h.z_grid[s]
        lab = labor_income_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)
        lab_tax = lab - τ_n * lab ^ (1-τ)  # HSV(2017)
        if ind_T[s]==1 && ind_Tu[s]==1
            add = lab - lab_tax - disutility_given_w_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid) + T + Tu
        elseif ind_T[s]==1 && ind_Tu[s]==0
            add = lab - lab_tax - disutility_given_w_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid) + T
        elseif ind_T[s]==0 && ind_Tu[s]==1
            add = lab - lab_tax - disutility_given_w_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid) + Tu
        else
            add = lab - lab_tax - disutility_given_w_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)
        end

        for i in eachindex(h.a_grid)
            x = _backwards_euler_x_helper(h.u, P′, i, s, next_values)
            a_next = h.a_grid[i]
            T_r = T_rev[i]
            if method == 1
                a_tmp[i, s]  = (x + a_next - add - T_r) / R
            elseif method == 2
                a_tmp[i, s]  = (x + a_next - add) / (Λ_ka*(R_k - δ) + (1.0 - Λ_ka)*R_b)             
            end
        end

        for i in eachindex(h.a_grid)
            a = h.a_grid[i] # R_2 a_2
            if covid == 0.0
                a_prime = interp1D(a, view(a_tmp, :, s), h.a_grid) #a_3(a_2,z_2)
                a_prime = clamp(a_prime, a_min, a_max)
                T_r = T_rev[i]
                if method == 1
                    x = R * a + add + T_r - a_prime
                elseif method == 2
                    x = (Λ_ka*(R_k - δ) + (1.0 - Λ_ka)*R_b) * a + add - a_prime
                end
                current_values.c_pol[i, s] = x + disutility_given_w_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)
                current_values.a_pol[i,s] = a_prime
                _extra_updates!(u, current_values, next_values, h, a_prime, P′, x, i, s)
            elseif covid == 1.0
                T_r = T_rev[i]
                x = Γ_c * c_ss[i, s] - disutility_given_w_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)
                a_prime = R * a + add + T_r - x
                if a_prime < a_min
                    a_prime = a_min
                    x = R * a + add + T_r - a_prime
                end
                current_values.c_pol[i, s] = x + disutility_given_w_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)
                current_values.a_pol[i,s] = a_prime
                _extra_updates!(u, current_values, next_values, h, a_prime, P′, x, i, s)                
            end
            current_values.n_pol[i, s] = labor_(h.v, g; w = w * z, τ_n = τ_n, n_ = n_ss[s], Γ_n = Γ_n, covid = covid)
        end
    end    
end