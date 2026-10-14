# ============================================================
# Technology methods in mutable structures
# ============================================================

output(t::CobbDouglasTechnologyMut; k, n) = t.A * k^t.α * n^(1- t.α)

mpk(t::CobbDouglasTechnologyMut; k, n) = t.α * t.A * (n / k)^(1 - t.α)

w_from_r(t::CobbDouglasTechnologyMut; r, τ_k) = (1 - t.α) * (t.A / t.μ)^(1/(1-t.α)) * (t.α / ((1 + τ_k) * (r + t.δ)))^(t.α / (1 - t.α))

k_from_r_n(t::CobbDouglasTechnologyMut; r, n, τ_k) = (t.α * t.A / (t.μ * (1 + τ_k) * (r + t.δ)))^(1/(1 - t.α)) * n 

τ_n_from_w(v::GHH, t::CobbDouglasTechnologyMut, g::TaxSys; n_aux, w, k) = (((((1 - t.α) * (t.A / t.μ)) / w)^(1 / t.α) * k) / n_aux)^((1+g.τ*v.ν)/v.ν)

function n_aux(h::Household, g::TaxSys; w)
    (; v, z_grid, Pss) = h
    n = zero(eltype(z_grid))
    for (p, z) in zip(Pss, z_grid)        
        n += p * z * ((1 - g.τ) * (w * z)^(1 - g.τ) / v.θ)^(v.ν / (1 + g.τ * v.ν))
    end
    return n
end

function taxes_HSV(h::Household, g::TaxSys; w, τ_n)
    (; v, z_grid, Pss) = h
    n_hsv = zero(eltype(z_grid))
    for (p, z) in zip(Pss, z_grid)        
        n_hsv += p * (w * z * labor(v, g; w = w * z, τ_n = τ_n))^(1 - g.τ)
    end
    T = w * labor_supply(h, g; w, τ_n) - τ_n * n_hsv
    return T
end

function agg_T(h::Household, g::TaxSys; T)
    (; Pss) = h
    (; ind_T) = g
    agg_T = zero(eltype(Pss))
    for i in eachindex(Pss)
        if ind_T[i]==1        
            agg_T += Pss[i] * T
        end
    end
    return agg_T
end

function agg_Tu(h::Household, g::TaxSys; Tu)
    (; Pss) = h
    (; ind_Tu) = g
    agg_Tu = zero(eltype(Pss))
    for i in eachindex(Pss)
        if ind_Tu[i]==1        
            agg_Tu += Pss[i] * Tu
        end
    end
    return agg_Tu
end

y(e) = output(e.t; e.k, e.n)

## Period 0

w_from_r_k(t::CobbDouglasTechnologyMut; r_k, τ_k) = (1 - t.α) * (t.A / t.μ)^(1/(1-t.α)) * (t.α / ((1 + τ_k) * r_k))^(t.α / (1 - t.α))

k_from_r_k_n(t::CobbDouglasTechnologyMut; r_k, n, τ_k) = (t.α * t.A / (t.μ * (1 + τ_k) * r_k))^(1/(1 - t.α)) * n 

function k_from_rk_hat_n(t::CobbDouglasTechnologyMut; r_k_hat, n, τ_k)
    r_k = r_k_hat + t.δ
    return (t.α * t.A / (t.μ * (1 + τ_k) * r_k))^(1/(1 - t.α)) * n 
end
#Solver finds r_k. Absent MIT shock r_k = r + t.δ, but at t=0 r_k =! r + t.δ.
# Define r_k_hat = r_k - t.δ. In t=0, r_k_hat =! r

# ============================================================
# Technology methods: Dealing with Pandemic Constraints
# ============================================================

output_(t::CobbDouglasTechnologyMut; k, n, d, covid) = covid == 1 ?  t.A * (d * k)^t.α * n^(1- t.α) : t.A * k^t.α * n^(1- t.α)

mpk_(t::CobbDouglasTechnologyMut; k, n, d, covid) = covid == 1 ? t.α * t.A * (n / (d * k))^(1 - t.α) : t.α * t.A * (n / k)^(1 - t.α)

w_from_r_(t::CobbDouglasTechnologyMut; r, τ_k, d, covid) = covid == 1 ? (1 - t.α) * (t.A / t.μ)^(1/(1-t.α)) * (d^(t.α /(1 - t.α))) * (t.α / ((1 + τ_k) * (r + t.δ)))^(t.α / (1 - t.α)) : (1 - t.α) * (t.A / t.μ)^(1/(1-t.α)) * (t.α / ((1 + τ_k) * (r + t.δ)))^(t.α / (1 - t.α))

k_from_r_n_(t::CobbDouglasTechnologyMut; r, n, τ_k, d, covid) = covid == 1 ? ((t.α * t.A * d^t.α ) / (t.μ * (1 + τ_k) * (r + t.δ)))^(1/(1 - t.α)) * n : (t.α * t.A / (t.μ * (1 + τ_k) * (r + t.δ)))^(1/(1 - t.α)) * n 

function agg_T_(g::TaxSys; T, z_dist)
    (; ind_T) = g
    agg_T = zero(eltype(z_dist))
    for i in eachindex(z_dist)
        if ind_T[i]==1        
            agg_T += z_dist[i] * T
        end
    end
    return agg_T
end

function agg_Tu_(g::TaxSys; Tu, z_dist)
    (; ind_Tu) = g
    agg_Tu = zero(eltype(z_dist))
    for i in eachindex(z_dist)
        if ind_Tu[i]==1        
            agg_Tu += z_dist[i] * Tu
        end
    end
    return agg_Tu
end

function taxes_HSV_(h::Household, g::TaxSys; w, τ_n, n_ss, Γ_n, covid, z_dist)
    (; v, z_grid) = h
    n_hsv = zero(eltype(z_grid))
    for (p, z, n_g) in zip(z_dist, z_grid, n_ss)        
        n_hsv += p * (w * z * labor_(v, g; w = w * z, τ_n = τ_n, n_ = n_g, Γ_n = Γ_n, covid = covid))^(1 - g.τ)
    end
    T = w * labor_supply_(h, g; w, τ_n, n_ss, Γ_n, covid, z_dist) - τ_n * n_hsv
    return T
end

## Period 0

w_from_r_k_(t::CobbDouglasTechnologyMut; r_k, τ_k, d, covid) = covid == 1 ? (1 - t.α) * (t.A / t.μ)^(1/(1-t.α)) * (d^(t.α /(1 - t.α))) * (t.α / ((1 + τ_k) * r_k))^(t.α / (1 - t.α)) : (1 - t.α) * (t.A / t.μ)^(1/(1-t.α)) * (t.α / ((1 + τ_k) * r_k))^(t.α / (1 - t.α))

k_from_r_k_n_(t::CobbDouglasTechnologyMut; r_k, n, τ_k, d, covid) = covid == 1 ? ((t.α * t.A * d^t.α ) / (t.μ * (1 + τ_k) * r_k))^(1/(1 - t.α)) * n : (t.α * t.A / (t.μ * (1 + τ_k) * r_k))^(1/(1 - t.α)) * n 

function k_from_rk_hat_n_(t::CobbDouglasTechnologyMut; r_k_hat, n, τ_k, d, covid)
    r_k = r_k_hat + t.δ
    k = covid == 1 ? ((t.α * t.A * d^t.α ) / (t.μ * (1 + τ_k) * r_k))^(1/(1 - t.α)) * n : (t.α * t.A / (t.μ * (1 + τ_k) * r_k))^(1/(1 - t.α)) * n 
    return k
end
