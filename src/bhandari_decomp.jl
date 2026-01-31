# ============================================================
# Decomposing Welfare into Aggregate Efficiency, Redistribution, and Insurance
# ============================================================

using SparseArrays
# linear index on (a,z) grid: i in 1:n_a (assets), s in 1:n_s (shock)
@inline linidx(i, s, n_a) = i + (s-1)*n_a

function build_transition_matrix(h; lower_index, lower_weight)
    P′ = h.Pprime
    n_a, n_s = size(lower_index)
    @assert size(lower_weight) == (n_a, n_s)
    @assert size(P′) == (n_s, n_s)

    # We know each (i,s) contributes to 2 destinations for each s′ → 2*n_s entries
    nnz_est = 2 * n_s * n_a * n_s
    rows = Vector{Int32}(undef, nnz_est)
    cols = Vector{Int32}(undef, nnz_est)
    vals = Vector{Float64}(undef, nnz_est)
    p = 0

    @inbounds for s in 1:n_s
        for i in 1:n_a
            ind    = lower_index[i, s]     # ∈ [1, n_a-1]
            weight = lower_weight[i, s]    # ∈ [0,1]
            col = linidx(i, s, n_a)        # source column

            for s′ in 1:n_s
                mass_fac = P′[s′, s]
                wL = mass_fac * weight
                wU = mass_fac * (1 - weight)

                # destination rows
                rL = linidx(ind,     s′, n_a)
                rU = linidx(ind + 1, s′, n_a)

                # add the two contributions
                p += 1; rows[p] = rL; cols[p] = col; vals[p] = wL
                p += 1; rows[p] = rU; cols[p] = col; vals[p] = wU
            end
        end
    end

    # Trim in case our estimate overshot
    rows = resize!(rows, p); cols = resize!(cols, p); vals = resize!(vals, p)

    T = sparse(rows, cols, vals, n_a*n_s, n_a*n_s)
    return T
end

# Given: pdf (n_a x n_s), pdf_next (n_a x n_s), and the arrays lower_index, lower_weight
h_trans = path_r[2].h #e_init.h
lower_index_trans = path_r[2].lower_index # e_init.ws.lower_index
lower_weight_trans = path_r[2].lower_weight # e_init.ws.lower_weight # 
pdf = path_r[2].pdf #e_init.ws.pdf
T = build_transition_matrix(h_trans; lower_index=lower_index_trans, lower_weight=lower_weight_trans)
vec_pdf     = vec(pdf)              # column-major vectorization
vec_pdfnext = T * vec_pdf
pdf_next_T  = reshape(vec_pdfnext, size(pdf))  # (n_a x n_s)

# Should match forward pass:
pdf_next_loop = similar(pdf)
forward_pdf!(pdf_next_loop; h=h_trans, pdf=pdf, lower_index=lower_index_trans, lower_weight=lower_weight_trans)
# forward_pdf!(path[i].pdf; h = h_path[i], path[i-1].pdf, path[i-1].lower_index, path[i-1].lower_weight)
maximum(abs, pdf_next_T .- pdf_next_loop)
maximum(abs, pdf_next_loop .- path_r[3].pdf)
pdf_next_T .- path_r[2].pdf

################################################
function build_all_transitions(path_r, e_init, h_path)
    T = length(path_r)
    Ts = Vector{SparseMatrixCSC{Float64,Int}}(undef, T)

    # @showprogress "Building transitions..." for t in 1:T
    for t in 1:T
        h_trans = h_path[t] #e_init.h
        lower_index_trans = t == 1 ? e_init.ws.lower_index : path_r[t - 1].lower_index # e_init.ws.lower_index
        lower_weight_trans = t == 1 ? e_init.ws.lower_weight : path_r[t - 1].lower_weight # e_init.ws.lower_weight # 
        Ts[t] = build_transition_matrix(h_trans; lower_index=lower_index_trans, lower_weight=lower_weight_trans)
    end
    return Ts
end

Ts =  build_all_transitions(path_r, e_init, h_path)
n_a, n_s = size(path_r[1].a_pol)
pdf0 = e_init.ws.pdf
pdfs = Vector{typeof(pdf0)}(undef, length(Ts) + 1)

v = vec(pdf0)
for t in 1:length(Ts)
    v = Ts[t] * v
    pdfs[t] = reshape(v, n_a, n_s)
end

## EXPECTED CONSUMPTION 
# Compute g_t = T_t' * vec(c_pol_t) without forming T_t explicitly:
function expected_consumption_all_states(Ts::AbstractVector{<:SparseMatrixCSC},
                                         c_pol_t::AbstractMatrix, t::Int)
    n_a, n_s = size(c_pol_t)
    b = vec(c_pol_t)               # c_t
    @inbounds for τ in t:-1:1      # multiply by T_s[τ]' ... T_s[1]'
        b = transpose(Ts[τ]) * b
    end
    return reshape(b, n_a, n_s)    
end

# Here I compute g_t for several horizons and return a Vector of matrices
function expected_consumption_all_states_many(Ts, c_pols; t_list = 1:length(Ts))
    n_a, n_s = size(c_pols[1])
    out = Vector{Matrix{eltype(c_pols[1])}}(undef, length(t_list))
    for (k, t) in pairs(t_list)
        out[k] = expected_consumption_all_states(Ts, c_pols[t], t)
    end
    return out
end

c_pols = [path_r[t].c_pol for t in 1:length(Ts)]
x_pols_B = [path_r[t].c_pol .- h.v.θ .* path_r[t].n_pol.^(1 + 1/h.v.ν) / (1 + 1/h.v.ν) for t in 1:length(Ts)]
x_pols_A = [e_init.ws.c_pol .- h.v.θ .* e_init.ws.n_pol.^(1 + 1/h.v.ν) / (1 + 1/h.v.ν) for t in 1:length(Ts)]
Ts_B =  build_all_transitions(path_r, e_init, h_path)

path_r_A = [ e_init.ws for _ in 1:length(path_r) ]
h_path_A = [e_init.h for _ in 1:length(path_r)]
Ts_A = build_all_transitions(path_r_A, e_init, h_path_A)

pdfs_A = Vector{typeof(pdf0)}(undef, length(Ts) + 1)
v = vec(pdf0)
for t in 1:length(Ts_A)
    v = Ts_A[t] * v
    pdfs_A[t] = reshape(v, n_a, n_s)
end

t_ = 300
# Expected c at horizon t for ALL starting (a0,z0):
g_t = expected_consumption_all_states(Ts, c_pols[t_], t_)
# g_t[i,s] = E_0[c_t | start at (a_i, z_s)]
v0 = vec(pdf0)                         # pdf over (a0,z0)
E_ct = dot(vec(g_t), v0)
E_ct - path_r[t_].C

function backchain_expected(Ts::AbstractVector{<:SparseMatrixCSC}, x_t::AbstractVector, t::Int)
    v = x_t
    @inbounds for τ in t:-1:1
        v = transpose(Ts[τ]) * v
    end
    return v
end

# Expected consumption at date t for ALL initial states (n_a×n_s matrix)
expected_c_all_states(Ts, c_pol_t::AbstractMatrix, t::Int) = reshape(backchain_expected(Ts, vec(c_pol_t), t), size(c_pol_t))

g_t2 = expected_c_all_states(Ts, c_pols[t_], t_)
v_t2 = vec(path_r[t_].pdf)
dot(vec(g_t2), v_t2) - path_r[t_].C

xx = x_pols[1][1,1]

(1 - h.u.β) * h.u.temporal(xx)
(1 - h.u.β) * log(xx)
ξ = Float64(get_inverse_ies(h.u))
β = h.u.β

# One period’s ϕ_{c,t} and Γ_{c,t}
function phi_and_gamma_t(gA_t::AbstractMatrix, gB_t::AbstractMatrix,
    pdf0::AbstractMatrix, β::Float64, ξ::Float64, t::Int)
    # c^Z_t(a0,z0) = sqrt(gA* gB)
    xZ = sqrt.(gA_t .* gB_t)
    ϕ  = (β^t) .* (1 - β) .* xZ.^(-1/ξ) .* xZ

    # Γ_{c,t} = log( E_{pdf0}[gB_t] ) - log( E_{pdf0}[gA_t] )
    EB = dot(vec(pdf0), vec(gB_t))
    EA = dot(vec(pdf0), vec(gA_t))
    # numerical safety
    ϵ = eps(Float64)
    Γ = log(max(EB, ϵ)) - log(max(EA, ϵ))
    return ϕ, Γ
end

function agg_efficiency(Ts_A, Ts_B, x_pols_A, x_pols_B, pdf0; β, ξ, t_list=1:length(Ts_B))
    n_a, n_s = size(x_pols_B[1])
    @assert size(pdf0) == (n_a, n_s)
    @assert all(size(c) == (n_a, n_s) for c in x_pols_A)
    @assert all(size(c) == (n_a, n_s) for c in x_pols_B)

    results = Vector{NamedTuple}(undef, length(t_list))
    total = 0.0
    for (k, t) in pairs(t_list)
        # E[x_t^J | a0,z0] for J=A,B
        gA_t = expected_c_all_states(Ts_A, x_pols_A[t], t)
        gB_t = expected_c_all_states(Ts_B, x_pols_B[t], t)

        ϕ, Γ = phi_and_gamma_t(gA_t, gB_t, pdf0, β, ξ, t)
        ϕ_= dot(vec(pdf0), vec(ϕ))   # ∫ pdf0 ϕ_{c,t}

        results[k] = (t=t, phi=ϕ, Gamma=Γ, phi_bar=ϕ_)
        total += ϕ_* Γ
    end
    return results, total
end

function redistribution(Ts_A, Ts_B, x_pols_A, x_pols_B, pdf0; β::Real, ξ::Real, t_list=1:length(Ts_B))
    n_a, n_s = size(x_pols_B[1])
    @assert size(pdf0) == (n_a, n_s)

    per_t = Vector{NamedTuple}(undef, length(t_list))
    total = 0.0
    ϵ = 1e-16  # numerical floor for logs/ratios

    for (k, t) in pairs(t_list)
        # Expected consumption for all initial states
        gA_t = expected_c_all_states(Ts_A, x_pols_A[t], t)
        gB_t = expected_c_all_states(Ts_B, x_pols_B[t], t)

        # φ_{x,t}(a0,z0) = β^t * ((1-β)*sqrt(gA*gB))^(1/(1-ξ))
        xZ = sqrt.(gA_t .* gB_t)
        ϕ  = (β^t) .* (1 - β) .* xZ.^(-1/ξ) .* xZ
        # Normalized weights over initial distribution pdf0
        EA = max(ϵ, dot(vec(pdf0), vec(gA_t)))
        EB = max(ϵ, dot(vec(pdf0), vec(gB_t)))
        wA = gA_t ./ EA
        wB = gB_t ./ EB

        # Δ = log wB - log wA  (safe logs)
        Δ  = log.(max.(wB, ϵ)) .- log.(max.(wA, ϵ))

        # Period contribution: ∫ pdf0 * ϕ * Δ
        term = dot(vec(pdf0), vec(ϕ .* Δ))
        total += term

        per_t[k] = (t=t, phi=ϕ, wA=wA, wB=wB, Delta=Δ, term=term)
    end
    return per_t, total
end

function insurance(Ts_A, Ts_B, x_pols_A, x_pols_B, pdf0; β::Real, ξ::Real, t_list=1:length(Ts_B))
    n_a, n_s = size(x_pols_B[1])
    @assert size(pdf0) == (n_a, n_s)
    ϵ = 1e-16  # floor for logs
    per_t = Vector{NamedTuple}(undef, length(t_list))
    total = 0.0

    for (k, t) in pairs(t_list)
        # --- Expectations of x_t for φ_{x,t} ---
        gA_t = expected_c_all_states(Ts_A, x_pols_A[t], t)
        gB_t = expected_c_all_states(Ts_B, x_pols_B[t], t)
        xZ = sqrt.(gA_t .* gB_t)
        ϕ  = (β^t) .* (1 - β) .* xZ.^(-1/ξ) .* xZ


        # --- Expectations of log x_t and (log x_t)^2 for variances ---
        xA_t = x_pols_A[t]; xB_t = x_pols_B[t]
        logxA = log.(max.(xA_t, ϵ))
        logxB = log.(max.(xB_t, ϵ))
        logxA2 = logxA .^ 2
        logxB2 = logxB .^ 2

        ElogA  = expected_c_all_states(Ts_A, logxA,  t)
        ElogB  = expected_c_all_states(Ts_B, logxB,  t)
        Elog2A = expected_c_all_states(Ts_A, logxA2, t)
        Elog2B = expected_c_all_states(Ts_B, logxB2, t)

        varA = Elog2A .- ElogA.^2
        varB = Elog2B .- ElogB.^2

        # Λ_{x,t} = -1/2 * (varB - varA)
        Λ = -0.5 .* (varB .- varA)

        # Period contribution: ∫ pdf0 * φ * (1/ξ) * Λ
        term = (1/ξ) * dot(vec(pdf0), vec(ϕ .* Λ))
        total += term

        per_t[k] = (t=t, phi=ϕ, varlogA=varA, varlogB=varB, Lambda=Λ, term=term)
    end

    return per_t, total
end

welfare_gains = sum((log.(path_r[1].v) .- log.( e_init.ws.v)).* e_init.ws.pdf)
t_list = 1:length(Ts_B)
results_eff, agg_eff = agg_efficiency(Ts_A, Ts_B, x_pols_A, x_pols_B, pdf0; β=β, ξ=ξ, t_list=t_list)
results_red, redist = redistribution(Ts_A, Ts_B, x_pols_A, x_pols_B, pdf0; β=β, ξ=ξ, t_list=t_list)
results_ins, insur = insurance(Ts_A, Ts_B, x_pols_A, x_pols_B, pdf0; β=β, ξ=ξ, t_list=t_list)

agg_eff / welfare_gains
redist / welfare_gains
insur / welfare_gains

(agg_eff + redist + insur)/ welfare_gains









ξ = Float64(get_inverse_ies(h.u))
β = h.u.β
##### SECOND STEP 
function cdf_vec(pdf::AbstractMatrix{<:Real})
    p = vec(pdf)                       # column-major (assets fastest)
    tot = sum(p)
    @assert tot > 0 "Distribution has zero total mass."
    p = p ./ tot                       # ensure it sums to 1
    return cumsum(p)
end

# Build the partition I_sorted from two distributions (A: no policy, B: policy at time t)
pdfA = e_init.ws.pdf 
pdfB = path_r[1].pdf
function quantile_partition(pdfA::AbstractMatrix, pdfB::AbstractMatrix; add_ends=false)
    @assert size(pdfA) == size(pdfB)
    FA = cdf_vec(pdfA)
    FB = cdf_vec(pdfB)

    I = vcat(FA, FB)
    if add_ends
        I = vcat(0.0, I, 1.0)
    end
    out = sort!(I)
    return out                 # this is Isorted
end

"""
Build the selector matrix T_j for policy j:
- Isorted :: Vector{Float64} partition points (keep duplicates)
- Fj      :: Vector{Float64} CDF of policy j in lexicographic order (length N = n_a*n_s)
Returns:
- Tj :: SparseMatrixCSC{Float64} of size (length(Isorted)) × N
       with exactly one 1 per row at k = searchsortedfirst(Fj, ι)
"""
function build_Tj(Isorted::AbstractVector{<:Real}, Fj::AbstractVector{<:Real})
    M = length(Isorted)
    N = length(Fj)
    # For each ι, find k = min { m : Fj[m] ≥ ι }
    ks = searchsortedfirst.(Ref(Fj), Isorted)  # vector of indices in 1..N

    # One-hot rows: row r selects column ks[r]
    rows = collect(1:M)
    cols = ks
    vals = ones(Float64, M)

    return sparse(rows, cols, vals, M, N)  # exactly one nonzero per row
end

"""
Compute vector m̃ of interval lengths from Isorted.
Normalization Isorted[1] must be 0, last must be 1.
"""
function interval_masses(Isorted::AbstractVector{<:Real})
    M = length(Isorted)
    #@assert first(Isorted) ≈ 0.0 "Isorted must start at 0"
    #@assert last(Isorted) ≈ 1.0 "Isorted must end at 1"
    mtilde = similar(Isorted, M)
    mtilde[1] = Isorted[1]        # =0.0
    @inbounds for i in 2:M
        mtilde[i] = Isorted[i] - Isorted[i-1]
    end
    return mtilde
end

# ---------- minimal helpers for this file ----------
# Backward-apply: g0 = H_1 … H_t x_t, with H_τ = (T_τ)' (T_τ maps t-1→t)
function backchain_expected(Ts::AbstractVector{<:SparseMatrixCSC}, x_t_vec::AbstractVector, t::Int)
    v = x_t_vec
    @inbounds for τ in t:-1:1
        v = transpose(Ts[τ]) * v
    end
    return v
end

# Build selector matrices T^A, T^B given Isorted and the two pdfs
function selectors_for_quantiles(Isorted, pdfA, pdfB)
    FA = cdf_vec(pdfA)
    FB = cdf_vec(pdfB)
    TA = build_Tj(Isorted, FA)
    TB = build_Tj(Isorted, FB)
    return TA, TB
end

# φ_x,t(ι) = β^t * ((1-β) * sqrt(ExA * ExB))^(1/(1-ξ))
@inline function phi_xt_quantile(ExA::AbstractVector, ExB::AbstractVector, pareto::AbstractVector; β::Real, ξ::Real, t::Integer)
    xZ = sqrt.(max.(ExA, 0.0) .* max.(ExB, 0.0))
    return (β^t) .* (1 - β) .* pareto .* xZ.^(-1/ξ) .* xZ
end


function welfare_decompose(Ts_A, Ts_B, x_pols_A, x_pols_B, pdfA, pdfB, pareto;
                             β::Real, ξ::Real, t_list = 1:length(Ts_B))
    n_a, n_s = size(x_pols_B[1])
    N = n_a*n_s
    @assert size(pdfA) == (n_a, n_s)
    @assert all(size(c) == (n_a, n_s) for c in x_pols_A)
    @assert all(size(c) == (n_a, n_s) for c in x_pols_B)
    #@assert all(size(p) == (n_a, n_s) for p in pdfB_list)
    @assert all(size(T) == (N,N) for T in Ts_A)
    @assert all(size(T) == (N,N) for T in Ts_B)

    # Build quantile partition ONCE (includes endpoints; keeps duplicates)
    Isorted = quantile_partition(pdfA, pdfB; add_ends=true)
    mtilde  = interval_masses(Isorted)
    # Build selectors once
    TA, TB = selectors_for_quantiles(Isorted, pdfA, pdfB)
    pareto_Isorted = TA * pareto

    #per_t = Vector{NamedTuple}(undef, length(t_list))
    agg_eff_total = 0.0
    redist_total  = 0.0
    ins_total     = 0.0
    ϵ = 1e-16
    ξf = float(ξ)

    for (_, t) in pairs(t_list)
        # 1) Quantile partition for this t (keep duplicates; include endpoints)
        #Isorted = quantile_partition(pdfA, pdfB_list[t]; add_ends=false)
        #mtilde  = interval_masses(Isorted)

        # 2) Build selectors T^A, T^B
        #TA, TB = selectors_for_quantiles(Isorted, pdfA, pdfB_list[t])

        # 3) Expectations on the quantile grid
        #    E[x_t^j | ι] = T^j * (H_1^j … H_t^j x_t^j)
        gA0 = backchain_expected(Ts_A, vec(x_pols_A[t]), t)    # size N
        gB0 = backchain_expected(Ts_B, vec(x_pols_B[t]), t)    # size N
        ExA = TA * gA0                                         # length |Isorted|
        ExB = TB * gB0

        # 4) φ and Γ
        phi = phi_xt_quantile(ExA, ExB, pareto_Isorted; β=β, ξ=ξ, t=t)         # vector |Isorted|
        denomA = max(ϵ, sum(mtilde .* ExA))
        denomB = max(ϵ, sum(mtilde .* ExB))
        Gamma = log(denomB) - log(denomA)                      # scalar

        # 5) Redistribution Δ = log wB - log wA, w^j = Ex^j / sum mtilde*Ex^j
        wA = ExA ./ denomA
        wB = ExB ./ denomB
        Delta = log.(max.(wB, ϵ)) .- log.(max.(wA, ϵ))         # |Isorted|

        # 6) Insurance Λ via conditional variances on quantiles
        logA  = log.(max.(x_pols_A[t], ϵ))
        logB  = log.(max.(x_pols_B[t], ϵ))
        logA2 = logA .^ 2
        logB2 = logB .^ 2

        ElogA0  = backchain_expected(Ts_A, vec(logA),  t)
        ElogB0  = backchain_expected(Ts_B, vec(logB),  t)
        Elog2A0 = backchain_expected(Ts_A, vec(logA2), t)
        Elog2B0 = backchain_expected(Ts_B, vec(logB2), t)

        ElogA   = TA * ElogA0
        ElogB   = TB * ElogB0
        Elog2A  = TA * Elog2A0
        Elog2B  = TB * Elog2B0

        varA = Elog2A .- ElogA.^2
        varB = Elog2B .- ElogB.^2
        Lambda = -0.5 .* (varB .- varA)                         # |Isorted|

        # 7) Period contributions (integrals over quantiles via mtilde)
        agg_eff_term = sum(mtilde .* (phi .* Gamma))
        redist_term  = sum(mtilde .* (phi .* Delta))
        ins_term     = (1/ξf) * sum(mtilde .* (phi .* Lambda))

        agg_eff_total += agg_eff_term
        redist_total  += redist_term
        ins_total     += ins_term
    end

    totals = (aggregate_efficiency = agg_eff_total,
              redistribution       = redist_total,
              insurance            = ins_total)
    #return per_t, totals
    return totals
end

path_A = [ e_init.ws for _ in 1:length(path_r) ]
h_path_A = [e_init.h for _ in 1:length(path_r)]
Ts_A = build_all_transitions(path_A, e_init, h_path_A)
Ts_B = build_all_transitions(path_r, e_init, h_path)
x_pols_A = [e_init.ws.c_pol .- h.v.θ .* e_init.ws.n_pol.^(1 + 1/h.v.ν) / (1 + 1/h.v.ν) for t in 1:length(Ts_A)]
x_pols_B = [path_r[t].c_pol .- h.v.θ .* path_r[t].n_pol.^(1 + 1/h.v.ν) / (1 + 1/h.v.ν) for t in 1:length(Ts_B)]
pdf0 = deepcopy(e_init.ws.pdf)
welfare_gains = sum((log.(path_r[1].v) .- log.( e_init.ws.v)).* e_init.ws.pdf)
welfare_gains_original = sum((path_r[1].v ./ e_init.ws.v .-1).* e_init.ws.pdf)

totals_A = welfare_decompose(Ts_A, Ts_B, x_pols_A, x_pols_B, pdf0, pdf0;
                                        β=β, ξ=ξ, t_list=1:100) #length(Ts_B)

totals_A.aggregate_efficiency
totals_A.redistribution
totals_A.insurance

## Percentages
totals_A.aggregate_efficiency / welfare_gains
totals_A.redistribution / welfare_gains
totals_A.insurance / welfare_gains

## Levels
totals_A.aggregate_efficiency / welfare_gains
totals_A.redistribution / welfare_gains
totals_A.insurance / welfare_gains


welfare_try = sum((log.(path_try[1].v) .- log.(e_init.ws.v)).* e_init.ws.pdf)
(totals_try.agg_efficiency + totals_try.redistribution + totals_try.insurance)/welfare_try
totals_A.aggregate_efficiency - totals_try.agg_efficiency

# PARETO WEIGHTS 
function pareto_weights_assets_posmass(α, a_grid::AbstractVector, pdf0::AbstractMatrix; tol::Real=0.0)
    n_a, n_s = size(pdf0)
    @assert length(a_grid) == n_a "a_grid length must match pdf0 rows (n_a)."

    # asset-only weights, repeated across shocks (assets-fast / column-major layout)
    w_a = exp.(float(α) .* float.(a_grid))              # n_a
    W   = repeat(reshape(w_a, n_a, 1), 1, n_s)          # n_a × n_s

    # mask positive-mass states
    M = pdf0 .> tol                                     # Bool n_a × n_s
    W_masked = W .* M                                   # zero out zero-mass states

    s = sum(W_masked)
    @assert s > 0 "No positive-mass states under the given tol=$(tol)."

    W_norm = W_masked ./ s
    return vec(W_norm)                                   # vectorized (assets fastest)
end

α     = -1.0
a_grid = e_init.h.a_grid
pdf0   = e_init.ws.pdf
w_a = pareto_weights_assets_posmass(α, a_grid, pdf0; tol=0.0)

totals_B = welfare_decompose(Ts_A, Ts_B, x_pols_A, x_pols_B, pdf0, pdf0, w_a;
                                        β=β, ξ=ξ, t_list=1:100) 

totals_B.aggregate_efficiency / welfare_gains
totals_B.redistribution / welfare_gains
totals_B.insurance / welfare_gains

using LinearAlgebra, SparseArrays

function compute_expectation(M::Int64, T::Int64, x::Matrix{Float64},
                             H::Vector{SparseMatrixCSC{Float64, Int64}})::Matrix{Float64}
    @assert size(x,1) == M && size(x,2) == T+1
    @assert length(H) == T || length(H) == T+1  # some codes store T+1; we only use 1..T
    E = Array{Float64}(undef, M, T+1)

    # two scratch buffers reused
    buf1 = similar(@view x[:,1])
    buf2 = similar(buf1)

    # t = 1 (horizon 0): no transitions, just copy the column
    @views copyto!(E[:,1], x[:,1])

    # t = 2..T+1 (horizons 1..T): apply H_{t-1} ... H_1
    for t in 2:T+1
        @views copyto!(buf1, x[:,t])         # start from the column (no copy of x if viewed)
        for τ in 1:(t-1)                      # exactly (t-1) multiplies
            mul!(buf2, H[τ], buf1)           # buf2 = H[τ] * buf1
            buf1, buf2 = buf2, buf1          # swap buffers (no alloc)
        end
        @views copyto!(E[:,t], buf1)
    end
    return E
end

function welfare_decomposition_deterministic(c_A::Matrix{Float64}, c_B::Matrix{Float64},
                                             H_A_ss::SparseMatrixCSC{Float64,Int64},
                                             H_B::Vector{SparseMatrixCSC{Float64,Int64}},
                                             μ_0::Vector{Float64},
                                             μ_A_ss::Vector{Float64}, μ_B_ss::Vector{Float64},
                                             T::Int, M::Int,
                                             β::Float64, ξ::Float64,pareto::Vector{Float64})

    # Compute the distribution of individuals in each t.
    μ_B = Array{Float64}(undef, M, T + 1)
    μ_B[:, 1] = μ_0
    for t = 2:T
        μ_B[:, t] = transpose(H_B[t-1]) * μ_B[:, t-1]
    end
    μ_B[:, T+1] = μ_B_ss

    agg_C_A = [dot(c_A[:, t], μ_A_ss) for t = 1:T+1]
    agg_C_B = [dot(c_B[:, t], μ_B[:, t]) for t = 1:T+1]
    agg_C_Z = sqrt.(agg_C_A .* agg_C_B)
    Γ_C = log.(agg_C_B) - log.(agg_C_A)


    expec_C_A = compute_expectation(M, T, c_A, H_A_ss)
    w_C_A = expec_C_A ./ transpose(agg_C_A)
    expec_C_B = compute_expectation(M, T, c_B, H_B)
    w_C_B = expec_C_B ./ transpose(agg_C_B)
    w_C_Z = sqrt.(w_C_A .* w_C_B)

    c_Z = transpose(agg_C_Z) .* w_C_Z

    Δ_C = log.(w_C_B) - log.(w_C_A)

    var_lnC_B = compute_expectation(M, T, (log.(c_B)) .^ 2, H_B) -
                compute_expectation(M, T, log.(c_B), H_B) .^ 2
    var_lnC_A = compute_expectation(M, T, (log.(c_A)) .^ 2, H_A_ss) -
                compute_expectation(M, T, log.(c_A), H_A_ss) .^ 2

    Λ_C = -(1 / 2) * (var_lnC_B - var_lnC_A)

    ϕ_C = transpose([β^t for t = 0:T]) .* pareto .* (1 - β) .* pareto .* c_Z.^(-1/ξ) .* c_Z #((c_Z) .^ (1 - σ))

    γ_C = 1/ξ #σ

    efficiency = sum(ϕ_C .* transpose(Γ_C), dims = 2)
    redistribution = sum(ϕ_C .* Δ_C, dims = 2)
    insurance = sum(γ_C .* ϕ_C .* Λ_C, dims = 2)

    return efficiency, redistribution, insurance
end

function welfare_decomposition_deterministic_safe(c_A::Matrix{Float64}, c_B::Matrix{Float64},
                                             H_A_ss::SparseMatrixCSC{Float64,Int64},
                                             H_B::Vector{SparseMatrixCSC{Float64,Int64}},
                                             μ_0::Vector{Float64},
                                             μ_A_ss::Vector{Float64}, μ_B_ss::Vector{Float64},
                                             T::Int, M::Int,
                                             β::Float64, ξ::Float64,pareto::Vector{Float64})

    # Compute the distribution of individuals in each t.
    μ_B = Array{Float64}(undef, M, T + 1)
    μ_B[:, 1] = μ_0
    for t = 2:T
        μ_B[:, t] = transpose(H_B[t-1]) * μ_B[:, t-1]
    end
    μ_B[:, T+1] = μ_B_ss

    ϵ = 1e-16  # or a bit larger if your scales are large

    # 1) Aggregates for Γ_C
    agg_C_A = [dot(c_A[:, t], μ_A_ss) for t = 1:T+1]
    agg_C_B = [dot(c_B[:, t], μ_B[:, t]) for t = 1:T+1]
    Γ_C     = log.(max.(agg_C_B, ϵ)) .- log.(max.(agg_C_A, ϵ))

    # 2) Expectations → weights (avoid 0 denominators)
    expec_C_A = compute_expectation(M, T, c_A, H_A_ss)
    expec_C_B = compute_expectation(M, T, c_B, H_B)

    denA = max.(transpose(agg_C_A), ϵ)   # 1×(T+1)
    denB = max.(transpose(agg_C_B), ϵ)
    w_C_A = expec_C_A ./ denA
    w_C_B = expec_C_B ./ denB

    # 3) Logs of weights (avoid log(0))
    Δ_C = log.(max.(w_C_B, ϵ)) .- log.(max.(w_C_A, ϵ))

    # 4) c_Z and φ (avoid 0^(-1/ξ))
    w_C_Z = sqrt.(w_C_A .* w_C_B)
    agg_C_Z = sqrt.(max.(agg_C_A, ϵ) .* max.(agg_C_B, ϵ))
    c_Z = transpose(agg_C_Z) .* w_C_Z                     # M×(T+1)
    cZ_safe = max.(c_Z, ϵ)
    ϕ_C = transpose([β^t for t = 0:T]) .* pareto .* (1 - β) .* pareto .* cZ_safe.^(-1/ξ) .* c_Z
    # (If you only want pareto once, drop one of the ". * pareto")

    # 5) Variance terms (avoid log(0))
    log_cB  = log.(max.(c_B, ϵ));   log2_cB = log_cB.^2
    log_cA  = log.(max.(c_A, ϵ));   log2_cA = log_cA.^2

    ElogB   = compute_expectation(M, T, log_cB,  H_B)
    ElogA   = compute_expectation(M, T, log_cA,  H_A_ss)
    Elog2B  = compute_expectation(M, T, log2_cB, H_B)
    Elog2A  = compute_expectation(M, T, log2_cA, H_A_ss)

    var_lnC_B = Elog2B .- ElogB.^2
    var_lnC_A = Elog2A .- ElogA.^2
    Λ_C = .-0.5 .* (var_lnC_B .- var_lnC_A)

    γ_C = 1/ξ
    efficiency    = sum(ϕ_C .* transpose(Γ_C), dims=2)
    redistribution= sum(ϕ_C .* Δ_C,                dims=2)
    insurance     = sum(γ_C .* ϕ_C .* Λ_C,         dims=2)
    return efficiency, redistribution, insurance
end

function stack_policies_matrix(pols::Vector{<:AbstractMatrix}, T::Int, c_ss::AbstractMatrix)
    n_a, n_s = size(pols[1])
    @assert size(c_ss) == (n_a, n_s)
    M = n_a * n_s
    @assert T <= length(pols) "T exceeds available policies length"

    C = Array{Float64}(undef, M, T+1)
    for t in 1:T
        C[:, t] = vec(pols[t])          # assets-fast order (matches vec/pdf)
    end
    C[:, T+1] = vec(c_ss)               # steady-state (or terminal) column
    return C
end

# If policy B also needs the final steady-state column but you don't have a separate
# c_B_ss matrix, you can reuse the last available policy matrix as an approximation:
function stack_policies_matrix_last(pols::Vector{<:AbstractMatrix}, T::Int)
    @assert T <= length(pols)
    return stack_policies_matrix(pols, T, pols[end])
end


# Dimensions
n_a, n_s = size(x_pols_A[1])
M = n_a * n_s
T = 200

# Steady-state consumption for A (no-policy, stationary)
c_A_ss = e_init.ws.c_pol

# If you have a policy steady state for B:
c_B_ss = path_r[end].c_pol

# Build c_A and c_B as M×(T+1) matrices
c_A = stack_policies_matrix(x_pols_A, T, c_A_ss)
c_B = stack_policies_matrix(x_pols_B, T, c_B_ss)  # or stack_policies_matrix_last(x_pols_B, T)

# Build H_A_ss and H_B in the types your function expects
H_A_ss = Ts_A[1]                 # T_ss is your stationary forward (N×N)
H_B     = [(Ts_B[t]) for t in 1:T]   # Ts_B[t] forward → backward CSC

# Initial distribution and steady states as vectors length M
μ_0    = vec(e_init.ws.pdf)         # initial (t=0) distribution
μ_A_ss = vec(e_init.ws.pdf)         # no-policy steady-state distribution
μ_B_ss = vec(path_r[end].pdf)       # policy steady-state (or your chosen terminal)

# Pareto weights vector over (a0,z0), length M (already normalized as you need)
pareto = ones(M)            # e.g., pareto_weights_assets_posmass(...)

# Now call deterministic decomposition:
efficiency, redistribution, insurance =
    welfare_decomposition_deterministic(c_A, c_B, H_A_ss, H_B, μ_0, μ_A_ss, μ_B_ss, T, M, β, ξ, pareto)

dot((efficiency + redistribution + insurance).*pareto, μ_0)
