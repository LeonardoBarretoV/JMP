# ============================================================
# Calibration Methods: Income Process, Unemployment, and Discretization
# ============================================================

# Income Process

function calibration(nP::Int64, nIID::Int64, ρ::Float64, σP::Float64, σIID::Float64; w0=1.0)

    # persistent component
    pP, y_vals = rouwenhorst(nP, ρ, σP, 0.0)
    y_vals = exp.(y_vals)

    # iid component
    if nIID > 1
        pIID, x_vals = rouwenhorst(nIID, 0.0, σIID, 0.0)
        x_vals = exp.(x_vals)
    else
        pIID = 1.0
        x_vals = 1.0
    end

    nZ = nP * nIID
    z_vals = Array{Float64,1}(undef, nZ)
    P = zeros(nZ, nZ)

    item = 1
    for i = 1:nP
        for j = 1:nIID
            z_vals[item] = y_vals[i] * x_vals[j]
            item += 1
        end
    end

    item = 1
    for i = 1:nP
        for j = 1:nIID
            item2 = 1
            for k = 1:nP
                for l = 1:nIID
                    P[item,item2] = pP[i,k] * pIID[j,l]
                    item2 += 1
                end
            end
            item += 1
        end
    end

    for p in eachrow(P)
        p .= p ./ sum(p)
    end

    # compute mean of ergodic
    invZ = real(inv(eigvecs(P))[nZ,:])
    invZ .= invZ ./ sum(invZ)
    meanZ = invZ' * z_vals
    # scale to get desired mean
    z_vals .= (w0 / meanZ) .* z_vals

    return P, z_vals
end

## Unemployment
function calibration_unemploy(nP::Int64, nIID::Int64, ρ::Float64, σP::Float64, σIID::Float64, Pee::Float64, Puu::Float64, z_u::Float64; w0=1.0)
    n_z = nP*nIID
    P_, z_vals_ = calibration(nP, nIID , ρ, σP, σIID)
    Pprime = copy(P_')
    Pprime_u = [Puu 1-Pee; 1-Puu Pee]

    BigP_prime = zeros(n_z+1, n_z+1)
    BigP_prime[1,1] = Pprime_u[1,1]
    BigP_prime[2:end,1] .= Pprime_u[2,1]/n_z
    BigP_prime[1,2:end] .= Pprime_u[1,2]
    BigP_prime[2:end, 2:end] = Pprime_u[2,2]*Pprime
    BigP = copy(BigP_prime')
    z_vals = [z_u; z_vals_]
    return BigP, z_vals
end

## Kaplan, Moll and Violante (2018) income process 
function calibration_kmv(ROOT::AbstractString)
    file_path = joinpath(ROOT, "kmv_process", "ymarkov_combined.txt")
    Pi_e = exp(readdlm(file_path))
    Pi_e ./= sum(Pi_e, dims=2)

    pi_e = ergodic(Pi_e)
    e_grid_short = exp.(readdlm(joinpath(ROOT, "kmv_process", "ygrid_combined.txt")))
    e_grid_short ./= sum(e_grid_short .* pi_e)
    return Pi_e, e_grid_short
end


function calibration_kmv_short(ROOT::AbstractString)
    Pi_e, e_grid_short = calibration_kmv(ROOT)
    pi_e = ergodic(Pi_e)
    #indices = 1:2:length(e_grid_short)  # Subsample indices
    indices = collect(2:2:length(e_grid_short))  # Subsample indices
    insert!(indices, 1, 1)
    e_grid_new = e_grid_short[indices]
        
    # Allocate new stationary distribution and transition matrix
    n_new = length(e_grid_new)
    Pi_new = zeros(n_new, n_new)

    for i in 1:n_new
        for j in 1:n_new
            # Identify original grid points contributing to the new grid points
            orig_i = indices[i]:min(indices[i]+1, length(e_grid_short))
            orig_j = indices[j]:min(indices[j]+1, length(e_grid_short))

            # Aggregate transitions with stationary distribution weights
            Pi_new[i, j] = sum(Pi_e[orig_i, orig_j] .* pi_e[orig_i]) / sum(pi_e[orig_i])
        end
    end

    # Ensure rows sum to 1
    Pi_new ./= sum(Pi_new, dims=2)
    return Pi_new, e_grid_new
end

function calibration_unemploy_kmv(ROOT::AbstractString, Pee::Float64, Puu::Float64, z_u::Float64; short = 0)
    P_, z_vals_ = short == 0 ? calibration_kmv(ROOT) : calibration_kmv_short(ROOT) 
    n_z = length(z_vals_)
    Pprime = copy(P_')
    Pprime_u = [Puu 1-Pee; 1-Puu Pee]

    BigP_prime = zeros(n_z+1, n_z+1)
    BigP_prime[1,1] = Pprime_u[1,1]
    #BigP_prime[2:end,1] .= Pprime_u[2,1]/n_z ##
    BigP_prime[2:end,1] = Pprime_u[2,1].*P_[1,:] ##
    BigP_prime[1,2:end] .= Pprime_u[1,2]
    BigP_prime[2:end, 2:end] = Pprime_u[2,2]*Pprime
    BigP = copy(BigP_prime')
    z_vals = [z_u; z_vals_]
    return BigP, z_vals
end    

function calibration_unemploy_kmv_(ROOT::AbstractString, Pee::Float64, Puu::Float64, z_u::Float64; short = 0)
    P_, z_vals_ = short == 0 ? calibration_kmv(ROOT) : calibration_kmv_short(ROOT)
    n_z = length(z_vals_)
    Pprime = copy(P_')
    Pprime_u = [Puu 1-Pee; 1-Puu Pee]
    unemp_rate = ergodic(Pprime_u')[1]
    #rescale = clamp(rescale, 0.0,1.0)
    if short == 1
        P_eu = zeros(n_z,1); P_eu[1:5] .= 0.09*0+0.05; P_eu[6:9] .= 0.08*0+0.05; P_eu[10:13] .= 0.04*0+0.01; P_eu[14:n_z] .= 0.02*0
    elseif short == 0
        P_eu = zeros(n_z,1); P_eu[1:10] .= 0.09*0 + 0.05; P_eu[11:18] .= 0.08*0 + 0.05; P_eu[19:26] .= 0.04; P_eu[26:n_z] .= 0.02
    end

    function compute_residual(rescale)
        BigP_prime = zeros(n_z+1, n_z+1)
        BigP_prime[1,1] = Pprime_u[1,1]
        #BigP_prime[2:end,1] .= Pprime_u[2,1]/n_z ##
        BigP_prime[2:end,1] = Pprime_u[2,1].*P_[1,:] ##
        BigP_prime[1,2:end] .= P_eu .* rescale
        for i=2:n_z+1
            BigP_prime[2:end, i] = (1.0 - P_eu[i-1] * rescale) * Pprime[:,i-1]
        end
        BigP = copy(BigP_prime')
        unemp_rate_2 = ergodic(BigP)[1]
        return unemp_rate - unemp_rate_2
    end
    rescale = fzero(compute_residual, 0.5)

    # **Recompute final values using optimized rescale**
    BigP_prime = zeros(n_z+1, n_z+1)
    BigP_prime[1,1] = Pprime_u[1,1]
    BigP_prime[2:end,1] = Pprime_u[2,1] .* P_[1,:]  
    BigP_prime[1,2:end] .= P_eu .* rescale
    for i = 2:n_z+1
        BigP_prime[2:end, i] = (1.0 - P_eu[i-1] * rescale) * Pprime[:,i-1]
    end
    BigP = copy(BigP_prime')

    z_vals = [z_u; z_vals_]
    return BigP, z_vals
end    

## McKay (2023) income process

function calibration_mckay(ROOT::AbstractString)
    z_grid_mc  = [0.256974677325607, 0.3069730037657587, 0.3666992639962444, 0.43804617528515916, 0.5232747117919538, 0.6250857545387845, 0.7467056819720169,
                  0.8919886134674335, 1.065538545862266, 1.2728552535045612, 1.5205085753731022, 2.9388587915700324, 24 ]# 16 24.846646635930057]
    file_path = joinpath(ROOT, "mckay_process", "TransitionMatrix.csv")
    df = CSV.read(file_path, DataFrame; header = false)  # ".." moves one level up
    P_mc = Matrix(df)
    return P_mc, z_grid_mc
end

function calibration_unemploy_mckay(ROOT::AbstractString, Pee::Float64, Puu::Float64, z_u::Float64)
    P_, z_vals_ = calibration_mckay(ROOT)
    n_z = length(z_vals_)
    Pprime = copy(P_')
    Pprime_u = [Puu 1-Pee; 1-Puu Pee]

    BigP_prime = zeros(n_z+1, n_z+1)
    BigP_prime[1,1] = Pprime_u[1,1]
    #BigP_prime[2:end,1] .= Pprime_u[2,1]/n_z ##
    BigP_prime[2:end,1] = Pprime_u[2,1].*P_[1,:] ##
    BigP_prime[1,2:end] .= Pprime_u[1,2]
    BigP_prime[2:end, 2:end] = Pprime_u[2,2]*Pprime
    BigP = copy(BigP_prime')
    z_vals = [z_u; z_vals_]
    return BigP, z_vals
end

function calibration_unemploy_mckay_(ROOT::AbstractString, Pee::Float64, Puu::Float64, z_u::Float64)
    P_, z_vals_ = calibration_mckay(ROOT)
    n_z = length(z_vals_)
    Pprime = copy(P_')
    Pprime_u = [Puu 1-Pee; 1-Puu Pee]
    unemp_rate = ergodic(Pprime_u')[1]
    P_eu = zeros(n_z,1); P_eu[1] = 0.17; P_eu[2] = 0.07; P_eu[3] = 0.045; P_eu[4:end-1] .= 0.025; P_eu[end] = 0.007;

    function compute_residual(rescale)
        BigP_prime = zeros(n_z+1, n_z+1)
        BigP_prime[1,1] = Pprime_u[1,1]
        BigP_prime[2:end,1] = Pprime_u[2,1].*P_[1,:] ## CHANGE THIS 
        BigP_prime[1,2:end] .= P_eu .* rescale
        for i=2:n_z+1
            BigP_prime[2:end, i] = (1.0 - P_eu[i-1] * rescale) * Pprime[:,i-1]
        end
        BigP = copy(BigP_prime')
        unemp_rate_2 = ergodic(BigP)[1]
        return unemp_rate - unemp_rate_2
    end
    # Find the optimal rescale that makes resid close to zero
    rescale = fzero(compute_residual, 0.5)

    # Recompute final values using optimized rescale
    BigP_prime = zeros(n_z+1, n_z+1)
    BigP_prime[1,1] = Pprime_u[1,1]
    BigP_prime[2:end,1] = Pprime_u[2,1] .* P_[1,:]
    BigP_prime[1,2:end] .= P_eu .* rescale
    for i = 2:n_z+1
        BigP_prime[2:end, i] = (1.0 - P_eu[i-1] * rescale) * Pprime[:,i-1]
    end
    BigP = copy(BigP_prime')
    z_vals = [z_u; z_vals_]
    return BigP, z_vals
end

function calibration_unemploy_mckay_init(ROOT::AbstractString, Pee::Float64, Puu::Float64, z_u::Float64; unemp_rate_target, Big_init_state)
    P_, z_vals_ = calibration_mckay(ROOT)
    n_z = length(z_vals_)
    Pprime = copy(P_')
    Pprime_u = [Puu 1-Pee; 1-Puu Pee]
    P_eu = zeros(n_z,1); P_eu[1] = 0.17; P_eu[2] = 0.07; P_eu[3] = 0.045; P_eu[4:end-1] .= 0.025; P_eu[end] = 0.007;

    function compute_residual(rescale; Big_init_state)
        BigP_prime = zeros(n_z+1, n_z+1)
        BigP_prime[1,1] = Pprime_u[1,1]
        BigP_prime[2:end,1] = Pprime_u[2,1].*P_[1,:]
        BigP_prime[1,2:end] .= P_eu .* rescale
        for i=2:n_z+1
            BigP_prime[2:end, i] = (1.0 - P_eu[i-1] * rescale) * Pprime[:,i-1]
        end
        BigP = copy(BigP_prime')
        updated_state = Big_init_state' * BigP
        unemp_rate_2 = updated_state[1]
        return unemp_rate_target - unemp_rate_2
    end
    # Find the optimal rescale that makes resid close to zero
    f = (rescale) -> compute_residual(rescale; Big_init_state)
    rescale_sol = fzero(f, 0.5)

    # **Recompute final values using optimized rescale**
    BigP_prime = zeros(n_z+1, n_z+1)
    BigP_prime[1,1] = Pprime_u[1,1]
    #BigP_prime[2:end,1] .= Pprime_u[2,1]/n_z # CHANGE THIS
    BigP_prime[2:end,1] = Pprime_u[2,1] .* P_[1,:]  
    BigP_prime[1,2:end] .= P_eu .* rescale_sol
    for i = 2:n_z+1
        BigP_prime[2:end, i] = (1.0 - P_eu[i-1] * rescale_sol) * Pprime[:,i-1]
    end
    BigP = copy(BigP_prime')
    z_vals = [z_u; z_vals_]
    return BigP, z_vals #, init_state_updated
end

# Discretization

function rouwenhorst(N::Integer, ρ::Real, σ::Real, μ::Real=0.0)
    σ_y = σ / sqrt(1 - ρ^2)
    p  = (1 + ρ) / 2
    Θ = [p 1 - p; 1 - p p]
    ψ = sqrt(N - 1) * σ_y
    m = μ / (1 - ρ)
    state_values, p = _rouwenhorst(p, p, m, ψ, N)

    return p, state_values
end

function _rouwenhorst(p::Real, q::Real, m::Real, Δ::Real, n::Integer)
    if n == 2
        return [m - Δ, m + Δ],  [p 1 - p; 1 - q q]
    else
        _, θ_nm1 = _rouwenhorst(p, q, m, Δ, n - 1)
        θN = p * [θ_nm1 zeros(n - 1, 1); zeros(1, n)] +
             (1 - p) * [zeros(n - 1, 1) θ_nm1; zeros(1, n)] +
             q * [zeros(1, n); zeros(n - 1, 1) θ_nm1] +
             (1 - q) * [zeros(1, n); θ_nm1 zeros(n - 1, 1)]
        θN[2:end - 1, :] ./= 2
        return range(m - Δ, stop=m + Δ, length=n), θN
    end
end

