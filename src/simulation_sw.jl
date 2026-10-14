# ============================================================
# [CSI] Monte-Carlo household panel and the Consumer-Sentiment "better off" shares
# ============================================================
# PROVENANCE -- this is a port of
#   Fourth Year\Aiyagari RPI\AiyagariOpen\src\simulation_interp1.jl   (2025-09-23)
# and the shares are computed with the estimator that actually produced the
# figure, which lives in MATLAB, NOT in the Julia file:
#   Fourth Year\Aiyagari RPI\AiyagariOpen\src\Matlab\Plots_histograms.m
#
#       u_gain(t,:)  = u_sim(t,:) ./ u_sim(t-4,:)
#       numerator    = sum(sum(u_gain(t:t+3,:) >= u_high, 'omitnan'))
#       denominator  = sum(sum(~isnan(u_gain(t:t+3,:))))
#       better       = round(100 * numerator / denominator)
#
# VERIFIED: with u_high in [1.0000001, 1.0005] and windows 5:8, 9:12, 13:16,
# 17:20, 21:24, 25:28, this reproduces the model column of CSI_answers_3.png
#   -> 30 43 44 29 27 28   EXACTLY, on sim_data_2.mat.
# The live Plots_histograms.m now uses u_high = 1.01, which gives 25 35 35 22 22 22;
# it was widened AFTER the figure numbers were taken.  The matching band survives
# as a commented-out line there: %pos = find(u_gain<=1.0002 & u_gain>=0.9998).
# The four share_better_* series in the Julia file use N or count(!=) as the
# denominator and reproduce NOTHING -- do not use them.
#
# The "same"/"worse" columns of that figure have no source in either language;
# they are computed here from the same band so they are at least defined.
#
# TWO UTILITY MEASURES are returned, because this model rations hours at
# n = s_t * ntilde_z and the original does not:
#   :v  stored value function v[a_idx, z_idx]      <- what the Fourth Year figure used
#   :x  flow GHH composite  c - theta*n^(1+1/nu)/(1+1/nu)  at the SAME (a_idx,z_idx)
# ------------------------------------------------------------
using Distributions, Random

"""
    simulate_panel_sw(e_init, path, h_path; N, T0, T_trans, seed)

Simulate N households through T0 steady-state quarters then T_trans transition
quarters.  Assets evolve by linear interpolation of the asset policy; utility is
read at a SNAPPED grid index (`searchsortedfirst`), matching the original.
Returns `(u_v, u_x, a_sim, z_sim)`, each T x N with T = T0 + T_trans.
"""
function simulate_panel_sw(e_init, path, h_path; N = 1_000_000, T0 = 10, T_trans = 40, seed = 9823)
    Random.seed!(seed)
    h = e_init.h
    a_grid = h.a_grid; Na = length(a_grid); Nz = length(h.z_grid)
    T = T0 + T_trans
    T_trans <= length(path) || error("T_trans = $(T_trans) exceeds the path length $(length(path))")
    v_ = h.v   # GHH preference struct, for the disutility term

    a_sim = zeros(Float64, T, N); z_sim = zeros(Int, T, N)
    u_v   = zeros(Float64, T, N); u_x   = zeros(Float64, T, N)

    pdf0 = vec(e_init.ws.pdf)
    draws = rand(Categorical(pdf0 ./ sum(pdf0)), N)
    @inbounds for n in 1:N
        z_idx, r = divrem(draws[n] - 1, Na)          # column-major: (a, z)
        a_sim[1, n] = a_grid[r + 1]; z_sim[1, n] = z_idx + 1
    end

    for t in 1:(T - 1)
        if t <= T0
            a_pol = e_init.ws.a_pol; c_pol = e_init.ws.c_pol
            n_pol = e_init.ws.n_pol; P = e_init.h.P; v = e_init.ws.v
        else
            j = t - T0
            a_pol = path[j].a_pol; c_pol = path[j].c_pol
            n_pol = path[j].n_pol; P = h_path[j].P; v = path[j].v
        end
        cumP = cumsum(P, dims = 2)     # one rand + a search per draw, not a Categorical per draw
        @inbounds for n in 1:N
            a_now = a_sim[t, n]; z_now = z_sim[t, n]
            a_sim[t+1, n] = interp1D(a_now, a_grid, view(a_pol, :, z_now))
            z_sim[t+1, n] = min(searchsortedfirst(view(cumP, z_now, :), rand()), Nz)
            ai = clamp(searchsortedfirst(a_grid, a_now), 1, Na)
            u_v[t, n] = v[ai, z_now]
            u_x[t, n] = c_pol[ai, z_now] - disutility(v_; n = n_pol[ai, z_now])
        end
    end
    let j = T_trans
        a_pol = path[j].a_pol; c_pol = path[j].c_pol; n_pol = path[j].n_pol; v = path[j].v
        @inbounds for n in 1:N
            ai = clamp(searchsortedfirst(a_grid, a_sim[T, n]), 1, Na); z = z_sim[T, n]
            u_v[T, n] = v[ai, z]
            u_x[T, n] = c_pol[ai, z] - disutility(v_; n = n_pol[ai, z])
        end
    end
    return (u_v = u_v, u_x = u_x, a_sim = a_sim, z_sim = z_sim)
end

"""
    csi_shares(u; wins, u_high)

The `Plots_histograms.m` estimator, transcribed exactly: ratio to four quarters
earlier, a symmetric band collapsed onto 1, NaN dropped from the denominator
(MATLAB `omitnan`).  Returns (better, same, worse), each a 6-vector of percentages,
UNROUNDED: the tail varies by ~0.2pp against 0.025pp of Monte-Carlo noise (5 seeds),
so rounding to integer would discard real signal.  Round at the display layer.
"""
function csi_shares(u; wins = [5, 9, 13, 17, 21, 25], u_high = 1.0002)
    u_low = 2.0 - u_high
    T, N = size(u)
    better = Float64[]; same = Float64[]; worse = Float64[]
    for w in wins
        (w + 3 <= T) || error("window $(w):$(w+3) exceeds T = $(T)")
        nb = 0; ns = 0; nw = 0; den = 0
        @inbounds for t in w:(w+3), n in 1:N
            g = u[t, n] / u[t-4, n]
            isnan(g) && continue                       # MATLAB omitnan
            den += 1
            if !isinf(g) && g <= u_high && g >= u_low; ns += 1      # collapsed onto 1
            elseif g >= u_high; nb += 1
            else; nw += 1 end
        end
        push!(better, 100 * nb / den)
        push!(same,   100 * ns / den)
        push!(worse,  100 * nw / den)
    end
    return (better = better, same = same, worse = worse)
end
