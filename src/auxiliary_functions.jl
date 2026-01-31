# ============================================================
# Auxiliary functions
# ============================================================
# helper function to compute the ergodic distribution from transition matrix P
function ergodic(P)
    for i in axes(P, 1)
        !(sum(@view P[i, :]) ≈ 1.0) && error("Matrix is not valid. Rows don't sum to one")
    end 
    # from the calibration function
    invariant_dist = real(inv(eigvecs(P))[size(P, 1), :])
    invariant_dist .= invariant_dist ./ sum(invariant_dist)
    return invariant_dist
end


# Linear interpolation function
function interp1D(x::Real, xvec, yvec)
    i = searchsortedfirst(xvec, x)
    ind = clamp(i, firstindex(xvec) + 1, lastindex(xvec))
    val = yvec[ind - 1] + (yvec[ind] - yvec[ind - 1]) / (xvec[ind] - xvec[ind - 1]) * (x - xvec[ind - 1])
    return val
end 


chebyshev(a, b) = @tullio (max) x = abs(a[i] - b[i])


logrange(start, stop, shifter, length) = (10^y - shifter for y in range(log10(start + shifter), log10(stop + shifter); length))


function grid(; start = 0.0, stop = 15.0, length = 500, scale = :log, shifter = 1 - start)  
    if scale == :linear 
        grid = collect(range(start, stop, length = length))
    elseif scale == :log 
        grid = collect(logrange(start, stop, shifter, length))
    else 
        @error "Scale neither linear nor log."
    end 
    return grid
end 

function calib_ss(h, t, g, τ_k_ss, Λ_t, Λ_tu, Λ_g, Λ_b, r_target, A)
    t_ = CobbDouglasTechnologyMut(t.α, A^((1 - t.α)), t.δ, t.μ)
    try
        @time e_init = stationary_equilibrium_init(
            h, t_, g, τ_k_ss, Λ_t, Λ_tu, Λ_g, Λ_b; 
            r_range = (0.0, 0.005), root_max = 1, verbose = true
        )
        resid = e_init.r - r_target
        println("resid = ", resid)
        println("A = ", A)
        return (; resid, A)
    catch e
        @warn "stationary_equilibrium_init failed at A = $A with error: $e"
        return (; resid = NaN, A)  # return a value that fails root-finding but doesn't crash
    end
end

function residual_beta!(F, x, h_base, t_base, g, τ_k_ss, Λ_t, Λ_tu, Λ_g, Λ_b, r_target)
    β = x[1]
    println("Trying β = $β")

    try
        v = h_base.v  # use existing θ from base
        u = EZ(ies = ies, ra = crra, β = β)
        h = Household(u = u, a_grid = h_base.a_grid, v = v, P = h_base.P, z_grid = h_base.z_grid)
        δ = (Λ_rk - Λ_f) / Λ_k - r_target
        τ_k_ss = Λ_f / ((r_target + δ) * Λ_k)
        A = 0.035 #0.05 #0.04 #0.03766691763295541
        α = Λ_rk
        μ = 1.0
        t = CobbDouglasTechnologyMut(α, A^((1 - α)), δ, μ)

        e = stationary_equilibrium_init(h, t, g, τ_k_ss, Λ_t, Λ_tu, Λ_g, Λ_b, Λ_bstar;
                                        r_range = (0.0, 0.01), root_max = 1, verbose = true)

        r_val = e.r
        F[1] = r_val - r_target
        println("Success: r = $r_val, residual = $(F[1])")
    catch err
        println("Error encountered: ", err)
        F[1] = 1e3  # penalty on failure
    end
end

# ============================================================
# Figures
# ============================================================

function generate_tables_detailed(e_init_ss, e_final_ss, transition_path; io = stdout)

    table1 = zeros(12, 1)
    
    pdf_marginal = sum(e_init_ss.ws.pdf, dims=2)[:, 1]
    a_grid = round.(e_init_ss.h.a_grid, digits = 10)
    qu_10 =  round.([quantile(a_grid, StatsBase.ProbabilityWeights(pdf_marginal), x) for x in 0.0:0.1:1.0], digits = 10)

    weight = similar(e_init_ss.ws.pdf)

    for (i, path) in enumerate([transition_path]) 
        w_gains = 100 .* (path.v[1] ./ e_init_ss.ws.v .- 1)

        _min = minimum(w_gains[e_init_ss.ws.pdf .> 0])
        _mean = sum(w_gains .* e_init_ss.ws.pdf)
     
        w_gains_rs = reshape(w_gains, 1, :)[1, :]
        wgains_summary = Float64[]

        for (qi, qii) in zip(qu_10, qu_10[2:end])
            if qii == qu_10[end]
                qii_ = qii + 10 * eps()
            else
                qii_ = qii
            end
            for i in axes(weight, 1), j in axes(weight, 2)
                weight[i, j] = (qi <= a_grid[i] < qii_) ? e_init_ss.ws.pdf[i, j] : 0.0
            end
            weight_rs = reshape(weight, 1, :)[1, :]
            push!(wgains_summary, mean(w_gains_rs, ProbabilityWeights(weight_rs)))
        end
    
        _poor = first(wgains_summary)
        _rich = last(wgains_summary)
        _middle = sum(wgains_summary[5:6])/2 
       table1[1:2, i] .= [_mean, _min]
       table1[3:end,i] .= [wgains_summary[j] for j=1:10]
    end 

    table2 = zeros(15, 2)
    for (i, e) in enumerate([e_init_ss, e_final_ss])    
        final_pdf_marginal = sum(e.ws.pdf, dims=2)[:, 1]
        qu = round.([quantile(a_grid, StatsBase.ProbabilityWeights(final_pdf_marginal), x) for x in 0.0:0.10:1.0], digits=10)
    
        total = sum(a_grid .* final_pdf_marginal)
        dist = [sum((a_grid .<= qi) .* (a_grid .* final_pdf_marginal)) * 100 / total for qi in qu] |> diff
        
        table2[:, i] .= [
            [e.b / y(e) * 100,
            e.r * 100,
            e.k / y(e),
            e.τ_n / y(e) * 100, 
            100 * sum(e.ws.pdf[1, :])];   # constrained households
            dist ]            
        
    end 

    out = [ 
        pretty_table(io, table2, header = ["Laissez Faire", "Portfolio"], 
            formatters = (v, i, j) -> i ∉ [2, 3, 4] ? round(Int, v) : round(v, digits = 2), 
            row_labels = ["b/y (%)", "r (%)", "k/y", "T/y (%)", "Shared constrained hh", "P10 Wealth Share", "P20 Wealth Share", "P30 Wealth Share" , "P40 Wealth Share", "P50 Wealth Share","P60 Wealth Share", "P70 Wealth Share", "P80 Wealth Share" , "P90 Wealth Share", "P99 Wealth Share"], title = "Steady State Statistics"
        ),
    
        pretty_table(io, table1, header = ["Debt 60%-80%"], formatters = ft_round(2), row_labels = ["mean", "min", "P10","P20","P30","P40","P50","P60","P70","P80","P90","P99"], title = "Welfare Gains at announcement")
    ]

    (io == String) && return out 
    return nothing
end 

function w_gains_path(e_init_ss, transition_path)
    pdf_marginal = sum(e_init_ss.ws.pdf, dims=2)[:, 1]
    a_grid = round.(e_init_ss.h.a_grid, digits = 10)
    qu_10 =  round.([quantile(a_grid, StatsBase.ProbabilityWeights(pdf_marginal), x) for x in 0.0:0.1:1.0], digits = 10)

    weight = similar(e_init_ss.ws.pdf)
    w_gains_path_1= zeros(length(qu_10)-1, length(transition_path))
    for (i, v) in enumerate(transition_path.v) 
        w_gains = 100 .* (v ./ e_init_ss.ws.v .- 1)

        _min = minimum(w_gains[e_init_ss.ws.pdf .> 0])
        _mean = sum(w_gains .* e_init_ss.ws.pdf)
     
        w_gains_rs = reshape(w_gains, 1, :)[1, :]
        wgains_summary = Float64[]

        for (qi, qii) in zip(qu_10, qu_10[2:end])
            if qii == qu_10[end]
                qii_ = qii + 10 * eps()
            else
                qii_ = qii
            end
            for i in axes(weight, 1), j in axes(weight, 2)
                weight[i, j] = (qi <= a_grid[i] < qii_) ? e_init_ss.ws.pdf[i, j] : 0.0
            end
            weight_rs = reshape(weight, 1, :)[1, :]
            push!(wgains_summary, mean(w_gains_rs, ProbabilityWeights(weight_rs)))
        end        
        w_gains_path_1[:,i] .= wgains_summary
    end 
    return w_gains_path_1
end 

function do_plots_gains(mat, T_fig; legend_pos = :rt, last_plot = false)
    fontsize_theme = Theme(fontsize = 15)
    set_theme!(fontsize_theme)
    x_jump = 10
    fig = Figure(resolution = (800, _φ * 800)) 

    ax = Axis(fig[1, 1], title = L"$\mathrm{1-10th~percentiles}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ylims!(ax,  (1.1*minimum(mat[:,1:T_fig]), 1.1*maximum(mat[:,1:T_fig])))
    ax.xticks = 0:x_jump:T_fig
    barplot!(ax, mat[1,1:T_fig], color=RGB(0, 0, 1), alpha=0.3, fontsize = 26)

    # Unemployment
    ax = Axis(fig[1, 2], title = L"$\mathrm{20-29th~percentiles}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ylims!(ax,  (1.1*minimum(mat[:,1:T_fig]), 1.1*maximum(mat[:,1:T_fig])))
    ax.xticks = 0:x_jump:T_fig
    barplot!(ax, mat[3,1:T_fig], color=RGB(0, 0, 1), alpha=0.3, fontsize = 26)


    ax = Axis(fig[1, 3], title = L"$\mathrm{40-49th~percentiles}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ylims!(ax,  (1.1*minimum(mat[:,1:T_fig]), 1.1*maximum(mat[:,1:T_fig])))
    ax.xticks = 0:x_jump:T_fig
    barplot!(ax, mat[5,1:T_fig], color=RGB(0, 0, 1), alpha=0.3, fontsize = 26)

        # Inflation
    ax = Axis(fig[2, 1], title = L"$\mathrm{70-79th~percentiles}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ylims!(ax,  (1.1*minimum(mat[:,1:T_fig]), 1.1*maximum(mat[:,1:T_fig])))
    ax.xticks = 0:x_jump:T_fig
    barplot!(ax, mat[8,1:T_fig], color=RGB(0, 0, 1), alpha=0.3, fontsize = 26)

    # Unemployment
    ax = Axis(fig[2, 2], title = L"$\mathrm{80-89th~percentiles}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ylims!(ax,  (1.1*minimum(mat[:,1:T_fig]), 1.1*maximum(mat[:,1:T_fig])))
    ax.xticks = 0:x_jump:T_fig
    barplot!(ax, mat[9,1:T_fig], color=RGB(0, 0, 1), alpha=0.3, fontsize = 26)


    ax = Axis(fig[2, 3], title = L"$\mathrm{90-99th~percentiles}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ylims!(ax,  (1.1*minimum(mat[:,1:T_fig]), 1.1*maximum(mat[:,1:T_fig])))
    ax.xticks = 0:x_jump:T_fig
    barplot!(ax, mat[10,1:T_fig], color=RGB(0, 0, 1), alpha=0.3, fontsize = 26)
    ################# 
    return fig
end

function do_plots_gains_long(mat, T_fig; legend_pos = :rt, last_plot = false)
    fontsize_theme = Theme(fontsize = 15)
    set_theme!(fontsize_theme)
    x_jump = 200
    fig = Figure(resolution = (800, _φ * 800)) 

    ax = Axis(fig[1, 1], title = L"$\mathrm{1-10th~percentiles}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ylims!(ax,  (1.1*minimum(mat[:,1:T_fig]), 1.1*maximum(mat[:,1:T_fig])))
    ax.xticks = 0:x_jump:T_fig
    barplot!(ax, mat[1,1:T_fig], color=RGB(0, 0, 1), alpha=0.3, fontsize = 26)

    # Unemployment
    ax = Axis(fig[1, 2], title = L"$\mathrm{20-29th~percentiles}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ylims!(ax,  (1.1*minimum(mat[:,1:T_fig]), 1.1*maximum(mat[:,1:T_fig])))
    ax.xticks = 0:x_jump:T_fig
    barplot!(ax, mat[3,1:T_fig], color=RGB(0, 0, 1), alpha=0.3, fontsize = 26)


    ax = Axis(fig[1, 3], title = L"$\mathrm{40-49th~percentiles}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ylims!(ax,  (1.1*minimum(mat[:,1:T_fig]), 1.1*maximum(mat[:,1:T_fig])))
    ax.xticks = 0:x_jump:T_fig
    barplot!(ax, mat[5,1:T_fig], color=RGB(0, 0, 1), alpha=0.3, fontsize = 26)

        # Inflation
    ax = Axis(fig[2, 1], title = L"$\mathrm{70-79th~percentiles}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ylims!(ax,  (1.1*minimum(mat[:,1:T_fig]), 1.1*maximum(mat[:,1:T_fig])))
    ax.xticks = 0:x_jump:T_fig
    barplot!(ax, mat[8,1:T_fig], color=RGB(0, 0, 1), alpha=0.3, fontsize = 26)

    # Unemployment
    ax = Axis(fig[2, 2], title = L"$\mathrm{80-89th~percentiles}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ylims!(ax,  (1.1*minimum(mat[:,1:T_fig]), 1.1*maximum(mat[:,1:T_fig])))
    ax.xticks = 0:x_jump:T_fig
    barplot!(ax, mat[9,1:T_fig], color=RGB(0, 0, 1), alpha=0.3, fontsize = 26)


    ax = Axis(fig[2, 3], title = L"$\mathrm{90-99th~percentiles}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ylims!(ax,  (1.1*minimum(mat[:,1:T_fig]), 1.1*maximum(mat[:,1:T_fig])))
    ax.xticks = 0:x_jump:T_fig
    barplot!(ax, mat[10,1:T_fig], color=RGB(0, 0, 1), alpha=0.3, fontsize = 26)
    ################# 
    return fig
end

function do_plots_slides(path, laissez_faire; legend_pos = :rt, last_plot = false)
    t = laissez_faire.t
    h = laissez_faire.h
    n = laissez_faire.n
    g = laissez_faire.g

    #ylist  = map(k -> output(t; k , n), path.k)
    #size = (400, 300)
    y0 = y(laissez_faire)

    path = path[1:100]

    #################
    #  SUBPLOT 1: DEBT 

    fontsize_theme = Theme(fontsize = 15)
    set_theme!(fontsize_theme)
    T_fig = 40
    x_jump = 10
    fig = Figure(resolution = (800, _φ * 800)) 


    # Inflation
    ax = Axis(fig[1, 1], title = L"$\mathrm{Inflation}~\pi_t $", ylabel = LaTeXString("%"), titlesize = 20)
    π_bar = 100 .*((1.0 .+ e_init.π_).^4 .-1)
    xlims!(ax,  (0, T_fig))
    #ylims!(ax, (0, 1.1 * maximum(path.k ./ y0)))
    ax.xticks = 0:x_jump:T_fig   
    lines!(ax, 100 .*((1.0 .+ path.π_).^4 .-1), color = :blue)
    scatter!(ax, [1.0], [π_bar], color = :black, markersize = 12, strokecolor = :black, strokewidth= 1, marker = :star5)


    # Unemployment
    ax = Axis(fig[1, 2], title = L"$\mathrm{Unemployment}~u_t $", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    un_bar = 100.0 * e_init.h.Pss[1]
    #ylims!(ax, (0, 1.1 * maximum(path.k ./ y0)))
    ax.xticks = 0:x_jump:T_fig   
    lines!(ax, 100 .* path.un, color = :blue)
    scatter!(ax, [1.0], [un_bar], color = :black, markersize = 12, strokecolor = :black, strokewidth= 1, marker = :star5)

    ax = Axis(fig[1, 3], title = L"$\mathrm{Debt}~q_{t}B_{t+1}/y_{t}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    #ylims!(ax,  (90,130))
    value_debt_bar = 100.0 * e_init.q * e_init.b / (4.0 * y0)
    ax.xticks = 0:x_jump:T_fig
    #value_debt = [q * b for (q, b) in zip([laissez_faire.q ; path.q[1:end-1]], path.b)]
    value_debt = [q * b for (q, b) in zip(path.q, [path.b[2:end]; path.b[end]])]
    lines!(ax, 100 .* (value_debt ./ path.y) ./4 , color = :blue)
    scatter!(ax, [1.0], [value_debt_bar], color = :black, markersize = 12, strokecolor = :black, strokewidth= 1, marker = :star5)

    #################
    #  SUBPLOT 2: DEBT ISSUANCES  

    ax = Axis(fig[2, 1], title =L"$\mathrm{Borrowing}~q_t \ell_t /y_{t}$", ylabel = LaTeXString("%"),  titlesize = 20)
    xlims!(ax,  (0, T_fig))
    #ylims!(ax,  (8,18))
    ax.xticks = 0:x_jump:T_fig
    l_bar = 100.0 * e_init.q *(e_init.b - (1.0 - g.ϑ)/(1.0 + e_init.π_) * e_init.b) / (4.0.* y0)
    lines!(ax, 100 .* (path.q .* path.l ./ (4.0.* path.y)) , color = :blue)
    scatter!(ax, [1.0], [l_bar], color = :black, markersize = 12, strokecolor = :black, strokewidth= 1, marker = :star5)


    #################
    #  SUBPLOT 3: LUMP SUM TRANSFERS

    ax = Axis(fig[2, 2], title = L"$\mathrm{Stimulus~checks}~T_t^s/y_{t}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig
    T_bar = 100.0 * e_init.lump_sum_T / (4.0.* y0)
    lines!(ax, 100 .* path.lump_sum_T ./ (4.0.* path.y), color = :blue, label = L"$T/y_0$")
    scatter!(ax, [1.0], [T_bar], color = :black, markersize = 12, strokecolor = :black, strokewidth= 1, marker = :star5)


    ax = Axis(fig[2, 3], title = L"$\mathrm{Unemp.~insurance}~T_t^u/y_{t}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig
    Tu_bar = 100.0 * e_init.lump_sum_Tu / (4.0.* y0)
    lines!(ax, 100 .* path.lump_sum_Tu ./ (4.0.* path.y), color = :blue, label = L"$T^u/y_0$")
    scatter!(ax, [1.0], [Tu_bar], color = :black, markersize = 12, strokecolor = :black, strokewidth= 1, marker = :star5)
    return fig
end

function do_plots_slides2(path, laissez_faire; legend_pos = :rt, last_plot = false)
    t = laissez_faire.t
    h = laissez_faire.h
    n = laissez_faire.n

    ylist  = map(k -> output(t; k , n), path.k)
    size = (400, 300)
    y0 = y(laissez_faire)

    path = path[1:100]
    #################
    #  B

    fontsize_theme = Theme(fontsize = 15)
    set_theme!(fontsize_theme)
    T_fig = 40
    x_jump = 10
    fig = Figure(resolution = (800, _φ * 800))
    # ######################
    # #  SUBPLOT 1
    #r_K = 100.0.*(path.r_k[1] .- t.δ)
    r_bar = 100.0 .*((1.0 .+ e_init.r).^4 .-1.0)
    r_b = 100 * ((1 + path.r[1])^4 - 1)
    r_k = 100 * ((1 + path.r_k[1])^4 - 1)
    ax = Axis(fig[1, 1], title = L"$\mathrm{Interest~rate}~r_t$", ylabel = LaTeXString("%"), titlesize = 20) 
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig   
    r_series = 100.0 .* ((1.0 .+ path.r).^4 .- 1.0)
    r_series[1] = NaN  # hide first value
    lines!(ax, r_series, color = :blue)
    scatter!(ax, [1.0], [r_bar], color = :black, markersize = 12, strokecolor = :black, strokewidth= 1, marker = :star5)
    label_str = "\$r^b_0 = $(Printf.@sprintf("%.1f", r_b))%\$"
    r_label = LaTeXString(label_str)
    label_str_k = "\$r^k_0 = $(Printf.@sprintf("%.1f", r_k))%\$"
    r_label_k = LaTeXString(label_str_k)
    text!(ax, r_label, position = (20, 0.95), textsize = 15, font = "CMU Serif")
    text!(ax, r_label_k, position = (20, 0.92), textsize = 15, font = "CMU Serif")

    # ######################
    # #  SUBPLOT 2

    ax = Axis(fig[1, 2], title = L"$\mathrm{Capital}~K_t/y_{0}$", titlesize = 20)
    xlims!(ax,  (0, T_fig))
    #ylims!(ax, (0, 1.1 * maximum(path.k ./ y0)))
    ylims!(ax, (0.999 * minimum(path.k ./ (4.0.*y0)), 1.001 * maximum(path.k ./ (4.0.*y0))))
    k_y_bar = e_init.k ./ (4.0 .* y0)
    ax.xticks = 0:x_jump:T_fig   
    lines!(ax, path.k ./ (4.0 .*y0), color = :blue)
    scatter!(ax, [1.0], [k_y_bar], color = :black, markersize = 12, strokecolor = :black, strokewidth= 1, marker = :star5)


    # ######################
    # # SUBPLOT 3 LABOR

    ax = Axis(fig[1, 3], title = L"$\mathrm{Labor}~N_t/N_{0}$", ylabel = LaTeXString("%"), titlesize = 20) 
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig   
    lines!(ax, 100 .*(path.n ./ laissez_faire.n .-1), color = :blue)

    # ######################
    # # SUBPLOT 4 WAGE
    ax = Axis(fig[2, 1], title = L"$\mathrm{Wage}~w_t/w_{0}$", ylabel = LaTeXString("%"), titlesize = 20) 
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig   
    lines!(ax, 100 .*( path.w ./ laissez_faire.w .-1), color = :blue)

    # ######################
    # # SUBPLOT 5 CONSUMPTION    
    ax = Axis(fig[2, 2], title = L"$\mathrm{Consumption}~C_t/C_{0}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    #ylims!(ax, (-5, 5))
    ax.xticks = 0:x_jump:T_fig
    c0 = laissez_faire.c
    lines!(ax, 100 .* (path.C ./ c0  .- 1), color = :blue)

    # ######################
    # # SUBPLOT 5 GDP  
    ax = Axis(fig[2, 3], title = L"$\mathrm{Output}~y_t/y_{0}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig
    lines!(ax, 100 .* (path.y ./ y0  .- 1), color = :blue)    
    return fig
end

function do_plots_slides3(path, laissez_faire; legend_pos = :rt, last_plot = false)
    t = laissez_faire.t
    h = laissez_faire.h
    n = laissez_faire.n


    y0 = y(laissez_faire)
    path = path[1:100]

    #################
    #  SUBPLOT 1: DEBT 

    fontsize_theme = Theme(fontsize = 15)
    set_theme!(fontsize_theme)
    T_fig = 40
    x_jump = 10
    fig = Figure(resolution = (800, _φ * 800)) 


    # Inflation
    lab_tax_bar = 100 * e_init.lab_tax / y0
    ax = Axis(fig[1, 1], title = L"$\mathrm{Net~labor~tax~revenue}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig
    lines!(ax, 100 .* path.lab_tax ./ path.y, color = :blue)
    scatter!(ax, [1.0], [lab_tax_bar], color = :black, markersize = 12, strokecolor = :black, strokewidth= 1, marker = :star5)


    # Unemployment
    ax = Axis(fig[1, 2], title = L"$\mathrm{Capital~tax~revenue}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig
    capital_tax_bar = 100* (e_init.τ_k * (e_init.r .+ t.δ)*e_init.k) / y0
    capital_taxes = (path.τ_k .* (path.r .+ t.δ).*path.k) ./ path.y
    capital_taxes[1] = path.τ_k[1] * (path.r_k[1])*path.k[1] / path.y[1]
    lines!(ax, 100 .* capital_taxes, color = :blue)
    scatter!(ax, [1.0], [capital_tax_bar], color = :black, markersize = 12, strokecolor = :black, strokewidth= 1, marker = :star5)


    ax = Axis(fig[1, 3], title =L"$\mathrm{Bond~price}~q$", ylabel = LaTeXString("%"),  titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig
    q_bar = e_init.q
    lines!(ax, path.q, color = :blue)
    scatter!(ax, [1.0], [q_bar], color = :black, markersize = 12, strokecolor = :black, strokewidth= 1, marker = :star5)
    
    # Labor tax
    ax = Axis(fig[2, 1], title = L"$\mathrm{Labor~tax}~\tau^n_t$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig   
    labor_tax_bar = 100 *(1.0 .- e_init.τ_n)
    lines!(ax, 100 .*(1.0 .- path.τ_n), color = :blue, label = L"$\mathrm{Labor}")
    scatter!(ax, [1.0], [labor_tax_bar], color = :black, markersize = 12, strokecolor = :black, strokewidth= 1, marker = :star5)
    
    # Capital tax
    ax = Axis(fig[2, 2], title = L"$\mathrm{Capital~tax}~\tau^k_t$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    cap_tax_bar =  100 *e_init.τ_k
    ax.xticks = 0:x_jump:T_fig   
    lines!(ax, 100 .*path.τ_k, color = :blue, label = L"$\mathrm{Capital}")
    scatter!(ax, [1.0], [cap_tax_bar], color = :black, markersize = 12, strokecolor = :black, strokewidth= 1, marker = :star5)
    ax = Axis(fig[2, 3], title = L"$\mathrm{Investment}~X_t / X_0$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig
    x0 = t.δ * laissez_faire.k
    #ylims!(ax,  (-0.5,3.1))
    lines!(ax, 100 .* (path.X ./ x0 .-1.0), color = :blue)
    #lines!(ax, 100 .* path.lump_sum_Tu ./ y0, color = :red, linestyle = :dash, label = L"$T^u/y_0$")
    #axislegend(position = legend_pos)

    return fig
end

function do_plots_slides_long(path, laissez_faire; legend_pos = :rt, last_plot = false)
    t = laissez_faire.t
    h = laissez_faire.h
    n = laissez_faire.n

    #ylist  = map(k -> output(t; k , n), path.k)
    #size = (400, 300)
    y0 = y(laissez_faire)
    #################
    #  SUBPLOT 1: DEBT 

    fontsize_theme = Theme(fontsize = 15)
    set_theme!(fontsize_theme)
    T_fig = length(path.r)#150
    x_jump = 200 #length(path.r) / 12
    fig = Figure(resolution = (800, _φ * 800))


    # Inflation
    ax = Axis(fig[1, 1], title = L"$\mathrm{Inflation}~\pi_t $", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    #ylims!(ax, (0, 1.1 * maximum(path.k ./ y0)))
    ax.xticks = 0:x_jump:T_fig   
    lines!(ax, 100 .* path.π_, color = :blue)

    # Unemployment
    ax = Axis(fig[1, 2], title = L"$\mathrm{Unemployment}~u_t $", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    #ylims!(ax, (0, 1.1 * maximum(path.k ./ y0)))
    ax.xticks = 0:x_jump:T_fig   
    lines!(ax, 100 .* path.un, color = :blue)


    ax = Axis(fig[1, 3], title = L"$\mathrm{Debt}~q_{t-1}B_t/y_{0}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    #ylims!(ax,  (90,130))
    ax.xticks = 0:x_jump:T_fig
    value_debt = [q * b for (q, b) in zip([laissez_faire.q ; path.q[1:end-1]], path.b)]
    lines!(ax, 100 .* value_debt ./ y0, color = :blue)


    #################
    #  SUBPLOT 2: DEBT ISSUANCES  

    ax = Axis(fig[2, 1], title =L"$\mathrm{Borrowing}~q_t \ell_t /y_{0}$", ylabel = LaTeXString("%"),  titlesize = 20)
    xlims!(ax,  (0, T_fig))
    #ylims!(ax,  (8,18))
    ax.xticks = 0:x_jump:T_fig
    lines!(ax, 100 .* (path.q .* path.l ./ y0), color = :blue)

    #################
    #  SUBPLOT 3: LUMP SUM TRANSFERS

    ax = Axis(fig[2, 2], title = L"$\mathrm{Stimulus~checks}~T_t^s/y_{0}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig
    #ylims!(ax,  (-0.5,3.1))
    lines!(ax, 100 .* path.lump_sum_T ./ y0, color = :blue, label = L"$T/y_0$")
    #lines!(ax, 100 .* path.lump_sum_Tu ./ y0, color = :red, linestyle = :dash, label = L"$T^u/y_0$")

    ax = Axis(fig[2, 3], title = L"$\mathrm{Unemp.~insurance}~T_t^u/y_{0}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig
    #ylims!(ax,  (-0.5,3.1))
    #lines!(ax, 100 .* path.lump_sum_T ./ y0, color = :blue, label = L"$T/y_0$")
    lines!(ax, 100 .* path.lump_sum_Tu ./ y0, color = :blue, label = L"$T^u/y_0$")
    return fig
end

function do_plots_slides_long2(path, laissez_faire; legend_pos = :rt, last_plot = false)
    t = laissez_faire.t
    h = laissez_faire.h
    n = laissez_faire.n

    ylist  = map(k -> output(t; k , n), path.k)
    size = (400, 300)
    y0 = y(laissez_faire)

    fontsize_theme = Theme(fontsize = 15)
    set_theme!(fontsize_theme)
    T_fig = length(path.r)#150
    x_jump = 200 #length(path.r) / 12
    fig = Figure(resolution = (800, _φ * 800))

    # ######################
    # #  SUBPLOT 1
    r_K = 100.0.*(path.r_k[1] .- t.δ)
    ax = Axis(fig[1, 1], title = L"$\mathrm{Interest~rate}~r_t$", ylabel = LaTeXString("%"), titlesize = 20) 
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig   
    lines!(ax, path.r .* 100, color = :blue)
    scatter!(ax, [1.0], [r_K], color = :black, markersize = 12, strokecolor = :black, strokewidth= 1, marker = :star5)


    # ######################
    # #  SUBPLOT 2

    ax = Axis(fig[1, 2], title = L"$\mathrm{Capital}~K_t/y_{0}$", titlesize = 20)
    xlims!(ax,  (0, T_fig))
    #ylims!(ax, (0, 1.1 * maximum(path.k ./ y0)))
    ylims!(ax, (0.999 * minimum(path.k ./ y0), 1.001 * maximum(path.k ./ y0)))
    ax.xticks = 0:x_jump:T_fig   
    lines!(ax, path.k ./ y0, color = :blue)

    # ######################
    # # SUBPLOT 3 LABOR

    ax = Axis(fig[1, 3], title = L"$\mathrm{Labor}~N_t/N_{0}$", ylabel = LaTeXString("%"), titlesize = 20) 
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig   
    lines!(ax, 100 .*(path.n ./ laissez_faire.n .-1), color = :blue)

    # ######################
    # # SUBPLOT 4 WAGE
    ax = Axis(fig[2, 1], title = L"$\mathrm{Wage}~w_t/w_{0}$", ylabel = LaTeXString("%"), titlesize = 20) 
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig   
    lines!(ax, 100 .*( path.w ./ laissez_faire.w .-1), color = :blue)

    # ######################
    # # SUBPLOT 5 CONSUMPTION    
    ax = Axis(fig[2, 2], title = L"$\mathrm{Consumption}~C_t/C_{0}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    #ylims!(ax, (-5, 5))
    ax.xticks = 0:x_jump:T_fig

    #agg_c = [y + (1 - t.δ) * k - kprime for (y, k, kprime) in zip(path.y, path.k, path.k[2:end])]
    c0 = y(laissez_faire) - t.δ * laissez_faire.k - laissez_faire.G

    #band!(ax, 1:length(agg_c), 100 .* (agg_c ./ c0  .- 1), [0.0 for _ in agg_c], color = (:blue, 0.2))
    lines!(ax, 100 .* (path.C ./ c0  .- 1), color = :blue)

    # ######################
    # # SUBPLOT 5 GDP  
    ax = Axis(fig[2, 3], title = L"$\mathrm{Output}~y_t/y_{0}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig
    #ylims!(ax,  (-0.5,3.1))
    #lines!(ax, 100 .* path.lump_sum_T ./ y0, color = :blue, label = L"$T/y_0$")
    lines!(ax, 100 .* (path.y ./ y0  .- 1), color = :blue) 
    return fig
end

function do_plots_slides_long3(path, laissez_faire; legend_pos = :rt, last_plot = false)
    t = laissez_faire.t
    h = laissez_faire.h
    n = laissez_faire.n
    y0 = y(laissez_faire)

    fontsize_theme = Theme(fontsize = 15)
    set_theme!(fontsize_theme)
    T_fig = length(path.r)#150
    x_jump = 200 #length(path.r) / 12
    fig = Figure(resolution = (800, _φ * 800))

    # Inflation
    ax = Axis(fig[1, 1], title = L"$\mathrm{Net~labor~tax~revenue}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig
    lines!(ax, 100 .* path.lab_tax ./ y0, color = :blue)

    # Unemployment
    ax = Axis(fig[1, 2], title = L"$\mathrm{Capital~tax~revenue}$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig
    capital_taxes = (path.τ_k .* (path.r .+ t.δ).*path.k) / y0
    capital_taxes[1] = path.τ_k[1] * (path.r_k[1])*path.k[1] / y0
    lines!(ax, 100 .* capital_taxes, color = :blue)

    ax = Axis(fig[1, 3], title =L"$\mathrm{Bond~price}~q$", ylabel = LaTeXString("%"),  titlesize = 20)
    xlims!(ax,  (0, T_fig))
    #ylims!(ax,  (8,18))
    ax.xticks = 0:x_jump:T_fig
    lines!(ax, path.q, color = :blue)
    
    # Labor tax
    ax = Axis(fig[2, 1], title = L"$\mathrm{Labor~tax}~\tau^n_t$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig   
    lines!(ax, 100 .*(1.0 .- path.τ_n), color = :red, linestyle = :dash, label = L"$\mathrm{Labor}")
    
    # Capital tax
    ax = Axis(fig[2, 2], title = L"$\mathrm{Capital~tax}~\tau^k_t$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig   
    lines!(ax, 100 .*path.τ_k, color = :blue, label = L"$\mathrm{Capital}")

    ax = Axis(fig[2, 3], title = L"$\mathrm{Investment}~X_t / X_0$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig
    x0 = t.δ * laissez_faire.k
    lines!(ax, 100 .* (path.X ./ x0 .-1.0), color = :blue)
    return fig
end

function do_plot_stdcon(path, laissez_faire; legend_pos = :rt, last_plot = false)
    t = laissez_faire.t
    h = laissez_faire.h
    n = laissez_faire.n

    y0 = y(laissez_faire)
    path = path[1:100]

    #################
    #  SUBPLOT 1: DEBT 

    fontsize_theme = Theme(fontsize = 15)
    set_theme!(fontsize_theme)
    T_fig = 40
    x_jump = 10
    fig = Figure(resolution = (800, _φ * 800)) 

    # #########################
    # # STD DEV
    ax = Axis(fig[1, 1],  title = L"$\sigma(\log c) - \sigma(\log c_0)$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    ax.xticks = 0:x_jump:T_fig  
    std_con = zeros(length(path))
    for i in eachindex(path)
        con = reshape(path[i].c_pol, 1, :)[1, :]
        allpdf = reshape(path[i].pdf, 1, :)[1, :]
        std_con[i] = 100 * StatsBase.std(log.(con), StatsBase.ProbabilityWeights(allpdf))
    end

    con_init = laissez_faire.ws.c_pol #consumption_alloc(h, g; R = 1 + laissez_faire.r, laissez_faire.w, laissez_faire.τ_n, laissez_faire.ws.a_pol)
    con_init_ = reshape(con_init,1, :)[1, :]
    allpdf_init = reshape(laissez_faire.ws.pdf, 1, :)[1, :]
    std_con_init = 100 * StatsBase.std(log.(con_init_), StatsBase.ProbabilityWeights(allpdf_init))

    band!(ax, 1:length(std_con), std_con .- std_con_init, [0 for _ in std_con], color = (:blue, 0.2))
    lines!(ax, std_con .- std_con_init, color = :blue)  
    return fig
end

function do_plot_stdcon_long(path, laissez_faire; legend_pos = :rt, last_plot = false)
    t = laissez_faire.t
    h = laissez_faire.h
    n = laissez_faire.n
    y0 = y(laissez_faire)

    #################
    #  SUBPLOT 1: DEBT 

    fontsize_theme = Theme(fontsize = 15)
    set_theme!(fontsize_theme)
    T_fig = length(path)
    x_jump = 200
    fig = Figure(resolution = (800, _φ * 800)) 

    # #########################
    # # STD DEV
    #
    ax = Axis(fig[1, 1],  title = L"$\sigma(\log c) - \sigma(\log c_0)$", ylabel = LaTeXString("%"), titlesize = 20)
    xlims!(ax,  (0, T_fig))
    #ylims!(ax, (0, 1.1 * maximum(path.k ./ y0)))
    ax.xticks = 0:x_jump:T_fig  
    std_con = zeros(length(path))
    for i in eachindex(path)
        con = reshape(path[i].c_pol, 1, :)[1, :]
        allpdf = reshape(path[i].pdf, 1, :)[1, :]
        std_con[i] = 100 * StatsBase.std(log.(con), StatsBase.ProbabilityWeights(allpdf))
    end

    con_init = laissez_faire.ws.c_pol #consumption_alloc(h, g; R = 1 + laissez_faire.r, laissez_faire.w, laissez_faire.τ_n, laissez_faire.ws.a_pol)
    con_init_ = reshape(con_init,1, :)[1, :]
    allpdf_init = reshape(laissez_faire.ws.pdf, 1, :)[1, :]
    std_con_init = 100 * StatsBase.std(log.(con_init_), StatsBase.ProbabilityWeights(allpdf_init))

    band!(ax, 1:length(std_con), std_con .- std_con_init, [0 for _ in std_con], color = (:blue, 0.2))
    lines!(ax, std_con .- std_con_init, color = :blue)  
    return fig
end