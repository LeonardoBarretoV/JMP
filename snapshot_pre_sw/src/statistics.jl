# ============================================================
# Statistics and Tables
# ============================================================
# Displays a statistical summary of the economy.

function generate_tables(e_init_ss, e_final_ss, transition_path; io = stdout)
    table1 = zeros(5, 1)
    
    pdf_marginal = sum(e_init_ss.ws.pdf, dims=2)[:, 1]
    a_grid = round.(e_init_ss.h.a_grid, digits = 10)
    qu_10 = round.([quantile(a_grid, StatsBase.ProbabilityWeights(pdf_marginal), x) for x in 0.0:0.1:1.0], digits = 10)

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

       table1[:, i] .= [_mean, _min, _poor, _middle, _rich]
    end 

    table2 = zeros(10, 2)
    for (i, e) in enumerate([e_init_ss, e_final_ss])    
        final_pdf_marginal = sum(e.ws.pdf, dims=2)[:, 1]
        qu = round.([quantile(a_grid, StatsBase.ProbabilityWeights(final_pdf_marginal), x) for x in 0.0:0.20:1.0], digits = 10)
    
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
            row_labels = ["b/y (%)", "r (%)", "k/y", "T/y (%)", "Shared constrained hh", "Q1 Wealth Share", "Q2 Wealth Share", "Q3 Wealth Share" , "Q4 Wealth Share", "Q5 Wealth Share"], title = "Steady State Statistics"
        ),
    
        pretty_table(io, table1, header = ["Debt 60%-80%"], formatters = ft_round(2), row_labels = ["mean", "min", "poor (≤10)", "middle (40-60)", "rich (>90)"], title = "Welfare Gains at announcement")
    ]

    (io == String) && return out 
    return nothing
end 

function generate_tables_detailed_period(e_init_ss, e_final_ss, transition_path, period; io = stdout)

    table1 = zeros(12, 1)
    
    pdf_marginal = sum(e_init_ss.ws.pdf, dims=2)[:, 1]
    a_grid = round.(e_init_ss.h.a_grid, digits = 10)
    qu_10 =  round.([quantile(a_grid, StatsBase.ProbabilityWeights(pdf_marginal), x) for x in 0.0:0.1:1.0], digits = 10)

    weight = similar(e_init_ss.ws.pdf)

    for (i, path) in enumerate([transition_path]) 
        w_gains = 100 .* (path.v[period] ./ e_init_ss.ws.v .- 1)

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
        #gains_p10 = zeros(10, 1)
        #for i=1:10
        #    gains_p10[i] = wgains_summary[i]
        #end

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