function generate_powergrid_fluctuations(pg::PowerGrid, fluc_node_idxs, fluctuations)
    @assert length(fluc_node_idxs) == length(fluctuations) "Give a time series for each node you want to simulate a fluctuation."
    nodes = deepcopy(pg.nodes) 
    nodes[fluc_node_idxs] .= fluctuations
    pg_fluc = PowerGrid(nodes, pg.lines)

    return pg_fluc
end

function calculate_performance_measures(solution::PowerGridSolution, frequency_symbol::Symbol = :x_1)
    ω_indices = findall(n -> frequency_symbol ∈ symbolsof(n), solution.powergrid.nodes)

    N = length(ω_indices)
    T = solution.dqsol.t[end]
    Δt = 0.01
    
    sol_ω = solution(0.0:Δt:T, ω_indices, frequency_symbol)
    ω_mean = mean(sol_ω, dims = 1)

    mean_norm = (1/T) * sum(abs2, ω_mean) * Δt |> sqrt             # √{1/T ∫ (1/N ∑ᵢωᵢ)² dt}
    sync_norm = (1/T) * sum(abs2,(sol_ω .- ω_mean)/N) * Δt |> sqrt # √{1/T ∫ [1/N ∑ᵢ(ωᵢ-1/N ∑ⱼωⱼ)]² dt}

    println("The mean norm is: ", mean_norm)
    println("The synchronization norm is: ", sync_norm)

    return mean_norm, sync_norm
end