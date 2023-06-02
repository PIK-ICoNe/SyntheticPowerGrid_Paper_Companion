function generate_powergrid_fluctuations(pg::PowerGrid, fluc_node_idxs, fluctuations)
    @assert length(fluc_node_idxs) == length(fluctuations) "Give a time series for each node you want to simulate a fluctuation."
    nodes = deepcopy(pg.nodes) 
    nodes[fluc_node_idxs] .= fluctuations
    pg_fluc = PowerGrid(nodes, pg.lines)

    return pg_fluc
end

function calculate_performance_measures(f, rocof, Δt = 0.01)
    t = collect(10.0:Δt:solution.dqsol.t[end])
    ω = f .* 2π
    f_max = findmax(f)[1]
    f_min = findmin(f)[1]

    RoCoF_max = findmax(RoCoFs)[1]
    RoCoF_min = findmin(RoCoFs)[1]

    println("The maximal frequency is: ", f_max)
    println("The minimal frequency is: ", f_min)

    println("The maximal RoCoF is: ", RoCoF_max)
    println("The minimal RoCoF is: ", RoCoF_min)

    N = length(ω_nodes)
    T = solution.dqsol.t[end]
    
    ω_mean = mean(f, dims = 1) .* 2π

    mean_norm = (1/T) * sum(abs2, ω_mean) * Δt |> sqrt             # √{1/T ∫ (1/N ∑ᵢωᵢ)² dt}
    sync_norm = (1/T) * sum(abs2,(ω .- ω_mean) / N) * Δt |> sqrt # √{1/T ∫ [1/N ∑ᵢ(ωᵢ-1/N ∑ⱼωⱼ)]² dt}

    println("The mean norm is: ", mean_norm)
    println("The synchronization norm is: ", sync_norm)

    return mean_norm, sync_norm, f_max, f_min, RoCoF_max, RoCoF_min
end