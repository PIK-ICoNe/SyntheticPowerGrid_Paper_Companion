function nodal_data(pg::PowerGrid)
    ω_indices = findall(n -> :x_1 ∈ symbolsof(n), pg.nodes)
    nodes = deepcopy(pg.nodes) 
    fluc_node_idxs = findall(typeof.(pg.nodes) .== PQAlgebraic) # Find all Load Buses in the grid
    P_set = map(i -> nodes[i].P, fluc_node_idxs) # Load their power set-points
    Q_set = map(i -> nodes[i].Q, fluc_node_idxs)

    return ω_indices, nodes, fluc_node_idxs, P_set, Q_set
end