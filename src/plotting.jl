function plot_fluc_results(solution::PowerGridSolution, fluc_node_idxs, ω_nodes; windowsize, t = solution.dqsol.t)
    num_nodes = length(solution.powergrid.nodes)
    s = solution.dqsol
    f_idx = findall(s -> occursin(Regex("x_1_"), String(s)), rhs(solution.powergrid).syms)

    plt_active_power = Plots.plot()
    for idx in fluc_node_idxs
        plt_active_power = Plots.plot!(t, transpose(solution(t, [idx], :p)), legend = false, ylabel = L"P[p.u.]", xlabel = L"t[s]")
    end

    plt_frequency = Plots.plot()
    for idx in ω_nodes
        plt_frequency = Plots.plot!(t, solution(t, idx, :x_1)./(2π), legend = false, ylabel = L"\Delta f[Hz]", xlabel = L"t[s]")
    end
    
    plt_rocof = Plots.plot()
    for i in eachindex(f_idx)
        RoCoF = s(t, Val{1}, idxs = f_idx[i]).u./(2π)
        RoCoF_run_mean = runmean(RoCoF, windowsize)
        plt_rocof = Plots.plot!(t, RoCoF_run_mean, legend = false, ylabel = L"RoCoF [Hz/s]", xlabel = L"t[s]")
    end

    return plt_active_power, plt_frequency, plt_rocof
end

"""
    function my_graph_plot(pg::PowerGrid, label_nodes = [])

Using GraphMakie to plot the power grid topology.
"""
function my_graph_plot(pg::PowerGrid, label_nodes = [])
    num_nodes = length(pg.nodes)
    load_idxs = findall(typeof.(pg.nodes) .== PQAlgebraic)

    node_color = fill(colorant"coral", num_nodes)
    node_color[load_idxs] .= colorant"teal"

    node_marker = fill(:circle, num_nodes)
    node_marker[load_idxs] .= :utriangle

    if label_nodes != []
        node_label = fill("", num_nodes)
        if length(label_nodes) == 1
            node_label[label_nodes] = string(label_nodes)
        else    
            map(n -> node_label[n] = string(n), label_nodes)
        end
        f, ax, p = graphplot(pg.graph, node_marker = node_marker, nlabels = node_label, node_color = node_color, node_size = 20)
        hidedecorations!(ax); hidespines!(ax)
    else
        f, ax, p = graphplot(pg.graph, node_marker = node_marker, node_color = node_color, node_size = 20)
        hidedecorations!(ax); hidespines!(ax)
    end
    return f
end