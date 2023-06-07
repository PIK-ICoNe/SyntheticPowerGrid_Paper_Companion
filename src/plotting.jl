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