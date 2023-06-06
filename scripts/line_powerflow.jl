using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using SyntheticPowerGrids
using Statistics
using PowerDynamics
using LaTeXStrings
using Plots
using EmbeddedGraphs
using Graphs
using DelimitedFiles
import SyntheticPowerGrids.validate_power_flow_on_lines
import SyntheticPowerGrids.get_effective_distances
default(grid = false, foreground_color_legend = nothing, bar_edges = false, lw = 3, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 12, tickfontsize = 10)

mean_len_km = 37.12856121212121
area_ge = 357111 # [km²]
a_ge = sqrt(area_ge)

##
file_path = joinpath(@__DIR__, "../data/powergrids/synthetic_power_grid_large.json")
pg = read_powergrid(file_path, Json) 
op = find_operationpoint(pg)
num_nodes = length(pg.nodes)
v_positions = readdlm(joinpath(@__DIR__, "../data/powergrids/power_grid_large_vertex_positions.txt"))

v_pos = Vector{Vector{Float64}}(undef, num_nodes)

for i in 1:num_nodes
    v_pos[i] = v_positions[i, :]
end

eg = EmbeddedGraph(pg.graph, v_pos)


##
L = get_effective_distances(eg; mean_len_km = mean_len_km, shortest_line_km = 0.06)

dist_nodes = EmbeddedGraphs.weights(eg)

dist_nodes_connected = vcat(dist_nodes...)
unconnected_idx = findall(iszero, dist_nodes_connected) # Unconnected nodes have a length of d = 0.0
deleteat!(dist_nodes_connected, unconnected_idx) 

d_mean = mean(dist_nodes_connected)

c_l = mean_len_km / d_mean

area = (2 * c_l)^2

maximum(dist_nodes_connected)

## Line Loading large grid
p_nodes = map(n -> pg.nodes[n].P, 1:num_nodes)
histogram(p_nodes, bins = 20)

mean(abs.(p_nodes)) # should be around 3.2

gen = findall(p_nodes .> 0.0)
con = findall(p_nodes .< 0.0)

sum(p_nodes[gen])
sum(p_nodes[con])

##
_, _, P, P_max = validate_power_flow_on_lines(op, :StaticLine)
P = vcat(values(P)...)
P_max = vcat(values(P_max)...)

line_loading = (P ./ P_max) .* 100

ll_mean = mean(line_loading)
ll_max = maximum(line_loading)

plt = histogram(line_loading, legend = false, xlabel = L"LL[%]", ylabel = L"p(LL)", normalized = true)
savefig(plt, "plots/line_loading_distribution_large_grid.pdf")

##
file_path = joinpath(@__DIR__, "../data/powergrids/synthetic_power_grid_peak_demand.json")
pg = read_powergrid(file_path, Json) 
op = find_operationpoint(pg)
num_nodes = length(pg.nodes)
v_positions = readdlm(joinpath(@__DIR__, "../data/powergrids/power_grid_peak_demand_vertex_positions.txt"))

v_pos = Vector{Vector{Float64}}(undef, num_nodes)

for i in 1:num_nodes
    v_pos[i] = v_positions[i, :]
end

eg = EmbeddedGraph(pg.graph, v_pos)

##
L = get_effective_distances(eg; mean_len_km = mean_len_km, shortest_line_km = 0.06)

dist_nodes = EmbeddedGraphs.weights(eg)

dist_nodes_connected = vcat(dist_nodes...)
unconnected_idx = findall(iszero, dist_nodes_connected) # Unconnected nodes have a length of d = 0.0
deleteat!(dist_nodes_connected, unconnected_idx) 

d_mean = mean(dist_nodes_connected)

c_l = mean_len_km / d_mean

area = (2 * c_l)^2
maximum(dist_nodes_connected)

##
p_nodes = map(n -> pg.nodes[n].P, 1:num_nodes)
histogram(p_nodes, bins = 20)

mean(abs.(p_nodes)) # should be around 3.2

gen = findall(p_nodes .> 0.0)
con = findall(p_nodes .< 0.0)

sum(p_nodes[gen])
sum(p_nodes[con])


##
_, _, P, P_max = validate_power_flow_on_lines(op, pg_struct)
P = vcat(values(P)...)
P_max = vcat(values(P_max)...)

line_loading = (P ./ P_max) .* 100

ll_mean = mean(line_loading)
ll_max = maximum(line_loading)

plt = histogram(line_loading, legend = false, xlabel = L"LL[%]", ylabel = L"p(LL)", normalized = true)
savefig(plt, "plots/line_loading_distribution_small_grid.pdf")