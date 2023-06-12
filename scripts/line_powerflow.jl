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
import SyntheticPowerGrids.get_ancillary_operationpoint
import SyntheticPowerGrids.get_initial_guess
default(grid = false, foreground_color_legend = nothing, bar_edges = false, lw = 3, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 12, tickfontsize = 10)

## LArge grid
file_path = joinpath(@__DIR__, "../data/powergrids/synthetic_power_grid_large.json")
pg = read_powergrid(file_path, Json) 
op = find_operationpoint(pg)
num_nodes = length(pg.nodes)

## Line Loading large grid
p_nodes = map(n -> pg.nodes[n].P, 1:num_nodes)
_, _, P, P_max = validate_power_flow_on_lines(op, :StaticLine)
P = vcat(values(P)...)
P_max = vcat(values(P_max)...)

line_loading = (P ./ P_max) .* 100

ll_mean = mean(line_loading)
ll_max = maximum(line_loading)
ll_large = line_loading

plt = histogram(line_loading, legend = false, xlabel = L"LL[\%]", ylabel = L"p(LL)", normalized = true)
savefig(plt, "plots/line_loading_distribution_large_grid.pdf")

## Peak demand grid
file_path = joinpath(@__DIR__, "../data/powergrids/synthetic_power_grid_peak_demand.json")
pg = read_powergrid(file_path, Json) 
op = find_operationpoint(pg)
num_nodes = length(pg.nodes)

##
_, _, P, P_max = validate_power_flow_on_lines(op, :StaticLine)
P = vcat(values(P)...)
P_max = vcat(values(P_max)...)

line_loading = (P ./ P_max) .* 100

ll_mean = mean(line_loading)
ll_max = maximum(line_loading)
ll_small = line_loading

## Normalized via area
file_path = joinpath(@__DIR__, "../data/powergrids/synthetic_power_grid_area_normalized.json")
pg = read_powergrid(file_path, Json) 
p_nodes = map(n -> pg.nodes[n].P, 1:num_nodes)
op = find_operationpoint(pg)
num_nodes = length(pg.nodes)

##
_, _, P, P_max = validate_power_flow_on_lines(op, :StaticLine)
P = vcat(values(P)...)
P_max = vcat(values(P_max)...)

line_loading = (P ./ P_max) .* 100

ll_mean = mean(line_loading)
ll_max = maximum(line_loading)
ll_small_normalization = line_loading
plt=histogram(ll_small_normalization, xaxis = L"LL[\%]", normalize = :probability, yaxis = L"p(LL)", legend = false)

savefig(plt, "plots/line_loading_area.pdf")

##
stephist(ll_large,  normalize = :probability, label = L"Large")
stephist!(ll_small,  normalize = :probability, label = L"Small")
plt = stephist!(ll_small_normalization, xaxis = L"LL[\%]", normalize = :probability, label = L"Norm")

savefig(plt, "plots/line_loading_comparison.pdf")