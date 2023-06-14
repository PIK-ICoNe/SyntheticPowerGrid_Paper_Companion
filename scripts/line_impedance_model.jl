using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using XLSX
using SyntheticPowerGrids
using DataFrames, Plots, CSV, Statistics
using LaTeXStrings, StatsPlots
using Colors
using EmbeddedGraphs, SyntheticNetworks
using Graphs
using KernelDensity
using StatsBase
import SyntheticPowerGrids.get_line_admittance_matrix
default(grid = false, foreground_color_legend = nothing, bar_edges = false, lw = 3, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 12, tickfontsize = 10)

## Parameters
# Chosen Per Unit System
P_base = 100 * 10^6
V_base = 380 * 10^3
Y_base = P_base / (V_base)^2 # Base Admittance

##
# Loading the SciGRID Dataset
path = string(@__DIR__) * "/../data/sci_grids/"
line_table = XLSX.readtable(joinpath(path, "links_de_power_151109.xlsx"), "links_de_power_151109")
df_lines = DataFrame(line_table) 

##
# Cables and Wires
df_380kV = filter(row -> row.voltage == 380000, df_lines) # Only look at lines in the 380kV level

cables = copy(df_380kV.cables)
cables_miss = findall(typeof.(cables) .== Missing)
deleteat!(cables, cables_miss)

wires = copy(df_380kV.wires)
wires_miss = findall(typeof.(wires) .== Missing)
deleteat!(wires, wires_miss)

##
# Plot Histograms
histogram(cables, legend = false, xlabel = "Number of Cables")
histogram(wires, legend = false, xlabel = "Number of Wires")

length(findall(cables .> 3)) / length(cables) * 100
length(findall(wires .< 4)) / length(wires) * 100

mean(cables)
mean(wires)
##
# R and X coefficients
k_c = (df_380kV.cables / 3)
k_w = (df_380kV.wires / 4)

unique(df_380kV.r_ohmkm .* (k_c .* k_w))
unique(df_380kV.x_ohmkm .* (k_c .* k_w))

##
# Mean Length [km] of German High Voltage Power grid lines using SciGRIDs
lengths = df_lines.length_m / 1000 # Conversion to kilometer
mean_len_km = mean(lengths) 

##
# Finding the extrema of the dataset
min_len = findmin(df_lines.length_m / 1000)[1]
max_len = findmax(df_lines.length_m / 1000)[1]

## Power Grid Generation
dist_nodes_vec = []
N = 100

for i in 1:1000
    g = generate_graph(RandomPowerGrid(N, 1, 1/5, 3/10, 1/3, 1/10, 0.0)) # Generate the Embedded Graph
    
    # Calculate the line lengths in [km]
    L_matrix = get_effective_distances(g, mean_len_km = mean_len_km, shortest_line_km = min_len)

    dest = dst.(edges(g.graph))
    source = src.(edges(g.graph))

    for i in eachindex(source) # Only look at nodes which are connected by a line
        push!(dist_nodes_vec, L_matrix[source[i], dest[i]])
    end
end
max_len = findmax(dist_nodes_vec)[1]

##
# Mean and standard deviation of the line lengths!
mean(dist_nodes_vec)
std(dist_nodes_vec)

mean(lengths)
std(lengths)

## 
# Plotting
c1 = colorant"coral"
c2 = colorant"teal"
num_bins = 100
##
p1 = histogram(lengths, label = "SciGRID", color = c1, lw = 0, xlims = [0.0, max_len], bins = num_bins, normalize = true, linecolor = :match)
p2 = histogram(dist_nodes_vec[1:end], color = c2, label = "Synthetic Networks", bins = num_bins, normalize = true, linecolor = :match)

plt = Plots.plot(p1, p2; layout = (2,1), size = (500, 500), xlabel = L"L [km]", xlims = [0.0, max_len])
savefig(plt, "plots/line_length.pdf")

##
# Shunt Capacitance of the lines
C_shunt_per_km = df_380kV.c_nfkm ./ ((df_380kV.wires ./ 4) .* (df_380kV.cables ./ 3)) # SciGRID Paper eq. (3)
unique!(C_shunt_per_km) # Capacitance per length [nF/km] -> Same as in the dena model 
