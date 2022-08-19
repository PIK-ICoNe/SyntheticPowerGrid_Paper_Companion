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
import SyntheticPowerGrids.get_line_admittance_matrix
default(grid = false, foreground_color_legend = nothing, bar_edges = false, lw = 3, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 12, tickfontsize = 10)

## Parameters
# Chosen Per Unit System
P_base = 400 * 10^6
V_base = 380 * 10^3

# Base Admittance
Y_base = P_base / (V_base)^2

##
# Loading the SciGrids Dataset
path = string(@__DIR__) * "/../data/sci_grids/"
line_table = XLSX.readtable(joinpath(path, "links_de_power_151109.xlsx"), "links_de_power_151109")

df_lines = DataFrame(line_table) 
df_380kV = filter(row -> row.voltage == 380000, df_lines) # Only look at lines in the 380kV level
lengths = df_380kV.length_m / 1000 # Conversion to kilometer

##
# Mean Length [km] of German High Voltage Power grid lines using Sci Grids
mean_len_km = get_mean_line_length()

##
# Finding most common line length
dens_scigrid = kde(lengths)
dens_max_scigrids = findmax(dens_scigrid.density)
peak_scigrids = dens_scigrid.x[dens_max_scigrids[2]]

##
# Finding the extrema of the dataset
min_len = findmin(df_380kV.length_m / 1000)[1]

## Power Grid Generation
Y_abs_vec = []
dist_nodes_vec = []

for i in 1:1000
    g = generate_graph(RandomPowerGrid(100, 1, 1/5, 3/10, 1/3, 1/10, 0.0)) # Generate the Embedded Graph
    
    # Calculate the line lengths in [km]
    L_matrix = get_geographic_distances(g, mean_len_km = mean_len_km, shortest_line_km = min_len)

    dest = dst.(edges(g.graph))
    source = src.(edges(g.graph))

    # Calculate the Admittance in [1 / Ω]
    Y, _ = get_line_admittance_matrix(g, L_matrix)
    Y_abs = vcat(abs.(Y)...)

    for i in eachindex(source) # Only look at nodes which are connected by a line
        push!(dist_nodes_vec, L_matrix[source[i], dest[i]])
        push!(Y_abs_vec, abs.(Y[source[i], dest[i]]))
    end
end
max_len = findmax(dist_nodes_vec)[1]

##
# Most common line length in the synthetic grids
dens_synthetic = kde(dist_nodes_vec / 1.0)
dens_max_synthetic = findmax(dens_synthetic.density)
peak_synthetic = dens_synthetic.x[dens_max_synthetic[2]]

## 
# Plotting
c1 = colorant"coral"
c2 = colorant"teal"

vline([peak_scigrids], color = c1, label = "", alpha = 0.6, linestyle = :dash)
vline!([peak_synthetic], color = c2, label = "", alpha = 0.6, linestyle = :dash)
StatsPlots.density!(dist_nodes_vec, color = c2, label = "Synthetic Networks", xlabel = L"L [km]", ylabel = L"p(L)")
plt = StatsPlots.density!(lengths, color = c1, label = "SciGrids Data", xlims = [min_len, max_len + 10], ylims = [0.0, dens_max_synthetic[1] + 0.001 ])
savefig(plt, "plots/line_length.pdf")

##
# Shunt Capacitance of the lines
C_shunt_per_km = df_380kV.c_nfkm ./ ((df_380kV.wires ./ 4) .* (df_380kV.cables ./ 3)) # SciGrid Paper eq. (3)
unique!(C_shunt_per_km) # Capacitance per length [nF/km] -> Same as in the dena model 

##
# Coupling constant / Admittance Magnitude
StatsPlots.density(Y_abs_vec ./ Y_base, xlabel = L"|Y_{ml}| \ [p.u.]", ylabel = L"p(|Y_{ml}|)", label = "Synthetic Grid", lw = 3)
plt = vline!([6.0], label = "Theoretical Physics Community", linestyle = :dash)
savefig(plt, "plots/coupling.pdf")

StatsPlots.density(Y_abs_vec, xlabel = L"|Y_{ml}| \ [1/Ω]", ylabel = L"p(|Y_{ml}|)", lw = 3, label = "Synthetic Networks")

# K ≈ 6 in basin stability lit.?, |Y_ml| = |K_ml|
# Coupling constant unit? Should be [p.u.]