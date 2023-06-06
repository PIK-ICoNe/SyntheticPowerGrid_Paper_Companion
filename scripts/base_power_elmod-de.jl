using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using XLSX
using DataFrames, Plots, Statistics
using LaTeXStrings, StatsPlots
using Colors
using Statistics
default(grid = false, foreground_color_legend = nothing, bar_edges = false, lw = 3, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 12, tickfontsize = 10)

##
# Analyzing the ELMOD-De Dataset
path = string(@__DIR__) * "/../data/elmod-de_v1/"
elmod_de = joinpath(path, "Data_Input.xlsx")

##
# Demand table
demand_table = XLSX.readtable(elmod_de, "H_demand", first_row = 1, column_labels = ["time_point", "Demand"])
df_demand = DataFrame(demand_table) 

P_tot = findmin(df_demand[:, :Demand])[1] # Taher et. al. Paper uses the off-peak scenario!
P_tot_on_peak = findmax(df_demand[:, :Demand])[1] # Taher et. al. Paper uses the off-peak scenario!

t_off_peak = findmin(df_demand[:, :Demand])[2]
off_peak_idx = String(df_demand.time_point[t_off_peak])

##
# Node Data
nodes_table = XLSX.readtable(elmod_de, "Node_data", first_row = 1, column_labels = ["Node_ID", "Country", "State", "Dena_zone", "Z6", "Load_share_off_peak", "Load_share_peak", "Run_of_river", "Photovoltaic", "Wind_onshore", "Wind_offshore", "Biomass", "Geothermal", "Coordinates", "Not_used"] )
df_nodes = DataFrame(nodes_table)
delete!(df_nodes, 1) # First row is parts of units of the data -> XLSX tables wont let me skip it properly, so I delete it here

num_nodes = eachindex(df_nodes.Node_ID)
node_names = String.(df_nodes.Node_ID)

P_i = df_nodes.Load_share_off_peak * P_tot # Nodal power demand!

##
# Renewable Energy data -> the off peak scenario is not included! Also ≈ 50% of the time series for the Renewable energy are missing! We should not include it in this study
pv_table = XLSX.readtable(elmod_de, "H_pv", first_row = 23, column_labels = ["time_point", "21", "22", "23", "24", "25", "26", "41", "42", "71", "72", "73", "74", "75", "76", "81", "82", "83", "84", "91", "92"])
df_pv = DataFrame(pv_table)

df_pv.time_point = String.(df_pv.time_point)

filter(row -> row.time_point == off_peak_idx, df_pv) 

## 
# Generating a new data frame to save all important variables
my_df = DataFrame(Node_ID = node_names, nodal_demand = P_i)

##
# Load the line data
lines_table = XLSX.readtable(elmod_de, "Grid_technical")
df_lines = DataFrame(lines_table)

# Add ID for the lines for make them searchable
line_id = collect(1:size(df_lines)[1])
df_lines[!, :line_id] = line_id 

mean(df_lines[:, :length]) # Sanity Check -> very close to the value in the SciGrid dataset! :-)
std(df_lines[:, :length]) 
findmin(df_lines[:, :length])[1]

# Find only lines in the 380kV level
df_380kV = filter(row -> row.voltage == 380, df_lines) 

lines_380kV = df_380kV.line_id

##
# Loading the grid topology
topology_table = XLSX.readtable(elmod_de, "Grid_topology")
df_topology = DataFrame(topology_table)

df_topology[:, :line_id] = line_id # Add same Line ID here as well

## 
# Find the connected nodes!
dfs = stack(df_topology, Not(:line_id))
dfs = filter(:value => val -> ! ismissing(val), dfs)
sort!(dfs, :line_id)

## 
# Save nodes with are connected by a 380kV line -> they are in the correct voltage level
nodes_380kV = []

for line in eachindex(lines_380kV)
    df = filter(row -> row.line_id == lines_380kV[line], dfs)
    append!(nodes_380kV, df.variable)
end

nodes_380kV = sort(String.(unique(nodes_380kV))) # Nodes will at least appear twice, String for filtering

##
# Availability for the different technologies during the off-peak scenario
ava_table = XLSX.readtable(elmod_de, "H_con")
df_ava = DataFrame(ava_table)

##
# Loading the plant data
plant_table = XLSX.readtable(elmod_de, "Plant_con", first_row = 3, column_labels = ["Plant_ID", "Country", "Node_ID", "State", "Dena_zone", "ZoDE_NE_6", "Technology", "Fuel", "Capacity", "Efficiency", "Emission", "Transport", "Not_used_1", "Name"])
df_plants = DataFrame(plant_table) 

capacity_nodes = zeros(num_nodes) # Installed nodal Capacity
availability_nodes = zeros(num_nodes) # Availability of the plant

for n in num_nodes # Run over all nodes
    df = filter(row -> row.Node_ID == node_names[n], df_plants) 
    
    if size(df)[1] != 0.0 # Nodes that don't have power plants connected to them
        mean_ava = 0.0
        for j in eachindex(df.Technology) # running over all sub power plants
            technology = df.Technology[j] # Technology of the plant
            mean_ava += df[j, :Capacity] .* df_ava[t_off_peak, technology] # Multiply the installed capacity times the availability of that technology during off-peak scenario
        end

        availability_nodes[n] = mean_ava
        capacity_nodes[n] = sum(df[:, :Capacity]) # Sum of all capacities at the node
    end
end

C_tot = sum(capacity_nodes) # Similar to the Taher et. al. Paper!
A_tot = sum(availability_nodes)

x = P_tot / A_tot
y = P_tot / C_tot # Same as in Taher paper

my_df.nodal_capacity = capacity_nodes
my_df.mean_nodal_availability = availability_nodes

##
# Nodal powers
ΔP_i = y * capacity_nodes - P_i
ΔP_i_availability = x * availability_nodes - P_i

my_df.delta_P = ΔP_i
my_df.delta_P_ava = ΔP_i_availability

##
# Filter the 380kV level
my_df.voltage_380kV = falses(num_nodes)

for n in num_nodes # Run over all nodes
    if df_nodes[n, :Node_ID] ∈ nodes_380kV
        my_df[n, :voltage_380kV] = true
    end
end

df_380kV = filter(row -> row.voltage_380kV == true, my_df) 

##
# Plotting
p1 = histogram(my_df.delta_P, xaxis = L"\Delta P [MW]", color = :red, lw = 0.0, c = colorant"coral", label = "Capacity", bins = 100)
p2 = histogram(my_df.delta_P_ava, xaxis = L"\Delta P_{ava} [MW]", lw = 0.0, c = colorant"teal", label = "Availability", bins = 100)

plt = Plots.plot(p1, p2; layout = (2,1), size = (500, 500), title = "All Voltage Levels")

##
#
p1 = histogram(df_380kV.delta_P, xaxis = L"\Delta P [MW]", color = :red, lw = 0.0, c = colorant"coral",  label = "Capacity", bins = 100)
p2 = histogram(df_380kV.delta_P_ava, xaxis = L"\Delta P_{ava} [MW]", lw = 0.0, c = colorant"teal", label = "Availability", bins = 100)

plt = Plots.plot(p1, p2; layout = (2,1), size = (500, 500), title = "380kV")

##
#
net_type = my_df.delta_P .>= 0.0
my_df.net_type = net_type

plt = @df my_df groupedhist(:delta_P_ava, group = :net_type, label = ["Net Consumer" "Net Generator"], c = [colorant"sandybrown" colorant"chocolate4"], normalize=true, xlabel = L"\Delta P [MW]", ylabel = L"p(\Delta P)")

savefig(plt, "plots/power_distribution_elmod_de.pdf")

##
# Calculating the Base power and P_0 for the bimodal distribution, P_base = 100MW
P_0_full = mean(abs.(my_df.delta_P)) # for the entire grid

P_0_380kV = mean(abs.(df_380kV.delta_P)) # 380kV capacity based
P_0_380kV_ava = mean(abs.(df_380kV.delta_P_ava)) # Based on availability data