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

##
# Analyzing the ELMOD-De Dataset
path = string(@__DIR__) * "/../data/elmod-de_v1/"
elmod_de = joinpath(path, "Data_Input.xlsx")

##
# Load the line data
lines_table = XLSX.readtable(elmod_de, "Grid_technical")
df_lines = DataFrame(lines_table)

# Add ID for the lines for make them searchable
line_id = collect(1:size(df_lines)[1])
df_lines[!, :line_id] = line_id 

# Find only lines in the 380kV level
df_380kV = filter(row -> row.voltage == 380, df_lines) 
mean(df_380kV[:, :length]) # Sanity Check -> very close to the value in the SciGrid dataset! :-)

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

nodes_380kV = String.(unique(nodes_380kV)) # Nodes will at least appear twice, String for filtering

##
# Calculating the mean hourly availability for the different technologies
ava_table = XLSX.readtable(elmod_de, "H_con")
df_ava = DataFrame(ava_table)

mean_availability_technology = Dict()

for i in eachindex(names(df_ava))
    mean_availability_technology[names(df_ava)[i]] = mean(df_ava[:, i])
end

##
# Loading the plant data
plant_table = XLSX.readtable(elmod_de, "Plant_con", first_row = 3, column_labels = ["Plant_ID", "Country", "Node_ID", "State", "Dena_zone", "ZoDE_NE_6", "Technology", "Fuel", "Capacity", "Efficiency", "Emission", "Transport", "Not_used_1", "Name"])
df_plants = DataFrame(plant_table) 

replace!(df_plants.Name, missing => "Unknown")

df_plants.Name = String.(df_plants.Name) # Turn name into strings for filtering
power_plants = unique(df_plants[:, :Name]) # Some nodes/buses have multiple sub power plants connected to them

capacity = Vector{Float64}(undef, length(power_plants)) # Installed Capacity
mean_availability_plant = zeros(length(power_plants))   # Mean availability of the plant

for plant in eachindex(power_plants) # Run over all power plants
    if df_plants.Node_ID[plant] ∈ nodes_380kV # Looking only for node in the 380kV level
        df = filter(row -> row.Name == power_plants[plant], df_plants) 
        mean_ava = 0.0

        for j in eachindex(df.Technology) # running over all sub power plants
            technology = df.Technology[j] # Technology of the plant
            mean_ava += df[j, :Capacity] .* mean_availability_technology[technology] # Multiply the installed capacity times the mean availability of that technology
        end

        mean_availability_plant[plant] = mean_ava
        capacity[plant] = sum(df[:, :Capacity]) # Sum of all capacities of the sub power plants 
    end
end

##
# P_base ??

## ToDo add pump storage files as well
mean(capacity)
mean(mean_availability_plant)

##
# Plotting
p1 = histogram(capacity, xaxis = "Capacity [MW]", color = :red, lw = 0.0, c = colorant"coral")
p2 = histogram(mean_availability_plant, xaxis = "Mean availability [MW]", lw = 0.0, c = colorant"teal")

plt = Plots.plot(p1, p2; layout = (2,1), size = (500, 500), legend = false)

##
# Demand table -> i don't think this is relevant for us...
demand_table = XLSX.readtable(elmod_de, "H_demand", first_row = 1, column_labels = ["time_point", "Demand"])
df_demand = DataFrame(demand_table) 

histogram(df_demand[:, :Demand])

