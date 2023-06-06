using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using Plots
using SyntheticPowerGrids
using PowerDynamics
using DelimitedFiles
using Random
using CairoMakie
seed = MersenneTwister(42) # Fixing the seed

## Generating a synthetic Power Grids consisting of droop controlled inverters
nodal_parameters = Dict(:τ_Q => 5.0, :K_P => 5, :K_Q => 0.1, :τ_P => 0.5)
nodal_dynamics = [(0.5, get_DroopControlledInverterApprox, nodal_parameters), (0.5, get_PQ, nothing)]

num_nodes = 100
pg_struct = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_struct)

file_path = joinpath(@__DIR__, "../data/powergrids/synthetic_power_grid_example.json")
write_powergrid(pg, file_path, Json)

## Plotting the grid structure
file_path = joinpath(@__DIR__, "../data/powergrids/synthetic_power_grid_example.json")
pg = read_powergrid(file_path, Json) 
f = my_graph_plot(pg) 
Makie.save(joinpath(@__DIR__, "../plots/grid_structure.svg"), f)

## Generate a test grid with the size of the ELMOD grid
parameters_third_order = Dict(:X => 1.0, :γ => 0.2, :α => 2.0) 
parameters_droop_controlled = Dict(:τ_Q => 5.0, :K_P => 5, :K_Q => 5.0, :τ_P => 1.0) 

nodal_dynamics = [(1/3, get_ThirdOrderMachineApprox, parameters_third_order), (1/3, get_DroopControlledInverterApprox, parameters_droop_controlled), (1/3, get_PQ, nothing)]

num_nodes = 438

pg_struct = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics, lines = :StaticLine, P0 = 3.2)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_struct)

file_path = joinpath(@__DIR__, "../data/powergrids/synthetic_power_grid_large.json")
write_powergrid(pg, file_path, Json)
writedlm(joinpath(@__DIR__, "../data/powergrids/power_grid_large_vertex_positions.txt"), pg_struct_new.embedded_graph.vertexpos)

## Generate a test grid with the peak demand 
num_nodes = 100

pg_struct = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics, lines = :StaticLine, P0 = 3.2)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_struct)

file_path = joinpath(@__DIR__, "../data/powergrids/synthetic_power_grid_peak_demand.json")
write_powergrid(pg, file_path, Json)
writedlm(joinpath(@__DIR__, "../data/powergrids/power_grid_peak_demand_vertex_positions.txt"), pg_struct_new.embedded_graph.vertexpos)

## Grid with power normalized via area
using Distributions
num_nodes = 100
P_base = 1.0e8 
P0 = 1.0

power_dist = MixtureModel(Normal[Normal(P0, P0/2), Normal(-P0, P0/2)]) # Distribution for the active power
P_vec = rand(power_dist, num_nodes)                                    # Power Generation / Consumption of the nodes
P_vec .-= sum(P_vec) / (num_nodes)                                     # Assure power balance

P_tot_on_peak = (86031 * 10^6/4) * (1 / P_base) # MW to p.u. 
consumers = findall(P_vec .< 0.0)
demand_synthetic = abs.(sum(P_vec[consumers]))

c_demand = P_tot_on_peak / demand_synthetic

P_re = c_demand * P_vec

mean(abs.(P_re))

##
pg_struct = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics, lines = :StaticLine , P_vec = P_re, maxiters = 10)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_struct)

file_path = joinpath(@__DIR__, "../data/powergrids/power_grid_area_normalized.json")
write_powergrid(pg, file_path, Json)
writedlm(joinpath(@__DIR__, "../data/powergrids/power_grid_area_normalized_vertex_positions.txt"), pg_struct_new.embedded_graph.vertexpos)

