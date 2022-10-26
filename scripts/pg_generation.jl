using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using Plots
using SyntheticPowerGrids
using PowerDynamics
using Random
using CairoMakie
seed = MersenneTwister(42) # Fixing the seed

##
# Generating a synthetic Power Grid consisting of droop controlled inverters
nodal_parameters = Dict(:τ_Q => 5.0, :K_P => 5, :K_Q => 0.1, :τ_P => 0.5)
nodal_dynamics = [(0.5, get_DroopControlledInverterApprox, nodal_parameters), (0.5, get_PQ, nothing)]
num_nodes = 100

pg_struct = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_struct)

##
# Saving the power grid for later usage
file_path = joinpath(@__DIR__, "../data/powergrids/synthetic_power_grid_example.json")
write_powergrid(pg, file_path, Json)

##
# Plotting the grid structure

f = my_graph_plot(pg)
Makie.save(joinpath(@__DIR__, "../plots/grid_structure.pdf"), f)
#read_powergrid(file_path, Json) # will work once the new pr is merged
