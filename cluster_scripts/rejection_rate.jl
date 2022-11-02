using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using SyntheticPowerGrids
using Statistics
using DelimitedFiles

##
num_runs = 100
num_nodes = collect(range(700, stop = 1300, length = 25))
rejection_vec = zeros(length(num_nodes), num_runs)

##
# Generating a synthetic Power Grid consisting of droop controlled inverters
parameters_third_order = Dict(:X => 1.0, :γ => 0.2, :α => 2.0) 
parameters_droop_controlled = Dict(:τ_Q => 5.0, :K_P => 5, :K_Q => 5.0, :τ_P => 1.0) 

nodal_dynamics = [(1/3, get_ThirdOrderMachineApprox, parameters_third_order), (1/3, get_DroopControlledInverterApprox, parameters_droop_controlled), (1/3, get_PQ, nothing)]

##
for n in eachindex(num_nodes)
    println(n/eachindex(num_nodes)[end])
    pg_struct = PGGeneration(num_nodes = Int64(num_nodes[n]), nodal_dynamics = nodal_dynamics, slack = true, lines = :StaticLine)
    for r in 1:num_runs
        _, _, _, rejections = generate_powergrid_dynamics(pg_struct)
        rejection_vec[n, r] = rejections.total
    end
end

##
#rejection_mean = map(x -> mean(rejection_vec[x, :]), eachindex(num_nodes))
#rejection_std = map(x -> std(rejection_vec[x, :]), eachindex(num_nodes))

##
# save Data
writedlm(joinpath(@__DIR__, "../data/rejection_rate.txt"), rejection_vec)
