using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using SyntheticPowerGrids
using Statistics
using DelimitedFiles

##
nodal_parameters = Dict(:X => 1.0, :γ => 0.2, :α => 2.0, :τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => [0.5, 1.0 , 5.0]) 

##
num_runs = 100
num_nodes = collect(range(100, stop = 700, length = 25))
rejection_vec = zeros(length(num_nodes), num_runs)

##
for n in eachindex(num_nodes)
    println(n/eachindex(num_nodes)[end])
    pg_struct = PGGeneration(num_nodes = num_nodes[n], nodal_parameters = nodal_parameters, loads = :PQAlgebraic, lines = :PiModelLine, generation_dynamics = :Mixed)
    for r in 1:num_runs
        _, _, rejections = random_PD_grid(pg_struct)
        rejection_vec[n, r] = rejections
    end
end

##
#rejection_mean = map(x -> mean(rejection_vec[x, :]), eachindex(num_nodes))
#rejection_std = map(x -> std(rejection_vec[x, :]), eachindex(num_nodes))

##
# save Data
writedlm(joinpath(@__DIR__, "../data/rejection_rate.txt"), rejection_vec)
