using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using SyntheticPowerGrids
using PowerDynamics
using OrdinaryDiffEq
using Plots
using PowerGridNoise
using Interpolations
using Statistics
using LaTeXStrings
using DelimitedFiles
default(grid = false, foreground_color_legend = nothing, bar_edges = false,  lw=1.5, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 15, tickfontsize = 10)

##
# Load solar time series
path = joinpath(@__DIR__, "../data/data_lambda0.1_diff0.001.txt")
Z = readdlm(path)
Z = f[:, 1]

nan_idxs = findall(typeof.(Z) .!== Float64) # find all NaNs
deleteat!(Z, nan_idxs)

Δt = 0.001
t = collect(range(0.0, length = length(Z), step = Δt))

Z_inter = linear_interpolation(t, Z)

plot(t, Z, legend = false, xlabel = L"t[s]", ylabel = L"Z")
plot!(t, Z_inter(t), legend = false)

p = 0.1
tspan = (0.0, 30.0)

##
# Generating a synthetic Power Grid consisting of droop controlled inverters
nodal_parameters_a = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :τ_P => 5.0) 
nodal_parameters_b = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :τ_P => 1.0) 
nodal_parameters_c = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :τ_P => 0.5)

nodal_dynamics = [(1/6, get_DroopControlledInverterApprox, nodal_parameters_a), (1/6, get_DroopControlledInverterApprox, nodal_parameters_b), (1/6, get_DroopControlledInverterApprox, nodal_parameters_c), (0.5, get_PQ, nothing)]
num_nodes = 20

a = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics, lines = :StaticLine)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(a)

##
# Accessing the node data from the grid
ω_indices = findall(n -> :x_1 ∈ symbolsof(n), pg.nodes)
nodes = deepcopy(pg.nodes) 
fluc_node_idxs = findall(typeof.(pg.nodes) .== PQAlgebraic) # Find all Load Buses in the grid
P_set = map(i -> nodes[i].P, fluc_node_idxs) # Load their power set-points
Q_set = map(i -> nodes[i].Q, fluc_node_idxs)

##
# Multi Node Fluctuations, completely correlated, exchange all PQAlgebraic with FluctuationNode
fluctuations_corr = map(f -> FluctuationNode(t -> P_set[f] * (1 + Z_inter(t)), t -> Q_set[f]), 1:length(fluc_node_idxs))
pg_solar_corr = generate_powergrid_fluctuations(pg, fluc_node_idxs, fluctuations_corr)

##
# Simulate a trajectory
ode = ODEProblem(rhs(pg_solar_corr), op.vec, tspan)
sol = solve(ode, Rodas4())

solution2 = PowerGridSolution(sol, pg_solar_corr)
plot(solution2, fluc_node_idxs, label = "Active Power",:p, lw = 3, ylabel = L"P[p.u.]", xlabel = L"t[s]", legend = false)

plt2 = plot(solution2, ω_indices, :x_1, legend = false, ylabel = L"ω[rad / s]", xlabel = L"t[s]")
savefig(plt2, "plots/solar_fluc/multi_node_solar_fluc_correlated.pdf")

calculate_performance_measures(solution2) # calculate performance measures

##
# Multi Node Fluctuations, completely uncorrelated
# Time series is not long enough...
t_end = floor(t[end] / num_nodes, digits=1)
tspan = (0.0, t_end)

fluctuations_uncorr = map(f -> FluctuationNode(t -> P_set[f] * (1 + Z_inter(t + t_end * f)), t -> Q_set[f]), 1:length(fluc_node_idxs))
pg_solar_uncorr = generate_powergrid_fluctuations(pg, fluc_node_idxs, fluctuations_uncorr)

##
# Simulate a trajectory
ode = ODEProblem(rhs(pg_solar_uncorr), op.vec, tspan)
sol = solve(ode, Rodas4())

solution3 = PowerGridSolution(sol, pg_solar_uncorr)
plot(solution3, fluc_node_idxs, label = "Active Power",:p, lw = 3, ylabel = L"P[p.u.]", xlabel = L"t[s]", legend = false)

plt3 = plot(solution3, ω_indices, :x_1, legend = false, ylabel = L"ω[rad / s]", xlabel = L"t[s]")
savefig(plt3, "plots/solar_fluc/multi_node_solar_fluc_uncorrelated.pdf")

calculate_performance_measures(solution3) # calculate performance measures
