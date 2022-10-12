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
default(grid = false, foreground_color_legend = nothing, bar_edges = false,  lw=1.5, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 15, tickfontsize = 10)

##
# Generating a synthetic Power Grid consisting of droop controlled inverters
nodal_parameters_a = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :τ_P => 5.0) 
nodal_parameters_b = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :τ_P => 1.0) 
nodal_parameters_c = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :τ_P => 0.5)

nodal_dynamics = [(1/6, get_DroopControlledInverterApprox, nodal_parameters_a), (1/6, get_DroopControlledInverterApprox, nodal_parameters_b), (1/6, get_DroopControlledInverterApprox, nodal_parameters_c), (0.5, get_PQ, nothing)]
num_nodes = 100

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
# Using an intermittent wind power fluctuation Langevin-type model to generate fluctuating time series
tspan = (0.0, 100.0)
D = 0.1 # Intermittence strength
p = 0.2 # Penetration parameter
x, t = wind_power_model(tspan, D = D)
x_inter = linear_interpolation(t, x) # Interpolate the time series

plot(t, x, idxs = 1, xlabel = "t[s]", ylabel = "x(t)", label = "Time series", lw = 3)
plot!(t, x_inter(t), idxs = 1,label = "Interpolated time series", line_style = :dash, xlabel = "t[s]", ylabel = "x(t)")

##
# Multi Node Fluctuations, completely correlated, exchange all PQAlgebraic with FluctuationNode
fluctuations_corr = map(f -> FluctuationNode(t -> P_set[f] + p * x_inter(t), t -> Q_set[f]), 1:length(fluc_node_idxs))
pg_wind_corr = generate_powergrid_fluctuations(pg, fluc_node_idxs, fluctuations_corr)

##
# Simulate a trajectory
ode = ODEProblem(rhs(pg_wind_corr), op.vec, tspan)
sol = solve(ode, Rodas4())

solution2 = PowerGridSolution(sol, pg_wind_corr)
plot(solution2, fluc_node_idxs, label = "Active Power",:p, lw = 3, ylabel = L"P[p.u.]", xlabel = L"t[s]", legend = false)

plt2 = plot(solution2, ω_indices, :x_1, legend = false, ylabel = L"ω[rad / s]", xlabel = L"t[s]")
savefig(plt2, "plots/wind_fluc/multi_node_fluc_correlated.pdf")

calculate_performance_measures(solution2) # calculate performance measures

##
# Multi Node Fluctuations , completely uncorrelated, exchange all PQAlgebraic with FluctuationNode
# Generate a time series for each node

flucs = [wind_power_model(tspan, D = D) for x in 1:length(fluc_node_idxs)]
x_inter = map(f -> linear_interpolation(flucs[f][2], flucs[f][1]), 1:length(fluc_node_idxs)) # Interpolate the time series
fluctuations_uncorr = map(f -> FluctuationNode(t -> P_set[f] + p * x_inter[f](t), t -> Q_set[f]), 1:length(fluc_node_idxs))

pg_wind_uncorr = pg_wind_corr = generate_powergrid_fluctuations(pg, fluc_node_idxs, fluctuations_uncorr)

##
# Simulate a trajectory

ode = ODEProblem(rhs(pg_wind_uncorr), op.vec, tspan)
sol = solve(ode, Rodas4())

solution3 = PowerGridSolution(sol, pg_wind_uncorr)
plot(solution3, fluc_node_idxs, label = "Active Power",:p, lw = 3, ylabel = L"P[p.u.]", xlabel = L"t[s]", legend = false)

plt3 = plot(solution3, ω_indices, :x_1, legend = false, ylabel = L"ω[rad / s]", xlabel = L"t[s]")
savefig(plt3, "plots/wind_fluc/multi_node_fluc_uncorrelated.pdf")

calculate_performance_measures(solution3) # calculate performance measures