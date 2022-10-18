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
# Using Data-Driven Load Profiles model to generate fluctuating time series
tspan = (0.0, 30.0)
p = 0.2 # Penetration parameter
P_fluc, t = load_profile_model(tspan)
P_mean = mean(P_fluc)
P_fluc = (P_fluc .- P_mean) ./ P_mean
P_fluc_inter = linear_interpolation(t, P_fluc) # Interpolate the time series

##
# Multi Node Fluctuations, completely correlated, exchange all PQAlgebraic with FluctuationNode
fluctuations_corr = map(f -> FluctuationNode(t -> P_set[f] + p * P_fluc_inter(t), t -> Q_set[f]), 1:length(fluc_node_idxs))
pg_demand_corr = generate_powergrid_fluctuations(pg, fluc_node_idxs, fluctuations_corr)

##
# Simulate a trajectory
ode = ODEProblem(rhs(pg_demand_corr), op.vec, tspan)
sol_corr = solve(ode, Rodas4())
pg_sol_corr_demand = PowerGridSolution(sol_corr, pg_demand_corr)

##
# Results
plt_corr_active_power, plt_corr_frequency, plt_corr_voltage, hist_corr_voltage, hist_corr_frequency  = plot_fluc_results(pg_sol_corr_demand, fluc_node_idxs, ω_indices)

savefig(plt_corr_active_power, "plots/demand_fluc/multi_node_demand_fluc_correlated_active_power.pdf")
savefig(plt_corr_frequency, "plots/demand_fluc/multi_node_demand_fluc_correlated_frequency.pdf")
savefig(plt_corr_voltage, "plots/demand_fluc/multi_node_demand_fluc_correlated_voltage.pdf")
savefig(hist_corr_voltage, "plots/demand_fluc/multi_node_demand_fluc_correlated_voltage_histogram.pdf")
savefig(hist_corr_frequency, "plots/demand_fluc/multi_node_demand_fluc_correlated_frequency_histogram.pdf")

calculate_performance_measures(pg_sol_corr_demand) # calculate performance measures

##
# Multi Node Fluctuations, completely uncorrelated, exchange all PQAlgebraic with FluctuationNode
# Generate a time series for each node

P_fluc_inter = linear_interpolation(t, P_fluc) # Interpolate the time series

flucs = [load_profile_model(tspan) for x in 1:length(fluc_node_idxs)]
t = map(f -> flucs[f][2], 1:length(fluc_node_idxs))
P_mean = map(f -> mean(flucs[f][1]), 1:length(fluc_node_idxs))
P_flucs = map(f -> (flucs[f][1] .- P_mean[f]) ./ P_mean[f], 1:length(fluc_node_idxs))

P_fluc_inter = map(f -> linear_interpolation(t[f], P_flucs[f]), 1:length(fluc_node_idxs)) # Interpolate the time series
fluctuations_uncorr = map(f -> FluctuationNode(t -> P_set[f] + p * P_fluc_inter[f](t), t -> Q_set[f]), 1:length(fluc_node_idxs))

##
pg_wind_uncorr = generate_powergrid_fluctuations(pg, fluc_node_idxs, fluctuations_uncorr)

##
# Simulate a trajectory
ode = ODEProblem(rhs(pg_wind_uncorr), op.vec, tspan)
sol_uncorr = solve(ode, Rodas4())
pg_sol_uncorr_demand = PowerGridSolution(sol_uncorr, pg_wind_uncorr)

##
# Results
plt_uncorr_active_power, plt_uncorr_frequency, plt_uncorr_voltage, hist_uncorr_voltage, hist_uncorr_frequency = plot_fluc_results(pg_sol_uncorr_demand, fluc_node_idxs, ω_indices)

savefig(plt_uncorr_active_power, "plots/demand_fluc/multi_node_demand_fluc_uncorrelated_active_power.pdf")
savefig(plt_uncorr_frequency, "plots/demand_fluc/multi_node_demand_fluc_uncorrelated_frequency.pdf")
savefig(plt_uncorr_voltage, "plots/demand_fluc/multi_node_demand_fluc_uncorrelated_voltage.pdf")
savefig(hist_uncorr_voltage, "plots/demand_fluc/multi_node_demand_fluc_uncorrelated_voltage_histogram.pdf")
savefig(hist_uncorr_frequency, "plots/demand_fluc/multi_node_demand_fluc_uncorrelated_frequency_histogram.pdf")

calculate_performance_measures(pg_sol_uncorr_demand) # calculate performance measures