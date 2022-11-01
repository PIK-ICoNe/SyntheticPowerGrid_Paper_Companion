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
# Loading a synthetic Power Grid consisting of droop controlled inverters
file_path = joinpath(@__DIR__, "../data/powergrids/synthetic_power_grid_example.json")
pg = read_powergrid(file_path, Json) 
op = find_operationpoint(pg)
ω_indices, nodes, fluc_node_idxs, P_set, Q_set = nodal_data(pg) # Accessing the node data from the grid

##
# Using Data-Driven Load Profiles model to generate fluctuating time series
tspan = (0.0, 100.0)
t_short = collect(0.0:0.01:100.0)

Δt = 10000.0
p = 0.2 # Penetration parameter
P_fluc, t = load_profile_model(tspan, Δt = Δt)
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
plt_corr_active_power, plt_corr_frequency, plt_corr_voltage = plot_fluc_results(pg_sol_corr_demand, fluc_node_idxs, ω_indices, t = t_short)
hist_corr_voltage, hist_corr_frequency = plot_histograms(pg_sol_corr_demand, ω_indices)
mean_norm, sync_norm = calculate_performance_measures(pg_sol_corr_demand) # calculate performance measures

savefig(plt_corr_active_power, "plots/demand_fluc/multi_node_demand_fluc_correlated_active_power.png")
savefig(plt_corr_frequency, "plots/demand_fluc/multi_node_demand_fluc_correlated_frequency.png")
savefig(plt_corr_voltage, "plots/demand_fluc/multi_node_demand_fluc_correlated_voltage.png")
savefig(hist_corr_voltage, "plots/demand_fluc/multi_node_demand_fluc_correlated_voltage_histogram.png")
savefig(hist_corr_frequency, "plots/demand_fluc/multi_node_demand_fluc_correlated_frequency_histogram.png")
writedlm("data/demand_fluctuations/performance_measures_demand_correlated.txt", [mean_norm, sync_norm])

##
# Multi Node Fluctuations, completely uncorrelated, exchange all PQAlgebraic with FluctuationNode
# Generate a time series for each node

P_fluc_inter = linear_interpolation(t, P_fluc) # Interpolate the time series

flucs = [load_profile_model(tspan, Δt = Δt) for x in 1:length(fluc_node_idxs)]
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
plt_uncorr_active_power, plt_uncorr_frequency, plt_uncorr_voltage,  = plot_fluc_results(pg_sol_uncorr_demand, fluc_node_idxs, ω_indices, t = t_short)
hist_uncorr_voltage, hist_uncorr_frequency = plot_histograms(pg_sol_uncorr_demand, ω_indices)
mean_norm, sync_norm = calculate_performance_measures(pg_sol_uncorr_demand) # calculate performance measures

savefig(plt_uncorr_active_power, "plots/demand_fluc/multi_node_demand_fluc_uncorrelated_active_power.png")
savefig(plt_uncorr_frequency, "plots/demand_fluc/multi_node_demand_fluc_uncorrelated_frequency.png")
savefig(plt_uncorr_voltage, "plots/demand_fluc/multi_node_demand_fluc_uncorrelated_voltage.png")
savefig(hist_uncorr_voltage, "plots/demand_fluc/multi_node_demand_fluc_uncorrelated_voltage_histogram.png")
savefig(hist_uncorr_frequency, "plots/demand_fluc/multi_node_demand_fluc_uncorrelated_frequency_histogram.png")
writedlm("data/demand_fluctuations/performance_measures_demand_uncorrelated.txt", [mean_norm, sync_norm])