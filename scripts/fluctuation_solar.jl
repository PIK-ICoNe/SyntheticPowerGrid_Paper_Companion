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
# Load all solar time series
file_paths = readdir("data/PVDataSet", join = true)
P_pv_arr = Vector{Vector{Float64}}(undef, length(file_paths)) 

iter = 1

for path in file_paths
    P_pv_temp = readdlm(path) # Load the time series
    P_pv_temp = P_pv_temp[:, 1] # Turn Matrix to an array

    nan_idxs = findall(typeof.(P_pv_temp) .!== Float64) # find all NaNs
    deleteat!(P_pv_temp, nan_idxs)

    P_pv_temp .-= mean(P_pv_temp) # Shifting the mean such that only the fluctuations are left

    P_pv_arr[iter] = P_pv_temp[191:end] # Throw away beginning of the time series!
    
    iter += 1
end

##
# Interpolate the time series
P_pv_corr = P_pv_arr[1]
Δt = 0.001
t = collect(range(0.0, length = length(P_pv_corr), step = Δt))
P_pv_inter = linear_interpolation(t, P_pv_corr)
tspan = (0.0, 100.0)

##
# Multi Node Fluctuations, completely correlated, exchange all PQAlgebraic with FluctuationNode
p = 0.2 
fluctuations_corr = map(f -> FluctuationNode(t -> P_set[f] + p * P_pv_inter(t), t -> Q_set[f]), 1:length(fluc_node_idxs))
pg_solar_corr = generate_powergrid_fluctuations(pg, fluc_node_idxs, fluctuations_corr)

##
# Simulate a trajectory
ode = ODEProblem(rhs(pg_solar_corr), op.vec, tspan)
sol_corr_solar = solve(ode, Rodas4())
pg_sol_corr_solar = PowerGridSolution(sol_corr_solar, pg_solar_corr)

##
# Results
plt_corr_active_power, plt_corr_frequency, plt_corr_voltage = plot_fluc_results(pg_sol_corr_solar, fluc_node_idxs, ω_indices, t = t_short)
hist_corr_voltage, hist_corr_frequency = plot_histograms(pg_sol_corr_solar, ω_indices)

savefig(plt_corr_active_power, "plots/solar_fluc/multi_node_solar_fluc_correlated_active_power.png")
savefig(plt_corr_frequency, "plots/solar_fluc/multi_node_solar_fluc_correlated_frequency.png")
savefig(plt_corr_voltage, "plots/solar_fluc/multi_node_solar_fluc_correlated_voltage.png")
savefig(hist_corr_voltage, "plots/solar_fluc/multi_node_solar_fluc_correlated_voltage_histogram.png")
savefig(hist_corr_frequency, "plots/solar_fluc/multi_node_solar_fluc_correlated_frequency_histogram.png")

mean_norm, sync_norm = calculate_performance_measures(pg_sol_corr_solar) # calculate performance measures
writedlm("data/solar_fluctuations/performance_measures_solar_correlated.txt", [mean_norm, sync_norm])

##
# Multi Node Fluctuations, completely uncorrelated
P_pv_inter = map(f -> linear_interpolation(t, P_pv_arr[f]), 1:length(fluc_node_idxs)) # Interpolate the time series
fluctuations_uncorr = map(f -> FluctuationNode(t -> P_set[f] + p * P_pv_inter[f](t), t -> Q_set[f]), 1:length(fluc_node_idxs))
pg_solar_uncorr = generate_powergrid_fluctuations(pg, fluc_node_idxs, fluctuations_uncorr)

##
# Simulate a trajectory
ode = ODEProblem(rhs(pg_solar_uncorr), op.vec, tspan)
sol_uncorr_solar = solve(ode, Rodas4())
pg_sol_uncorr_solar = PowerGridSolution(sol_uncorr_solar, pg_solar_uncorr)

##
# Results
plt_uncorr_active_power, plt_uncorr_frequency, plt_uncorr_voltage = plot_fluc_results(pg_sol_uncorr_solar, fluc_node_idxs, ω_indices, t = t_short)
hist_uncorr_voltage, hist_uncorr_frequency = plot_histograms(pg_sol_uncorr_solar, ω_indices)

savefig(plt_uncorr_active_power, "plots/solar_fluc/multi_node_solar_fluc_uncorrelated_active_power.png")
savefig(plt_uncorr_frequency, "plots/solar_fluc/multi_node_solar_fluc_uncorrelated_frequency.png")
savefig(plt_uncorr_voltage, "plots/solar_fluc/multi_node_solar_fluc_uncorrelated_voltage.png")
savefig(hist_uncorr_voltage, "plots/solar_fluc/multi_node_solar_fluc_uncorrelated_voltage_histogram.png")
savefig(hist_uncorr_frequency, "plots/solar_fluc/multi_node_solar_fluc_uncorrelated_frequency_histogram.png")

mean_norm, sync_norm = calculate_performance_measures(pg_sol_uncorr_solar) # calculate performance measures
writedlm("data/solar_fluctuations/performance_measures_solar_uncorrelated.txt", [mean_norm, sync_norm])