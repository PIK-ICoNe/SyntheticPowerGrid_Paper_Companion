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
# Generating a synthetic Power Grid consisting of droop controlled inverters
nodal_parameters = Dict(:τ_Q => 5.0, :K_P => 5, :K_Q => 0.1, :τ_P => 0.5)

nodal_dynamics = [(0.5, get_DroopControlledInverterApprox, nodal_parameters), (0.5, get_PQ, nothing)]
num_nodes = 20

a = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(a)

##
# Load solar time series
path = joinpath(@__DIR__, "../data/data_lambda0.1_diff0.001.txt")
P_pv = readdlm(path)
P_pv = P_pv[:, 1]

nan_idxs = findall(typeof.(P_pv) .!== Float64) # find all NaNs
deleteat!(P_pv, nan_idxs)

P_pv .-= mean(P_pv) # Shifting the mean such that only the fluctuations are left

##
# Interpolate the time series
Δt = 0.001
t = collect(range(0.0, length = length(P_pv), step = Δt))
P_pv_inter = linear_interpolation(t, P_pv)
tspan = (0.0, 30.0)

##
# Accessing the node data from the grid
ω_indices = findall(n -> :x_1 ∈ symbolsof(n), pg.nodes)
nodes = deepcopy(pg.nodes) 
fluc_node_idxs = findall(typeof.(pg.nodes) .== PQAlgebraic) # Find all Load Buses in the grid
P_set = map(i -> nodes[i].P, fluc_node_idxs) # Load their power set-points
Q_set = map(i -> nodes[i].Q, fluc_node_idxs)

##
# Multi Node Fluctuations, completely correlated, exchange all PQAlgebraic with FluctuationNode
fluctuations_corr = map(f -> FluctuationNode(t -> P_set[f] + P_pv_inter(t), t -> Q_set[f]), 1:length(fluc_node_idxs))
pg_solar_corr = generate_powergrid_fluctuations(pg, fluc_node_idxs, fluctuations_corr)

##
# Simulate a trajectory
ode = ODEProblem(rhs(pg_solar_corr), op.vec, tspan)
sol_corr_solar = solve(ode, Rodas4())
pg_sol_corr_solar = PowerGridSolution(sol_corr_solar, pg_solar_corr)

##
# Results
plt_corr_active_power, plt_corr_frequency, plt_corr_voltage, hist_corr_voltage, hist_corr_frequency = plot_fluc_results(pg_sol_corr_solar, fluc_node_idxs, ω_indices)

savefig(plt_corr_active_power, "plots/solar_fluc/multi_node_solar_fluc_correlated_active_power.png")
savefig(plt_corr_frequency, "plots/solar_fluc/multi_node_solar_fluc_correlated_frequency.png")
savefig(plt_corr_voltage, "plots/solar_fluc/multi_node_solar_fluc_correlated_voltage.png")
savefig(hist_corr_voltage, "plots/solar_fluc/multi_node_solar_fluc_correlated_voltage_histogram.png")
savefig(hist_corr_frequency, "plots/solar_fluc/multi_node_solar_fluc_correlated_frequency_histogram.png")

mean_norm, sync_norm = calculate_performance_measures(pg_sol_corr_solar) # calculate performance measures
writedlm("data/solar_fluctuations/performance_measures_solar_correlated.txt", [mean_norm, sync_norm])

##
# Multi Node Fluctuations, completely uncorrelated
# Time series is not long enough...
t_end = floor(t[end] / num_nodes, digits=1)
tspan = (0.0, t_end)

fluctuations_uncorr = map(f -> FluctuationNode(t -> P_set[f] + P_pv_inter(t + t_end * f), t -> Q_set[f]), 1:length(fluc_node_idxs))
pg_solar_uncorr = generate_powergrid_fluctuations(pg, fluc_node_idxs, fluctuations_uncorr)

##
# Simulate a trajectory
ode = ODEProblem(rhs(pg_solar_uncorr), op.vec, tspan)
sol_uncorr_solar = solve(ode, Rodas4())
pg_sol_uncorr_solar = PowerGridSolution(sol_uncorr_solar, pg_solar_uncorr)

##
# Results
plt_uncorr_active_power, plt_uncorr_frequency, plt_uncorr_voltage, hist_uncorr_voltage, hist_uncorr_frequency = plot_fluc_results(pg_sol_uncorr_solar, fluc_node_idxs, ω_indices)

savefig(plt_uncorr_active_power, "plots/solar_fluc/multi_node_solar_fluc_uncorrelated_active_power.png")
savefig(plt_uncorr_frequency, "plots/solar_fluc/multi_node_solar_fluc_uncorrelated_frequency.png")
savefig(plt_uncorr_voltage, "plots/solar_fluc/multi_node_solar_fluc_uncorrelated_voltage.png")
savefig(hist_uncorr_voltage, "plots/solar_fluc/multi_node_solar_fluc_uncorrelated_voltage_histogram.png")
savefig(hist_uncorr_frequency, "plots/solar_fluc/multi_node_solar_fluc_uncorrelated_frequency_histogram.png")

mean_norm, sync_norm = calculate_performance_measures(pg_sol_uncorr_solar) # calculate performance measures
writedlm("data/solar_fluctuations/performance_measures_solar_uncorrelated.txt", [mean_norm, sync_norm])