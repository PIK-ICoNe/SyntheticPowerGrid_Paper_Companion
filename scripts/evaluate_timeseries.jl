using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using PowerDynamics
using OrdinaryDiffEq
using PowerGridNoise  
using Interpolations
using Statistics
using Plots
using DelimitedFiles
using LaTeXStrings

default(grid = false, foreground_color_legend = nothing, bar_edges = false,  lw=3, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 18, labelfontsize = 18, tickfontsize = 15)
## 
step_size = 0.01
tspan_plt = (100.0, 200.0)
plt_start = Int(tspan_plt[1] / step_size)
plt_end = Int(tspan_plt[2] / step_size) 
t_plt = collect(0.0:step_size:tspan_plt[2]-tspan_plt[1])

## Load the correlated data
f_d = readdlm("data/demand_fluctuations/corr/frequencies.txt")
f_w = readdlm("data/wind_fluctuations/corr/frequencies.txt")
f_s = readdlm("data/solar_fluctuations/corr/frequencies.txt")

p_d = readdlm("data/demand_fluctuations/corr/powers.txt")
p_w = readdlm("data/wind_fluctuations/corr/powers.txt")
p_s = readdlm("data/solar_fluctuations/corr/powers.txt")

rocof_d = readdlm("data/demand_fluctuations/corr/rocof.txt")
rocof_w = readdlm("data/wind_fluctuations/corr/rocof.txt")
rocof_s = readdlm("data/solar_fluctuations/corr/rocof.txt")

## Performance Measures
calculate_performance_measures(f_d, rocof_d; T = 10000.0, Δt = step_size)
calculate_performance_measures(f_w, rocof_w; T = 10000.0, Δt = step_size)
calculate_performance_measures(f_s, rocof_s; T = 1000.0, Δt = 0.001)

## Plot power time series
plt_p_d = plot(t_plt, p_d[plt_start:plt_end, :], legend = false, ylabel = L"P[p.u.]", xlabel = L"t[s]")
savefig(plt_p_d, "plots/demand_fluc/multi_node_demand_fluc_correlated_active_power.png")

plt_p_w = plot(t_plt, p_w[plt_start:plt_end, :], legend = false, ylabel = L"P[p.u.]", xlabel = L"t[s]")
savefig(plt_p_w, "plots/wind_fluc/multi_node_wind_fluc_correlated_active_power.png")

plt_p_s = plot(t_plt, p_s[plt_start:plt_end, :], legend = false, ylabel = L"P[p.u.]", xlabel = L"t[s]")
savefig(plt_p_s, "plots/solar_fluc/multi_node_solar_fluc_correlated_active_power.png")

## Plot frequency time series
plt_f_d = plot(t_plt, f_d[plt_start:plt_end, :], legend = false, ylabel = L"\Delta f [Hz]", xlabel = L"t[s]")
savefig(plt_f_d, "plots/demand_fluc/multi_node_demand_fluc_correlated_frequency.png")

plt_f_w = plot(t_plt, f_w[plt_start:plt_end, :], legend = false, ylabel = L"\Delta f [Hz]", xlabel = L"t[s]")
savefig(plt_f_w, "plots/wind_fluc/multi_node_wind_fluc_correlated_frequency.png")

plt_f_s = plot(t_plt, f_s[plt_start:plt_end, :], legend = false, ylabel = L"\Delta f [Hz]", xlabel = L"t[s]")
savefig(plt_f_s, "plots/solar_fluc/multi_node_solar_fluc_correlated_frequency.png")

## Load the uncorrelated data
f_d = readdlm("data/demand_fluctuations/uncorr/frequencies.txt")
f_w = readdlm("data/wind_fluctuations/uncorr/frequencies.txt")
f_s = readdlm("data/solar_fluctuations/uncorr/frequencies.txt")

p_d = readdlm("data/demand_fluctuations/uncorr/powers.txt")
p_w = readdlm("data/wind_fluctuations/uncorr/powers.txt")
p_s = readdlm("data/solar_fluctuations/uncorr/powers.txt")

rocof_d = readdlm("data/demand_fluctuations/uncorr/rocof.txt")
rocof_w = readdlm("data/wind_fluctuations/uncorr/rocof.txt")
rocof_s = readdlm("data/solar_fluctuations/uncorr/rocof.txt")

## Performance Measures
calculate_performance_measures(f_d, rocof_d; T = 10000.0, Δt = step_size)
calculate_performance_measures(f_w, rocof_w; T = 10000.0, Δt = step_size)
calculate_performance_measures(f_s, rocof_s; T = 1000.0, Δt = 0.001)

## Plot power time series
plt_p_d = plot(t_plt, p_d[plt_start:plt_end, :], legend = false, ylabel = L"P[p.u.]", xlabel = L"t[s]")
savefig(plt_p_d, "plots/demand_fluc/multi_node_demand_fluc_uncorrelated_active_power.png")

plt_p_w = plot(t_plt, p_w[plt_start:plt_end, :], legend = false, ylabel = L"P[p.u.]", xlabel = L"t[s]")
savefig(plt_p_w, "plots/wind_fluc/multi_node_wind_fluc_uncorrelated_active_power.png")

plt_p_s = plot(t_plt, p_s[plt_start:plt_end, :], legend = false, ylabel = L"P[p.u.]", xlabel = L"t[s]")
savefig(plt_p_s, "plots/solar_fluc/multi_node_solar_fluc_uncorrelated_active_power.png")

## Plot frequency time series
plt_f_d = plot(t_plt, f_d[plt_start:plt_end, :], legend = false, ylabel = L"\Delta f [Hz]", xlabel = L"t[s]")
savefig(plt_f_d, "plots/demand_fluc/multi_node_demand_fluc_uncorrelated_frequency.png")

plt_f_w = plot(t_plt, f_w[plt_start:plt_end, :], legend = false, ylabel = L"\Delta f [Hz]", xlabel = L"t[s]")
savefig(plt_f_w, "plots/wind_fluc/multi_node_wind_fluc_uncorrelated_frequency.png")

plt_f_s = plot(t_plt, f_s[plt_start:plt_end, :], legend = false, ylabel = L"\Delta f [Hz]", xlabel = L"t[s]")
savefig(plt_f_s, "plots/solar_fluc/multi_node_solar_fluc_uncorrelated_frequency.png")