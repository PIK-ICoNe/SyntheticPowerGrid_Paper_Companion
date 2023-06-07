Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using PowerDynamics
using OrdinaryDiffEq
using PowerGridNoise  
using Interpolations
using Statistics
using DelimitedFiles

## 
step_size = 0.01
tspan = (0.0, 1000.0)
t = collect(tspan[1]:step_size:tspan[2])

tspan_plt = (100.0, 200.0)
plt_start = Int(tspan_plt[1] / step_size)
plt_end = Int(tspan_plt[2] / step_size) 
t_plt = collect(0.0:step_size:tspan_plt[2]-tspan_plt[1])

## Load the data
f_d_uc = readdlm("data/demand_fluctuations/uncorr/frequencies.txt")
f_w_uc = readdlm("data/wind_fluctuations/uncorr/frequencies.txt")
f_s_uc = readdlm("data/solar_fluctuations/uncorr/frequencies.txt")

f_d_c = readdlm("data/demand_fluctuations/corr/frequencies.txt")
f_w_c = readdlm("data/wind_fluctuations/corr/frequencies.txt")
f_s_c = readdlm("data/solar_fluctuations/corr/frequencies.txt")

p_d_uc = readdlm("data/demand_fluctuations/uncorr/powers.txt")
p_w_uc = readdlm("data/wind_fluctuations/uncorr/powers.txt")
p_s_uc = readdlm("data/solar_fluctuations/uncorr/powers.txt")

p_d_c = readdlm("data/demand_fluctuations/corr/powers.txt")
p_w_c = readdlm("data/wind_fluctuations/corr/powers.txt")
p_s_c = readdlm("data/solar_fluctuations/corr/powers.txt")

rocof_d_uc = readdlm("data/demand_fluctuations/uncorr/rocof.txt")
rocof_w_uc = readdlm("data/wind_fluctuations/uncorr/rocof.txt")
rocof_s_uc = readdlm("data/solar_fluctuations/uncorr/rocof.txt")

rocof_d_c = readdlm("data/demand_fluctuations/corr/rocof.txt")
rocof_w_c = readdlm("data/wind_fluctuations/corr/rocof.txt")
rocof_s_c = readdlm("data/solar_fluctuations/corr/rocof.txt")

## Performance Measures
calculate_performance_measures(f_d_uc, rocof_d_uc; T = 1000.0, Δt = 0.01)
calculate_performance_measures(f_d_c, rocof_d_c; T = 1000.0, Δt = 0.01)

calculate_performance_measures(f_w_uc, rocof_w_uc; T = 1000.0, Δt = 0.01)
calculate_performance_measures(f_w_c, rocof_w_c; T = 1000.0, Δt = 0.01)

calculate_performance_measures(f_s_uc, rocof_s_uc; T = 1000.0, Δt = 0.01)
calculate_performance_measures(f_s_c, rocof_s_c; T = 1000.0, Δt = 0.01)

## Plot power time series
plt_p_d_c = plot(t_plt, p_d_c[plt_start:plt_end, :], legend = false, ylabel = L"P[p.u.]", xlabel = L"t[s]")
savefig(plt_p_d_c, "plots/demand_fluc/multi_node_demand_fluc_correlated_active_power.png")

plt_p_d_uc = plot(t_plt, p_d_uc[plt_start:plt_end, :], legend = false, ylabel = L"P[p.u.]", xlabel = L"t[s]")
savefig(plt_p_d_uc, "plots/demand_fluc/multi_node_demand_fluc_uncorrelated_active_power.png")

plt_p_w_c = plot(t_plt, p_w_c[plt_start:plt_end, :], legend = false, ylabel = L"P[p.u.]", xlabel = L"t[s]")
savefig(plt_p_w_c, "plots/wind_fluc/multi_node_wind_fluc_correlated_active_power.png")

plt_p_w_uc = plot(t_plt, p_w_uc[plt_start:plt_end, :], legend = false, ylabel = L"P[p.u.]", xlabel = L"t[s]")
savefig(plt_p_w_uc, "plots/wind_fluc/multi_node_wind_fluc_uncorrelated_active_power.png")

plt_p_s_c = plot(t_plt, p_s_c[plt_start:plt_end, :], legend = false, ylabel = L"P[p.u.]", xlabel = L"t[s]")
savefig(plt_p_s_c, "plots/solar_fluc/multi_node_solar_fluc_correlated_active_power.png")

plt_p_s_uc = plot(t_plt, p_s_uc[plt_start:plt_end, :], legend = false, ylabel = L"P[p.u.]", xlabel = L"t[s]")
savefig(plt_p_s_uc, "plots/solar_fluc/multi_node_solar_fluc_uncorrelated_active_power.png")

## Plot frequency time series
plt_f_d_c = plot(t_plt, f_d_c[plt_start:plt_end, :], legend = false, ylabel = L"\Delta f [Hz]", xlabel = L"t[s]")
savefig(plt_f_d_c, "plots/demand_fluc/multi_node_demand_fluc_correlated_frequency.png")

plt_f_d_uc = plot(t_plt, f_d_uc[plt_start:plt_end, :], legend = false, ylabel = L"\Delta f [Hz]", xlabel = L"t[s]")
savefig(plt_f_d_uc, "plots/demand_fluc/multi_node_demand_fluc_uncorrelated_frequency.png")

plt_f_w_c = plot(t_plt, f_w_c[plt_start:plt_end, :], legend = false, ylabel = L"\Delta f [Hz]", xlabel = L"t[s]")
savefig(plt_f_w_c, "plots/wind_fluc/multi_node_wind_fluc_correlated_frequency.png")

plt_f_w_uc = plot(t_plt, f_w_uc[plt_start:plt_end, :], legend = false, ylabel = L"\Delta f [Hz]", xlabel = L"t[s]")
savefig(plt_f_w_uc, "plots/wind_fluc/multi_node_wind_fluc_uncorrelated_frequency.png")

plt_f_s_c = plot(t_plt, f_s_c[plt_start:plt_end, :], legend = false, ylabel = L"\Delta f [Hz]", xlabel = L"t[s]")
savefig(plt_f_s_c, "plots/solar_fluc/multi_node_solar_fluc_correlated_frequency.png")

plt_f_s_uc = plot(t_plt, f_s_uc[plt_start:plt_end, :], legend = false, ylabel = L"\Delta f [Hz]", xlabel = L"t[s]")
savefig(plt_f_s_uc, "plots/solar_fluc/multi_node_solar_fluc_uncorrelated_frequency.png")
