using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using Plots
using Statistics
using LaTeXStrings
using Distributions
using DelimitedFiles
using KernelDensity
using StatsBase

default(grid = false, foreground_color_legend = nothing, bar_edges = false,  lw=3, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 18, labelfontsize = 18, tickfontsize = 15)

##
step_size = 0.01

f_demand = readdlm("data/demand_fluctuations/uncorr/frequencies.txt")
f_wind = readdlm("data/wind_fluctuations/uncorr/frequencies.txt")
f_solar = readdlm("data/solar_fluctuations/uncorr/frequencies.txt")

f_d_mean = map(n -> mean(f_demand[n, :]), 1:size(f_demand)[1])
f_w_mean = map(n -> mean(f_wind[n, :]), 1:size(f_wind)[1])
f_s_mean = map(n -> mean(f_solar[n, :]), 1:size(f_solar)[1])

## Statistics
mean(f_demand)
mean(f_wind)
mean(f_solar)

std(f_demand)
std(f_wind)
std(f_solar)

## Auto correlation
lags = collect(0:1:size(f_demand)[1]-1) 
dts = lags * step_size / 60

f_d_ac = autocor(f_d_mean, lags)
f_w_ac = autocor(f_w_mean, lags)
f_s_ac = autocor(f_s_mean, lags)

plot(dts, f_w_ac, label = L"f_w")
plot!(dts, f_d_ac, label = L"f_d")
plot!(dts, f_s_ac, ylabel = L"c(\Delta t)", xlabel = L"\Delta t [min]", label = L"f_s")

## Probability Density Function
pdf_demand = kde(vcat(f_demand...))
pdf_wind = kde(vcat(f_wind...))
pdf_solar = kde(vcat(f_solar...))

plot(pdf_demand.x, pdf_demand.density ./ findmax(pdf_demand.density)[1], label = L"f_d")
plot!(pdf_wind.x, pdf_wind.density ./ findmax(pdf_wind.density)[1], label = L"f_w")
plot!(pdf_solar.x, pdf_solar.density ./ findmax(pdf_solar.density)[1], label = L"f_s", xaxis = L"\Delta f [Hz]")

## Frequency Increment Statistics
Θ = 0.2
steps_Θ = Int(0.2 / step_size)

ΔΘf = []

for s in 1:(length(f_w_mean) - steps_Θ)
    for n in 1:size(f_wind)[2]
        append!(ΔΘf, f_wind[s, n] - f_wind[s + steps_Θ, n])
    end
end

ΔΘf = ΔΘf .|> Float64

pdf_increments = kde(ΔΘf)
σ_ΔΘf = std(ΔΘf) 
d_normal = fit(Normal{Float64}, ΔΘf)
pdf_normal = map(i -> pdf(d_normal, pdf_increments.x[i]), 1:length(pdf_increments.x))

plot(pdf_increments.x / σ_ΔΘf, pdf_normal)
plot!(pdf_increments.x / σ_ΔΘf, pdf_increments.density, legend = false, xlabel = L"\Delta_{\Theta} f / \sigma_{\Delta_{\Theta} f}", ylabel = L"p(\Delta_{\Theta} f)", yaxis = :log, ylims = [10^-5, 200])
