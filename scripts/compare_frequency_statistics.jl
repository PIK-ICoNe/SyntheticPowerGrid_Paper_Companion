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
using Colors
using Plots.Measures

default(grid = false, foreground_color_legend = nothing, bar_edges = false,  lw=3, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 18, labelfontsize = 18, tickfontsize = 15)

##
step_size = 0.01

f_demand = readdlm("data/demand_fluctuations/uncorr/frequencies.txt")
f_wind = readdlm("data/wind_fluctuations/uncorr/frequencies.txt")
f_solar = readdlm("data/solar_fluctuations/uncorr/frequencies.txt")

f_d_mean = mean(f_demand, dims = 2)
f_w_mean = mean(f_wind, dims = 2)
f_s_mean = mean(f_solar, dims = 2)

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

plot(dts, f_s_ac, ylabel = L"c(\Delta t)", xlabel = L"\Delta t [min]", label = L"f_s", c = colorant"goldenrod2")
plot!(dts, f_w_ac, label = L"f_w", c = colorant"teal", ls = :dashdot)
plt_ac = plot!(dts, f_d_ac, label = L"f_d", c = colorant"coral", ls = :dot)
savefig(plt_ac, "plots/autocorrelation.pdf")

## Probability Density Function
pdf_demand = kde(vcat(f_demand...))

d_d_n = fit(Normal{Float64}, vcat(f_demand...))
pdf_d_n = map(i -> pdf(d_d_n, pdf_demand.x[i]), 1:length(pdf_demand.x))

pdf_wind = kde(vcat(f_wind...))
d_w_n = fit(Normal{Float64}, vcat(f_wind...))
pdf_w_n = map(i -> pdf(d_w_n, pdf_wind.x[i]), 1:length(pdf_wind.x))

pdf_solar = kde(vcat(f_solar...))
d_s_n = fit(Normal{Float64}, vcat(f_solar...))
pdf_s_n = map(i -> pdf(d_s_n, pdf_solar.x[i]), 1:length(pdf_solar.x))

plot(pdf_solar.x, pdf_solar.density ./ findmax(pdf_solar.density)[1], c = colorant"goldenrod3", label = L"f_s", yaxis = L"PDF")
plt_pdf_s = plot!(pdf_solar.x, pdf_s_n ./ findmax(pdf_s_n)[1], label = "", c = colorant"goldenrod3", ls = :dashdot, alpha = 0.7)

plot(pdf_wind.x, pdf_wind.density ./ findmax(pdf_wind.density)[1], label = L"f_w", c = colorant"teal")
plt_pdf_w = plot!(pdf_wind.x, pdf_w_n ./ findmax(pdf_w_n)[1], label = "", c = colorant"teal", ls = :dashdot, alpha = 0.7)

plot(pdf_demand.x, pdf_demand.density ./ findmax(pdf_demand.density)[1], label = L"f_d", c = colorant"coral")
plt_pdf_d = plot!(pdf_demand.x, pdf_d_n ./ findmax(pdf_d_n)[1], label = "", c = colorant"coral", ls = :dashdot, alpha = 0.7)

plt_pdf = Plots.plot(plt_pdf_s, plt_pdf_w, plt_pdf_d; layout=(1, 3), size=(1500, 500), xaxis = L"\Delta f [Hz]", left_margin=10mm, bottom_margin = 10mm)
savefig(plt_pdf, "plots/pdf_frequency.pdf")

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

plot(pdf_increments.x / σ_ΔΘf, pdf_normal, c = colorant"chocolate3", label = L"Normal", ls = :dash)
plt_increments = plot!(pdf_increments.x / σ_ΔΘf, pdf_increments.density, c = colorant"teal", label = L"f_w", xlabel = L"\Delta_{\Theta} f / \sigma_{\Delta_{\Theta} f}", ylabel = L"p(\Delta_{\Theta} f)", yaxis = :log, ylims = [10^-5, 200])

savefig(plt_increments, "plots/frequency_increments.pdf")
