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
f_demand = readdlm("data/demand_fluctuations/uncorr/frequencies.txt")
f_wind = readdlm("data/wind_fluctuations/uncorr/frequencies.txt")
f_solar = readdlm("data/solar_fluctuations/uncorr/frequencies.txt")

## Statistics
mean(f_demand)
mean(f_wind)
mean(f_solar)

std(f_demand)
std(f_wind)
std(f_solar)

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

## R2
r2_demand = round(R_squared(pdf_demand.density, pdf_d_n), digits = 3)
r2_wind = round(R_squared(pdf_wind.density, pdf_w_n), digits = 3)
r2_solar = round(R_squared(pdf_solar.density, pdf_s_n), digits = 3)

## Plotting
x_lims = [-0.23, 0.23]

plot(pdf_solar.x, pdf_solar.density ./ sum(pdf_solar.density), c = colorant"goldenrod3", label = L"f_{solar}", yaxis = L"p(Δf)")
annotate!(-0.05, 0.002, text(L"R^2= " * latexstring(r2_solar) , :black, :right, 15))
plt_pdf_s = plot!(pdf_solar.x, pdf_s_n ./ sum(pdf_s_n), label = "", c = colorant"goldenrod3", ls = :dashdot, alpha = 0.7, xlims = x_lims)

plot(pdf_wind.x, pdf_wind.density ./ sum(pdf_wind.density), label = L"f_{wind}", c = colorant"teal")
annotate!(-0.05, 0.0014, text(L"R^2= " * latexstring(r2_wind) , :black, :right, 15))
plt_pdf_w = plot!(pdf_wind.x, pdf_w_n ./ sum(pdf_w_n), label = "", c = colorant"teal", ls = :dashdot, alpha = 0.7, xlims = x_lims)

plot(pdf_demand.x, pdf_demand.density ./ sum(pdf_demand.density), label = L"f_{demand}", c = colorant"coral")
annotate!(-0.05, 0.0013, text(L"R^2= " * latexstring(r2_demand) , :black, :right, 15))
plt_pdf_d = plot!(pdf_demand.x, pdf_d_n ./ sum(pdf_d_n), label = "", c = colorant"coral", ls = :dashdot, alpha = 0.7, xlims = x_lims)

plt_pdf = Plots.plot(plt_pdf_w, plt_pdf_d, plt_pdf_s; layout=(1, 3), size=(1500, 400), xaxis = L"\Delta f [Hz]", left_margin=10mm, bottom_margin = 10mm)
savefig(plt_pdf, "plots/pdf_frequency.pdf")
