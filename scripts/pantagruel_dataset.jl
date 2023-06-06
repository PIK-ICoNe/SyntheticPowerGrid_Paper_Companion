using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using CSV
using DataFrames, Plots, Statistics
using LaTeXStrings, StatsPlots
using Colors
using Statistics
default(grid = false, foreground_color_legend = nothing, bar_edges = false, lw = 3, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 12, tickfontsize = 10)

##
# Analyzing the ELMOD-De Dataset
path = string(@__DIR__) * "/../data/pantagruel/buses.csv"
file = CSV.File(open(path))