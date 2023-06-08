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

## Analyzing the pantagruel dataset
## Loads
path = string(@__DIR__) * "/../data/pantagruel/buses.csv"
df_buses = DataFrame(CSV.File(open(path)))

df_buses_380kV = filter(row -> row."Voltage [KV]" == 380, df_buses) 

demand_MW = df_buses_380kV[:, "Demand [MW]"]

plt = histogram(demand_MW, xaxis = L"Demand \, \, [MW]", legend = false, normalize = :probability)
savefig(plt, "plots/pantagruel_demand.pdf")

## Lines
nodes_380kV = df_buses_380kV.ID
path = string(@__DIR__) * "/../data/pantagruel/lines.csv"
df_lines = DataFrame(CSV.File(open(path)))

line_id = collect(1:size(df_lines)[1])
df_lines[!, :line_id] = line_id 

##
lines_380kV = []

for node in nodes_380kV
    df = filter(row -> row.ID1 == node, df_lines)
    df = filter(row -> row.ID2 == node, df_lines)
    append!(lines_380kV, df[:, "line_id"])
end
lines_380kV = unique(lines_380kV)

##
df_lines_380kv = df_lines[lines_380kV, :]

x_per_km = 265 * 10^-3 # Ω/km for 380kV level: https://storage.googleapis.com/plos-corpus-prod/10.1371/journal.pone.0213550/1/pone.0213550.s002.pdf?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=wombat-sa%40plos-prod.iam.gserviceaccount.com%2F20230606%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20230606T112614Z&X-Goog-Expires=86400&X-Goog-SignedHeaders=host&X-Goog-Signature=70d54da2103307b5e52316dda214e78ac2971b7b5d020012952f29df83e9d74a90e182f34332199ca8c1e10f6465d2b71229165e20f3a9b654fa4564c3be8dfc2188c142f96b68259be9a1e5c250f400958d7f5684c35718df3a4e9ec6076a112001db222344615020ee0fbca5320c27a84d27bf9a08e3e2ecf7f88fcdce8fb891a5b83d5ae09833878aefd97944c7b46472754ad8e1496b784061f1508aad1ba78cf64cc75116b459b8a054074af4abb4521f02f4efed286585f2d7a5824ccf8ddb82ff6a947666f444156282ee99e43b23fa2697b44574a7d55c50ec17b1bef9c455517852c795782b02da45d4de56c563034bc675f0ebbb2a0c86d2765b00

X_ohm = df_lines_380kv[:, "Reactance [Ohm]"]

lengths = X_ohm * 10^3 / x_per_km # I assume that their reactance is in m Ω and not in Ω. Otherwise I can not explain this
minimum(lengths)
mean(lengths) 
std(lengths)

plt = histogram(lengths, xaxis = L"l [km]", legend = false, normalize = :probability)
savefig(plt, "plots/pantagruel_lines.pdf")
