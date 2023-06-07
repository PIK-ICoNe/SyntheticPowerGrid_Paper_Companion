module SyntheticPowerGridsPaper
    using SyntheticPowerGrids
    using SyntheticNetworks
    using XLSX
    using PowerDynamics
    using DataFrames
    using Statistics
    using Plots
    using Colors
    using LaTeXStrings
    using GraphMakie
    using CairoMakie
    
    include("lines.jl")
    export get_mean_line_length

    include("pg_fluctuations.jl")
    export calculate_performance_measures, generate_powergrid_fluctuations

    include("plotting.jl")
    export my_graph_plot

    include("utils.jl")
    export nodal_data
end
