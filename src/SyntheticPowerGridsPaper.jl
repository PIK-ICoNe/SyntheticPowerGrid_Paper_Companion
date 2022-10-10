module SyntheticPowerGridsPaper
    using SyntheticPowerGrids
    using SyntheticNetworks
    using XLSX
    using PowerDynamics
    using DataFrames
    using Statistics
    
    include("lines.jl")
    export get_mean_line_length

    include("pg_fluctuations.jl")
    export calculate_performance_measures, generate_powergrid_fluctuations
end
