module SyntheticPowerGridsPaper
    using SyntheticPowerGrids
    using SyntheticNetworks
    using XLSX
    using DataFrames
    using Statistics
    
    include("lines.jl")
    export get_mean_line_length
end
