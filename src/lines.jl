
function get_mean_line_length()
    path = string(@__DIR__) * "/../data/sci_grids/"
    line_table = XLSX.readtable(joinpath(path, "links_de_power_151109.xlsx"), "links_de_power_151109")

    df_lines = DataFrame(line_table) # Read in Sci Grid Data Set
    df_380kV = filter(row -> row.voltage == 380000, df_lines) # We only need lines with a voltage of 380kV
    mean(df_380kV.length_m) / 1000 # Take the mean length, conversion to km
end