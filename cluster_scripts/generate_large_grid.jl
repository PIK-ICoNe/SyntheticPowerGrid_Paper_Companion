using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using SyntheticPowerGrids
using PowerDynamics
#import SyntheticPowerGrids.validate_power_flow_on_lines

##
nodal_parameters = Dict(:η => 3 * 10^(-3), :α => 5.0,  :κ => π/2)
nodal_dynamics = [(1.0, get_dVOCapprox, nodal_parameters)]

num_nodes = 10

pg_struct = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_struct)

## Saving the power grid for later usage
file_path = joinpath(@__DIR__, "../data/powergrids/synthetic_power_grid_large.json")
write_powergrid(pg, file_path, Json)

##
#p_nodes = map(n -> pg.nodes[n].P, 1:num_nodes)

#histogram(p_nodes, bins = 20)

##
#_, _, P, P_max = validate_power_flow_on_lines(op, pg_struct)

##
#P = vcat(values(P)...)
#P_max = vcat(values(P_max)...)

#line_loading = (P ./ P_max) .* 100

#histogram(line_loading)