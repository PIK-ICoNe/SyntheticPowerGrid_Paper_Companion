using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using SyntheticPowerGrids
using PowerDynamics
using OrdinaryDiffEq
using PowerGridNoise
using Interpolations
using Statistics
using DelimitedFiles

## Loading a synthetic Power Grid consisting of droop controlled inverters
file_path = joinpath(@__DIR__, "../data/powergrids/synthetic_power_grid_example.json")
pg = read_powergrid(file_path, Json) 
op = find_operationpoint(pg)
ω_nodes, nodes, fluc_node_idxs, P_set, Q_set, f_idx= nodal_data(pg) # Accessing the node data from the grid

## Using an intermittent wind power fluctuation Langevin-type model to generate fluctuating time series
tspan = (0.0, 5000.0)
step_size = 0.01
t_stops = collect(tspan[1]:step_size:tspan[end])

Δt = 10000.0
D = 0.1 # Intermittence strength
x, t = wind_power_model(tspan, D = D, Δt = Δt)
x_inter = linear_interpolation(t, x) # Interpolate the time series

## Multi Node Fluctuations, completely correlated, exchange all PQAlgebraic with FluctuationNode
fluctuations = map(f -> FluctuationNode(t -> P_set[f] + x_inter(t), t -> Q_set[f]), 1:length(fluc_node_idxs))
pg = generate_powergrid_fluctuations(pg, fluc_node_idxs, fluctuations)

ode = ODEProblem(rhs(pg), op.vec, tspan)
sol = solve(ode, Rodas4())
pg_sol = PowerGridSolution(sol, pg)

f = sol(t_stops, idxs = f_idx).u./(2π)
p = pg_sol(t_stops, fluc_node_idxs, :p)
rocof = sol(t_stops, Val{1}, idxs = f_idx).u./(2π)

writedlm("data/wind_fluctuations/corr/frequencies.txt", f)
writedlm("data/wind_fluctuations/corr/powers.txt", transpose(p))
writedlm("data/wind_fluctuations/corr/rocof.txt", rocof)

## Multi Node Fluctuations , completely uncorrelated, exchange all PQAlgebraic with FluctuationNode / Generate a time series for each node
flucs = [wind_power_model(tspan, D = D, Δt = Δt) for x in 1:length(fluc_node_idxs)]
x_inter = map(f -> linear_interpolation(flucs[f][2], flucs[f][1]), 1:length(fluc_node_idxs)) # Interpolate the time series
fluctuations = map(f -> FluctuationNode(t -> P_set[f] + x_inter[f](t), t -> Q_set[f]), 1:length(fluc_node_idxs))
pg = generate_powergrid_fluctuations(pg, fluc_node_idxs, fluctuations)

ode = ODEProblem(rhs(pg), op.vec, tspan)
sol = solve(ode, Rodas4())
pg_sol = PowerGridSolution(sol, pg)

f = sol(t_stops, idxs = f_idx).u./(2π)
p = pg_sol(t_stops, fluc_node_idxs, :p)
rocof = sol(t_stops, Val{1}, idxs = f_idx).u./(2π)

writedlm("data/wind_fluctuations/uncorr/frequencies.txt", f)
writedlm("data/wind_fluctuations/uncorr/powers.txt", transpose(p))
writedlm("data/wind_fluctuations/uncorr/rocof.txt", rocof)
