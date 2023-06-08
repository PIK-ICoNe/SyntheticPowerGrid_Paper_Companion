using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using PowerDynamics
using OrdinaryDiffEq
using PowerGridNoise  
using Interpolations
using Statistics
using DelimitedFiles

# This file generates the frequency time series for the 
## Loading a synthetic Power Grid consisting of droop controlled inverters
file_path = joinpath(@__DIR__, "../data/powergrids/synthetic_power_grid_example.json")
pg = read_powergrid(file_path, Json) 
op = find_operationpoint(pg)
ω_nodes, nodes, fluc_node_idxs, P_set, Q_set, f_idx = nodal_data(pg) # Accessing the node data from the grid

## 
tspan = (0.0, 5000.0)
step_size = 0.01
t_stops = collect(tspan[1]:step_size:tspan[end])
Δt = 10000.0

## Correlated
P_fluc, t = load_profile_model(tspan, Δt = Δt)
P_mean = mean(P_fluc)
P_fluc = (P_fluc .- P_mean) ./ P_mean
P_fluc_inter = linear_interpolation(t, P_fluc) # Interpolate the time series

fluctuations = map(f -> FluctuationNode(t -> P_set[f] + P_fluc_inter(t), t -> Q_set[f]), 1:length(fluc_node_idxs))
pg = generate_powergrid_fluctuations(pg, fluc_node_idxs, fluctuations)

ode = ODEProblem(rhs(pg), op.vec, tspan)
sol = solve(ode, Rodas4())
pg_sol = PowerGridSolution(sol, pg)

f = sol(t_stops, idxs = f_idx).u./(2π)
p = pg_sol(t_stops, fluc_node_idxs, :p)
rocof = sol(t_stops, Val{1}, idxs = f_idx).u./(2π)
writedlm("data/demand_fluctuations/corr/frequencies.txt", f)
writedlm("data/demand_fluctuations/corr/powers.txt", transpose(p))
writedlm("data/demand_fluctuations/corr/rocof.txt", rocof)

## Uncorrelated
flucs = [load_profile_model(tspan, Δt = Δt) for x in 1:length(fluc_node_idxs)]
t = map(f -> flucs[f][2], 1:length(fluc_node_idxs))
P_mean = map(f -> mean(flucs[f][1]), 1:length(fluc_node_idxs))
P_flucs = map(f -> (flucs[f][1] .- P_mean[f]) ./ P_mean[f], 1:length(fluc_node_idxs))

P_fluc_inter = map(f -> linear_interpolation(t[f], P_flucs[f]), 1:length(fluc_node_idxs)) # Interpolate the time series
fluctuations = map(f -> FluctuationNode(t -> P_set[f] + P_fluc_inter[f](t), t -> Q_set[f]), 1:length(fluc_node_idxs))

pg = generate_powergrid_fluctuations(pg, fluc_node_idxs, fluctuations)

ode = ODEProblem(rhs(pg), op.vec, tspan)
sol = solve(ode, Rodas4())
pg_sol = PowerGridSolution(sol, pg)

f = sol(t_stops, idxs = f_idx).u./(2π)
p = pg_sol(t_stops, fluc_node_idxs, :p)
rocof = sol(t_stops, Val{1}, idxs = f_idx).u./(2π)
writedlm("data/demand_fluctuations/uncorr/frequencies.txt", f)
writedlm("data/demand_fluctuations/uncorr/powers.txt", transpose(p))
writedlm("data/demand_fluctuations/uncorr/rocof.txt", rocof)
