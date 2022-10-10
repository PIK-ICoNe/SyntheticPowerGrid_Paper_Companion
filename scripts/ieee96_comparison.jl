using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using SyntheticPowerGrids
using PowerDynamics
using OrdinaryDiffEq
using Plots
using PowerGridNoise
using Interpolations
using Statistics
using LaTeXStrings
default(grid = false, foreground_color_legend = nothing, bar_edges = false,  lw=1.5, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 15, tickfontsize = 10)

##
# Define your own function for a node type
# Here we use simply the swing equation!

function get_swingLVS(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters::Dict)
    H = nodal_parameters[:H] # Inertia Constant
    Ω = nodal_parameters[:Ω] # Rated Frequency
    D = nodal_parameters[:D] # Damping Coefficient
    Γ = nodal_parameters[:Γ] # Voltage stability Coefficient

    SwingEqLVS(H = H, P = P_set, D = D, Ω = Ω, Γ = Γ, V = V_set)
end

num_nodes = 73
parameters_swing = Dict(:D => 2.0, :H => 3, :Γ => 10, :Ω => 2π * 50)
nodal_dynamics = [(33/num_nodes, get_swingLVS, parameters_swing), (1-33/num_nodes, get_PQ, nothing)]

a = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics, lines = :StaticLine)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(a)

##
# Single node fluctuation
# Accessing the node data from the grid
ω_indices = findall(n -> :ω ∈ symbolsof(n), pg.nodes)
nodes = deepcopy(pg.nodes) 
pq_node_idxs = findall(typeof.(pg.nodes) .== PQAlgebraic) # Find all Load Buses in the grid
P_set = map(i -> nodes[i].P, pq_node_idxs) # Load their power set-points
Q_set = map(i -> nodes[i].Q, pq_node_idxs)

##
# Using an intermittent wind power fluctuation Langevin-type model to generate fluctuating time series
tspan = (0.0, 100.0)
D = 0.1 # Intermittence strength
p = 0.2 # Penetration parameter
x, t = wind_power_model(tspan, D = D)
x_inter = linear_interpolation(t, x) # Interpolate the time series

##
fluctuation = FluctuationNode(t -> P_set[1] + p * x_inter(t), t -> Q_set[1])
pg_fluc = generate_powergrid_fluctuations(pg, [pq_node_idxs[1]], [fluctuation])

##
# Simulate a trajectory
ode = ODEProblem(rhs(pg_fluc), op.vec, tspan)
sol = solve(ode, Rodas4())

solution = PowerGridSolution(sol, pg_fluc)
plot(solution, pq_node_idxs, label = "Active Power",:p, lw = 3, ylabel = L"P[p.u.]", xlabel = L"t[s]", legend = false)

plt2 = plot(solution, ω_indices, :ω, legend = false, ylabel = L"ω[rad / s]", xlabel = L"t[s]")
savefig(plt2, "plots/single_node_fluc_73_nodes.pdf")

calculate_performance_measures(solution, :ω) # calculate performance measures
