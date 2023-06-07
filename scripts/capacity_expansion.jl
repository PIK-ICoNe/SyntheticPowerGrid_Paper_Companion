using SyntheticPowerGrids

##
# In this example we will use the probabilistic capacity expansion algorithm that has been introduced in the paper. 
num_nodes = 100

nodal_parameters = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => 0.5)
nodal_dynamics =  [(1/2, get_DroopControlledInverterApprox, nodal_parameters), (1/2, get_PQ, nothing)]

## Here we trigger the algorithm by using `probabilistic_capacity_expansion = true`.
pg_struct = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics, probabilistic_capacity_expansion = true, P0 = 2.6) 

pg, op, pg_struct_updated, rejections = generate_powergrid_dynamics(pg_struct)

# Here we can see that the number of cables in the transmission lines has not been increased!
unique(pg_struct_updated.cables_vec)