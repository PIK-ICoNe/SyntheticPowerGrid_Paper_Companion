using SyntheticPowerGrids
#using SyntheticNetworks
#using Graphs
##
# In this example we will use the probabilistic capacity expansion algorithm that has been introduced in the paper. 
# For this we will generate our own topology. 

num_nodes = 100

#own_graph = generate_graph(RandomPowerGrid(num_nodes, [1, 1/5, 3/10, 1/3, 1/10, 0.0]...)) # Generate embedded graph
#e = edges(own_graph.graph)
#cables_vec = 3 * ones(Int, length(e))

#L = get_effective_distances(own_graph; mean_len_km = 37.12856121212121, shortest_line_km = 0.06)
#Y, Y_shunt = get_line_admittance_matrix(own_graph; L_matrix = L, cables_vec = cables_vec, num_nodes = num_nodes) 

# We use again nodal and edge parameters to hand over the admittances:
#edge_parameters = Dict(:Y => Y, :Y_shunt => Y_shunt) 
nodal_parameters = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => 0.5)
nodal_dynamics =  [(1/2, get_DroopControlledInverterApprox, nodal_parameters), (1/2, get_PQ, nothing)]

##
# Here we trigger the algorithm by using `probabilistic_capacity_expansion = true`.
pg_struct = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics, probabilistic_capacity_expansion = true, P0 = 2) #, validators = false, embedded_graph = own_graph, edge_parameters = edge_parameters, cables_vec = cables_vec, P0 = 2)
pg, op, pg_struct_updated, rejections = generate_powergrid_dynamics(pg_struct)

# Here we can see that the number of cables in the transmission lines has not been increased!
unique(pg_struct_updated.cables_vec)