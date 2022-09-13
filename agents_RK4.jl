###Introduction
So this is just our first test of using the agents.jl thang. What we'll try to do is, we have create a flock, with a single leader. The leader is still counted as part of the flock, but doesn't care about COM and just moves around randomly.Other members of the flock try to minimise their distance to the COM of the flock. The idea is that the leader will have greater weighting in teh COM calculation as well, and so what should hopefully happen is the flock will move as one, following the leader. 


###Create the ODE function that determines the movement of the agent, i.e, account for the force and EOM
#Note that we only need to specify kn, since all arrays are passed by pointer in Julia 
function move_gradient(agent, model, k_pos, k_vel, kn)
	#In general, look through all possible agents
	r = 2*sqrt((spacesize(model)[1])^2 + (spacesize(model)[2])^2) #Calculate the maximum possible distance any neighbour could be. Note that this is for continuous space. For other stuff like grids, use manhattan or chebyshev 
	neighbours = nearby_ids(agent, model, r) #neigbours is now an iterable of all other agents j\neq i
	
	#Iterate and add the contribution to acceleration of all neighbouring agents
	 x = k_pos[1]
         y = k_pos[2]
         vx = k_vel[1]
         vy = k_vel[2]
         kn[1] = vx
         kn[2] = vy

	for neighbour in neighbours
		#Generate the set of equations that determines the change in position and movement
		x_j = neighbour.pos[1]
		y_j = neighbour.pos[2]
        	kn[3] += -(x-x_j)/((x-x_j)^2+(y-y_j)^2)^1.5
        	kn[4] += -(y-y_j)/((x-x_j)^2+(y-y_j)^2)^1.5
	end
end



###Create the agent
mutable subtype agent <: AgentType
	id::Int
	pos::NTuple{2, Float64}
	vel::NTuple{2, Float64}
	acc::NTuple{2, Float64}
	shepherd::Bool #So the idea is that for the shepherd, they won't care, and they just move around freely
end
	


###Create the initialisation function
using Random #for reproducibility
function initialise()
	#Create the space
	space = ContinuousSpace(extent::NTuple{2, <:Real})
	#Create the properties of the model
	
	#Create the rng
	rng = 
	#Create the model
	
	#Populate the model
	


###Create the agent step function. This will call upon the force or acceleration function. I'm assuming that this function will be applied to each agent in turn
#We will use RK4 to simulate each step
t = 0.0 #Actually, pretty sure that we should set t and dt parameters in the model initialisation properties
dt = 0.001
function agent_step!(agent, model)
	#As according to RK4, generate the 4 slopes of movement, kn. Same as in UQ/PHYS4070, the kn[1] and kn[2] (1 indexed for Julia) represent changes to the position x and y, while kn[3] and kn[4] represent changes to the velocity x and y components  
	k1 = [0.0, 0.0, 0.0, 0.0]
	k2 = [0.0, 0.0, 0.0, 0.0]
	k3 = [0.0, 0.0, 0.0, 0.0]
	k4 = [0.0, 0.0, 0.0, 0.0]
	
	#Update the slopes, note that technically we should be using the vectorised dot operators, but Julia seems to allow us to be lazy when working with vectors
	#Now, why have we separated the position and velocity as two different vectors unlike PHYS4070? Because the pos is intrinsically a 2D vector for Julia Agents. 
	move_gradient(agent, model, agent.pos, agent.vel, k1)
	move_gradient(agent, model, agent.pos .+ dt/2*k1[1:2], agent.vel .+ dt/2*k1[3:4],  k2)
	move_gradient(agent, model, agent.pos .+ dt/2*k2[1:2], agent.vel .+ dt/2*k2[3:4], k3)
	move_gradient(agent, model, agent.pos .+ dt*k3[1:2], agent.vel .+ dt*k3[3:4], k4);
	
	
	#Calculate the new agent positions and velocities based off RK4
	new_agent_pos = Tuple(agent.pos .+ dt/6*(k1[1:2] .+ 2*k2[1:2] .+ 2*k3[1:2] .+ k4[1:2]))
	new_agent_vel = Tuple(agent.vel .+ dt/6*(k1[3:4] .+ 2*k2[3:4] .+ 2*k3[3:4] .+ k4[3:4]))
	

	#Update the agent position and velocities
	move_agent!(agent, new_agent_pos, model)
	agent.vel = new_agent_vel
	 	
end
	
###Initialise the model
model = initialise()


###Test the model has been initialised and works


###Animate

	
