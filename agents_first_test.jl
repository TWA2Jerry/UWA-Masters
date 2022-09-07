###Introduction
So this is just our first test of using the agents.jl thang. What we'll try to do is, we have create a flock, with a single leader. The leader is still counted as part of the flock, but doesn't care about COM and just moves around randomly.Other members of the flock try to minimise their distance to the COM of the flock. The idea is that the leader will have greater weighting in teh COM calculation as well, and so what should hopefully happen is the flock will move as one, following the leader. 

#Create the agent
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
	#Create the model
	rng = 
	#Populate the model


#Create the agent and model step functions

###Create the step function. This will call upon the force or acceleration function
#We will use RK4 to simulate each step
delta_t = 0.01
function agent_step()
	


