###Introduction
#In this template, we generate an agent, but now it moves with some velocity and acceleration. However, the acceleration will be randomised at each time step and the velocity updated accordingly
#One thing that we're looking out for in particular is if we can update position, velocity and acceleration, which will be N tuples, with vectors generated from our ODE or movement gradient function. Actually no, in accordance with RK4, or even Euler integration, acceleration isn't an quantity we need to keep track of.

###Preliminaries
using Agents
using Random
###Create the movement gradient function that gives the rate of change in position and velocity, I'll use RK4, cause why not. Actually, no RK4 since we don't have a consistent acceleration function
function move_gradient(agent, k_pos, k_vel, kn)
	print("move_gradient called\n") 
	x = k_pos[1]
         y = k_pos[2]
         vx = k_vel[1]
         vy = k_vel[2]
	
	#Generate the acceleration, and set the value of kn
	kn = (vx, vy, rand(-5.0:0.001:5.0), rand(-5.0:0.001:5.0))
	print("kn set to a value of ", kn, "\n")
	return kn
end

###Create the agent
mutable struct random_agent <: AbstractAgent
	id::Int
	pos::NTuple{2, Float64}
	vel::NTuple{2, Float64}
	group::Int
end
	
print("Agent template created")

###Create the initialisation function
using Random #for reproducibility
function initialise(; seed = 123)
	#Create the space
	space = ContinuousSpace((100.0, 100.0); periodic = true)
	
	#Create the properties of the model
	properties = Dict(:t => 0.0, :dt => 0.1)
	
	#Create the rng
	rng = Random.MersenneTwister(seed)
	
	print("Before model")

	#Create the model
	model = ABM(
		random_agent, space; 
		properties, rng, scheduler = Schedulers.fastest
	)	

	#Populate the model
	agent1 = random_agent(1, (1.0, 1.0), (-1.0, -1.0), 1)
	agent2 = random_agent(2, (rand(-5.0:0.001:5.0), rand(-5.0:0.001:5.0)), (rand(-5.0:0.001:5.0), rand(-5.0:0.001:5.0)), 2)	

	add_agent!(agent1, model)
	add_agent!(agent2, model)
	return model
end  



###Create the agent step function. This will call upon the force or acceleration function. I'm assuming that this function will be applied to each agent in turn
function agent_step!(agent, model)		
	#Update the agent position and velocities
	if(agent.group == 1)
	print("Step!\n", agent.group)
	dt = model.dt
	k1 = (0.0, 0.0, 0.0, 0.0)
        k2 = (0.0, 0.0, 0.0, 0.0)
        k3 = (0.0, 0.0, 0.0, 0.0)
        k4 = (0.0, 0.0, 0.0, 0.0)

	#Update the slopes, note that technically we should be using the vectorised dot operators, but Julia seems to allow us to be lazy when working with vectors
        #Now, why have we separated the position and velocity as two different vectors unlike PHYS4070? Because the pos is intrinsically a 2D vector for Julia Agents.
        k1 = move_gradient(agent, agent.pos, agent.vel, k1)
        #move_gradient(agent, model, agent.pos .+ dt/2*k1[1:2], agent.vel .+ dt/2*k1[3:4],  k2)
        #move_gradient(agent, model, agent.pos .+ dt/2*k2[1:2], agent.vel .+ dt/2*k2[3:4], k3)
        #move_gradient(agent, model, agent.pos .+ dt*k3[1:2], agent.vel .+ dt*k3[3:4], k4);
	
	#Update the agent position and velocity
	new_agent_pos = agent.pos .+ dt .* (k1[1:2])
        new_agent_vel = agent.vel .+ dt .* (k1[3:4])
	agent.vel = new_agent_vel
	print(k1, "\n")
	print(new_agent_pos, new_agent_vel, "\n")
	move_agent!(agent, new_agent_pos, model)	
	end
end


###Create the model_step function
function model_step!(model)
	model.t += model.dt
end
	
###Initialise the model
model = initialise()


###Test the model has been initialised and works
using InteractiveDynamics
using CairoMakie # choosing a plotting backend

figure, _ = abmplot(model)
figure # returning the figure displays it
save("two_group.png", figure)

###Animate
model = initialise();
abmvideo(
    "two_group.mp4", model, agent_step!, model_step!;
    framerate = 4, frames = 20,
    title = "Random movement"
)
	
