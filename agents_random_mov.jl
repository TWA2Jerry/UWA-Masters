###Introduction
#In this template, we simply generate an agent on the model in continuous space, and move it to arbitrary positions each step
#This template has been shown to work, see random_mov.mp4

###Preliminaries
using Agents

###Create the agent
mutable struct random_agent <: AbstractAgent
	id::Int
	pos::NTuple{2, Float64}
end
	
print("Agent template created")

###Create the initialisation function
using Random #for reproducibility
function initialise(; seed = 123)
	#Create the space
	space = ContinuousSpace((100.0, 100.0), periodic=true)
	
	#Create the properties of the model
	properties = Dict(:t => 0.0, :dt => 0.01)
	
	#Create the rng
	rng = Random.MersenneTwister(seed)
	
	print("Before model")

	#Create the model
	model = ABM(
		random_agent, space; 
		properties, rng, scheduler = Schedulers.fastest
	)	

	#Populate the model
	agent1 = random_agent(1, (1.0, 1.0))
	add_agent!(agent1, (1.0, 1.0), model)
	return model
end  



###Create the agent step function. This will call upon the force or acceleration function. I'm assuming that this function will be applied to each agent in turn
function agent_step!(agent, model)		
	#Update the agent position and velocities
	move_agent!(agent, model)	
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
save("random_mov.png", figure)

###Animate
model = initialize();
abmvideo(
    "random_mov.mp4", model, agent_step!, model_step!;
    framerate = 4, frames = 20,
    title = "Random movement"
)
	
