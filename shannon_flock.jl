###Introduction
#In this template, we generate an agent, but now it moves with some velocity and acceleration. However, the acceleration will be randomised at each time step and the velocity updated accordingly
#One thing that we're looking out for in particular is if we can update position, velocity and acceleration, which will be N tuples, with vectors generated from our ODE or movement gradient function. Actually no, in accordance with RK4, or even Euler integration, acceleration isn't an quantity we need to keep track of.

###Preliminaries
using Agents
using Random
using VoronoiCells



###Defining a norm function, cause why not
function norm (v)
	sum_of_squares = 0.0
	for i in 1:length(v)
		sum_of_squares += (v[i])^2
	end
	
	return sqrt(sum_of_squares)
end



###Function that calculates the voronoi area of a given agent a given position
function voronoi_area(pts)
        Cells = voronoicells(pts)
        Areas = voronoiarea[Cells]
        return  Areas
end



###Function that determines the gradient of movement
function move_gradient(agent, model,  kn, q, m)
	#Calculate the unit vector in the current direction of motion
	unit_v = agent.vel ./ norm(agent.vel)
	agent_speed = norm(agent.vel)
	vix = unit_v[1]
	viy = unit_v[2]
	positions = []
	neighbours = nearby_agents(agent, model)
	for neighbour in neighbours
		pushfirst!(positions, neighbour.position)	
	end		

	#Iterate through all the possible places the agent can move, keeping track of which one minimises area assuming static neighbour positions, though we make sure that if none of the moves optimises the current area, don't move at all
	pushfirst!(positions, agent.pos)
	min_area = agent.A #The agent's current DOD area
	deleteat!(positions, 1)
	min_direction = [0.0, 0.0]
	for i in 0:q #For every direction
		conflict = 0
		direction_of_move = [cos(i*pi/q)*vix - sin(i*pi/q)*viy, sin(i*pi/q)*vix + cos(i*pi/q)*viy]
		for j in 1:m #For every position up to m
			new_agent_pos = agent.pos .+ j .* direction_of_move .* agent_speed
		
			#Check first if there are no other agents in the potential position, note that we don't need to keep updating nearest neighbours since we assume the neighbours of a given agent are static
			for neighbour_position in positions
				if norm(new_agent_pos .- neighbour_position) < 2
					conflict = 1
					break
				end			
			end
			
			if (conflict == 1)		
				break
			end
		end
		
		if (conflict == 1) #Look at another direction if there's a agent in at or close to the potential position
                	continue
                end

		#If there are no other agents in the potential position (no conflicts), go ahead and evaluate the new DOD
                pushfirst!(positions, new_agent_pos)
                new_areas = voronoi_area(positions)
                new_agent_area = new_areas[1]
                if (new_area < min_area)
                	min_area = new_area
                        min_direction = direction_of_move
                end


	end

	#It really doesn't have to be like this, but at least just for the simple SHH model of Dr.Algar, we can simply return a velocity
	kn[1] += (direction_of_move .* agent_speed)[1]
	kn[2] += (direction_of_move .* agent_speed)[2]
end



###Create the agent
mutable struct bird <: AbstractAgent
	id::Int
	pos::NTuple{2, Float64}
	vel::NTuple{2, Float64}
	A::Float64 #The area of the agent's DOD
end
	


print("Agent template created")



###Create the initialisation function
using Random #for reproducibility
function initialise(; seed = 123, no_birds = 100)
	#Create the space
	space = ContinuousSpace((100.0, 100.0); periodic = true)
	
	#Create the properties of the model
	properties = Dict(:t => 0.0, :dt => 1.0)
	
	#Create the rng
	rng = Random.MersenneTwister(seed)
	
	print("Before model")

	#Create the model
	model = ABM(
		bird, space; 
		properties, rng, scheduler = Schedulers.fastest
	)	


	#Generate random initial positions for each bird, then calculate the DoDs
	initial_positions = []
	for i in 1:no_birds
		rand_position = Tuple(100*rand(Float64, 2))
		pushfirst!(initial_positions, rand_position)
	end

	#Calculate the DoDs based off the initial positions
	initial_dods = voronoi_area(initial_positions)
	
	#Now make the agents with their respective DoDs and add to the model
	for i in 1:no_birds
		agent = bird(i, initial_positions[i], Tuple(rand(Float64, 2)), initial_dods[i])
		add_agent!(agent, initial_positions[i], model)	
	end	

	return model
end  



###Create the agent step function. This will call upon the force or acceleration function. I'm assuming that this function will be applied to each agent in turn
function agent_step!(agent, model)		
	#Update the agent position and velocities, but only if it is a 
	#print("Step!\n", agent.planet)
	dt = model.dt
	k1 = [0.0, 0.0, 0.0, 0.0]
        #k2 = [0.0, 0.0, 0.0, 0.0]
        #k3 = [0.0, 0.0, 0.0, 0.0]
        #k4 = [0.0, 0.0, 0.0, 0.0]

	#Update the slopes, note that technically we should be using the vectorised dot operators, but Julia seems to allow us to be lazy when working with vectors
        #Now, why have we separated the position and velocity as two different vectors unlike PHYS4070? Because the pos is intrinsically a 2D vector for Julia Agents.
        move_gradient(agent, model, k1, 8, 1)
        #move_gradient(agent, model, agent.pos .+ dt/2*k1[1:2], agent.vel .+ dt/2*k1[3:4],  k2)
        #move_gradient(agent, model, agent.pos .+ dt/2*k2[1:2], agent.vel .+ dt/2*k2[3:4], k3)
        #move_gradient(agent, model, agent.pos .+ dt*k3[1:2], agent.vel .+ dt*k3[3:4], k4);
	
	#Update the agent position and velocity
	new_agent_pos = Tuple(agent.pos .+ dt .* k1[1:2])
        new_agent_vel = Tuple(agent.vel .+ dt .* k1[3:4])
	agent.vel = new_agent_vel
	print(new_agent_pos, "\n")
	#print(k1, "\n")
	#print(new_agent_pos, new_agent_vel, "\n")
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
save("moon_planet.png", figure)
	


###Animate
model = initialise();
period_time = 2*pi*19.0^1.5
dt = 0.001
num_moves = floor(Int64,period_time/dt)
abmvideo(
    "moon_planet.mp4", model, agent_step!, model_step!;
    framerate = 4, frames = 500,
    title = "Moon orbiting planet"
)
	

