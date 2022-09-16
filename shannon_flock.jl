###Introduction
#In this template, we generate an agent, but now it moves with some velocity and acceleration. However, the acceleration will be randomised at each time step and the velocity updated accordingly
#One thing that we're looking out for in particular is if we can update position, velocity and acceleration, which will be N tuples, with vectors generated from our ODE or movement gradient function. Actually no, in accordance with RK4, or even Euler integration, acceleration isn't an quantity we need to keep track of.

###Preliminaries
using Agents
using Random
using VoronoiCells

###Function that calculates the voronoi area of a given agent a given position
function voronoi_area(pts)
        Cells = voronoicells(pts)
        Areas = voronoiarea[Cells]
        return  Areas[1]
end

###Function that determines the gradient of movement
function move_gradient(agent, model, agent_speed)
	#Calculate the unit vector in the current direction of motion
	unit_v = agent.vel ./ norm(agent.vel)
	vix = unit_v[1]
	viy = unit_v[2]
	positions = []
	neighbours = nearby_agents(agent, model)
	for neighbour in neighbours
		pushfirst!(positions, neighbour.position)	
	end		

	#Iterate through all the possible places the agent can move, keeping track of which one minimises area assuming static neighbour positions
	min_area = 100000000.0
	min_direction = [0.0]
	for i in range 0:7
		direction_of_move = [cos(i*pi/8)*vix - sin(i*pi/8)*viy, sin(i*pi/8)*vix + cos(i*pi/8)*viy]
		new_agent_pos = agent.pos .+ direction_of_move .* agent_speed
		
		#Check first if there are no other agents in the potential position
		conflict = 0
		for neighbour_position in positions
			if norm(new_agent_pos .- neighbour_position) > 2
				conflict = 1
				break
			end			
		end
		
		if conflict = 1
			continue
		end
		
		#If there are no other agents in the potential position, go ahead and evaluate the new DOD
		pushfirst!(positions, new_agent_pos)
		new_areas = voronoi_area(positions)
		new_agent_area = new_areas[1]
		if new_area < min_area
			min_area = new_area
			min_direction = direction_of_move
		end	
	end

	return Tuple(agent.pos .+ mind_direction .* agent_speed .* model.dt)
end


###Create the agent
mutable struct celestial_object <: AbstractAgent
	id::Int
	pos::NTuple{2, Float64}
	vel::NTuple{2, Float64}
	planet::Bool
end
	
print("Agent template created")

###Create the initialisation function
using Random #for reproducibility
function initialise(; seed = 123)
	#Create the space
	space = ContinuousSpace((100.0, 100.0); periodic = true)
	
	#Create the properties of the model
	properties = Dict(:t => 0.0, :dt => 0.01)
	
	#Create the rng
	rng = Random.MersenneTwister(seed)
	
	print("Before model")

	#Create the model
	model = ABM(
		celestial_object, space; 
		properties, rng, scheduler = Schedulers.fastest
	)	

	#Calculate, in planetary units such that m_planet = 1, and G = 1, the velocity required for the moon to achieve a circular orbit around the planet 
	r = 19.0
	moon_speed = sqrt(1.0/r)
	moon_period = 2*pi*r^1.5
	vi_x = 0.0
	vi_y = moon_speed

	#Populate the model
	agent1 = celestial_object(1, (69.0, 50.0), (vi_x, vi_y), 0) #The moon
	agent2 = celestial_object(2, (50.0, 50.0), (0.0, 0.0), 1) #The planet	

	add_agent!(agent1, (69.0, 50.0), model)
	add_agent!(agent2, (50.0, 50.0), model)
	return model
end  



###Create the agent step function. This will call upon the force or acceleration function. I'm assuming that this function will be applied to each agent in turn
function agent_step!(agent, model)		
	#Update the agent position and velocities, but only if it is a 
	if(agent.planet == 0)
	#print("Step!\n", agent.planet)
	dt = model.dt
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
	
	#Update the agent position and velocity
	new_agent_pos = Tuple(agent.pos .+ dt/6*(k1[1:2] .+ 2*k2[1:2] .+ 2*k3[1:2] .+ k4[1:2]))
        new_agent_vel = Tuple(agent.vel .+ dt/6*(k1[3:4] .+ 2*k2[3:4] .+ 2*k3[3:4] .+ k4[3:4]))
	agent.vel = new_agent_vel
	print(new_agent_pos, "\n")
	#print(k1, "\n")
	#print(new_agent_pos, new_agent_vel, "\n")
	move_agent!(agent, new_agent_pos, model)	
	end
	
	#If the agent is the planet, do nothing. This is just to check that we are indeed, doing nothing
	#if(agent.planet == 1) 
	#print("Agent 2 position is ", agent.pos, "\n")
	#end
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
	

