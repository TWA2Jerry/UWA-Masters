###Introduction
#In this template, we generate an agent, but now it moves with some velocity and acceleration. However, the acceleration will be randomised at each time step and the velocity updated accordingly
#One thing that we're looking out for in particular is if we can update position, velocity and acceleration, which will be N tuples, with vectors generated from our ODE or movement gradient function. Actually no, in accordance with RK4, or even Euler integration, acceleration isn't an quantity we need to keep track of.

###Preliminaries
using Agents
using Random
out_file = open("moon_proj_dat_jul.txt", "w")

###Function for the norm of a vector, cause why not
function norm(v)
	sum_of_squares = 0.0
	for i in 1:length(v)
		sum_of_squares += v[i]^2
	end
	return sqrt(sum_of_squares)
end


###Create the movement gradient function that gives the rate of change in position and velocity, using RK4
function move_gradient(agent, model, k_pos, k_vel, kn, n) #Note that we have to specify k_pos, k_vel since it's not the actual agent pos and vel, but the agent pos and vel given the addition of a kn
	#Determine the neighbours of the agent. In this case, we want the moon, and we look through all of possible space
	r = sqrt((spacesize(model)[1])^2 + (spacesize(model)[2])^2) #Calculate the maximum possible distance any neighbour could be. Note that this is for continuous space. For other stuff like grids, use manhattan or chebyshev
        #neighbours = nearby_agents(agent, model, r) #neigbours is now an iterable of all other agents j\neq i
	allagent_iterable = allagents(model)
	x = k_pos[1]
        y = k_pos[2]
        vx = k_vel[1]
        vy = k_vel[2]
	kn[1] = vx
	kn[2] = vy
	#print("For kn = ", kn, "the number of neighbours detected is ", neighbours, "\n")
	no_neighbours_detected = 0
	for neighbour in allagent_iterable
        	#Generate the set of equations that determines the change in position and movement
		if(neighbour.id == agent.id)
			continue
		end
		no_neighbours_detected += 1
		#print("For n = ", n, "with the ", no_neighbours_detected, "-th neighbour", "Neighbour position at $(neighbour.pos)\n")
		x_j = neighbour.pos[1]
		y_j = neighbour.pos[2]
                kn[3] += -1.0*(x-x_j)/((x-x_j)^2+(y-y_j)^2)^1.5
                kn[4] += -1.0*(y-y_j)/((x-x_j)^2+(y-y_j)^2)^1.5
        end

	
	#= x_j = 50.0
	y_j = 50.0
	kn[3] = -1.0*(x-x_j)/((x-x_j)^2+(y-y_j)^2)^1.5
        kn[4] = -1.0*(y-y_j)/((x-x_j)^2+(y-y_j)^2)^1.5
	=#
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
	#space = ContinuousSpace((100.0, 100.0); periodic = true, spacing = )

	
	#Create the properties of the model
	properties = Dict(:t => 0.0, :dt => 0.001)
	
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
	moon_period = 2*pi*r^(1.5)
	vi_x = 0.0
	vi_y = moon_speed

	#Populate the model
	agent1 = celestial_object(1, (69.0, 50.0), (vi_x, vi_y), 0) #The moon
	agent2 = celestial_object(2, (50.0, 50.0), (0.0, 0.0), 1) #The planet	

	add_agent!(agent1, (69.0, 50.0), model)
	add_agent!(agent2, (50.0, 50.0), model)
	r_mp = norm(agent1.pos .- agent2.pos)
	write(out_file, "$(agent1.pos[1]) $(agent1.pos[2])) $((agent1.pos .- agent2.pos)[1]) $((agent1.pos .- agent2.pos)[2]) $r_mp\n")
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
        move_gradient(agent, model, agent.pos, agent.vel, k1, 1)
	#print("The value of k1 after the ODE function is ", k1, "\n")
	move_gradient(agent, model, agent.pos .+ dt/2*k1[1:2], agent.vel .+ dt/2*k1[3:4],  k2, 2)
        #print("The value of k2 after the ODE function is ", k2, "\n")
	move_gradient(agent, model, agent.pos .+ dt/2*k2[1:2], agent.vel .+ dt/2*k2[3:4], k3, 3)
        move_gradient(agent, model, agent.pos .+ dt*k3[1:2], agent.vel .+ dt*k3[3:4], k4, 4);
	
	#Update the agent position and velocity
	new_agent_pos = Tuple(agent.pos .+ dt/6 .*(k1[1:2] .+ 2*k2[1:2] .+ 2*k3[1:2] .+ k4[1:2]))
        new_agent_vel = Tuple(agent.vel .+ dt/6 .*(k1[3:4] .+ 2*k2[3:4] .+ 2*k3[3:4] .+ k4[3:4]))
	agent.vel = new_agent_vel
	#print(new_agent_pos, "\n")
	#print(k1, "\n")
	#print(new_agent_pos, new_agent_vel, "\n")
	r_mp = norm(new_agent_pos .- (50.0, 50.0))
	write(out_file,"$(new_agent_pos[1]) $(new_agent_pos[2]) $(new_agent_pos[1]-50.0) $(new_agent_pos[2]-50.0)  $r_mp\n") 
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

moon_period = 2*pi*19.0^1.5
print("Moon period calculated to be ", moon_period, "\n")
rounded_moon_period = round(moon_period, digits = 3)
print("Rounded moon period calculated to be ", rounded_moon_period, "\n")
N = Int(floor(0.5*rounded_moon_period/0.001))
print("N calculated to be ", N, "\n")
for i in 1:260108
	step!(model, agent_step!)
end
###Test the model has been initialised and works
#= using InteractiveDynamics
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
    framerate = 10, frames = 1000,
    title = "Moon orbiting planet"
)
=#	
