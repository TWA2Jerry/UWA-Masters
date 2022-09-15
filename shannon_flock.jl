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
	return 	Areas[1]
end


###



###Function that determines the gradient of movement
function movement_gradient(agent, model, agent_speed)
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
		new_area = voronoi_area(positions)
		if new_area < min_area
			min_area = new_area
			min_direction = direction_of_move
		end	
	end
	
end

###Create the function that takes a given agent, model and then find the vertices that comprise the voronoi cell
function voronoi_vertices(agent, model)

end


###Create the movement gradient function that gives the rate of change in position and velocity, I'll use RK4, cause why not. Actually, no RK4 since we don't have a consistent acceleration function
function move_gradient(agent, model, k_pos, k_vel, kn)
        #Determine the neighbours of the agent. In this case, we want the moon, and we look through all of possible space

