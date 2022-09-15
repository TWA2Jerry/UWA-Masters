###Introduction
#In this template, we generate an agent, but now it moves with some velocity and acceleration. However, the acceleration will be randomised at each time step and the velocity updated accordingly
#One thing that we're looking out for in particular is if we can update position, velocity and acceleration, which will be N tuples, with vectors generated from our ODE or movement gradient function. Actually no, in accordance with RK4, or even Euler integration, acceleration isn't an quantity we need to keep track of.

###Preliminaries
using Agents
using Random
using VoronoiCells

###Function that calculates the voronoi area of a given agent a given position


###Function that generates the half planes assoiated with a given agent to its neighbours



###Create the function that takes a given agent, model and then find the vertices that comprise the voronoi cell
function voronoi_vertices(agent, model)

end


###Create the movement gradient function that gives the rate of change in position and velocity, I'll use RK4, cause why not. Actually, no RK4 since we don't have a consistent acceleration function
function move_gradient(agent, model, k_pos, k_vel, kn)
        #Determine the neighbours of the agent. In this case, we want the moon, and we look through all of possible space

