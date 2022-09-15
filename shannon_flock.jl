###Introduction
#In this template, we generate an agent, but now it moves with some velocity and acceleration. However, the acceleration will be randomised at each time step and the velocity updated accordingly
#One thing that we're looking out for in particular is if we can update position, velocity and acceleration, which will be N tuples, with vectors generated from our ODE or movement gradient function. Actually no, in accordance with RK4, or even Euler integration, acceleration isn't an quantity we need to keep track of.

###Preliminaries
using Agents
using Random
using VoronoiCells


###Create the function returns the area of a given Voronoi cell given the vertices which define the cell
function voronoi_area(vertices)
        num_points = length(vertices)
        A = 0.0
        for i in range 1:num_points
                j = (i+1)%num_points
                xi = vertices[i][1]
                yi = vertices[i][2]
                xj = vertices[j][1]
                yj = vertices[j][2]
                A += 0.5 * (yi + yj)* (xi - xj)
        end

        return A
end



###Function that generates the half planes assoiated with a given agent to its neighbours



###Create the function that takes a given agent, model and then find the vertices that comprise the voronoi cell
function voronoi_vertices(agent, model)

end


###Create the movement gradient function that gives the rate of change in position and velocity, I'll use RK4, cause why not. Actually, no RK4 since we don't have a consistent acceleration function
function move_gradient(agent, model, k_pos, k_vel, kn)
        #Determine the neighbours of the agent. In this case, we want the moon, and we look through all of possible space

