function voronoi_area_brute(agent, model, rho)
	ri = agent.pos
	half_planes = [] #A vector containing all the half planes
	#For each neighbouring agent	
	for neighbour in all_agents_iterable
                if(neighbour.id == agent.id)
                        continue
                end
                #Calculate the line/half plane for that neighbour
		r_ji = neighbour.pos .- ri
                half_plane_point = 0.5 .* r_ji .+ ri

                #Calculate the appropriate vector pq which lies parallel to the line in a direction such that the inner region is to the left of the vector
                v_jix = -1.0 * (0.5 * r_ji[2])
                v_jiy = 0.5 * r_ji[1] #Hopefully you can see that this is literally just v = [-sin(\theta), \cos(\theta)]
                pq = [v_jix, v_jiy]
                angle = atan(v_jiy, v_jix)
                is_box = 0 #This is just to differentiate between the box and actual line segments later
                half_plane = [angle, pq, half_plane_point, is_box]
                push!(half_planes, half_plane)

	end
#For each angle
	for n in 0:thetaN
		#Set the default shortest distance to be rho 
			
		#Generate the line corresponding to the angle from the agent. Let's call this the angle line 

		#For each neighbouring agent
			#Calculate the intersect of the angle line with the neighbour half plane
			#Calculate the distance from agent to that intersect
				#If the distancce is less than the current shortest distance in that angle, set it to this thang. 
		#Add to the current running total of the voronoi area, d\theta*shortest_distance
end
