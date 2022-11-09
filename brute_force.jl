function voronoi_area_brute(ri, neighbour_positions, rho, thetaN)
	total_area = 0.0
	#Generate and record the angle theta_ij for all j, as well as dij in a vector or tuple
	theta_dij = []
	vertices = []
	for rj in neighbour_positions
		rij = rj .- ri
		dij = norm(rij)
		thetaij = atan(rij[2], rij[1])
		half_plane_point = 0.5 .* r_ji .+ ri
                #Calculate the appropriate vector pq which lies parallel to the line in a direction such that the inner region is to the left of the vector
                v_jix = -1.0 * (0.5 * r_ji[2])
                v_jiy = 0.5 * r_ji[1] #Hopefully you can see that this is literally just v = [-sin(\theta), \cos(\theta)]
                pq = [v_jix, v_jiy]
                angle = atan(v_jiy, v_jix)
                is_box = 0 #This is just to differentiate between the box and actual line segments later
                half_plane = [angle, pq, half_plane_point, is_box]

		push!(theta_dij, [thetaij, dij, half_plane])
	end
	
	argminD = zeros(Int64, thetaN)

	#For each angle
	for n in 1:(thetaN)
		theta = 0.0 + n*2*pi/thetaN
		#Set the default shortest distance to be rho 
		dtheta = rho
		lowest_j = 0
		#For each neighbour, calculate the distance from the agent to the fence with the neighbour at the particular angle
		for j in 1:length(theta_dij)
			dp = theta_dij[j]
			thetaij = dp[1]
			dij = dp[2]
			dfj = dij/2 * sec(thetaij - theta)
		#If the distancce is less than the current shortest distance in that angle, set it to this thang. 		
			if(dfj < 0)
				continue
			end

			if (dfj < dtheta)
				dtheta = dfj
				lowest_j = j
			end
		end
		#Add to the current running total of the voronoi area, d\theta*shortest_distance
		argminD[n] = lowest_j
	end

	#Go through the vector of lowest indices
	current_lowest_index = argminD[1]
	for i in 2:length(argminD)
		#Each time you hit a new index, calculate the intersect between the line of that index and the previous line
		if (argminD[i] != current_lowest_index)
			if(current_lowest_index == 0)
				
			else if(argminD[i] == 0)

			else

			end
		end
		#Finally, iff the first and last indices are different, then calculate their intersect as well
	end
	return total_area
end
