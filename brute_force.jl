function voronoi_area_brute(ri, neighbour_positions, rho, thetaN)
	total_area = 0.0
	#Generate and record the angle theta_ij for all j, as well as dij in a vector or tuple
	theta_dij = []
	for rj in neighbour_positions
		rij = rj .- ri
		dij = norm(rij)
		thetaij = atan(rij[2], rij[1])
		push!(theta_dij, [thetaij, dij])
	end

	#For each angle
	for n in 0:(thetaN-1)
		theta = 0.0 + n*2*pi/thetaN
		#Set the default shortest distance to be rho 
		dtheta = rho	
		#For each neighbour, calculate the distance from the agent to the fence with the neighbour at the particular angle
		for dp in theta_dij
			thetaij = dp[1]
			dij = dp[2]
			dfj = dij/2 * sec(thetaij - theta)
		#If the distancce is less than the current shortest distance in that angle, set it to this thang. 		
			if(dfj < 0)
				continue
			end

			if (dfj < dtheta)
				dtheta = dfj
			end
		end
		#Add to the current running total of the voronoi area, d\theta*shortest_distance
		total_area += dtheta * 2*pi/thetaN	
	end

	return total_area
end
