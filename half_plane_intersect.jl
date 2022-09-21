using DataStructures

function intercept_points (ri, neighbouring_points)
	#ri represents the position of our agent i, neighbouring points should be a vector containing the positions of the neighbouring agents (the positions should also be represented as vectors)

	#Look at each of the neighbours of the agent, and generate the half planes
	half_planes = [] #The vector that will contain the half plane structures, which will be vectors comprised of the point and vector defining the half plane
	for point in neighbouring_points	
		#Use the half-way point as the point p
		r_ji = point .- ri
		half_plane_point = 0.5 .* v_ji .+ ri

		#Calculate the appropriate vector pq which lies parallel to the line in a direction such that the inner region is to the left of the vector
		v_jix = -1.0 .* (0.5 .* r_ji[2])
		v_jiy = 0.5 .* r_ji[1] #Hopefully you can see that this is literally just v = [-sin(\theta), \cos(\theta)]
		pq = [v_jix, v_jiy]
	end

	#deque for the half planes/lines, I mean, technically you could just use Julia vectors with pushfirst and whatnot, but eh
	

	//Put the first half plane in, and calculate its intercepts with the circle. Do this using the fact that y = mx+c, and the equation of the circle is x^ + y^2 = \rho ^2

	//Iterate through the other half planes
		//For all previous intersections, if the intersection is outside the about-to-be added half-plane, remove that corresponding half-plane from the front of the queue. Of course, if it was an intersection with the circle, don't remove the corresponding half plane just yet

		//Remove any half planes from the back of the queue, again, don't do it if 

		//So yeah, we should keep a dequeue also of the intercepts, and go through the dequeue of intersects instead, but keeping a track of whether the intersect is in fact an intersect with the circle. If it is with the circle, do not remove the latest line from the lines/half planes dequeue. 

		//Add the half plane for this iteration, and also add its front end intersect with the circle

end
