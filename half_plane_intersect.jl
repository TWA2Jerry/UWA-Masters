function intercept_points 

	//Look at each of the neighbours of the agent, and generate the half planes
		//Use the half-way point as the point p

		//Calculate the appropriate vector pq which lies parallel to the line in a direction such that the inner region is to the left of the vector

	//deque for the half planes/lines

	//Put the first half plane in, and calculate its intercepts with the circle. Do this using the fact that y = mx+c, and the equation of the circle is x^ + y^2 = \rho ^2

	//Iterate through the other half planes
		//For all previous intersections, if the intersection is outside the about-to-be added half-plane, remove that corresponding half-plane from the front of the queue. Of course, if it was an intersection with the circle, don't remove the corresponding half plane just yet

		//Remove any half planes from the back of the queue, again, don't do it if 

		//So yeah, we should keep a dequeue also of the intercepts, and go through the dequeue of intersects instead, but keeping a track of whether the intersect is in fact an intersect with the circle. If it is with the circle, do not remove the latest line from the lines/half planes dequeue. 

		//Add the half plane for this iteration, and also add its front end intersect with the circle

end
