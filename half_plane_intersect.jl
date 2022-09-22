using DataStructures
eps = 0.0000000001


###Function for calculating the intersection between two points
function inter(h1, h2)
        #h1 and h2 represent the half planes we want to calculate the line intersections for
        m1 = h1[2][2]/h1[2][1]
        m2 = h2[2][2]/h2[2][1]
        c1 = [2] -
        xint = (c2-c1)/(m1-m2)
end



###Function for calculating a cross product
function cross_product(v1, v2)
	return v1[1] * v2[2] - v1[2]*v2[1]
end



###Function for calculating whether or not a point lies within a half plane
function outside(half_plane, point)
	return cross(half_plane[2], point - half_plane[3]) < -eps
end




###Function for calculating whether or not a point lies within a half plane
function outside(half_plane, point)
	return cross(half_plane[2], point - half_plane[3]) < -eps
end



###Function for generating the set of vertices defining the voronoi cell
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
		angle = atan(v_ijy, v_ijx)
                half_plane = [angle, pq, half_plane_point]
		push!(half_planes, half_plane)
	end

	#deque for the half planes/lines, I mean, technically you could just use Julia vectors with pushfirst and whatnot, but eh
	dq = []	

	#Add in the bounding box lines, and sort the vector of half planes according to their angles
	bottom_side = [0.0, [50.0, 0.0], [50.0, 0.0]]
	right_side = [pi/2, [0.0, 50.0], [100.0, 50.0]]
	top_side = [pi, [-50.0, 0.0], [50.0, 100.0]]
	left_side = [-pi/2, [0.0, -50.0], [0.0, 50.0]]

	push!(half_planes, bottom_side)
	push!(half_planes, right_side)
	push!(half_planes, top_side)
	push!(half_planes, left_side)

	sort(half_planes)


	#Iterate through the half planes, adding them one at a time, cutting out any redundant half planes as we go
	len = 0
	for i in 1:length(half_planes)	
		#For all previous intersections, if the intersection is outside the about-to-be added half-plane, remove that corresponding half-plane from the front of the queue. Of course, if it was an intersection with the circle, don't remove the corresponding half plane just yet
		while(len > 1 && outside(half_planes[i], inter(dq[len], dq[len-1])))
			pop!(dq)
			len -= 1
		end
		
		#Remove any half planes from the back of the queue, again, don't do it if 
		while(len > 1 && outside(half_planes[i], inter(dq[1], dq[2])))
			popfirst!(dq)
			lne -= 1
		end
		
		
		#Remove any half planes from the back of the queue, again, don't do it if 
		while(len > 1 && outside(half_planes[i], inter(dq[1], dq[2])))
			popfirst!(dq)
			lne -= 1
		end
		
		#Add the new half plane
		push!(dq, half_planes[i])
		len += 1
	end

	#Do a final cleanup
	while(len > 2 && outside(dq[0], inter(dq[len], dq[len-1])))
		pop(dq)
		len -= 1
	end

	while(len > 2 && outside(dq[len], inter(dq[1], dq[2])))
		popfirst!(dq)
		len -= 1
	end

	#Report empty intersection (if number of edges is less than 3?)
	if(len < 3)
		print("Yo this intersection be empty")
		return -1
	end


	#Having found the voronoi cell with the bounded box method, we now account for the fact that we have a bounding circle and not a box, and so get rid of the box line segments first
	for i in 1:length(dq)
		if(dq[i][4])
			deleteat!(dq, i)
		end
	end

	#Now, go through and start calculating the intersects between the non-redundant lines, but if there is no valid intersect, then use the circle
	

end
