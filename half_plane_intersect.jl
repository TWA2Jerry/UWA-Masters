using DataStructures
eps = 0.0000000001

function norm (v)
        sum_of_squares = 0.0
        for i in 1:length(v)
                sum_of_squares += (v[i])^2
        end

        return sqrt(sum_of_squares)
end


###Function for calculating the intersection between two points
function inter(h1, h2)
        #h1 and h2 represent the half planes we want to calculate the line intersections for
        m1 = h1[2][2]/h1[2][1]
        m2 = h2[2][2]/h2[2][1]
        c1 = h1[3][2] - m1*h1[3][1]
	c2 = h2[3][2] - m2*h2[3][1]
        xint = (c2-c1)/(m1-m2)
	yint = m1 * xint
	return [xint, yint]
end



###Function for calculating a cross product
function cross(v1, v2)
	return v1[1] * v2[2] - v1[2]*v2[1]
end



###Function for calculating whether or not a point lies within a half plane, returning 1 if it does lie outside
function outside(half_plane, point)
	return cross(half_plane[2], point - half_plane[3]) < -eps
end



###Function for generating the set of vertices defining the voronoi cell
function intercept_points (ri, neighbouring_points, rho)
	#ri represents the position of our agent i, neighbouring points should be a vector containing the positions of the neighbouring agents (the positions should also be represented as vectors)

	#Look at each of the neighbours of the agent, and generate the half planes
	half_planes = [] #The vector that will contain the half plane structures, which will be vectors comprised of the point and vector defining the half plane
	for point in neighbouring_points	
		#Use the half-way point as the point p
		r_ji = point .- ri
		half_plane_point = 0.5 .* r_ji .+ ri

		#Calculate the appropriate vector pq which lies parallel to the line in a direction such that the inner region is to the left of the vector
		v_jix = -1.0 * (0.5 * r_ji[2])
		v_jiy = 0.5 * r_ji[1] #Hopefully you can see that this is literally just v = [-sin(\theta), \cos(\theta)]
		pq = [v_jix, v_jiy]
		angle = atan(v_ijy, v_ijx)
		is_box = 0 #This is just to differentiate between the box and actual line segments later
                half_plane = [angle, pq, half_plane_point, is_box]
		push!(half_planes, half_plane)
	end

	#deque for the half planes/lines, I mean, technically you could just use Julia vectors with pushfirst and whatnot, but eh
	dq = []	

	#Add in the bounding box lines, and sort the vector of half planes according to their angles, note that the 1 at the end of the vector defining the half plane is simply to characterise them as box bounds so we can delete them later
	bottom_side = [0.0, [50.0, 0.0], [50.0, 0.0], 1]
	right_side = [pi/2, [0.0, 50.0], [100.0, 50.0], 1]
	top_side = [pi, [-50.0, 0.0], [50.0, 100.0], 1]
	left_side = [-pi/2, [0.0, -50.0], [0.0, 50.0], 1]

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
	
		#Check for parallel half planes
		if(len > 0 && norm(cross(half_planes[i][2], dq[len][2]))< eps) #Check if parallel by if the cross product is less than eps. Note that norm also works on scalars (returns abs val)
			if(half_planes[i][2] .* dq[len][2] < 0.0)
				print("Uh, Houston, we may have a anti-parallel pair")
				return -1
				if (out(half_planes[i], dq[len][3])) #Check if the last line in the 
				pop!(dq)
				len -= 1
			else continue
			end
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
	vertices = []
	dql = length(dq)
	for i in 1:length(dq)
		#Calculate the intersect between two thangs, and make sure they be valid
		intersect_i = inter(dq[i], dq[(i+1)%length(dq)])
		vhalf_int = intersect_i .- vertices[length(vertices)] #This is the vector from the last intersect to the new potential intersect
		v_proper = vhalf_int .* dq[i][2]
		if(v_proper < 0)
			#Calculate the appropriate intersect of the half plane dq[i] with the circle
			m = dq[i][2][2]/dq[i][2][1]
			c = dq[i][3][2] - m*dq[i][3][1]			
			
			x1 = (-(2*m*c) - sqrt((2*m*c)^2 - 4*(m^2+1)*(c^2-rho^2)))/(2*(m^2+1))
			y1 = m*x1 + c

			x2 = (-(2*m*c) - sqrt((2*m*c)^2 - 4*(m^2+1)*(c^2-rho^2)))/(2*(m^2+1))
			y2 = m*x2 + c
			
			#Okay, we should really check if the solutions aren't imaginary, but eh
			vhalf_int1 = [x1, y1] .- vertices[length(vertices)] #This is the vector from the last vertex to the intersect of the base edge with the circle
			circle_intersect_i = vhalf_int1 .* dq[i][2] < 0? [x2, y2] : [x1, y1] #This is to see we of the intersects is right, by testing if the first needs us to move "backwards" from the last vertex 
			push!(vertices, [circle_intersect_i,1])
			

			#Calculate the appropriate intersect of the half plane dq[(i+1)%length(dq)] with the circle		
                        m = dq[(i+1)%dql][2][2]/dq[(i+1)%dql][2][1]
                        c = dq[(i+1)%dql][3][2] - m*dq[(i+1)%dql][3][1]

                        x1 = (-(2*m*c) - sqrt((2*m*c)^2 - 4*(m^2+1)*(c^2-rho^2)))/(2*(m^2+1))
                        y1 = m*x1 + c

                        x2 = (-(2*m*c) - sqrt((2*m*c)^2 - 4*(m^2+1)*(c^2-rho^2)))/(2*(m^2+1))
                        y2 = m*x2 + c

                        #Okay, we should really check if the solutions aren't imaginary, but eh
                        vhalf_int2 = [x1, y1] .- vertices[0] #This is the vector from the last vertex to the intersect of the base edge with the circle
                        circle_intersect_ip1 = vhalf_int2 .* dq[(i+1)%dql][2] < 0? [x1, y1] : [x2, y2] #This is to see we of the intersects is right, by testing if the first needs us to move "backwards" from the last vertex
                        push!(vertices, [circle_intersect_ip1, 1])

	
			#Add these intersects to the list of edges, but label them as being circle edges
		else
			#Just add the intersect already calculated
			push!(vertices, [intersect_i, 0])
		end
	end

end
