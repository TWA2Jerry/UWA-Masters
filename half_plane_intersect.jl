eps = 0.0000000001
inf = 1000000000000.0
function norm(v)
        sum_of_squares = 0.0
        for i in 1:length(v)
                sum_of_squares += (v[i])^2
        end

        return sqrt(sum_of_squares)
end

function dot(v1, v2)
	if(length(v1) != length(v2))
	   print("Bruh these vectors ain't got the same length no?\n")
	   return -1
	end

	sum = 0.0
	for i in length(v1)
		sum += v1[i]*v2[i]
	end
	return sum
end

###Function for calculating the intersection between two points
function inter(h1, h2)
        #h1 and h2 represent the half planes we want to calculate the line intersections for
	#print("Calculating the intersection for $(h1[2]) and $(h2[2])\n")
	#m1 = h1[2][2]/h1[2][1]
	#m2 = h2[2][2]/h2[2][1]
	#m1 = 0.0
	#m2 = 0.0
	if(h1[2][1] == 0.0)
		#print("Infinite gradient detected for m1\n")
		m1 = inf
	else 
		m1 = h1[2][2]/h1[2][1]
	end
	#print("m1 found to be $m1\n")
	if(h2[2][1] == 0.0)
		m2 = inf
		#print("Infinite gradient detected for m2\n")
	else
		m2 = h2[2][2]/h2[2][1]
	end	
	if(abs(m1 - m2) < abs(eps))
		print("Parallel planes yo\n")
		exit()
		return -1
	end  
	##print("m1 - m2 found to be $(m1-m2)\n")
	c1 = h1[3][2] - m1*h1[3][1]
	c2 = h2[3][2] - m2*h2[3][1]
        xint = (c2-c1)/(m1-m2)
	yint = -1
	if(abs(m1 - inf) < abs(eps))
		yint = m2 * xint + c2
	else 
		yint = m1 * xint + c1
	end
	#print("Intersect calculated as $([xint, yint])\n")
	return [xint, yint]
end



###Function for calculating the magnitude of the cross product between two vectors v1 and v2
function cross(v1, v2)
	return v1[1] * v2[2] - v1[2]*v2[1]
end



###Function for calculating whether or not a point lies within a half plane, returning 1 if it does lie outside
function outside(half_plane, point)
	#What we're doing here is, we calculate the vector from the point on the half plane fence to the point that we're considering (which is the intersection). Call this vector b. The vector of the half plane is a. Now, according to the right hand rule for the calculation of the cross product, the point (intersection) will be to the right of the half-plane if the cross product points in the negative z direction, that is, ax*by-ay*bx < 0.
	return cross(half_plane[2], point .- half_plane[3]) < -eps
end



###Function for generating the set of vertices defining the voronoi cell
function voronoi_cell(ri, neighbouring_points, rho)
	#ri represents the position of our agent i for whom we wish to calculate the voronoi cell, neighbouring points should be a vector containing the positions of the neighbouring agents (the positions should also be represented as vectors)

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
		angle = atan(v_jiy, v_jix)
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

	sort!(half_planes)
	#print("After sorting, half_planes is given by \n")
	#=for i in 1:length(half_planes)
		print("$(half_planes[i])\n")
	end
	=#
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
			len -= 1
		end
	
		#Check for parallel half planes
		if(len > 0 && norm(cross(half_planes[i][2], dq[len][2]))< eps) #Check if parallel by if the cross product is less than eps. Note that norm also works on scalars (returns abs val)
			if(dot(half_planes[i][2], dq[len][2]) < 0.0)
				#print("Uh, Houston, we may have a anti-parallel pair")
				return -1
			end
			if (outside(half_planes[i], dq[len][3])) #Check if the last line in the dq is outside the half plane we're about to add 
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
	while(len > 2 && outside(dq[1], inter(dq[len], dq[len-1])))
		pop!(dq)
		len -= 1
	end

	while(len > 2 && outside(dq[len], inter(dq[1], dq[2])))
		popfirst!(dq)
		len -= 1
	end

	#Report empty intersection (if number of edges is less than 3?)
	if(len < 3)
		print("Yo this intersection be empty")
		exit()
		return -1
	end

	#print("dq processing complete, the deqeue is given by $dq")
	#Having found the voronoi cell with the bounded box method, we now account for the fact that we have a bounding circle and not a box, and so get rid of the box line segments first
		
	i = 1
	while (i <= length(dq))
		if(dq[i][4]==1 || norm(dq[i][3] .- ri) > rho)
			deleteat!(dq, i)
			continue
		end
		i += 1
	end
	
	#=
	print("We have now removed all bounding box and redundant half-planes; for the agent position of $ri, the remaining half planes are (given by their vectors) \n")
	for half_plane in dq
		#print("The half plane is $half_plane\n")
		print("The half plane angle is $half_plane\n")
	end
	=#

	#Now, go through and start calculating the intersects between the non-redundant lines, but if there is no valid intersect, then use the circle
	vertices = []
	dql = length(dq)
	for i in 1:length(dq)
		if(dql == 1) #To handle the case if there's only one fence for an agent
			#print("Single fence agent detected\n")
			m = dq[i][2][2]/dq[i][2][1]
                        #print("Gradient for i is $m \n")
                        c = dq[i][3][2] - m*dq[i][3][1]
                        #print("c for i is $c\n")

                        a = 1+m^2
                        b = -2*ri[1]+2*m*c-2*m*ri[2]
                        d = ri[1]^2 + c^2 - 2*ri[2]*c + ri[2]^2 - rho^2
                        x1 = (-(b) - sqrt((b)^2 - 4*(a)*(d)))/(2*(a))
                        #print("x1 calculated to be $x1\n")
                        y1 = m*x1 + c

                        x2 = (-(b) + sqrt((b)^2 - 4*(a)*(d)))/(2*(a))
                        y2 = m*x2 + c

			push!(vertices, [[x1, y1], 1, i])
			push!(vertices, [[x2, y2], 1, i])
			break
		end
		#Calculate the intersect between two thangs, and make sure they be valid
		#print("Iteration of $i\n")
		v_proper = -1
		outside_circle = 0
		intersect_i = inter(dq[i], dq[(i)%length(dq)+1])
		#print("The intersect is $intersect_i\n")
		if(intersect_i != -1) #Only consider looking at whether or not the intersect is "in front" if the planes aren't parallel
			#print("Passed non-parallel condition\n")
			if(norm(ri .- intersect_i) <= rho)
				#print("Passed within circle condition\n")
				if(i == 1)
					#print("Passed i = 1 condition\n")
					v_proper = 1.0
				else 
					vhalf_int = intersect_i .- vertices[length(vertices)][1] #This is the vector from the last intersect to the new potential intersect
					v_proper = dot(vhalf_int, dq[i][2])
					if(norm(intersect_i .- vertices[length(vertices)][1]) < eps)
						v_proper = -1
					end
					#print("Dot product of old->new intersect with new plane is $v_proper\n")
				end
			else
				print("Intersect was calculated to be outside the circle.\n")
				outside_circle = 1
			end
		end
		if(v_proper < 0.0)
			#print("Still no valid intersect detected\n")
			#Calculate the appropriate intersect of the half plane dq[i] with the circle
			m = dq[i][2][2]/dq[i][2][1]
			#print("Gradient for i is $m \n")
			c = dq[i][3][2] - m*dq[i][3][1]	
			#print("c for i is $c\n")
			
			a = 1+m^2
			b = -2*ri[1]+2*m*c-2*m*ri[2]
			d = ri[1]^2 + c^2 - 2*ri[2]*c + ri[2]^2 - rho^2
			x1 = (-(b) - sqrt((b)^2 - 4*(a)*(d)))/(2*(a))
			#print("x1 calculated to be $x1\n")
			y1 = m*x1 + c

			x2 = (-(b) + sqrt((b)^2 - 4*(a)*(d)))/(2*(a))
			y2 = m*x2 + c
			
			#Okay, we should really check if the solutions aren't imaginary, but eh
			#REDUNDANT vhalf_int1 = [x1, y1] .- vertices[length(vertices)] #This is the vector from the last vertex to the intersect of the base edge with the circle
			a1_a2 = [x1, y1] .- [x2, y2] #Calculate the vector from the second to first intersect with the circle
			#print("The value of a1_a2 is $a1_a2 and the value of the half plane is $(dq[i][2])\n")
			#print("The dot product of a1_a2 with the vector of the half segment is $(a1_a2 .*  dq[i][2])\n")
			circle_intersect_i = dot(a1_a2, dq[i][2]) >=  0.0 ? [x1, y1] : [x2, y2] #This is to see we of the intersects is right, by testing if the first needs us to move "backwards" from the last vertex 
			
			if(i > 1)       
				vec_i = circle_intersect_i .- vertices[length(vertices)][1]
				direction_of_vec = dot(vec_i, dq[i][2])
				#print("The direction of the circle intersect from the last vertex is $direction_of_vec\n")
                        end
			push!(vertices, [circle_intersect_i,1, i])

			#Calculate the appropriate intersect of the half plane dq[(i+1)%length(dq)] with the circle		
                        m = dq[(i)%dql+1][2][2]/dq[(i)%dql+1][2][1]
			c = dq[(i)%dql+1][3][2] - m*dq[(i)%dql+1][3][1]

			a = 1+m^2
                        b = -2*ri[1]+2*m*c-2*m*ri[2]
                        d = ri[1]^2 + c^2 - 2*ri[2]*c + ri[2]^2 - rho^2
                        x1 = (-(b) - sqrt((b)^2 - 4*(a)*(d)))/(2*(a))
                        #print("x1 calculated for i+1 to be $x2\n")
                        y1 = m*x1 + c

                        x2 = (-(b) + sqrt((b)^2 - 4*(a)*(d)))/(2*(a))
                        y2 = m*x2 + c


                        #Okay, we should really check if the solutions aren't imaginary, but eh
                        #REDUNDANT vhalf_int2 = [x1, y1] .- vertices[0] #This is the vector from the last vertex to the intersect of the base edge with the circle
			b1_b2 = [x1, y1] .- [x2, y2] #Calculation of the vector from the second intersect to first intersect 
			circle_intersect_ip1 = dot(b1_b2, dq[(i)%dql+1][2]) <= 0 ? [x1, y1] : [x2, y2] #This is to see we of the intersects is right, by testing if the first needs us to move "backwards" from the last vertex
                        push!(vertices, [circle_intersect_ip1, 1, i%dql+1])
				
			i_vec = circle_intersect_i .- ri
			ip1_vec = circle_intersect_ip1 .- ri
			angle_i = atan(i_vec[2], i_vec[1])
			angle_ip1 = atan(ip1_vec[2], ip1_vec[1])
			if(outside_circle == 1)
				bro = intersect_i .- ri
				bruh = atan(bro[2], bro[1])
				#print("For an original intersect outside the circle, it had an angle of $bruh. ")
				print("The original intersect outside the circle for half planes $i and $(i%dql+1) was $intersect_i")
			end
			#print("Due to an invalid intersect, intersects with circle calculated instead. Intersects were (by angle, i followed by i+1) $angle_i and $angle_ip1\n")
			print("Due to an invalid intersect between half planes $i and $(i%dql+1), intersects with circle calculated instead. Intersects were (i and i+1): $circle_intersect_i, $circle_intersect_ip1\n")

			#Add these intersects to the list of edges, but label them as being circle edges
		else
			#Just add the intersect already calculated
			push!(vertices, [intersect_i, 0, i])
		end
	end
	if(length(vertices) == 2 && norm(vertices[1][1] .- vertices[2][1]) < eps)
		print("Single intersect calculated (area of 0). The dq which generated this was $dq\n")
	end
	return vertices

end

#Square test case
#neighbouring_positions = [[60.0, 50.0], [50.0, 60.0], [40.0, 50.0], [50.0, 40.0]]

#Triangle test case
#agent_pos = [50.0, 50.0]
#neighbouring_positions = [[55.0, 55.0], [45.0, 55.0], [50.0, 45.0]]
#

#Weird triangle test case
#=
neighbouring_positions = [[55.0, 55.0], [45.0, 55.0], [55.0, 45.0]]
=#

#Only two neighbours, parallel planes
#neighbouring_positions = [[45.0, 55.0], (55.0, 45.0)]

#=
agent_pos = [50.0, 50.0]
rho = 10.0
vertices = voronoi_cell(agent_pos, neighbouring_positions, rho)
print("\n\nNow displaying the vertices\n")
for vertex in vertices
	print(vertex[1])
end
=#
