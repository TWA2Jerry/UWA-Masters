###Function for generating the set of vertices defining the voronoi cell
function voronoi_cell_bounded(ri, neighbouring_points, rho, temp_half_planes = [], vel = [0.0,0.0], relic = [atan(0.0), [0.0, 0.0], [0.0,0.0], 0])
	#ri represents the position of our agent i for whom we wish to calculate the voronoi cell, neighbouring points should be a vector containing the positions of the neighbouring agents (the positions should also be represented as vectors)
###This is the section for deriving the original voronoi cell
	#Look at each of the neighbours of the agent, and generate the half planes
	#print(outside([1, [1,2], [3,4]], [5,6]))
	half_planes = [] #The vector that will contain the half plane structures, which will be vectors comprised of the point and vector defining the half plane
	for point in neighbouring_points	
		#Use the half-way point as the point p
		r_ji = point .- ri
		half_plane_point = 0.5 .* r_ji .+ ri

		#Calculate the appropriate vector pq which lies parallel to the line in a direction such that the inner region is to the left of the vector
		v_jix = -1.0 * (0.5 * r_ji[2])
		v_jiy = 0.5 * r_ji[1] #Hopefully you can see that this is literally just v = [-sin(\theta), \cos(\theta)]
		pq = [v_jix, v_jiy]
		#print("$pq\n")
		angle = atan(v_jiy, v_jix)
		is_box = 0 #This is just to differentiate between the box and actual line segments later
		half_plane = [angle, pq, Tuple(half_plane_point), is_box]
		push!(half_planes, half_plane)
	end


	#deque for the half planes/lines, I mean, technically you could just use Julia vectors with pushfirst and whatnot, but eh
	
	dq = []	

	#Add in the bounding box lines, and sort the vector of half planes according to their angles, note that the 1 at the end of the vector defining the half plane is simply to characterise them as box bounds so we can delete them later
	
	bottom_side = [0.0, [50.0, 0.0], Tuple([0.0, -1000.0]), 1]
	right_side = [pi/2, [0.0, 50.0], Tuple([1000.0, 0.0]), 1]
	top_side = [pi, [-50.0, 0.0], Tuple([0.0, 1000.0]), 1]
	left_side = [-pi/2, [0.0, -50.0], Tuple([-1000.0, 0.0]), 1]
	push!(half_planes, bottom_side)
	push!(half_planes, right_side)
	push!(half_planes, top_side)
	push!(half_planes, left_side)
	
	#Add the half plane that bounds the area to the area in front of the agent
        fw_point = ri
        fw_x = -1.0*(-vel[2])
        fw_y = -vel[1] #you might be confused about the negative sign, remember that this is meant to be the vector of the half plane, which is the vector fo the velocity rotated 90 degrees clockwise. 
        fw_pq = [fw_x, fw_y]
        angle = atan(fw_y, fw_x)
        fw_is_box = 2
        fw_half_plane = [angle, fw_pq, Tuple(ri), fw_is_box]
        push!(half_planes, fw_half_plane)

	#For the relic version of stemler vision, we also need to retain the relic half plane as a bounding half plane for all sampled positions
	#print("About to push the relic, which is $relic\n")
	#push!(half_planes, relic)

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
		print("Yo this intersection be empty, the position considered was $ri\n")
		exit()
		return -1
	end
	#print("dq processing complete, the deqeue is given by $dq")
	





###
###This is the section where we account for the circle of vision
	#Having found the voronoi cell with the bounded box method, we now account for the fact that we have a bounding circle and not a box, and so get rid of the box line segments first
	#print("\n\n\nCommencing bounded DOD calculation\n")	
	i = 1
	while (i <= length(dq))
		if(dq[i][4]==1 || norm(dq[i][3] .- ri) >= rho)
			deleteat!(dq, i)
			continue
		end
		i += 1
	end
	
	#=		
	print("We have now removed all bounding box and redundant half-planes; for the agent position of $ri, the remaining half planes are (given by their vectors) \n")
	for half_plane in dq
		print("The half plane is $half_plane\n")
	end
	=#	
	
	#print("In the voronoi function, a was modified to a value of $a\n")
	replace_vector(temp_half_planes, dq)
	#Now, go through and start calculating the intersects between the non-redundant lines, but if there is no valid intersect, then use the circle
	vq = []
	newdq = []
	dql = length(dq)
	len = 0
	vlen = 0
	for i in 1:dql
		#print("Single fence agent detected\n")
		m = 0.0
		if(dq[i][4] == 1)
			print("Bounding box fence detected\n")
			exit()
		end
		if(abs(dq[i][2][1]) > eps) 
			m = dq[i][2][2]/dq[i][2][1]
		else 
			m = inf
		end
                #print("Gradient for i is $m \n")
		x1 = 0.0
		y1 = 0.0
		x2 = 0.0
		y2 = 0.0

		if(abs(m-inf) < eps)
			x1 = dq[i][3][1]
			x2 = dq[i][3][1]
			if(rho^2 - (x1 - ri[1])^2 <0)
				print("Circle intercept negative. This was for the infinite gradient case, with half plane $(dq[i]).\n")
				exit()
			end
			y1 = -sqrt(rho^2 - (x1 - ri[1])^2) + ri[2]
			y2 = sqrt(rho^2 - (x2 - ri[1])^2) + ri[2]
		else
			c = dq[i][3][2] - m*dq[i][3][1]
                	#print("c for $i is $c\n")
                	a = 1+m^2
                	b = -2*ri[1]+2*m*c-2*m*ri[2]
                	d = ri[1]^2 + c^2 - 2*ri[2]*c + ri[2]^2 - rho^2
			if((b)^2 - 4*(a)*(d) < 0)
				print("Circle intercept negative. This was for the normal case for the half plane $(dq[i]). The value of the discriminant was $((b)^2 - 4*(a)*(d))\n")
				exit()
			end
			x1 = (-(b) - sqrt((b)^2 - 4*(a)*(d)))/(2*(a))
                	y1 = m*x1 + c

                	x2 = (-(b) + sqrt((b)^2 - 4*(a)*(d)))/(2*(a))
                	y2 = m*x2 + c
		end
		#print("x1, y1, x2, y2 calculated to be $x1, $y1, $x2, $y2\n")
		
		a1_a2 = [x1, y1] .- [x2, y2] #Calculate the vector from the second to first intersect with the circle
                f_circle_intersect_i = dot(a1_a2, dq[i][2]) >=  0.0 ? [x1, y1] : [x2, y2] #This is to see we of the intersects is right, by testing if the first needs us to move "backwards" from the last vertex 
		#print("f_circle was selected to be $f_circle_intersect_i because a1_a2 was $a1_a2, the dequeue vector was $(dq[i][2]) resulting in a dot product of $(dot(a1_a2, dq[i][2]))\n")
		b_circle_intersect_i = dot(a1_a2, dq[i][2]) <=  0.0 ? [x1, y1] : [x2, y2]

		while(vlen >= 1 && outside(dq[i], vq[vlen][1]))
			if(vq[vlen][3] != 0)
				if(len <= 0)
                                        print("Tried to delete a half plane from dequeue when there wasn't one\n")
                                        exit()
                                end
				pop!(newdq)
				len -= 1
			end
			#print("Popping from the back of the newdq. The back is $(vq[vlen])\n")
			pop!(vq)
                        vlen -= 1
                end

                #Remove any half planes from the back of the queue, again, don't do it if 
		while(vlen >= 1 && outside(dq[i], vq[1][1]))
			if(vq[1][2] != 0)
				if(len <= 0)
					print("Tried to delete a half plane from dequeue when there wasn't one\n")
					exit()
				end
				popfirst!(newdq)
				len -= 1
			end
			#print("Popping from the front of the dequeue. The front is $(vq[1])\n")
			popfirst!(vq)
                        vlen -= 1
                end

		#Check for parallel half planes is no longer needed as far as I'm aware, because we...okay maybe it is...because if you cut out none, you'll still end up trying to compute parallel plane intersects
		if(len > 0 && norm(cross(dq[i][2], newdq[len][2]))< eps && outside(newdq[len], dq[i][3])) #Check if parallel by if the cross product is less than eps. Note that norm also works on scalars (returns abs val)
                        continue
                end

		invalid_half_plane = 0
		if (len >= 1)
			#Determine the intersect of hp_i with hp_(i-1)
			intersect_i = inter(dq[i], newdq[len])
			is_outside = 0
			invalid = 0
			if(intersect_i == -1 || norm(intersect_i .- ri) > rho)
				is_outside = 1	
			end

			forward_after_inter = intersect_i .+ 1.0 .* dq[i][2]
			if(intersect_i == -1 || outside(newdq[len], forward_after_inter))
				invalid = 1
			end

			if(is_outside == 0 && invalid == 0)
				push!(vq, [intersect_i, i-1, i])
				vlen += 1
				#print("Normal intersect pushed for i = $i. Intersect was $intersect_i\n")
			elseif(!outside(newdq[len], b_circle_intersect_i))
				if(dq[i][4] == 2)
					push!(vq, [b_circle_intersect_i, 0, -1])
				else
					push!(vq, [b_circle_intersect_i, 0, i])
				end
				vlen += 1
				#print("Circle intersect pushed for i = $i. Intersect was $b_circle_intersect_i\n")
			else 
				invalid_half_plane = 1
			end	

		else #Note that we've changed the format to account for what to label the forwards half plane is if it was the artifical bounding half plane
			if(dq[i][4] == 2)
                                        push!(vq, [b_circle_intersect_i, 0, -1])
                        else
                                        push!(vq, [b_circle_intersect_i, 0, i])                                
			end

			vlen += 1
			#print("Circle intersect pushed for i = $i. Intersect was $b_circle_intersect_i\n")
		end
		
		if(invalid_half_plane == 1)
			continue
		end

		#Add the foward intersect
		if(dq[i][4] == 2)
			push!(vq, [f_circle_intersect_i, -1, 0])
                else
                       	push!(vq, [f_circle_intersect_i, i, 0])
                end

		vlen += 1
		#print("Forward intersect pushed for i = $i. Intersect was $f_circle_intersect_i\n")
		#Add the new half plane
                push!(newdq, dq[i])
                len += 1
		#print("Half plane $i added. The dequeues vq and dq are now \n")
		#print("$vq\n")
		#print("$newdq\n")
        end

	#print("Commencing final cleanup\n")
        #Do a final cleanup
	while(vlen >= 2 && outside(newdq[1], vq[vlen][1]))

		if(vq[vlen][3] != 0)
			#print("Popping from the back of the newdq because the intersect was $(vq[vlen]). The back is $(newdq[len])\n")
			if(len <= 0)
                                print("Tried to delete a half plane from dequeue when there wasn't one\n")
                        	exit()
                        end

			pop!(newdq)
                	len -= 1
		end
		pop!(vq)
		vlen -= 1

        end

	while(vlen >= 2 && outside(newdq[len], vq[1][1]))
		if(vq[1][2] != 0)
			#print("Popping from the front of newdq. The front is $(vq[1])\n")
			if(len <= 0)
                                print("Tried to delete a half plane from dequeue when there wasn't one\n")
                        	exit()
                        end
			popfirst!(newdq)
                	len -= 1
		end
                popfirst!(vq)
                vlen -= 1

        end
		
	#=	
	print("After cleanup, the final half planes were\n")
	for hp in newdq
		print("$hp\n")
	end
	=#

	#Finally, look at the link between the first and last half-planes, if it's valid, add it, if it's not, then the circle intersects would've already been added. 
	#print("Commencing link between first and last planes\n")
	if (len > 1)
                        #Determine the intersect of hp_i with hp_(i-1)
                        intersect_last = inter(newdq[len], newdq[1])
                        is_outside = 0
                        invalid = 0
                        if(norm(intersect_last .- ri) > rho)
                                is_outside = 1
                        end

			forward_after_inter = intersect_last .+ 1.0 .* newdq[1][2]
                        if(outside(newdq[len], forward_after_inter))
                                invalid = 1
                        end

                        if(is_outside == 0 && invalid == 0)
                                push!(vq, [intersect_last, len, 1])
                                vlen += 1
                        end
	end

	if(length(vq) == 2 && norm(vq[1][1] .- vq[2][1]) < eps)
		#print("Single intersect calculated (area of 0). The dq which generated this was $dq\n")
	end

	
	return vq

end
