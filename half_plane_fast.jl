include("intersect_check.jl")
const eps::Float64 = 0.0000000001
const inf::Float64 = 1000000000000.0
function norm(v::Tuple{Float64, Float64})
        sum_of_squares::Float64 = 0.0
        for i in 1:length(v)
                sum_of_squares += (v[i])^2
        end

        return sqrt(sum_of_squares)
end

function norm(v::Float64)
        sum_of_squares::Float64 = 0.0
        for i in 1:length(v)
                sum_of_squares += (v[i])^2
        end

        return sqrt(sum_of_squares)
end

function norm(v::Vector{Float64})
        sum_of_squares::Float64 = 0.0
        for i in 1:length(v)
                sum_of_squares += (v[i])^2
        end

        return sqrt(sum_of_squares)
end

function dot(v1::Tuple{Float64, Float64}, v2::Tuple{Float64, Float64})
	if(length(v1) != length(v2))
	   print("Bruh these vectors ain't got the same length no?\n")
	   return -1
	end

	sum = 0.0
	for i in 1:length(v1)
		sum += v1[i]*v2[i]
	end
	return sum
end

###Function for calculating the intersection between two points
function inter(h1::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}, h2::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}, eps::Float64, inf::Float64)
        #h1 and h2 represent the half planes we want to calculate the line intersections for
	#print("Calculating the intersection for $(h1[2]) and $(h2[2])\n")
	#m1 = h1[2][2]/h1[2][1]
	#m2 = h2[2][2]/h2[2][1]
	m1::Float64 = 0.0
	m2::Float64 = 0.0
	if(abs(abs(h1[1]) - pi/2) < 0.000001)
		#print("Infinite gradient detected for m1\n")
		m1 = inf
	else 
		m1 = h1[2][2]/h1[2][1]
	end
	#print("m1 found to be $m1\n")
	if(abs(abs(h2[1]) - pi/2) < 0.000001)
		m2 = inf
		#print("Infinite gradient detected for m2\n")
	else
		m2 = h2[2][2]/h2[2][1]
	end	
	if(abs(m1 - m2) < abs(eps))
		print("Parallel planes yo\n")
		#exit()
		return ((-1.0, -1.0), -1)
	end  
	##print("m1 - m2 found to be $(m1-m2)\n")
	c1::Float64 = h1[3][2] - m1*h1[3][1]
	c2::Float64 = h2[3][2] - m2*h2[3][1]
        xint::Float64 = (c2-c1)/(m1-m2)
	yint::Float64 = -1.0
	if(abs(m1 - inf) < abs(eps))
		yint = m2 * xint + c2
	else 
		yint = m1 * xint + c1
	end
	#print("Intersect calculated as $([xint, yint])\n")
	return ((xint, yint), 1)
end



###Function for calculating the magnitude of the cross product between two vectors v1 and v2
function cross(v1::Tuple{Float64, Float64}, v2::Tuple{Float64, Float64})
	return v1[1] * v2[2] - v1[2]*v2[1]
end



###Function for calculating whether or not a point lies within a half plane, returning 1 if it does lie outside

function outside(half_plane::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}, point::Tuple{Float64, Float64}, eps::Float64, inf::Float64)
        #What we're doing here is, we calculate the vector from the point on the half plane fence to the point that we're considering (which is the intersection). Call this vector b. The vector of the half plane is a. Now, according to the right hand rule for the calculation of the cross product, the point (intersection) will be to the right of the half-plane if the cross product points in the negative z direction, that is, ax*by-ay*bx < 0.
        return cross(half_plane[2], point .- half_plane[3]) < -eps
end

#print(outside([1, [1,2], [3,4]], [5,6]))



###Function for generating the set of vertices defining the voronoi cell
function voronoi_cell(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}, ri::Tuple{Float64, Float64}, neighbouring_points::Vector{Tuple{Tuple{Float64, Float64}, Int64}}, rho::Float64,eps::Float64, inf::Float64 ,temp_half_planes::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}} = [], vel = (0.0,0.0), relic::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64} = (atan(0.0), (0.0, 0.0), (0.0,0.0), 0))
	#ri represents the position of our agent i for whom we wish to calculate the voronoi cell, neighbouring points should be a vector containing the positions of the neighbouring agents (the positions should also be represented as vectors)
###This is the section for deriving the original voronoi cell
	#Look at each of the neighbours of the agent, and generate the half planes
	#print(outside([1, [1,2], [3,4]], [5,6]))
	
	#deque for the half planes/lines, I mean, technically you could just use Julia vectors with pushfirst and whatnot, but eh
	dq::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}} = []
	point::Tuple{Float64, Float64} = (0.0, 0.0) 
	id::Int64 = 0
	for point_tup in neighbouring_points	
		#Use the half-way point as the point p
		point = point_tup[1]
		id = point_tup[2]
		r_ji = point .- ri
		half_plane_point = 0.5 .* r_ji .+ ri

		#Calculate the appropriate vector pq which lies parallel to the line in a direction such that the inner region is to the left of the vector
		v_jix = -1.0 * (0.5 * r_ji[2])
		v_jiy = 0.5 * r_ji[1] #Hopefully you can see that this is literally just v = [-sin(\theta), \cos(\theta)]
		pq = (v_jix, v_jiy)
		#print("$pq\n")
		angle = atan(v_jiy, v_jix)
		is_box = id #This is just to differentiate between the box and actual line segments later
		half_plane = (angle, pq, half_plane_point, is_box)
		push!(dq, half_plane)
	end

	
	#= This stuff is now going to be done using the voronoi_cell_bounded function in half_plane_bounded.jl
	#Add the half plane that bounds the area to the area in front of the agent
        fw_point = ri
        fw_x = -1.0*(-vel[2])
        fw_y = -vel[1] #you might be confused about the negative sign, remember that this is meant to be the vector of the half plane, which is the vector fo the velocity rotated 90 degrees clockwise. 
        fw_pq = [fw_x, fw_y]
        angle = atan(fw_y, fw_x)
        fw_is_box = 2
        fw_half_plane = (angle, fw_pq, Tuple(ri), fw_is_box)
	push!(dq, fw_half_plane)

	#For the relic version of stemler vision, we also need to retain the relic half plane as a bounding half plane for all sampled positions
	#print("About to push the relic, which is $relic\n")
	#push!(half_planes, relic)
	=#

###
###This is the section where we account for the circle of vision
	#Having found the voronoi cell with the bounded box method, we now account for the fact that we have a bounding circle and not a box, and so get rid of the box line segments first
	#print("\n\n\nCommencing bounded DOD calculation\n")	
	i = 1
	while (i <= length(dq))
		if(dq[i][4]== -1 || norm(dq[i][3] .- ri) >= rho || abs(norm(dq[i][3] .- ri)-rho) < 10^(-12))
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

	sort!(dq)
	
	#print("In the voronoi function, a was modified to a value of $a\n")
	replace_vector(temp_half_planes, dq)
	#Now, go through and start calculating the intersects between the non-redundant lines, but if there is no valid intersect, then use the circle
	vq::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = []
	newdq::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}} = []
	dql::Int64 = length(dq)
	len::Int64 = 0
	vlen::Int64 = 0
	for i in 1:dql
		#print("Single fence agent detected\n")
		m::Float64 = 0.0
		if(dq[i][4] == -1)
			print("Bounding box fence detected\n")
			AgentsIO.save_checkpoint("simulation_save.jld2", model)
			exit()
		end
		if(abs(abs(dq[i][1]) - pi/2) > 0.00001) 
			m = dq[i][2][2]/dq[i][2][1]
		else 
			m = inf
		end
                #print("Gradient for i is $m \n")
		x1::Float64 = 0.0
		y1::Float64 = 0.0
		x2::Float64 = 0.0
		y2::Float64 = 0.0

		if(abs(m-inf) < eps)
			x1 = dq[i][3][1]
			x2 = dq[i][3][1]
			if(rho^2 - (x1 - ri[1])^2 <0)
				print("Circle intercept negative. This was for the infinite gradient case, with half plane $(dq[i]). The value of ri was $ri\n" )
				AgentsIO.save_checkpoint("simulation_save.jld2", model)
				exit()
			end
			y1 = -sqrt(rho^2 - (x1 - ri[1])^2) + ri[2]
			y2 = sqrt(rho^2 - (x2 - ri[1])^2) + ri[2]
		else
			c::Float64 = dq[i][3][2] - m*dq[i][3][1]
                	#print("c for $i is $c\n")
                	a::Float64 = 1+m^2
                	b::Float64 = -2*ri[1]+2*m*c-2*m*ri[2]
                	d::Float64 = ri[1]^2 + c^2 - 2*ri[2]*c + ri[2]^2 - rho^2
			discriminant::Float64 = (b)^2 - 4*(a)*(d)
			if((b)^2 - 4*(a)*(d) < 0)
				if(discriminant/b^2 < eps*10)
					discriminant = 0.0
				else
					print("Circle intercept negative. This was for the normal case for the half plane $(dq[i]). The value of the discriminant was $((b)^2 - 4*(a)*(d))\n")
					AgentsIO.save_checkpoint("simulation_save.jld2", model)
					exit()
				end
			end
			x1 = (-(b) - sqrt(discriminant))/(2*(a))
                	y1 = m*x1 + c

                	x2 = (-(b) + sqrt(discriminant))/(2*(a))
                	y2 = m*x2 + c
		end
		#print("x1, y1, x2, y2 calculated to be $x1, $y1, $x2, $y2\n")
		
		a1_a2 = (x1, y1) .- (x2, y2) #Calculate the vector from the second to first intersect with the circle
                f_circle_intersect_i::Tuple{Float64, Float64} = dot(a1_a2, dq[i][2]) >=  0.0 ? (x1, y1) : (x2, y2) #This is to see we of the intersects is right, by testing if the first needs us to move "backwards" from the last vertex 
		#print("f_circle was selected to be $f_circle_intersect_i because a1_a2 was $a1_a2, the dequeue vector was $(dq[i][2]) resulting in a dot product of $(dot(a1_a2, dq[i][2]))\n")
		b_circle_intersect_i::Tuple{Float64, Float64} = dot(a1_a2, dq[i][2]) <=  0.0 ? (x1, y1) : (x2, y2)


		while(vlen >= 1 && outside(dq[i], vq[vlen][1], eps, inf))
			if(vq[vlen][3] != 0)
				if(len <= 0)
                                        print("Tried to delete a half plane from dequeue when there wasn't one\n")
                                        AgentsIO.save_checkpoint("simulation_save.jld2", model)
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
		while(vlen >= 1 && outside(dq[i], vq[1][1], eps, inf))
			if(vq[1][2] != 0)
			if(len <= 0)
					print("Tried to delete a half plane from dequeue when there wasn't one\n")
					AgentsIO.save_checkpoint("simulation_save.jld2", model)
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
		if(len > 0 && norm(cross(dq[i][2], newdq[len][2]))< eps && outside(newdq[len], dq[i][3], eps, inf)) #Check if parallel by if the cross product is less than eps. Note that norm also works on scalars (returns abs val)
                        continue
                end

		invalid_half_plane::Int64 = 0
		if (len >= 1)
			#Determine the intersect of hp_i with hp_(i-1)
			#print("\nThe time for collecting an intersection was ")
			intersect_info::Tuple{Tuple{Float64, Float64}, Int64} =  inter(dq[i], newdq[len], eps, inf)
			intersect_i::Tuple{Float64, Float64} = intersect_info[1]
			is_outside::Int64 = 0
			invalid::Int64 = 0
			if(intersect_info[2] == -1 || norm(intersect_i .- ri) > rho)
				is_outside = 1	
			end

			forward_after_inter = intersect_i .+ 1.0 .* dq[i][2]
			if(intersect_info[2] == -1 || outside(newdq[len], forward_after_inter, eps, inf))
				invalid = 1
			end

			if(is_outside == 0 && invalid == 0)
				push!(vq, (intersect_i, dq[i-1][4], dq[i][4]))
				vlen += 1
				#print("Normal intersect pushed for i = $i. Intersect was $intersect_i\n")
			elseif(!outside(newdq[len], b_circle_intersect_i, eps, inf))
				if(dq[i][4] == -1)
					push!(vq, (b_circle_intersect_i, 0, -1))
				else
					push!(vq, (b_circle_intersect_i, 0, dq[i][4]))
				end
				vlen += 1
				#print("Circle intersect pushed for i = $i. Intersect was $b_circle_intersect_i\n")
			else 
				invalid_half_plane = 1
			end	

		else #Note that we've changed the format to account for what to label the forwards half plane is if it was the artifical bounding half plane
			if(dq[i][4] == 2)
                                        push!(vq, (b_circle_intersect_i, 0, -1))
                        else
                                        push!(vq, (b_circle_intersect_i, 0, dq[i][4]))                                
			end

			vlen += 1
			#print("Circle intersect pushed for i = $i. Intersect was $b_circle_intersect_i\n")
		end
		
		if(invalid_half_plane == 1)
			continue
		end

		#Add the foward intersect
		if(dq[i][4] == 2)
			push!(vq, (f_circle_intersect_i, -1, 0))
                else
                       	push!(vq, (f_circle_intersect_i, dq[i][4], 0))
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
	while(vlen >= 2 && outside(newdq[1], vq[vlen][1], eps, inf))

		if(vq[vlen][3] != 0)
			#print("Popping from the back of the newdq because the intersect was $(vq[vlen]). The back is $(newdq[len])\n")
			if(len <= 0)
                                print("Tried to delete a half plane from dequeue when there wasn't one\n")
                        	AgentsIO.save_checkpoint("simulation_save.jld2", model)
				exit()
                        end

			pop!(newdq)
                	len -= 1
		end
		pop!(vq)
		vlen -= 1

        end

	while(vlen >= 2 && outside(newdq[len], vq[1][1], eps, inf))
		if(vq[1][2] != 0)
			#print("Popping from the front of newdq. The front is $(vq[1])\n")
			if(len <= 0)
                                print("Tried to delete a half plane from dequeue when there wasn't one\n")
                        	AgentsIO.save_checkpoint("simulation_save.jld2", model)
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
                        intersect_last_info::Tuple{Tuple{Float64, Float64}, Int64} = inter(newdq[len], newdq[1], eps, inf)
                        intersect_last::Tuple{Float64, Float64} = intersect_last_info[1]
			is_outside = 0
                        invalid = 0
			if(norm(intersect_last .- ri) > rho)
                                is_outside = 1
                        end

			forward_after_inter = intersect_last .+ 1.0 .* newdq[1][2]
                        if(intersect_last_info[2] == -1 || outside(newdq[len], forward_after_inter, eps, inf))
                                invalid = 1
                        end

                        if(is_outside == 0 && invalid == 0)
                                push!(vq, (intersect_last, newdq[len][4], 1))
                                vlen += 1
                        end
	end

	if(length(vq) == 2 && norm(vq[1][1] .- vq[2][1]) < eps)
		#print("Single intersect calculated (area of 0). The dq which generated this was $dq\n")
	end

	###Check that the cell is indeed the intersection of all half planes
	flag::Int32 = intersect_check(vq, dq)	
	
	return vq

end

function replace_vector(a, b)
	for i in 1:length(a)
		pop!(a)
	end

	for i in 1:length(b)
		push!(a, b[i])
	end
end

#Square test case
#neighbouring_positions = [[60.0, 50.0], [50.0, 60.0], [40.0, 50.0], [50.0, 40.0]]

#Triangle test case
#agent_pos = [50.0, 50.0]
#neighbouring_positions = [[55.0, 55.0], [45.0, 55.0], [50.0, 45.0]]
#

#Weird triangle test case
#neighbouring_positions = [[55.0, 55.0], [45.0, 55.0], [55.0, 45.0]]


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
