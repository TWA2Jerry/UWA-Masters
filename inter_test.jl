const eps::Float64 = 0.0000000001
const inf::Float64 = 1000000000000.0

using Random

function inter(h1, h2, eps::Float64, inf::Float64)
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
                return -1
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
        return [xint, yint]
end

function inter_opt1(h1::Tuple{Float64, Vector{Float64}, Tuple{Float64, Float64}, Int64}, h2::Tuple{Float64, Vector{Float64}, Tuple{Float64, Float64}, Int64}, eps::Float64, inf::Float64)
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
                return -1
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
        return [xint, yint]
end

function inter_opt2(h1::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}, h2::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}, eps::Float64, inf::Float64)
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
                return -1
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
        return [xint, yint]
end


function run_normal(dq)
	for i in 1:100000000
                inter(dq[1], dq[2], eps, inf)
        end
	return 
end

function run_tupled(dq)
	for i in 1:100000000
                inter_opt1(dq[1], dq[2], eps, inf)
        end
	return 
end	

function run_opt2(dq)
        for i in 1:100000000
                inter_opt2(dq[1], dq[2], eps, inf)
        end
        return 
end  

function run_test()
	ri = (0.0, 0.0)
	neighbouring_points = []
	for i in 1:2
		neighbour_point = Tuple(100*rand(Float64, 2))
		push!(neighbouring_points, neighbour_point)	
	end

	dq = []
	dq_tupled = []	
	dq_opt2 = []
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
                half_plane_tupled = (angle, pq, Tuple(half_plane_point), is_box)
		half_plane_opt2 = (angle, Tuple(pq), Tuple(half_plane_point), is_box)	
		push!(dq, half_plane)
		push!(dq_tupled, half_plane_tupled)
		push!(dq_opt2, half_plane_opt2)
end #end for the for loop over neighbouring points
	
	for i in 1:10
                print("The time for calculating an intersection in the usual unoptimised manner is ")
                @time inter(dq[1], dq[2], eps, inf)
        end

	for i in 1:10
                print("The time for calculating an intersection in the tupled manner is ")
                @time inter_opt1(dq_tupled[1], dq_tupled[2], eps, inf)
        end
	
	print("The time for running lots of the normal is ")
	@time run_normal(dq)

	print("The time for running lots of the tupled is ")
	@time run_tupled(dq_tupled)	
	
	print("The time for running lots of opt2 is ")
	@time run_opt2(dq_opt2)
end #end for the function run test

run_test()
