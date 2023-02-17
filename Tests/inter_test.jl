const eps = 0.0000000001
const inf = 1000000000000.0

using Random

function inter(h1, h2)
        #h1 and h2 represent the half planes we want to calculate the line intersections for
        #print("Calculating the intersection for $(h1[2]) and $(h2[2])\n")
        #m1 = h1[2][2]/h1[2][1]
        #m2 = h2[2][2]/h2[2][1]
        #m1 = 0.0
        #m2 = 0.0
        if(abs(h1[2][1]) < eps)
                #print("Infinite gradient detected for m1\n")
                m1 = ((h1[2][1] < 0 && h1[2][2] < 0) || (h1[2][1] > 0 && h1[2][1] > 0)) ? inf : -inf
        else
                m1 = h1[2][2]/h1[2][1]
        end
        #print("m1 found to be $m1\n")
        if(abs(h2[2][1]) < eps)
                m2 = ((h2[2][1] < 0 && h2[2][2] < 0) || (h2[2][1] > 0 && h2[2][1] > 0)) ? inf : -inf
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

function main()
	#Warm up the intersect function first
	half_plane_a = [pi/2, [0.0, 5.0], [5.0, 0.0], 0]
	half_plane_b = [pi, [5.0, 0.0], [0.0, 5.0], 0]
	a = @time inter(half_plane_a, half_plane_b)
	print("The intersect of the warm up half planes is $(a)\n")

	#Generate three positions, put them into an array. 
	initial_pos = []
	for i in 1:3
		position = 100*(rand(Float64, 2)) .+ (50.0, 50.0)
		push!(initial_pos, position)
	end
		
	#Using the first position as a reference, generate the half planes according to our formula in the actual simulations
	half_planes_holder = []
	for j in 2:3
		ri = initial_pos[1]
		#Use the half-way point as the point p
		r_ji = initial_pos[j] .- ri
                half_plane_point = 0.5 .* r_ji .+ ri

                #Calculate the appropriate vector pq which lies parallel to the line in a direction such that the inner region is to the left of the vector
                v_jix = -1.0 * (0.5 * r_ji[2])
                v_jiy = 0.5 * r_ji[1] #Hopefully you can see that this is literally just v = [-sin(\theta), \cos(\theta)]
                pq = [v_jix, v_jiy]
                #print("$pq\n")
                angle = atan(v_jiy, v_jix)
                is_box = 0 #This is just to differentiate between the box and actual line segments later
                half_plane = [angle, pq, Tuple(half_plane_point), is_box]	
		push!(half_planes_holder, half_plane)
	end

	print("The time for intersect to be calculated was ") 
	@time inter(half_planes_holder[1], half_planes_holder[2])
	#Send the half planes to the intercept function for testing!
end

main()

half_plane_a = [pi/2, [0.0, 5.0], [5.0, 0.0], 0]
        half_plane_b = [pi, [5.0, 0.0], [0.0, 5.0], 0]
        a = @time inter(half_plane_a, half_plane_b)
        print("The intersect of the warm up half planes is $(a)\n")

        #Generate three positions, put them into an array.
        initial_pos = []
        for i in 1:3
                position = 100*(rand(Float64, 2)) .+ (50.0, 50.0)
                push!(initial_pos, position)
        end

        #Using the first position as a reference, generate the half planes according to our formula in the actual simulations
        half_planes_holder = []
        for j in 2:3
                ri = initial_pos[1]
                #Use the half-way point as the point p
                r_ji = initial_pos[j] .- ri
                half_plane_point = 0.5 .* r_ji .+ ri

                #Calculate the appropriate vector pq which lies parallel to the line in a direction such that the inner region is to the left of the vector
                v_jix = -1.0 * (0.5 * r_ji[2])
                v_jiy = 0.5 * r_ji[1] #Hopefully you can see that this is literally just v = [-sin(\theta), \cos(\theta)]
                pq = [v_jix, v_jiy]
                #print("$pq\n")
                angle = atan(v_jiy, v_jix)
                is_box = 0 #This is just to differentiate between the box and actual line segments later
                half_plane = [angle, pq, Tuple(half_plane_point), is_box]
                push!(half_planes_holder, half_plane)
        end

        print("The time for intersect to be calculated was ")
        @time inter(half_planes_holder[1], half_planes_holder[2])

