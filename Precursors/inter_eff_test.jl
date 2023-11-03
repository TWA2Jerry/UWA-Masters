include("../some_math_functions.jl")
const ri = (0.0, 0.0)

function generate_hp()
	point::Tuple{Float64, Float64} = 2 .* Tuple(rand(Float64, 2)) .- (1.0, 1.0)
	r_ji = point .- ri
                half_plane_point = 0.5 .* r_ji .+ ri

                #Calculate the appropriate vector pq which lies parallel to the line in a direction such that the inner region is to the left of the vector
                v_jix = -1.0 * (0.5 * r_ji[2])
                v_jiy = 0.5 * r_ji[1] #Hopefully you can see that this is literally just v = [-sin(\theta), \cos(\theta)]
                pq = (v_jix, v_jiy)
                #print("$pq\n")
                angle = atan(v_jiy, v_jix)
                is_box = rand(Int64) #This is just to differentiate between the box and actual line segments later
                half_plane = (angle, pq, half_plane_point, is_box)
end

#half_plane_1 = (pi/2, (0.0, 1.0), (0.0, 0.0), 1)
#half_plane_2 = (pi/1.0, (1.0, 0.0), (0.0, 1.0), 2)

half_plane_1 = generate_hp()
half_plane_2 = generate_hp()

function test_inter_eff()
	for i in 1:1e8
		bro = inter(half_plane_1, half_plane_2,  eps, inf)
	end
end


@time test_inter_eff()
