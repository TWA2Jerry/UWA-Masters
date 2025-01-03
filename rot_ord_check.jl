function rot_ord_check(origin::Tuple{Float64, Float64}, points::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}})
	all_good::Int64 = 1
	flicked::Int64 = 0
	vec_angles::Vector{Float64} = []
	for i::Int64 in 1:length(points)
		vec_to_point::Tuple{Float64, Float64} = (points[i][1][1] - origin[1], points[i][1][2] - origin[2])
		angle_of_vec::Float64 = atan(vec_to_point[2], vec_to_point[1])
		push!(vec_angles, angle_of_vec)
		if(i != 1 && vec_angles[i] < vec_angles[i-1] && norm((points[i][1][1] - points[i-1][1][1], points[i][1][2]-points[i-1][1][2])) > 10^(-7) && (-points[i][1][1] + points[i-1][1][1])*(points[i][1][2]+points[i-1][1][2]) > 0.1)
			if(flicked == 0)
				flicked = 1
			else
				all_good = 0
			end
		end
	end
	return all_good
end
