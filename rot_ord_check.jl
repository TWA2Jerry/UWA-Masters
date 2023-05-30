function rot_ord_check(origin, points)
	all_good = 1
	flicked = 0
	vec_angles = []
	for i in 1:length(points)
		vec_to_point = [points[i][1][1] - origin[1], points[i][1][2] - origin[2]]
		angle_of_vec = atan(vec_to_point[2], vec_to_point[1])
		push!(vec_angles, angle_of_vec)
		if(i != 1 && vec_angles[i] < vec_angles[i-1] && norm((points[i][1][1] - points[i-1][1][1], points[i][1][2]-points[i-1][1][2])) > 10^(-9))
			if(flicked == 0)
				flicked = 1
			else
				all_good = 0
			end
		end
	end
	return all_good
end
