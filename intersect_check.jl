function intersect_check(cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}}, half_planes::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}})
	flag = 0
	#print("Intersect checking\n")
	for i in 1:length(cell)
		vertex_tup = cell[i]
		vertex = vertex_tup[1]
		for j in 1:length(half_planes)
			half_plane = half_planes[j]
			if(outside(half_plane, vertex, inf, eps) == 1 && cell[i][2] != half_plane[4] && cell[i][3] != half_plane[4]) 
				print("yo we got a violation. This is for cell vertex $(vertex) and half plane $(half_plane)\n")
				exit()
				flag = 1
			end
		end 
	end
	return flag
end
