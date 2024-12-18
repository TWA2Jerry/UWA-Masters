include("./half_plane_fast.jl")
include("voronoi_area_file.jl")

function convex_cross(v1, v2)
        return v1[1] * v2[2] - v1[2]*v2[1]
end

###Function for updating the convex hull, returns the points of the convex hull so the area of the convex hull can be calculated
function update_convex_hull(model)
	points = []
	all_agents_iterable = allagents(model)
	for agent in all_agents_iterable
		push!(points, [agent.pos[1], agent.pos[2], agent.id])
		#print("The point added was $([agent.pos[1], agent.pos[2], Int64(agent.id)]), the agent id is $(agent.id)\n")
	end
	
	for i in 1:nagents(model)
		convex_hull_point[i] = 0
	end

	CH = convex_hull(points)
	CH_for_area = []
	
	for point in CH
		convex_hull_point[Int64(point[3])] = 1
		push!(CH_for_area, [[point[1], point[2]], 1, 1]) #Note that the points used for the CH_for_area consists of the actual point itself and 1,1 because the voronoi area calcualtor we're going to use requires a 2 and 3 index, checking them to see if the points are circular. 
	end

	return CH_for_area
end


###Function that calculates the convex hull of a set of points P
function convex_hull(P)
	#Create the vector H which will contain the points of the convex hull
	H = []

	#k will represent the number of points currently in the Convex Hull
	k = 0

	n = length(P)

	if(n <= 3) 
		return P
	end
	
	#Sort the list of points, pretty sure that Julia automatically checks second coordinate if first are equal, and sorts in ascending order
	sort!(P) 
	
	for i in 1:n
		while (k >= 2 && convex_cross(H[k] .- H[k-1], P[i] .- H[k-1]) <= 0)
			k -= 1
			pop!(H)
		end
		push!(H, P[i])
		k += 1
	end
	t = k+1
	for i in reverse(1:n)
		if(i < 2)
			break
		end
		while (k >= t && convex_cross(H[k] .- H[k-1], P[i-1] .- H[k-1]) <= 0)
                        k -= 1
                        pop!(H)
                end
                push!(H, P[i-1])
                k += 1
	end

	return H
end


#=
#Test set, simple square
points = [[40, 40], [60, 40], [60, 60], [40, 60]]
CH = convex_hull(points)
print("The points are given by")
for point in CH
	print(point)
end
=#


function convex_hull_in_area(model, top_left, bottom_right)
	relevant_points = []
	for i in 1:nagents(model)
		if (model[i].pos[1] > top_left[1] && model[i].pos[1] < bottom_right[1] && model[i].pos[2] < top_left[2] && model[i].pos[2] > bottom_right[2])
			push!(relevant_points, (model[i].pos[1], model[i].pos[2], i))
		end
	end

	CH = convex_hull(relevant_points)
        
	#=
	CH_for_area = []
	for point in CH
                convex_hull_point[Int64(point[3])] = 1
        	push!(CH_for_area, [[point[1], point[2]], 1, 1])
	end
	=#
	return CH	
end

function convex_in_area(model, top_left, bottom_right)
	CH  = convex_hull_in_area(model, top_left, bottom_right)
	CH_for_area = []
        for point in CH
                convex_hull_point[Int64(point[3])] = 1
                push!(CH_for_area, [[point[1], point[2]], 1, 1])
        end
	convex_hull_area = voronoi_area(model, -1, CH_for_area, rho)
end

function convex_density(model, top_left, bottom_right)
	convex_area = convex_in_area(model, top_left, bottom_right)
	no_agents_in_area = 0
	for i in 1:nagents(model)
                if (model[i].pos[1] > top_left[1] && model[i].pos[1] < bottom_right[1] && model[i].pos[2] < top_left[2] && model[i].pos[2] > bottom_right[2])
                	no_agents_in_area += 1
		end
        end
	density = convex_area/no_agents_in_area
	print("The number of agents in area was $no_agents_in_area\n")	
	return density
end
