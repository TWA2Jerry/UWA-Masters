function cross(v1, v2)
        return v1[1] * v2[2] - v1[2]*v2[1]
end

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
	sort(P) 
	
	for i in 1:n
		while (k >= 2 && cross(H[k-1] .- H[k-2], P[i] .- H[k-2]) <= 0)
			k -= 1
			pop!(H)
		end
		push!(H, P[i])
		k += 1
	end
	t = k+1
	for i in reverse(1:n)
		while (k >= t && cross(H[k-1] .- H[k-2], P[i] .- H[k-2]) <= 0)
                        k -= 1
                        pop!(H)
                end
                push!(H, P[i])
                k += 1
	end

	return H
end

