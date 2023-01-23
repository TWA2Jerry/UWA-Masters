include("./half_plane_alt.jl")
#=
function cross(v1, v2)
        return v1[1] * v2[2] - v1[2]*v2[1]
end
=#

###Function that calculates the convex hull of a set of points P
function convex_hull(P)
	#Create the vector H which will contain the points of the convex hull
	H = []

	#k will represent the number of points currently in the Convex Hull
	k = 0

	n = length(P)
	print("The length of P is $n\n")
	if(n <= 3) 
		return P
	end
	
	#Sort the list of points, pretty sure that Julia automatically checks second coordinate if first are equal, and sorts in ascending order
	sort!(P) 
	
	print("After sorting, P, the list of all points, is given by\n")
	print(P)
	print("Now commencing calculation\n")
	
	for i in 1:n
		while (k >= 2 && cross(H[k] .- H[k-1], P[i] .- H[k-1]) <= 0)
			k -= 1
			pop!(H)
			print("Popped point. ")
		end
		push!(H, P[i])
		print("Added point $(P[i]), Hull is now $H\n")
		k += 1
	end
	t = k+1 #As far as I'm aware, t is a marker between the construction of the lower and upper hull so that we don't remove points from the already constructed lower hull.
	for i in reverse(1:n)
		if(i < 2)
			break
		end
		while (k >= t && cross(H[k] .- H[k-1], P[i-1] .- H[k-1]) <= 0)
                        k -= 1
                        pop!(H)
			print("Popped point. ")
                end
                push!(H, P[i-1])
		print("Added point $(P[i-1]), Hull is now $H\n")
                k += 1
	end

	#Note that as far as I'm aware, the last point of the Hull vector H will be the same as the first one, I don't think this is a problem though. 
	print("Convex hull was calculated to be $H\n")
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
