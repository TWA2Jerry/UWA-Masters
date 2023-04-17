using Random
using InteractiveDynamics
using CairoMakie

points = [Tuple(rand(Float64, 2)) for i in 1:10]

figure  = Makie.scatter(points)
for i in 1:length(points)
	text!(points[i], text=  "$i", align = (:center, :top))
end



save("./Makie_test.png", figure)
