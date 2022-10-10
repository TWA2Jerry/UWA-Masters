using VoronoiCells
using GeometryBasics
using Plots
using Random

rng = Random.MersenneTwister(123)

rect = Rectangle(Point2(0, 0), Point2(100, 100))
points = [Point2(100*rand(Float64, 2)) for _ in 1:10]
tess = voronoicells(points, rect)

print(tess.Cells[1])
scatter(points, markersize = 6, label = "generators")
annotate!([(points[n][1] + 0.02, points[n][2] + 0.03, Plots.text(n)) for n in 1:10])
display(plot!(tess, legend=:topleft))
savefig("voronoi_pack_test.png")


