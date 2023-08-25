using VoronoiCells
using GeometryBasics
no_birds::Int64 =100
const rect_bound::Float64 = 1000.0
tracked_agent::Int64 = 1
R::Float64 = 100.0
rect = Rectangle(Point2(0,0), Point2(Int64(rect_bound), Int64(rect_bound)))
pack_positions = Vector{Point2{Float64}}(undef, no_birds)
	for i in 1:no_birds
		angle_per_bird = 2*pi/(no_birds-1)
                initial_pos = i != tracked_agent ? (R*cos(angle_per_bird*(i-1)), R*sin(angle_per_bird*(i-1))) .+ (300.0, 300.0) : (300.0, 300.0)
                rand_vel = (-R*sin(angle_per_bird*(i-1)), R*cos(angle_per_bird*(i-1)))
		pack_positions[i] = initial_pos
	end
#pack_positions[1] = [0.0, 0.0]
#pack_positions[2] = [1.0, 0.0]
 init_tess = voronoicells(pack_positions, rect)
        init_tess_areas = voronoiarea(init_tess)
print("Done\n")
