include("./half_plane_alt.jl")
include("./convex_hull.jl")

bro = [40,40]
bruh = [60,40]
bro2 = [60, 60]
bruh2 = [40, 60]
points = [bro, bruh, bro2, bruh2]

cross(bro, bruh)
convex_hull(points)
