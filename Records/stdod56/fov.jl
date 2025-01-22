using CairoMakie
include("../../prog.h")
include("../../draw_circle_part.jl")

model = initialise()

o = model[1].pos
r = 100.0
theta1 = 0.0
theta2 = Float64(pi)

sl1 = draw_straight_line(o, r, theta1)
sl2 = draw_straight_line(o, r, theta2)
cseg = circle_seg(o, r, theta1, theta2)

sl1p::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
sl2p::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
csegp::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)

for i in 1:length(sl1[1])
       push!(sl1p, (sl1[1][i], sl1[2][i]))
       end

for i in 1:length(sl2[1])
       push!(sl2p, (sl2[1][i], sl2[2][i]))
       end

for i in 1:length(cseg[1])
       push!(csegp, (cseg[1][i], cseg[2][i]))
       end

append!(sl1p, csegp)
append!(sl1p, sl2p) 

poly(sl1p, color = :red, alpha = 0.2)

