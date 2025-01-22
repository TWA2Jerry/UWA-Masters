function circle_seg(o::Tuple{Float64, Float64}, r::Float64, theta_1::Float64, theta_2::Float64)
	theta::Vector{Float64} = LinRange(theta_1, theta_2, 1000)
	return o[1] .+  r*cos.(theta), o[2] .+  r*sin.(theta)
end

function draw_circle_seg(o::Tuple{Float64, Float64}, r::Float64, theta_1::Float64, theta_2::Float64)
	print("Draw circle called\n")
	circle_points::Tuple{Vector{Float64}, Vector{Float64}} = circle_seg(o, r, theta_1, theta_2)
	Makie.lines!(circle_points[1], circle_points[2], lw = 0.5, color = :black)
end	

function draw_straight_line(o::Tuple{Float64, Float64}, r::Float64, theta::Float64)
	rp::Vector{Float64} = LinRange(0.0, r, 1000)
	return o[1] .+  rp .* cos(theta), o[2] .+  rp .* sin(theta)
end
