function circle_seg(o, r, theta_1, theta_2)
	theta = LinRange(theta_1, theta_2, 500)
	return o[1] .+  r*cos.(theta), o[2] .+  r*sin.(theta)
end

function draw_circle_seg(o, r, theta_1, theta_2)
	print("Draw circle called\n")
	circle_points = circle_seg(o, r, theta_1, theta_2)
	Plots.plot!(circle_points[1], circle_points[2], lw = 0.5)
end	

