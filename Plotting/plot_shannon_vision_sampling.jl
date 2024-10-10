using Plots
using Plots.PlotMeasures
using CairoMakie
using LaTeXStrings

points = Vector{Tuple{Float64, Float64}}(undef, 0)
colours = []
	
push!(points, (0, 0))
push!(colours, :black)

for angle in -4:4
	for radius in 1:3
		push!(points, (cos(angle * 2*pi/8)*radius, sin(angle*2*pi/8)*radius))
		if(abs(angle) < 2) push!(colours, :green)
		else push!(colours, :orange) end
		#push!(colours, :green)
	end	
end


fig = Figure()
ax = Axis(fig[1,1]; 
	aspect = 1,
	xticklabelsize = 30,
	yticklabelsize = 30, 
	xgridvisible = false, 
    ygridvisible = false,	
)

Makie.scatter!(ax, 
	points,
	marker = :circle,
	markersize = 40,
	color = colours 
)

Makie.arrows!(ax, 
	[0],
	[0],
	[0.7],
	[0],
	linewidth = 5,
	arrowsize = 20,
)

Makie.text!(ax, 
	(0.15, 0.15),
	text = L"\hat{v}_i",
	align = (:left, :bottom),
	fontsize = 50
)

Makie.text!(ax, 
    (-0.15, -0.1),
    text = L"i",
    align = (:right, :top),
    fontsize = 50
)

colsize!(fig.layout, 1, Aspect(1, 1.0))
resize_to_layout!(fig)	
hidedecorations!(ax)

