using CairoMakie
using Plots

#=
function plot_model(model; marker_arg = :circle, markersize_arg = 30)
	model_positions = Vector{Tuple{Float64, Float64}}(undef, 0)
	model_rotations = Vector{Float64}(undef, 0)
	for i in 1:nagents(model)
		push!(model_positions, model[i].pos)
		push!(model_rotations, atan(model[i].vel[2], model[i].vel[1]))
	end

	fig, ax, thing = Makie.scatter(model_positions, 
		axis = (; 
			limits = (0, rect_bound, 0, rect_bound), 
			aspect = 1,
			xticklabelsize = 30,
			yticklabelsize = 30,
			xgridvisible = false, 
			ygridvisible = false
		),
		marker = marker_arg,
		markersize = markersize_arg,
		rotations = model_rotations,
		color = :black,
	)
	
	return fig, ax, thing	
end
=#

function plot_model(model; marker_arg = :circle, markersize_arg = 30)
    model_positions = Vector{Tuple{Float64, Float64}}(undef, 0)
    model_rotations = Vector{Float64}(undef, 0)
    for i in 1:nagents(model)
        push!(model_positions, model[i].pos)
        push!(model_rotations, atan(model[i].vel[2], model[i].vel[1]))
    end

	
	fig = Figure(size = (100, 100), figure_padding = (20, 50, 20, 50))
	
	ax = Axis(fig[1,1],
		limits = (0, rect_bound, 0, rect_bound), 
        aspect = 1,
        xticklabelsize = 30,
        yticklabelsize = 30,
        xgridvisible = false, 
        ygridvisible = false,
	)

	Makie.scatter!(ax, model_positions, 
		size = (100, 100),
		marker = marker_arg,
        markersize = markersize_arg,
        rotations = model_rotations,
        color = :black,
	)
	
	colsize!(fig.layout, 1, Aspect(1, 1.0))
    resize_to_layout!(fig)


	return fig, ax
end

