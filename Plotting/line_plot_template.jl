using CairoMakie
using Plots

function give_line_plot(data_points; label_arg = "", xlabel_arg = "", ylabel_arg = "")
	fig, ax, thing = Makie.lines(data_points,
		linewidth = 10,
		axis = (;
			xlabelsize = 50,
			ylabelsize = 50,
			xlabel = xlabel_arg,
			ylabel = ylabel_arg,
			xticklabelsize = 30,
			yticklabelsize = 30,
			xgridvisible = false,
			ygridvisible = false
		),
	
		label = label_arg,
	)

	axislegend(ax,
		labelsize = 30
	)
	return fig, ax, thing
end

function give_line_plot!(data_points; label_arg = "", xlabel_arg = "", ylabel_arg = "")
    Makie.lines!(data_points,
        linewidth = 10,
		label = label_arg
    )

    return
end

