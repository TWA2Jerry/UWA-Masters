using CairoMakie
using Plots

function give_line_plot(data_points; label_arg = "", xlabel_arg = "", ylabel_arg = "", linewidth_arg = 10, labelvisible_arg = true, limits_arg = (nothing, nothing), figure_padding_arg = (:auto, :auto, :auto, :auto))
	fig, ax, thing = Makie.lines(data_points,
		linewidth = linewidth_arg,
		axis = (;
			limits = limits_arg,
			#width = 900,
			xlabelsize = 50,
			ylabelsize = 50,
			xlabel = xlabel_arg,
			ylabel = ylabel_arg,
			xticklabelsize = 30,
			yticklabelsize = 30,
			xgridvisible = false,
			ygridvisible = false,
			#xlabelpadding = 2.0
		),

		label = label_arg,
		#figure_padding = (:auto, :auto, :auto, :auto)
	)

	if(labelvisible_arg == true)
	axislegend(ax,
		labelsize = 30
	)
	end
	return fig, ax, thing
end

function give_line_plot!(data_points; label_arg = "", xlabel_arg = "", ylabel_arg = "")
    Makie.lines!(data_points,
        linewidth = 10,
		label = label_arg
    )

    return
end

