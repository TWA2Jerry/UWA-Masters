using CairoMakie
using LaTeXStrings
using CSV
using DataFrames
include("../../prog.h")

f = Figure(size= (2400, 300))
main_figures = f[1,1] = GridLayout()

offset = 5/4

ax1 = Axis(f[1,1],
	width = 1000,
	xlabelsize = 50,
	xlabel = "Time Step",
	ylabelsize = 50, 
	ylabel = L"\Phi_R",
	xticklabelsize = 30,
	yticklabelsize = 30,
	yticks = ([0.6, 0.3, 0.0], ["0.6", "0.3", "0.0"]),
	limits = ((0, 10001), nothing),
	xgridvisible = false, 
	ygridvisible = false	
)

ax2 = Axis(f[1,1],
    width = 1000,
    ylabelsize = 50,
    ylabel = L"\Phi_S",
    yticklabelsize = 30,
    yaxisposition = :right,
	limits = ((0, 10001), (0.0, 1.0)),
	xgridvisible = false,
	ygridvisible = false,
	xticklabelsvisible = false
)

#rowsize!(main_figures, 2, Fixed(400))
#hidedecorations!(main_axis)
#hidespines!(main_axis)
hidexdecorations!(ax2, grid=false, ticks = false, ticklabels = false)
###Process rot_o data 
df1 = DataFrame(CSV.File("./stdod79_22000_mdata")) #Rot o data
Makie.lines!(ax1, 
		df1[:, :step], 
		df1[:, :rot_o_alt],
		color = :red,
		linewidth = 5)
#errorbars!(main_axis, df1[:, :target_area_arg], df1[:, :rot_o_alt_mean_mean], df1[:, :rot_o_alt_mean_std])

Makie.lines!(ax2, 
		df1[:, :step], 
		df1[:, :mean_speed],
		color = :blue,
		linewidth = 20)
#errorbars!(main_axis, df2[:, :target_area_arg], df2[:, :rot_o_alt_mean_mean], df2[:, :rot_o_alt_mean_std])
hidespines!(ax1, :t)
hidespines!(ax2, :t)
#Makie.xlims!(10*sqrt(12), 10000*sqrt(12))


#axislegend()
resize_to_layout!(f)
