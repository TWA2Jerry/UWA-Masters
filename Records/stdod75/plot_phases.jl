using CairoMakie
using LaTeXStrings
using CSV
using DataFrames
include("../../prog.h")

f = Figure(size= (2000, 300))
main_figures = f[1,1] = GridLayout()

offset = 5/4

main_axis = Axis(f[1,1],
	width = 800,
	xscale = log10,
	xlabelsize = 50,
	xlabel = "tDOD",
	ylabelsize = 50, 
	ylabel = "Rotational Order",
	xticklabelsize = 30,
	yticklabelsize = 30,
	xticks = ([1*sqrt(12), 10*sqrt(12),  100*sqrt(12), 1000*sqrt(12)], [L"1\times\sqrt{12}", L"10\times\sqrt{12}", L"100\times\sqrt{12}", L"1000\times\sqrt{12}"]),
	limits = ((1*sqrt(12)/offset, 20000*sqrt(12)/offset), nothing),
	xgridvisible = false, 
	ygridvisible = false	
)
#rowsize!(main_figures, 2, Fixed(400))
#hidedecorations!(main_axis)
#hidespines!(main_axis)

###Process rot_o data 
df1 = DataFrame(CSV.File("../stdod66/stdod66_combined_mdata"))
df2 = DataFrame(CSV.File("../stdod71/stdod71_combined_mdata"))
df3 = DataFrame(CSV.File("../stdod73/stdod73_combined_mdata"))

Makie.lines!(main_axis, 
		df1[:, :target_area_arg], 
		df1[:, :rot_o_alt_mean_mean],
		color = :red,
		linewidth = 3)
errorbars!(main_axis, df1[:, :target_area_arg], df1[:, :rot_o_alt_mean_mean], df1[:, :rot_o_alt_mean_std])

Makie.lines!(main_axis, 
		df2[:, :target_area_arg], 
		df2[:, :rot_o_alt_mean_mean],
		color = :blue,
		linewidth = 3)
errorbars!(main_axis, df2[:, :target_area_arg], df2[:, :rot_o_alt_mean_mean], df2[:, :rot_o_alt_mean_std])
hidespines!(main_axis, :t, :r)
#Makie.xlims!(10*sqrt(12), 10000*sqrt(12))

Makie.lines!(main_axis, 
	df3[:, :target_area_arg], 
	df3[:, :rot_o_alt_mean_mean],
	color = :grey,
	linewidth = 3)
errorbars!(main_axis, df3[:, :target_area_arg], df3[:, :rot_o_alt_mean_mean], df3[:, :rot_o_alt_mean_std])

resize_to_layout!(f)
