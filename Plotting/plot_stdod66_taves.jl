using Plots
using Plots.PlotMeasures
#using CairoMakie
using LaTeXStrings
using CSV
using DataFrames
#ms_taves_file = open("../stdod66_taves.txt", "r")
#tdods_ms_lines = readlines(ms_taves_file)


tdods = []
ms = []

df = DataFrame(CSV.File("../Records/stdod66/stdod66_taves"))

#=
for i in 1:length(tdods_ms_lines)
	line = tdods_ms_lines[i]
	split_line = parse.(Float64, split(line, " "))
	push!(tdods, split_line[1])
	push!(ms, split_line[4])	
end


for i in 1:length(tdods)
	print("$(tdods[i]) $(ms[i])\n")
end
=#

xticks_vec = [10*sqrt(12), 100*sqrt(12), 1000*sqrt(12)]
xticks_string_vec = [L"$10\times\sqrt{12}$", L"$100\times\sqrt{12}$", L"$1000\times\sqrt{12}$"]
fig = Plots.plot(df[:, :target_area_arg], df[:, :rot_o_alt_mean_mean], 
	yerror = df[:, :rot_o_alt_mean_std],
	plot_title = L"Rotational Order($\Phi_R$) vs tDOD", 
	plot_titlefontsize = 50,
	size = (1800, 1200),
	linewidth = 10,
	legend=false,
	grid = false,
	xguide = "tDOD", xlims = (10*sqrt(12), 20000), xscale = :log10, xticks = (xticks_vec, xticks_string_vec), xguidefontsize = 60, xtickfontsize= 50, 
	yguide = L"$\Phi_R$", yguidefontsize=60, ytickfontsize= 50,
	#xticks = ([10*sqrt(12), 100*sqrt(12), 1000*sqrt(12)], [10*sqrt(12), 100*sqrt(12), 1000*sqrt(12)]))
	left_margin = 10mm,
	bottom_margin = 10mm,
	dpi = 400
)

#=display(Plots.plot(tdods, ms, plot_title = L"Mean speed ($\Phi_S$) vs tDOD", plot_titlefontsize = 40,
        xaxis = (xlabel = "tDOD")
)
)=#
savefig("../Records/stdod66_rot_o_alt_tave.pdf")
#display(fig)
