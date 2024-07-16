using Plots
using Plots.PlotMeasures
using CairoMakie
using LaTeXStrings

fig = Plots.plot(mdf[:, 2], plot_title = L"Mean Speed($\Phi_S$) vs Time", 
	plot_titlefontsize = 50,
	size = (1800, 1200),
	linewidth = 10,
	legend=false,
	grid = false,
	xguide = "Time Step", xguidefontsize = 60, xtickfontsize= 50, 
	yguide = L"$\Phi_S$", yguidefontsize=60, ytickfontsize= 50,
	#xticks = ([10*sqrt(12), 100*sqrt(12), 1000*sqrt(12)], [10*sqrt(12), 100*sqrt(12), 1000*sqrt(12)]))
	left_margin = 10mm,
	bottom_margin = 10mm,
	dpi = 400
)

fig2 = Plots.plot(mdf[:, 3], plot_title = L"Rotational Order($\Phi_R$) vs Time",
    plot_titlefontsize = 50,
    size = (1800, 1200),
    linewidth = 10,
    legend=false,
    grid = false,
    xguide = "Time Step", xguidefontsize = 60, xtickfontsize= 50,
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
savefig("./stdod52_rot_o_alt_tave2.pdf")
#display(fig)
