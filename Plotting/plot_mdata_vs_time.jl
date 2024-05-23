using Plots.PlotMeasures
using LaTeXStrings
fig = Plots.plot(mdf[:,5],  labelsize=15, xtickfontsize = 15, ytickfontsize=15, legends=false, ylabel = L"$\Phi_{R*}$", xlabel = "Time step", xguidefontsize = 30, yguidefontsize=30, left_margin = 5mm, bottom_margin=5mm, right_margin=5mm)
