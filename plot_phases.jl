using CairoMakie


f = Figure()
main_figures = f[1,1] = GridLayout()
main_axis = Axis(f[1,1])
upper_figs = main_figures[1,1] = GridLayout()
lower_figs = main_figures[2,1]

upper_ax1 = Axis(upper_figs[1,1]; limits = (0, 100, 0, 100), aspect = 1)
upper_ax2 = Axis(upper_figs[1,2]; limits = (0, 1000, 0, 1000), aspect =1)
upper_ax3 = Axis(upper_figs[1,3]; limits = (0, 10000, 0, 10000), aspect = 1)

lower_ax1 = Axis(lower_figs)

#=
hidedecorations!(upper_ax1)
hidespines!(upper_ax1)
hidedecorations!(upper_ax2)
hidespines!(upper_ax2)
hidedecorations!(upper_ax3)
hidespines!(upper_ax3)
=#

hidedecorations!(lower_ax1)
hidespines!(lower_ax1)
