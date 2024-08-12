using Plots
using Plots.PlotMeasures
using CairoMakie
using LaTeXStrings

rot_o_l_bias_file = open("Records/mv82/rot_o_v_no_bird.txt", "r")
l_bias_rot_o_lines = readlines(rot_o_l_bias_file)

bias_rot_o_array = Vector{Tuple{Float64, Float64}}(undef, 0)

size_array = []
rot_o_array = []
std_array = []

for line in l_bias_rot_o_lines
        split_line = parse.(Float64, split(line, " "))
        push!(bias_rot_o_array, (split_line[1], split_line[2]))
		push!(size_array, split_line[1])
		push!(rot_o_array, split_line[2])
		#push!(std_array, split_line[3])
end


#Makie.plot(bias_rot_o_array)
#fig = Plots.plot(bias_rot_o_array)
#fig = Makie.scatter(bias_rot_o_array)
fig = Plots.plot(size_array, rot_o_array, plot_title = L"Rotational Order($\Phi_{R*}$) vs Size",
    plot_titlefontsize = 50,
    size = (1800, 1200),
    left_margin = 10mm,
    bottom_margin = 10mm,
    dpi = 400,
    linewidth = 10,
    legend=false,
    grid = false,
    xguide = "Cluster Size", xguidefontsize = 60, xtickfontsize= 50,
    yguide = L"$\Phi_{R*}$", yguidefontsize=60, ytickfontsize= 50
)

