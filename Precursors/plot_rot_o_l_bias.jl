using Plots
using CairoMakie

rot_o_l_bias_file = open("Records/mv60/mv60.txt", "r")
l_bias_rot_o_lines = readlines(rot_o_l_bias_file)

bias_rot_o_array = Vector{Tuple{Float64, Float64}}(undef, 0)

for line in l_bias_rot_o_lines
	split_line = parse.(Float64, split(line, " "))
	push!(bias_rot_o_array, (split_line[1], split_line[2]))
end

#Makie.plot(bias_rot_o_array)
fig = Plots.plot(bias_rot_o_array)
