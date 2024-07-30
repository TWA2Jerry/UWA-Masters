using Plots
using CairoMakie

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

pop!(size_array)
pop!(rot_o_array)

#Makie.plot(bias_rot_o_array)
#fig = Plots.plot(bias_rot_o_array)
#fig = Makie.scatter(bias_rot_o_array)
fig = Plots.plot(size_array, rot_o_array,
	xlabel = "Cluster size",	
	ylabel = "Average group rot o"
) 
