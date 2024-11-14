using Plots
using CairoMakie
include("line_plot_template.jl")

data_file = open("../Records/mv118/mv118_prop_one_cluster_size.txt", "r")
data_lines = readlines(data_file)

data_array = Vector{Tuple{Float64, Float64}}(undef, 0)

for line in data_lines
        split_line = parse.(Float64, split(line, " "))
		push!(data_array, (split_line[1], split_line[2]))	
end

fig, ax = give_line_plot(data_array, labelvisible_arg = false, xlabel_arg = "Cluster Size", ylabel_arg = "Stability")
