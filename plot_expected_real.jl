include("some_math_functions.jl")
include("global_vars.jl")
using Plots
thingo = open("agent_vals.txt", "r")

lines = readlines(thingo)
parsed_lines = []
for i in 1:length(lines)
        parsed_line = parse.(Float64, split(lines[i], " "))
        push!(parsed_lines, parsed_line)
end

differences_to_plot = []
actuals = []
expecteds = []
no_steps = 5000
for block_number in 0:no_steps-2
	differences = []
	for agent_num in 1:no_birds
		projected_area = parsed_lines[block_number * no_birds + agent_num][5] 
		actual_area = parsed_lines[(block_number+1) * no_birds + agent_num][3]  
		push!(actuals, actual_area)
		push!(expecteds, projected_area)
		push!(differences, projected_area - actual_area)
	end
	mean_difference = mean(differences)
	push!(differences_to_plot, mean_difference)
end

#figure = Plots.plot(differences_to_plot)
#display(figure)

figure = Plots.plot(expecteds, actuals, seriestype=:scatter, label="Expecteds vs actual")
display(figure)

