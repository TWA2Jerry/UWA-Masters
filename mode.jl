using StatsBase

thingo = open("agent_vals.txt", "r")

lines =  readlines(thingo)
parsed_lines = []
for i in 1:length(lines)
	parsed_line = parse.(Float64, split(lines[i], " "))
	#print("$(parsed_line[4])\n")
	push!(parsed_lines, parsed_line)
end

distances = [parsed_lines[i][4] for i in 1:length(parsed_lines)]


#print(length(lines))
#print(lines[1][4])
#typeof(lines[1])
print("The most often occurring was $(mode(distances))\n")
