using DelimitedFiles

bro_file = open("bro.txt", "r")
lines = readlines(bro_file)

a = readdlm("bro.txt", ' ', Float64, '\n')
print("a is $a\n")
print("a[1] is $(a[1,3])\n")

print(lines[1])
print(lines[1][1])
print(lines[1][2])
averages_array = []

for i in 1:length(lines[1])
	push!(averages_array, [])
end

for line in lines
	split_line = parse.(Float64, split(line, " "))
	for i in 1:length(split_line)
		push!(averages_array[i], Float64(split_line[i]))
	end
end

print(averages_array)
