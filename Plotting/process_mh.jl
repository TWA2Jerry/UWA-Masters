###Just reading in teh mean happinesses for a given simulation
mh_file = open("../adf_file.txt", "r")
wall_file = open("mh_wall.txt", "w")

###Get ready to store and process the mean happinesses
mh_array = []
mh_lines = readlines(mh_file)
for line in mh_lines
	split_line = parse.(Float64, split(line, " "))
        push!(mh_array, split_line[3])
end

close(mh_file)

bin_width = 0.005
function bin_number(x)
	return Int64(floor(x/bin_width))
end

function rounded(x)
	bin_width * (bin_number(x) + 0.5)
end

max_bin::Int64 = 0
for i in 500000:length(mh_array)
	global max_bin = bin_number(mh_array[i]) > max_bin ? bin_number(mh_array[i]) : max_bin	
end

frequency = zeros(max_bin+1) #This creates a vector of zeroes representing the frequency of each bin of unhappiness 
for i in 1:length(mh_array)
	frequency[bin_number(mh_array[i])+1] += 1
end

###Now, see where the "wall" or severe drop in unhappiness is. 
max_change=0
max_change_bin = -1
for i in 2:length(frequency)
	global max_change
	global max_change_bin
	if(max_change < abs(frequency[i]-frequency[i-1]))
		max_change = abs(frequency[i]-frequency[i-1])
		max_change_bin = (i-1)*bin_width + 0.5*bin_width
	end
end

write(wall_file, "$max_change_bin\n")
close(wall_file)
