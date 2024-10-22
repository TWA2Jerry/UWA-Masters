function calculate_entropy(values)
	#Just for unhappiness
	no_bins = 100
	bins::Array{Int64} = zeros(Int64, no_bins)
	bin_width = 1.0/no_bins
	no_unhappiness = 0
	for unhappiness in values
		bin = Int64(trunc(unhappiness/bin_width))+1
		bins[bin] += 1
		no_unhappiness += 1
	end	
	#print(no_unhappiness)
	entropy = 0.0
	for i in 1:no_bins
		P_i = bins[i]/no_unhappiness 
		if(isnan(P_i)) print("Bro\n") end
		if(P_i < 0.0) print("Dude\n") end
		entropy += P_i < eps ? 0 : -P_i * log(P_i)
	end

	return entropy
end
