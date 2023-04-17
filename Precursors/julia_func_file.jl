function norm(v)
	sum_of_squares = 0.0
	for i in 1:length(v)
		sum_of_squares += v[i]
	end
	return sqrt(sum_of_squares)
end
