function polarisation(model)
	sum_of_vels::Tuple{Float64, Float64} = (0.0, 0.0)
	for i::Int64 in 1:no_birds
		sum_of_vels = sum_of_vels .+ 1/no_birds .* model[i].vel	
	end

	polarisation_order::Float64 = norm(sum_of_vels)
	return polarisation_order
end
