function plot_ave_across_df(df, column)
        dude = combine(groupby(df, column), :rot_o => mean)
        data = Vector{Tuple{Float64, Float64}}(undef, 0)
        for i in 1:(size(df, 1)-1)
                push!(data, (dude[i, 1], dude[i, 2]))
        end
        fig = Makie.lines(data)
end 

function plot_filtered_ave_across_df(df, column)
	dude = combine(groupby(df, :left_bias_arg), :rot_o => mean)
        data = Vector{Tuple{Float64, Float64}}(undef, 0)
        for i in 1:size(df, 1)
                push!(data, (dude[i, 1], dude[i, 2]))
        end
        fig = Makie.lines(data)
end

function return_filtered_data(df)
	bro = filter(df -> df.step >= 55000)
	return bro
end

###This is what's used for finding the mean rot_o for long term cluster size vs rot_o simulations. 
function return_filtered_rot_o_mean(df)
	fdf = filter(row -> row.step >= 55000, df)
	com_df = combine(groupby(fdf, :no_bird), :rot_o => mean)
	return com_df[1, 1]
end
