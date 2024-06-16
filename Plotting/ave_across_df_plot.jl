function plot_ave_across_df(df, column)
        dude = combine(groupby(df, :left_bias_arg), :rot_o => mean)
        data = Vector{Tuple{Float64, Float64}}(undef, 0)
        for i in 1:size(df, 1)
                push!(data, (dude[i, 1], dude[i, 2]))
        end
        fig = Makie.lines(data)
end 

