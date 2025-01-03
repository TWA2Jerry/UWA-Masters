Makie.hist(adf[200000:500100, 3], bins = range(0.0, 1.0, 100), normalization = :probability, axis =(xticklabelsize = 30, yticklabelsize=30, xlabelsize = 40, ylabelsize=40, ylabel = "Relative Frequency", xlabel = "Unhappiness"))

function plot_happiness_hist(filtered_data)
	fig, ax = Makie.hist(filtered_data[:, 3], bins = range(0.0, 1.0, 100), normalization = :probability, axis =(xticklabelsize = 30, yticklabelsize=30, xlabelsize = 40, ylabelsize=40, ylabel = "Relative Frequency", xlabel = "Unhappiness"))
	return fig, ax
end
