using LaTeXStrings
fig, ax = Makie.lines(mdf[:, 6], axis =(xticklabelsize = 30, yticklabelsize=30, xlabelsize = 40, ylabelsize=40, ylabel = L"\textrm{Rotational order } (\Phi_{R*})", xlabel = "Step", xgridvisible = false, ygridvisible = false))
