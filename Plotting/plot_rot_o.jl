using LaTeXStrings
fig, ax = Makie.lines(mdf[:, 5], axis =(xticklabelsize = 30, yticklabelsize=30, xlabelsize = 40, ylabelsize=40, ylabel = L"\textrm{Rotational order } (\Phi_R)", xlabel = "Step", xgridvisible = false, ygridvisible = false))
