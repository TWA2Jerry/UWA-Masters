using CairoMakie
using CSV
using DataFrames

include("../../Plotting/line_plot_template.jl")

data::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
mdf = DataFrame(CSV.File("mv119_ms_rot_o"))

for i in 1:nrow(mdf)
	push!(data, (mdf[i, :step], mdf[i, :rot_o]))
end 

fig, ax = give_line_plot(data, xlabel_arg = "Step", ylabel_arg = L"\Phi_{R*}", labelvisible_arg = false, limits_arg = (nothing, (0, 1.0)), linewidth_arg = 5, figure_padding_arg = (1,1,1,15))

###Plotting of rot_o_alt from stdod model for comparison
stdodmdf  = DataFrame(CSV.File("../stdod63/stdod63_ms_rot_o_dat"))
#lines!(stdodmdf[:, :rot_o], color=:red)
lines!([(0.0,  0.0711), (75000, 0.0711)], color = :red)


