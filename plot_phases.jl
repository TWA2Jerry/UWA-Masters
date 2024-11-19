using CairoMakie
using LaTeXStrings
include("prog.h")

f = Figure(size= (400, 300))
main_figures = f[1,1] = GridLayout()
main_axis = Axis(f[1,1])
upper_figs = main_figures[1,1] = GridLayout()
lower_figs = main_figures[2,1]

upper_ax1 = Axis(upper_figs[1,1]; limits = (0, 800, 0, 800), aspect = 1)
upper_ax2 = Axis(upper_figs[1,2]; limits = (0, 1000, 0, 1000), aspect =1)
upper_ax3 = Axis(upper_figs[1,3]; limits = (0, 10000, 0, 10000), aspect = 1)

lower_ax1 = Axis(lower_figs, 
	xscale = log10,
	#xticks = ([10, 10*sqrt(12), 100, 100*sqrt(12), 1000, 1000*sqrt(12), 10000], ["10", L"10\times\sqrt{12}", "100", L"100\times\sqrt{12}", "1000", L"1000\times\sqrt{12}", "10000"]),
	xticks = ([10*sqrt(12), 100, 100*sqrt(12), 1000, 1000*sqrt(12), 10000, 6000*sqrt(12)], [L"10\times\sqrt{12}", "100", L"100\times\sqrt{12}", "1000", L"1000\times\sqrt{12}", "10000", L"6000*\sqrt{12}"]),
	#limits = ((10, 10000*sqrt(12)), nothing),
	limits = ((10*sqrt(12), 20000*sqrt(12)), nothing),
	xgridvisible = false, 
	ygridvisible = false	
)
#rowsize!(main_figures, 2, Fixed(400))
hidedecorations!(upper_ax1)
hidespines!(upper_ax1)
hidedecorations!(upper_ax2)
hidespines!(upper_ax2)
hidedecorations!(upper_ax3)
hidespines!(upper_ax3)

hidedecorations!(main_axis)
hidespines!(main_axis)

#hidedecorations!(lower_ax1)
#hidespines!(lower_ax1)

crystal_model  = AgentsIO.load_checkpoint("./Records/stdod61/stdod61_crystal")
Makie.scatter!(upper_ax1, [crystal_model[i].pos for i in 1:nagents(crystal_model)], marker = arrow_marker, markersize = 3, rotations = [atan(crystal_model[i].vel[2], crystal_model[i].vel[1]) for i in 1:nagents(crystal_model)])

liquid_model = AgentsIO.load_checkpoint("./Records/stdod61/stdod61_liquid")
Makie.scatter!(upper_ax2, [liquid_model[i].pos for i in 1:nagents(liquid_model)], marker = arrow_marker, markersize = 3, rotations = [atan(liquid_model[i].vel[2], liquid_model[i].vel[1]) for i in 1:nagents(liquid_model)])

gas_model = AgentsIO.load_checkpoint("./Records/stdod61/stdod61_gaseous")
Makie.scatter!(upper_ax3, [gas_model[i].pos for i in 1:nagents(liquid_model)], marker = arrow_marker, markersize = 3, rotations = [atan(gas_model[i].vel[2], gas_model[i].vel[1]) for i in 1:nagents(gas_model)])

###Process rot_o data 
rot_o_alt_taves_file = open("Plotting/rot_order_alt_tave.txt", "r")
tdods_rot_o_alt_lines = readlines(rot_o_alt_taves_file)
tdods = Vector{Float64}(undef, 0)
rot_o_alts = Vector{Float64}(undef, 0)
for i in 2:length(tdods_rot_o_alt_lines)
    line = tdods_rot_o_alt_lines[i]
    split_line = parse.(Float64, split(line, " "))
    push!(tdods, split_line[1])
    push!(rot_o_alts, split_line[2])
end


Makie.lines!(lower_ax1, tdods, rot_o_alts)
hidespines!(lower_ax1, :t, :r)
#Makie.xlims!(10*sqrt(12), 10000*sqrt(12))

Makie.bracket!(lower_ax1,  log10(10*sqrt(12)), 0.1, log10(100*sqrt(12)), 0.1, text = "Crystal", style = :curly)
Makie.bracket!(lower_ax1,  log10(100*sqrt(12)), 0.1, log10(2000*sqrt(12)), 0.1, text = "Liquid", style = :curly)
Makie.bracket!(lower_ax1, log10(2000*sqrt(12)), 0.1, log10(10000*sqrt(12)), 0.1, text = "Gas", style = :curly)

#Makie.bracket!(upper_ax1,100, 50, 100*sqrt(12), 50, text = "Crystal", style = :curly)
