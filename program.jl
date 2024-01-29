const no_simulations::Int64 = 1
const no_steps::Int64 = 1000

###Define IO. files
compac_frac_file = open("compaction_frac.txt", "w")
mean_a_file = open("mean_area.txt", "w")
rot_o_file = open("rot_order.txt", "w")
rot_o_alt_file = open("rot_order_alt.txt", "w")
mean_speed_file = open("mean_speed.txt", "w")
rot_alt_target_ave_file = open("rot_order_alt_tave.txt", "w")
pos_vels_file = open("pos_vels.txt", "w")

###Load the functions
include("order_parameters.jl")
include("io_file.jl")
include("main.jl")

###Run ABM
#=
for i in 1:no_simulations
	run_ABM(i, 0.0)
end
=#

adata = [happiness, :true_A, :perimeter_squared, :no_neighbours, :rot_o_alt, :rot_o_alt_corr, agent_regularity]
mdata = [mean_radial_distance, rot_o_alt, random_happiness, mean_no_moves, polarisation, random_radius, mean_happiness, rot_o, mean_no_neighbours, no_collabs]

model = initialise(target_area_arg = 500*sqrt(12), simulation_number_arg = 1, no_bird = no_birds)
adf, mdf = @time run!(model, agent_step!, model_step!, no_steps; adata, mdata)

do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)
do_more_io_stuff(adf, mdf)
