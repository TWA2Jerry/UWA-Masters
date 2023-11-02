const no_simulations::Int64 = 1
const no_steps::Int64 = 5000

###Define IO. files
compac_frac_file = open("compaction_frac.txt", "w")
mean_a_file = open("mean_area.txt", "w")
rot_o_file = open("rot_order.txt", "w")
rot_o_alt_file = open("rot_order_alt.txt", "w")
mean_speed_file = open("mean_speed.txt", "w")
rot_alt_target_ave_file = open("rot_order_alt_tave.txt", "w")
pos_vels_file = open("pos_vels.txt", "w")

###Load the functions
include("io_file.jl")
include("main.jl")

###Run ABM
for i in 1:no_simulations
	run_ABM(i, 1400.0)
end

do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)
#do_more_io_stuff(adf, mdf)
