const no_simulations::Int64 = 15
const no_steps::Int64 = 30000

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

