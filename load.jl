###Define the IO files for loading
compac_frac_file = open("compaction_frac.txt", "a")
mean_a_file = open("mean_area.txt", "a")
rot_o_file = open("rot_order.txt", "a")
rot_o_alt_file = open("rot_order_alt.txt", "a")
mean_speed_file = open("mean_speed.txt", "a")

###Load the function, variable and struct definitions
include("io_file.jl")
include("main_functions.jl")

###

