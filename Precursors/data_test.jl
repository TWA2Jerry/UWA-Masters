###Define IO files
compac_frac_file = open("compaction_frac.txt", "w")
mean_a_file = open("mean_area.txt", "w")
rot_o_file = open("rot_order.txt", "w")
rot_o_alt_file = open("rot_order_alt.txt", "w")
mean_speed_file = open("mean_speed.txt", "w")

###Load the functions
include("main.jl")
include("io_file.jl")

using Statistics

###Run ABM
model = initialise(1000*sqrt(12), 1, no_birds)
adata = [(:A, Statistics.mean)]
adf, _ = run!(model, agent_step!, model_step!, no_steps)
write(compac_frac_file, "\n")
write(mean_a_file, "\n")
write(rot_o_file, "\n")
write(rot_o_alt_file, "\n")
write(mean_speed_file, "\n")
do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)

display(adf)
print(adf[1,1])




