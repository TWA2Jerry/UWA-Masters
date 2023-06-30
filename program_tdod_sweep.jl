###Define IO. files
rot_alt_target_ave_file = open("rot_order_alt_tave.txt", "w") #This is the file that'll allow us to track DOD vs target DOD
compac_frac_file = open("compaction_frac.txt", "w")
mean_a_file = open("mean_area.txt", "w")
rot_o_file = open("rot_order.txt", "w")
rot_o_alt_file = open("rot_order_alt.txt", "w")
mean_speed_file = open("mean_speed.txt", "w")

###Load the functions
include("io_file.jl")
include("main.jl")

###Run ABM
for target_DOD in [0.0, sqrt(12), 2*sqrt(12), 10*sqrt(12), 20*sqrt(12),sqrt(12)*100, 200*sqrt(12), sqrt(12)*1000,  sqrt(12)*1000*2, 15700.0] 
	#truncate the file between each target DOD 
	global compac_frac_file = open("compaction_frac.txt", "w")
	global mean_a_file = open("mean_area.txt", "w")
	global rot_o_file = open("rot_order.txt", "w")
	global rot_o_alt_file = open("rot_order_alt.txt", "w")
	global mean_speed_file = open("mean_speed.txt", "w")

	truncate(compac_frac_file, 0)
        truncate(mean_a_file, 0)
        truncate(rot_o_file, 0)
        truncate(rot_o_alt_file, 0)
        truncate(mean_speed_file, 0)
	for i in 1:no_simulations
		run_ABM(i, target_DOD)
	end
	#=close(compac_frac_file)
	compac_frac_file = open("compaction_frac.txt", "r")
	bro = readlines(compac_frac_file)
	print(bro[1])
	exit()=#
	do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)
	
	###Do the stuff specific for target_DOD plotting
        rot_o_alt_ave_file = open("rot_o_alt_ave.txt", "r")
	rot_o_alt_simulation_mean_lines = readlines(rot_o_alt_ave_file)

	mean_rot_o_alt::Float64 = 0.0
	for line in rot_o_alt_simulation_mean_lines
		rot_o_alt_simulation_mean_splits = parse.(Float64, split(line, " "))
		mean_rot_o_alt += rot_o_alt_simulation_mean_splits[2]/(no_steps+1)
	end

	write(rot_alt_target_ave_file, "$target_DOD $mean_rot_o_alt\n")
	
end
