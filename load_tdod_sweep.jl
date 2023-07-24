###Define some constants
const no_steps::Int64 = 5000
const no_simulations::Int64 = 7

###Define the IO files for loading
rot_alt_target_ave_file = open("rot_order_alt_tave.txt", "a") #This is the file that'll allow us to track DOD vs target DOD
mean_speed_tave_file = open("mean_speed_tave.txt", "a")
compac_frac_file = open("compaction_frac.txt", "a")
mean_a_file = open("mean_area.txt", "a")
rot_o_file = open("rot_order.txt", "a")
rot_o_alt_file = open("rot_order_alt.txt", "a")
mean_speed_file = open("mean_speed.txt", "a")

###Load the function, variable and struct definitions
include("io_file.jl")
include("main.jl")

###Do the loading
loaded_model = AgentsIO.load_checkpoint("simulation_save.jld2")
print("Model loaded\n")

#=
###Simulate the model for however many more steps as necessary
step!(loaded_model, agent_step!, model_step!, no_steps-loaded_model.n)
write(compac_frac_file, "\n")
write(mean_a_file, "\n")
write(rot_o_file, "\n")
write(rot_o_alt_file, "\n")
write(mean_speed_file, "\n")


###Simulate for however many more simulations you need to do
for i in loaded_model.simulation_number+1:no_simulations
	run_ABM(i, loaded_model.target_area)
end

do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)

###Do the stuff specific for target_DOD plotting
rot_o_alt_ave_file = open("rot_o_alt_ave.txt", "r")
	rot_o_alt_simulation_mean_lines = readlines(rot_o_alt_ave_file)

	mean_rot_o_alt::Float64 = 0.0
        for i in Int64(no_steps/2):no_steps+1
                line = rot_o_alt_simulation_mean_lines[i]
                rot_o_alt_simulation_mean_splits = parse.(Float64, split(line, " "))
                global mean_rot_o_alt += rot_o_alt_simulation_mean_splits[2]/((no_steps+1)-Int64(no_steps/2)+1)
        end

        #write(rot_alt_target_ave_file, "$(loaded_model.target_area) $mean_rot_o_alt\n")

        ##This is for keeping track of mean speed vs tdod
        mean_speed_ave_file = open("mean_speed_ave.txt", "r")
        mean_speed_mean_lines = readlines(mean_speed_ave_file)

        mean_speed::Float64 = 0.0
        for i in Int64(no_steps/2):no_steps+1
                line = mean_speed_mean_lines[i]
                mean_speed_simulation_mean_splits = parse.(Float64, split(line, " "))
                global mean_speed += mean_speed_simulation_mean_splits[2]/((no_steps+1)-Int64(no_steps/2)+1)
        end

        write(mean_speed_tave_file, "$(loaded_model.target_area) $mean_speed\n")
=#

###And now, resume running the simulation as usual for the rest of the target areas
for target_DOD in [450.0, 200*sqrt(12), 1000.0, 500*sqrt(12)]
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
        for i in Int64(no_steps/2):no_steps+1
                line = rot_o_alt_simulation_mean_lines[i]
                rot_o_alt_simulation_mean_splits = parse.(Float64, split(line, " "))
                mean_rot_o_alt += rot_o_alt_simulation_mean_splits[2]/((no_steps+1)-Int64(no_steps/2)+1)
        end

        write(rot_alt_target_ave_file, "$target_DOD $mean_rot_o_alt\n")

        ##This is for keeping track of mean speed vs tdod
        mean_speed_ave_file = open("mean_speed_ave.txt", "r")
        mean_speed_mean_lines = readlines(mean_speed_ave_file)

        mean_speed::Float64 = 0.0
        for i in Int64(no_steps/2):no_steps+1
                line = mean_speed_mean_lines[i]
                mean_speed_simulation_mean_splits = parse.(Float64, split(line, " "))
                mean_speed += mean_speed_simulation_mean_splits[2]/((no_steps+1)-Int64(no_steps/2)+1)
        end

        write(mean_speed_tave_file, "$target_DOD $mean_speed\n")


end

