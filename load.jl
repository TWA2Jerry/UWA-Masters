###Define the IO files for loading
compac_frac_file = open("compaction_frac.txt", "a")
mean_a_file = open("mean_area.txt", "a")
rot_o_file = open("rot_order.txt", "a")
rot_o_alt_file = open("rot_order_alt.txt", "a")
mean_speed_file = open("mean_speed.txt", "a")

###Load the function, variable and struct definitions
include("io_file.jl")
include("main_functions.jl")

###Do the loading
loaded_model = AgentsIO.load_checkpoint("simulation_save.jld2")
print("Model loaded\n")

###Simulate the model for however many more steps as necessary
step!(loaded_model, agent_step!, model_step!, no_steps-loaded_model.n)
write(compac_frac_file, "\n")
write(mean_a_file, "\n")
write(rot_o_file, "\n")
write(rot_o_alt_file, "\n")
write(mean_speed_file, "\n")

###
Simulate for however many more simulations you need to do
for i in model.simulation_number+1:no_simulations
	run_ABM(i)
end

do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)
