include("io_file.jl")

function run_ABM()
        global compac_frac_file
        global mean_a_file
        global rot_o_file
        global rot_o_alt_file
        global mean_speed_file
for i in 1:no_simulations
        model = initialise(1000.0*sqrt(12), i)
        #figure, _ = abmplot(model)
        #save("./Simulation_Images/shannon_flock_n_=_$(0).png", figure)
        step!(model, agent_step!, model_step!, no_steps)
        write(compac_frac_file, "\n")
        write(mean_a_file, "\n")
        write(rot_o_file, "\n")
        write(rot_o_alt_file, "\n")
        write(mean_speed_file, "\n")
end

        do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)

end #This should be the end of the function or running the ABM
