#include("prog.h")
###Run ABM
#=
for i in 1:no_simulations
        run_ABM(i, 0.0)
end
=#

adata = [happiness, :true_A, :perimeter_squared, :no_neighbours, :rot_o_alt, :rot_o_alt_corr, agent_regularity]
mdata = [mean_radial_distance, rot_o_alt, random_happiness, mean_no_moves, polarisation, random_radius, mean_happiness, rot_o, mean_no_neighbours, no_collabs, num_in_bin]

target_dods = [100.0]

areas_of_interest = [sqrt(12), 50*sqrt(12), 100*sqrt(12), 200*sqrt(12), 1000*sqrt(12), 2000*sqrt(12), 18000.0, 23000.0, 10000*pi*rho^2]
lower_upper_areas = []
for upper_area in areas_of_interest
	for lower_area in areas_of_interest
		if(lower_area > upper_area) continue end
		push!(lower_upper_areas, (lower_area, upper_area))
	end
end

parameters = Dict(
        :simulation_number_arg => [i for i::Int64 in 1:no_simulations],
        :target_area_arg => target_dods,
	:area_args => lower_upper_areas,	
)

#= Original program runner.
model = initialise(target_area_arg = 1000.0*sqrt(12), simulation_number_arg = 1, no_bird = no_birds)
adf, mdf = @time run!(model, agent_step!, model_step!, no_steps; adata, mdata)
=#

###New thingo for running, just because there's never reason you wouldn't use this general method of running possibly multiple params
adf, mdf  = paramscan(parameters, initialise; adata, mdata, agent_step!, model_step!, n = no_steps, parallel = true)

#do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)
#do_more_io_stuff(adf, mdf)

