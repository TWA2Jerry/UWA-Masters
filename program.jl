const no_simulations::Int64 = 1
const no_steps::Int64 = 5000
include("prog.h")

###Run ABM
#=
for i in 1:no_simulations
	run_ABM(i, 0.0)
end
=#

target_dods = [2000*sqrt(12)]
alphaps = [1.0, 1000.0]

adata = [:delta_dod_var, :collab_p, :selfish_p]
mdata = [rot_o_alt, polarisation, rot_o, no_collabs]

parameters = Dict(
        :seed => [i for i::Int64 in 1:no_simulations],
        :target_area_arg => target_dods,
	:alpha_p_arg => alphaps 
)

model = initialise(target_area_arg = 3000.0*sqrt(12), simulation_number_arg = 1, no_bird = no_birds, alpha_p_arg = 1.0, r_arg = rho)
adf, mdf = @time run!(model, agent_step!, model_step!, no_steps; adata, mdata)

###New thingo for running, just because there's never reason you wouldn't use this general method of running possibly multiple params
#adf, mdf  = paramscan(parameters, initialise; adata, mdata, agent_step!, model_step!, n = no_steps)

#do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)
#do_more_io_stuff(adf, mdf)
