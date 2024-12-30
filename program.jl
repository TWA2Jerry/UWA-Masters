const no_simulations::Int64 = 20
const no_steps::Int64 = 5000


include("prog.h")

adata = [(normalised_true_area, StatsBase.mean)]
mdata = [mean_speed, rot_o_alt, rot_o, polarisation]

target_dods = [2000.0*sqrt(12)]
q_args = [8]
qp_args = [1]
m_args = [100]


parameters = Dict(
        :simulation_number_arg => [i for i::Int64 in 1:no_simulations],
        :target_area_arg => target_dods,
		:q_arg => q_args,
		:qp_arg => qp_args,
		:m_arg => m_args,
        #:left_bias_arg => left_biases
) 

#model = initialise(target_area_arg = 1000.0*sqrt(12), simulation_number_arg = 1, no_bird = no_birds)
#adf, mdf = @time run!(model, agent_step!, model_step!, no_steps; adata, mdata)

###New thingo for running, just because there's never reason you wouldn't use this general method of running possibly multiple params
adf, mdf  = paramscan(parameters, initialise; adata, mdata, agent_step!, model_step!, n = no_steps)

#do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)
#do_more_io_stuff(adf, mdf)
