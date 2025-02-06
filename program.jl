const no_simulations::Int64 = 1
const no_steps::Int64 = 5000


include("prog.h")

adata = [(normalised_true_area, StatsBase.mean)]
mdata = [mean_speed, rot_o_alt, rot_o, polarisation]

target_dods = [1000*sqrt(12)]
q_args = [8]
qp_args = [1]
m_args = [100]
fov_args = [180.0] #The FOV parameter (total angular vision in degrees). Please note that since Voronoi cells are only convex, we can only handle an FOV such that 0 < FOV < 180, or FOV = 360 degrees.

parameters = Dict(
        :simulation_number_arg => [i for i::Int64 in 1:no_simulations],
        :target_area_arg => target_dods,
		:q_arg => q_args,
		:qp_arg => qp_args,
		:m_arg => m_args,
		:fov_arg => fov_args,
        #:left_bias_arg => left_biases
) 

#Simulation marker details. The best settings are marker = :circle, markersize = 10, or marker = arrow_marker, markersize = 4
simulation_marker = :circle 
simulation_markersize = 10



#model = initialise(target_area_arg = 1000.0*sqrt(12), simulation_number_arg = 1, no_bird = no_birds)
#adf, mdf = @time run!(model, agent_step!, model_step!, no_steps; adata, mdata)

###New thingo for running, just because there's never reason you wouldn't use this general method of running possibly multiple params
adf, mdf  = paramscan(parameters, initialise; adata, mdata, agent_step!, model_step!, n = no_steps)

#do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)
#do_more_io_stuff(adf, mdf)


