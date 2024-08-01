const no_simulations::Int64 = 1
const no_steps::Int64 = 10000

###Define IO. files
#=compac_frac_file = open("compaction_frac.txt", "w")
mean_a_file = open("mean_area.txt", "w")
rot_o_file = open("rot_order.txt", "w")
rot_o_alt_file = open("rot_order_alt.txt", "w")
mean_speed_file = open("mean_speed.txt", "w")
rot_alt_target_ave_file = open("rot_order_alt_tave.txt", "w")
pos_vels_file = open("pos_vels.txt", "w")

###Load the functions
include("order_parameters.jl")
include("io_file.jl")
include("main.jl")
=#
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

model = initialise(target_area_arg = 4000.0*sqrt(12), simulation_number_arg = 1, no_bird = no_birds, alpha_p_arg = 1.0)
adf, mdf = @time run!(model, agent_step!, model_step!, no_steps; adata, mdata)

###New thingo for running, just because there's never reason you wouldn't use this general method of running possibly multiple params
#adf, mdf  = paramscan(parameters, initialise; adata, mdata, agent_step!, model_step!, n = no_steps)

#do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)
#do_more_io_stuff(adf, mdf)
