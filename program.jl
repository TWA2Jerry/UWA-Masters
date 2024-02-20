const no_simulations::Int64 = 1
const no_steps::Int64 = 5000

###Define IO. files
compac_frac_file = open("compaction_frac.txt", "w")
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

###Run ABM
#=
for i in 1:no_simulations
	run_ABM(i, 0.0)
end
=#

target_dodsrange = LinRange{Float64}(0.0, 3*1000*sqrt(12), 100)
target_dods = [element for element in target_dodsrange]
target_dodsrange = LinRange{Float64}(4*1000*sqrt(12), pi*rho^2, 50)
for tdod in target_dodsrange
	push!(target_dods, tdod)
end

adata = [happiness, :true_A]
mdata = [mean_radial_distance, rot_o_alt, random_happiness, mean_no_moves, polarisation, random_radius, mean_happiness, rot_o]

target_dods = [1000.0*sqrt(12)]

parameters = Dict(
        :simulation_number_arg => [i for i::Int64 in 1:no_simulations],
        :target_area_arg => target_dods,
        #:left_bias_arg => left_biases
) 


###New thingo for running, just because there's never reason you wouldn't use this general method of running possibly multiple params
adf, mdf  = paramscan(parameters, initialise; adata, mdata, agent_step!, model_step!, n = no_steps, parallel = true)

#do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)
#do_more_io_stuff(adf, mdf)
