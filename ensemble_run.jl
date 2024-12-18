##Define the no_simulations and steps
const no_steps = 5000
const no_simulations::Int64 = 20

##Include the IO files for the previous order parameters that we wanted. 
compac_frac_file = open("compaction_frac.txt", "w")
mean_a_file = open("mean_area.txt", "w")
rot_o_file = open("rot_order.txt", "w")
rot_o_alt_file = open("rot_order_alt.txt", "w")
mean_speed_file = open("mean_speed.txt", "w")
rot_alt_target_ave_file = open("rot_order_alt_tave.txt", "w")
no_moves_file = open("no_moves.txt", "w")
r_tave_file = open("radius_tave.txt", "w")
pos_vels_file = open("pos_vels.txt", "w")

##Include the main functions
include("io_file.jl")
include("order_parameters.jl")
include("main.jl")

using Statistics

function rot_o_alt(model)
	agents_iterable = allagents(model)
	return rot_ord_alt(agents_iterable)	
end

function rot_o(model)
        agents_iterable = allagents(model)
        return rot_ord(agents_iterable)
end

adata = [happiness, :tdodr, :nospots]
#adata = [(happiness, Statistics.mean)]
mdata = [mean_radial_distance, rot_o_alt, random_happiness, mean_no_moves, polarisation, random_radius, mean_happiness, rot_o]



#Define the parameters we want to scan over
target_dods = [0.0, 1*sqrt(12),2*sqrt(12), 5*sqrt(12), 10*sqrt(12), 20*sqrt(12), 1000*sqrt(12)]
target_dods = LinRange{Float64}(0.0, pi*rho^2, 500)
target_dods = [element for element in target_dods]
#left_biases::Vector{Float64} = LinRange(0.5, 1.0, 11)

parameters = Dict(
	:seed => [i for i::Int64 in 1:no_simulations],
	:target_area_arg => target_dods,
	#:left_bias_arg => left_biases
)

a = [i for i::Int64 in 1:no_simulations]
print(typeof(a))

##Run the ABM using paramscan, and with changing the seed
adf, mdf  = paramscan(parameters, initialise; adata, mdata, agent_step!, model_step!, n = no_steps)
print(adf[5, 2])
print(mdf[no_steps, 3])

do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)
do_more_io_stuff(adf, mdf)
write_rtave(mdf, target_dods)
