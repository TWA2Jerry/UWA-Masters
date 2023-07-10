##Define the no_simulations and steps
const no_steps = 5000
const no_simulations = 1

##Include the IO files for the previous order parameters that we wanted. 
compac_frac_file = open("compaction_frac.txt", "w")
mean_a_file = open("mean_area.txt", "w")
rot_o_file = open("rot_order.txt", "w")
rot_o_alt_file = open("rot_order_alt.txt", "w")
mean_speed_file = open("mean_speed.txt", "w")
rot_alt_target_ave_file = open("rot_order_alt_tave.txt", "w")

##Include the main functions
include("io_file.jl")
include("main.jl")

using Statistics

#Define the data we want
function happiness(agent)
 	return abs((agent.A - 1000*sqrt(12))/(pi/2*rho^2-1000*sqrt(12)))
end

function mean_radial_distance(model)
	com::Tuple{Float64, Float64} = (0.0, 0.0)
	n::Int64 = nagents(model)
	for i in 1:n
		com =  (com .+ 1/n .* (model[i].pos))	
	end

	mrd::Float64 = 0.0
	for i in 1:n
		mrd += 1/n*norm(model[i].pos .- com)	
	end
	
	return mrd
end

function rot_o_alt(model)
	agents_iterable = allagents(model)
	return rot_ord_alt(agents_iterable)	
end

adata = [(happiness, Statistics.mean)]
mdata = [mean_radial_distance, rot_o_alt]

#Define the parameters we want to scan over
parameters = Dict(
	:target_area_arg => [1000*sqrt(12)]
)

##Run the ABM using paramscan, and with changing the seed
adf, mdf  = paramscan(parameters, initialise; adata, mdata, agent_step!, model_step!, n = no_steps)

print(adf[5, 2])
print(mdf[no_steps, 3])

rot_o_alt_ave_file = open("ensemble_rot_o_alt.txt", "w")
for i in 1:no_steps+1
	write(rot_o_alt_ave_file, "$i $(mdf[i, 3])\n")	
end

close(rot_o_alt_ave_file)
#do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)
