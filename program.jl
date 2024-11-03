include("prog.h")
###Run ABM
#=
for i in 1:no_simulations
        run_ABM(i, 0.0)
end
=#

no_hp_behind_ignored::Vector{Int32} = Vector{Int32}(undef, 0)

adata = [happiness, :true_A, :perimeter_squared, :no_neighbours, :rot_o_alt, :rot_o_alt_corr, :distance, agent_regularity, :best_A, :A, :direction]
mdata = [mean_speed, rot_o_alt, polarisation, rot_o, ave_group_rot_o, max_group_rot_o, max_rot_o_group_size, no_groups]

target_dods = [100.0]
#=
a = range(0.65, 0.7, 2)
left_biases = Vector{Float64}(undef, 0)
for i in a
	push!(left_biases, i)
end
=#

left_biases = [0.5, 0.75]
lower_upper_areas = [(1*sqrt(12), 2000*sqrt(12))]

no_birds_vec = [5]

parameters = Dict(
        :simulation_number_arg => [i for i::Int64 in 1:no_simulations],
        :area_args => [(1*sqrt(12), 2000*sqrt(12))],
	:no_bird => no_birds_vec
	#:target_area_arg => target_dods,
	#:left_bias_arg => left_biases	
)

#model = initialise(target_area_arg = 100.0, simulation_number_arg = 1, no_bird = no_birds, area_args = (1*sqrt(12), 2500*sqrt(12)), left_bias_arg = 0.5)
#adf, mdf = @time run!(model, agent_step!, model_step!, no_steps; adata, mdata)


###New thingo for running, just because there's never reason you wouldn't use this general method of running possibly multiple params
_, mdf  = paramscan(parameters, initialise; mdata, agent_step!, model_step!, n = no_steps)

#do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)
#do_more_io_stuff(adf, mdf)

