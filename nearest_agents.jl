include("some_math_functions.jl")
 
function find_nearest_agents(agent::bird, model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister})
	distance_id_v::Vector{Tuple{Float64, Int64}} = []
	for i in 1:nagents(model)
		if(model[i].id == agent.id) continue end
		distance::Float64 = norm(model[i].pos .- agent.pos)
		push!(distance_id_v, (distance, model[i].id))
	end	
	sort!(distance_id_v)
	ids::Vector{Int64} = []
	for i in 1:length(distance_id_v)
		push!(ids, distance_id_v[i][2])
	end
	return ids
end	
