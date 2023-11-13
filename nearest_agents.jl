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

function find_turn(agent::bird, model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister})
	left_turns::Int32 = 0
	agent_angle::Float64 = atan(agent.vel[2], agent.vel[1])
	for i in 1:nagents(model)
		if(model[i].id == agent.id) continue end
		neighbour = model[i]
		neighbour_vector::Tuple{Float64, Float64} = neighbour.pos .- agent.pos
		rnv::Tuple{Float64, Float64} = rotate_vector(-agent_angle, neighbour_vector)
		neighbour_angle = atan(rnv[2], rnv[1])
		if(abs(neighbour_angle) > pi/2) 
			continue
		elseif neighbour_angle < 0.0
			left_turns -= 1
		else
			left_turns += 1
		end	
	end	

	return left_turns >= 0 ? 1 : -1
end

function find_cop(agent::bird, model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister})
	agent_angle::Float64 = atan(agent.vel[2], agent.vel[1])
        positions::Vector{Tuple{Float64, Float64}} = []
	for i in 1:nagents(model)
                if(model[i].id == agent.id) continue end
                neighbour = model[i]
                neighbour_vector::Tuple{Float64, Float64} = neighbour.pos .- agent.pos
                rnv::Tuple{Float64, Float64} = rotate_vector(-agent_angle, neighbour_vector)
                neighbour_angle = atan(rnv[2], rnv[1])
                if(abs(neighbour_angle) > pi/2) 
                        continue
        	else
			push!(positions, neighbour.pos)
		end
	end     
	return center_of_pos(positions)
end	
