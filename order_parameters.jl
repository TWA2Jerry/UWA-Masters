include("agent_definition.jl")
include("rot_ord.jl")
include("some_math_functions.jl")
include("voronoi_area_file.jl")
#Define the data we want

function happiness(agent::bird)
        return abs((agent.A - agent.tdod)/max(pi/2*rho^2-agent.tdod, abs(0.0-agent.tdod)))
end

function center_of_mass(model)
        com::Tuple{Float64, Float64} = (0.0, 0.0)
        n::Int64 = nagents(model)
        for i in 1:n
                com =  (com .+ 1/n .* (model[i].pos))
        end
	return com
end

function radial_distance(agent, com)
        return norm(agent.pos .- com)
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

function polarisation(model)
	sum_of_vels::Tuple{Float64, Float64} = (0.0, 0.0)
	for i::Int64 in 1:no_birds
		sum_of_vels = sum_of_vels .+ 1/no_birds .* model[i].vel	
	end

	polarisation_order::Float64 = norm(sum_of_vels)
	return polarisation_order
end

function random_happiness(model)
	return happiness(model[model.tracked_agent]) 
end

function mean_no_moves(model)
	return 1-model.no_moves/nagents(model)
end

function random_radius(model)
        return radial_distance(model[model.tracked_agent], center_of_mass(model))
end

function mean_happiness(model)
	meanhappiness::Float64 = 0.0
	n::Int32 = nagents(model)
	for i in 1:nagents(model)
		meanhappiness += happiness(model[i])/n
	end
	return meanhappiness
end

function thingo(model)
 	
		kn::Vector{Float64} = [0.0, 0.0, 0.0, 0.0]
                q::Int64 = 8
                m::Int64 = 100
                next_pos_area = move_gradient_alt(model[i], model, kn, q, m, rho, model.target_area)
                area_next_step::Float64 = next_pos_area[2]
                next_pos::Tuple{Float64, Float64} = next_pos_area[1]
                distance_to_next_pos::Float64 = norm(next_pos .- model[i].pos)
end

function no_neighbours(cell)
	neighbour_count::Int32 = 0
	for i in 1:length(cell)
		if(cell[i][3] != 0) neighbour_count += 1 end
	end
	return neighbour_count
end

function mean_no_neighbours(model)
	n::Int32 = nagents(model)
	ave_no_neighbours::Float64 = 0.0
	all_agents_iterable = allagents(model)	
	for agent_i in all_agents_iterable

                neighbour_positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = []
                for agent_j in all_agents_iterable
                        if(agent_i.id == agent_j.id)
                                continue
                        end
                        push!(neighbour_positions, (agent_j.pos, agent_j.id))
                end
                ri::Tuple{Float64, Float64} = agent_i.pos
                vix::Float64 = agent_i.vel[1]
                viy::Float64 = agent_i.vel[2]
                relic_x::Float64 = -1.0*(-viy)
                relic_y::Float64 = -vix
                relic_pq::Tuple{Float64, Float64} = (relic_x, relic_y)
                relic_angle::Float64 = atan(relic_y, relic_x)
                relic_is_box::Int64 = 2
                relic_half_plane::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64} = (relic_angle, relic_pq, agent_i.pos, relic_is_box)
		temp_hp::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}}= []
                #print("The time for calculating a cell was\n")
                new_cell_i::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = voronoi_cell(model, ri, neighbour_positions, rho, eps, inf, temp_hp, agent_i.vel, [relic_half_plane])
		neighbour_count = no_neighbours(new_cell_i)
		ave_no_neighbours += neighbour_count/n
        end	
	return ave_no_neighbours
end

function regularity_metric(cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}}, area::Float64)
	perimeter_squared::Float64 = 0.0
	v::Int32 = length(cell)
	for i in 1:v
		perimeter_squared += distance(cell[i][1], cell[i%v+1][1])
	end
	return area/perimeter_squared
end
