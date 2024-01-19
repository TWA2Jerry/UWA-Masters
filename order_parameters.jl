include("agent_definition.jl")
include("rot_ord.jl")
include("some_math_functions.jl")
include("polygon_area_file.jl")
#Define the data we want

function happiness(agent::bird)
        return abs((agent.A - agent.tdod)/max(pi/2*rho^2-agent.tdod, abs(0.0-agent.tdod)))
end

function center_of_mass(positions::Vector{Tuple{Float64, Float64}})
	com::Tuple{Float64, Float64} = (0.0, 0.0)
        n::Int64 = length(positions)
        for i in 1:n
                com =  (com .+ 1/n .* (positions[i]))
        end
        return com
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

function no_neighbours(cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}})
	neighbour_count::Int32 = 0
	for i in 1:length(cell)
		if(cell[i][3] != 0) neighbour_count += 1 end
	end
	return neighbour_count
end

function neighbours(cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}})
	neighbour_set::Vector{Int64} = Vector{Int32}(undef, 0) #We could make this a set, but the vec should work since every vertex can be associated with one unique neighbour
	for i in 1:length(cell)
                if(cell[i][3] != 0) push!(neighbour_set, cell[i][3]) end
        end
	return neighbour_set
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


function sides_squared(vertices::Vector{Tuple{Float64, Float64}})
         perimeter_squared::Float64 = 0.0
         v::Int32 = length(vertices)
         for i in 1:v
                 perimeter_squared += distance(vertices[i], vertices[i%v+1])^2
         end
        return perimeter_squared
end


function cell_sides_squared(cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}})
	perimeter_squared::Float64 = 0.0
        v::Int32 = length(cell)
        for i in 1:v
                perimeter_squared += distance(cell[i][1], cell[i%v+1][1])^2
        	if(cell[i][2] == 0 || cell[i%v+1][2] == 0)
			perimeter_squared = 0.0 #This is used because I don't think it makes sense to calculate regularities for cells with circles
			break
		end
	end
	return perimeter_squared
end

function regularity_metric(vertices::Vector{Tuple{Float64, Float64}})
	A::Float64 = polygon_area(vertices)
	perimeter_squared::Float64 = sides_squared(vertices)
	return A/perimeter_squared
end

function regularity_metric(cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}}, area::Float64)
	perimeter_squared::Float64 = 0.0
	v::Int32 = length(cell)
	for i in 1:v
		perimeter_squared += distance(cell[i][1], cell[i%v+1][1])^2
	end
	return area/perimeter_squared
end

function agent_regularity(agent::bird)
	return abs(agent.true_A/agent.sides_squared - regularities[agent.no_neighbours])/(regularities[agent.no_neighbours])
end

function return_regularities()
        regularities = Dict{Int32, Float64}([])
	r::Float64 = 1.0
        for n in 3:20
                vertexes::Vector{Tuple{Float64, Float64}} = Array{Tuple{Float64, Float64}}(undef, 0)
                for i in 0:n
                        coord::Tuple{Float64, Float64} = (r*cos(2*pi*i/n), r*sin(2*pi*i/n))
                        push!(vertexes, coord)
                end

        	regularities[n] = regularity_metric(vertexes)
	end
        return regularities
end



regularities::Dict{Int32, Float64} = return_regularities()
print("Regularities calculated\n")

function regularity_dev(agent::bird)
	no_sides::Int32 = agent.no_neighbours
	return (agent_regularity(agent)-regularities[no_sides])/(regularities[no_sides])
end

function agent_neighbour_correlation(agent::bird, neighbours::Vector{Int64}, model)
	neighbour_rot_o_alts::Vector{Float64} = Vector{Float64}(undef, 0)
	for nid in neighbours
		push!(neighbour_rot_o_alts, model[nid].rot_o_alt)
	end
	d = abs(agent.rot_o_alt - mean(neighbour_rot_o_alts))
	correlation::Float64 = 1-d
	return correlation
end

function model_correlation(model)
	
	#Go through all agents
		#Calculate their voronoi cell and therefore who their neighbours are

		#Calculate how similar their rot_o_alt parameter is to that of their neighbours
end
