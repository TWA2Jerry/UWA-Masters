include("give_agent_cell.jl")
include("some_math_functions.jl")

###Function takes in the agent objects. It calculates the center of the group (average of all agent positions). It then calculates the cross product of r_ig and v_i, the relative position of the agent to the group center and its velocity. 

function rot_o_generic(r_com::Tuple{Float64, Float64}, velocity::Tuple{Float64, Float64})
	#rcom is the agent position relative to com, velocity is the velocity of the agent (direction they face, not including speed)
	rot_o_raw::Float64 = norm(r_com) > eps ? cross((r_com) ./ norm(r_com), velocity) : 0.0
	return rot_o_raw
end 

function rot_o_alt(model)
        agents_iterable = allagents(model)
        return rot_ord_alt(agents_iterable)
end

function rot_o(model)
        agents_iterable = allagents(model)
        return rot_ord(agents_iterable)
end


function rot_ord(agents)
	num_agents = length(agents)

	##Calculate the group center
	r_g::Tuple{Float64, Float64} = (0.0, 0.0)
	for agent in agents
		r_g = r_g .+ agent.pos
	end

	r_g = r_g ./ length(agents)

	##Iterate through the agents again and find the absolute value of the sum of the cross products r_ig \times v_i
	num_neg = 0
	num_pos = 0
	rot_order::Float64 = 0.0
	for agent in agents
		rot_order += cross((agent.pos .- r_g) ./ norm(agent.pos .- r_g), agent.vel .* agent.speed)
		if(cross((agent.pos .- r_g) ./ norm(agent.pos .- r_g), agent.vel .* agent.speed) < 0.0)
			num_neg += 1
		else
			num_pos += 1
		end
	end
		
	rot_order = abs(rot_order)/length(agents)

	#print("The number of negative was $num_neg, the number of pos was $num_pos\n")
	return rot_order
end

function rot_ord_alt(agents)
        num_agents = length(agents)
	#print("Alt rotational order function here. Number of agents calculated to be $num_agents\n")
        ##Calculate the group center
        r_g = (0.0, 0.0)
        for agent in agents
                r_g = r_g .+ agent.pos
        end

        r_g = r_g ./ num_agents
	#print("Alt rotational order function here. r_g calculated to be $r_g\n")
        ##Iterate through the agents again and find the absolute value of the sum of the cross products r_ig \times v_i
        rot_order = 0.0
        for agent in agents
		rot_order += abs(cross((agent.pos .- r_g) ./ norm(agent.pos .- r_g), agent.vel .* agent.speed))
        	int_rot = abs(cross((agent.pos .- r_g) ./ norm(agent.pos .- r_g), agent.vel .* agent.speed))
		int_rad = dot((agent.pos .- r_g) ./ norm(agent.pos .- r_g),agent.vel .* agent.speed) 
		if(agent.id == 1)
			#print("Cross product calculated to be $int_rot, while the radial component was $int_rad. Total was $(int_rot^2 + int_rad^2)\n")
		end
	end

        rot_order = rot_order/num_agents
	#print("Rotational order calculated to be $rot_order\n")

        return rot_order
end

function dominant_rotation(model)
        num_agents = nagents(model)
	agents = allagents(model)

        ##Calculate the group center
        r_g::Tuple{Float64, Float64} = (0.0, 0.0)
        for agent in agents
                r_g = r_g .+ agent.pos
        end

        r_g = r_g ./ length(agents)

        ##Iterate through the agents again and find the absolute value of the sum of the cross products r_ig \times v_i
        num_neg = 0
        num_pos = 0
        rot_order::Float64 = 0.0
        for agent in agents
                rot_order += cross((agent.pos .- r_g) ./ norm(agent.pos .- r_g), agent.vel .* agent.speed)
                if(cross((agent.pos .- r_g) ./ norm(agent.pos .- r_g), agent.vel .* agent.speed) < 0.0)
                        num_neg += 1
                else
                        num_pos += 1
                end
        end

        rot_order = abs(rot_order)/length(agents)

	dominant_rotation_direction = num_pos > num_neg ? +1 : -1

        #print("The number of negative was $num_neg, the number of pos was $num_pos\n")
        return dominant_rotation_direction
end

###
function find_rep(a, rep)
	while(a != rep[a])
		a = rep[a]
	end
	return a
end

###Function that takes two integers, representing agents, and a vector/array that holds the reps of the agents, 
function merge(agent_a, agent_b, rep, size)
	a = agent_a
	b = agent_b
	rep_a = find_rep(a, rep)
	rep_b = find_rep(b, rep)
	
	if(rep_a == rep_b)
		return ##This is for if the two agents are already part of the same group
	end	

	if(size[rep_a] > size[rep_b])
		size[rep_a] += size[rep_b]
		rep[rep_b] = rep_a
	else
		size[rep_b] += size[rep_a]	
		rep[rep_a] = rep_b
	end
	
	return 
end

###
function group_ids(adj, rep, size)
	n = length(rep)
	for i in 1:n	
		rep[i] = i
		size[i] = 1
	end

	for i in 1:n
		for neighbour in adj[i]
			merge(i, neighbour, rep, size)
		end
	end

	return 
end

###Function that constructs the adjacency list rep of the voronoi graph
function construct_voronoi_adj_list(model, adj)
	n::Int64 = nagents(model)
        for i in 1:n
                adj[i] = Vector{Int64}(undef, 0)
                cell = give_agent_cell(model[i], model)
                agent_neighbours = neighbours(cell)
                for neighbour in agent_neighbours
                        push!(adj[i], neighbour)
                end
        end
	return
end

###Function that constructs the groups, and finds the rot_o of each group
function ave_group_rot_o(model)
	n::Int64 = nagents(model)
	adj::Array{Vector} = Array{Vector}(undef, n)
	construct_voronoi_adj_list(model, adj)
	
	rep::Vector{Int64} = Vector{Int64}(undef, n)
	size::Vector{Int64} = Vector{Int64}(undef, n)

	group_ids(adj, rep, size)

	groups = Set()
	for i in 1:n
		group_i = find_rep(i, rep)
		push!(groups, group_i)
	end

	no_groups = length(groups)

	group_dict = Dict(group => Vector{Int64}(undef,0) for group in groups)
	for i in 1:n
		rep_i = find_rep(i, rep)
		push!(group_dict[rep_i], i) 
	end

	ave_rot_o::Float64 = 0.0
	for group in groups
		positions::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
		com::Tuple{Float64, Float64} = (0.0, 0.0)
		for i in group_dict[group]
			push!(positions, model[i].pos)
		end
		com = center_of_mass(positions)
			
		group_rot_o::Float64 = 0.0
		for i in group_dict[group]
			group_rot_o += model[i].speed * rot_o_generic(model[i].pos .- com, model[i].vel)/size[group]
		end	
		print("rot_ord thang here. Rot o for group of size $(size[group]) is $group_rot_o\n")
		ave_rot_o += abs(group_rot_o)/no_groups
	end
	
	print("Number of groups detected was $no_groups and the average rot_o was $ave_rot_o\n")

	#=
	for group in groups
		print("$(size[group])\n")
	end
	=#
	return ave_rot_o
end

###Function that constructs the groups, and finds the rot_o of each group
size_rot_o_array = []
function record_group_rot_o(model)
        n::Int64 = nagents(model)
        adj::Array{Vector} = Array{Vector}(undef, n)
        construct_voronoi_adj_list(model, adj)

        rep::Vector{Int64} = Vector{Int64}(undef, n)
        size::Vector{Int64} = Vector{Int64}(undef, n)

        group_ids(adj, rep, size)

        groups = Set()
        for i in 1:n
                group_i = find_rep(i, rep)
                push!(groups, group_i)
        end

        no_groups = length(groups)

        group_dict = Dict(group => Vector{Int64}(undef,0) for group in groups)
        for i in 1:n
                rep_i = find_rep(i, rep)
                push!(group_dict[rep_i], i)
        end

        ave_rot_o::Float64 = 0.0
        for group in groups
                positions::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
                com::Tuple{Float64, Float64} = (0.0, 0.0)
                for i in group_dict[group]
                        push!(positions, model[i].pos)
                end
                com = center_of_mass(positions)

                group_rot_o::Float64 = 0.0
                for i in group_dict[group]
                        group_rot_o += rot_o_generic(model[i].pos .- com, model[i].vel)/size[group]
                end
                push!(size_rot_o_array, (size[group], abs(group_rot_o)))
		#print("Rot o for group of size $(size[group]) is $group_rot_o\n")
        end

        #print("Number of groups detected was $no_groups and the average rot_o was $ave_rot_o\n")
	
	#=
        for group in groups
                print("$(size[group])\n")
        end
	=#

        return 
end

