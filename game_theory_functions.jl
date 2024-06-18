include("order_parameters.jl")
include("global_vars.jl")
function cl(sl::Int32, r::Float64, L::Float64)
	return sl*r/L
end

function pl(rl, cl, alpha)
	return rl-alpha*cl
end

function pl_quick(agent_l::bird, model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}; alpha::Float64 = 1.0, r::Float64 = rho)
	l::Int64 = agent_l.id
	rl::Float64 = rl_quick(l, r, model)
	return rl - alpha*cl(agent_l.collaborator, r, rect_bound)
end

function pl_selfish_quick(agent_l::bird, model; alpha::Float64 = 1.0, alpha_p::Float64 = 0.5, r::Float64 = rho)
	l::Int64 = agent_l.id
        #rl::Float64 = rl_quick(l, r, model)
	#return rl + agent_l.delta_dod_var
	delta_max = max(agent_l.tdod, abs(pi*rho^2-agent_l.tdod))
	return alpha_p*agent_l.delta_dod_var/delta_max
end

function wlm(pl::Float64, pm::Float64, beta::Float64 = 1.0)
	return 1/(1+exp((pl-pm)/beta))
end

function change_strat(agent_l::bird, model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}; alpha::Float64 = 1.0, alpha_p::Float64 = 0.5, beta::Float64 = 1.0, r::Float64 = 0.5*rho)
	##Calculate the profit the agent would have under both the selfish and collab strats
	pl::Float64 = pl_quick(agent_l, model, alpha=alpha, r= 0.5*rho)
	pl_selfish::Float64 = pl_selfish_quick(agent_l, model, alpha_p=alpha_p, r= r)
	neighbour_pos_vec::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
	for i in 1:no_birds
		push!(neighbour_pos_vec, (model[i].pos, i)) 
	end	
	neighbour_vec::Vector{Int64} = neighbours_l_r(agent_l.id, rho, neighbour_pos_vec)
	alphap::Float64 = alpha	
	wl::Float64 = 0.0
	if(length(neighbour_vec) > 0) 
		m::Int64 = rand(neighbour_vec) 	
		pm::Float64 = 0.0
		if(agent_l.collaborator == 0) 
			wl = wlm(pl_selfish, pl, 1.0)
		else
			wl = wlm(pl, pl_selfish, 1.0)
		end
	end
	prob::Float64 = rand(Float64)

	if(model.n == no_steps-1)
		print("Agent $(agent_l.id) here. Selfish profit was $pl_selfish, the change in dod is $(agent_l.delta_dod_var). Alpha_p is $alpha_p\n")
	end
	
	return prob < wl ? (1, agent_l.collaborator == 1 ? 0 : 1) : (0, 1) #This line says, if the probability sampled is less than wl, we need to change (1) and change to the opp. strat of the agents' current strat, otherwise do nothing. 
end

function translate_positions(positions::Vector{Tuple{Float64, Float64}}, xoffset::Float64 = 0.0, yoffset::Float64 = 0.0) 
	translated_positions::Vector{Tuple{Float64, Float64}} = [positions[i] .+ (xoffset, yoffset) for i in 1:length(positions)]
	return translated_positions
end

function add_translated_positions(original_positions::Vector{Tuple{Float64, Float64}}, all_positions::Vector{Tuple{Float64, Float64}}, xoffset::Float64 = 0.0, yoffset::Float64 = 0.0)
	temp_positions::Vector{Tuple{Float64, Float64}} = translate_positions(original_positions, xoffset, yoffset)
	for trans_position in temp_positions
		push!(all_positions, trans_position)
	end

	return
end

#Function that takes in a vector of neighbour positions and adds to that vec their translations for periodic boundaries 
function translate_periodic_quick(positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}}, period::Float64 = rect_bound)
	only_positions::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
	all_translated_positions::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
	index_to_identity::Vector{Int64} = Vector{Int64}(undef, 0)	

	for i in 1:length(positions)
		push!(only_positions, positions[i][1])
		push!(index_to_identity, positions[i][2])
		push!(all_translated_positions, positions[i][1])
	end
		
	add_translated_positions(only_positions, all_translated_positions, period, 0.0)
	add_translated_positions(only_positions, all_translated_positions, period, period)
	add_translated_positions(only_positions, all_translated_positions, 0.0, period)
	add_translated_positions(only_positions, all_translated_positions, -period, period)
	add_translated_positions(only_positions, all_translated_positions, -period, 0.0)
	add_translated_positions(only_positions, all_translated_positions, -period, -period)
	add_translated_positions(only_positions, all_translated_positions, 0.0, -period)
	add_translated_positions(only_positions, all_translated_positions, period, -period)	

	#print("The length of positions is $(length(positions)), while the length of only_positions is $(length(only_positions))\n")

	no_neighbours::Int64 = Int64(length(positions))	
	for i in (no_neighbours+1):length(all_translated_positions)
		index = (i-1)%(no_neighbours)+1
		#print("Index is $index, i is $i\n")
		push!(positions, (all_translated_positions[i], index_to_identity[index]))	
	end
	#print("Length of positions after is $(length(positions))\n")
	return
end

function translate_periodic_quick(positions::Vector{Tuple{Float64, Float64}}, period::Float64 = rect_bound)
        only_positions::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
        all_translated_positions::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
        index_to_identity::Vector{Int64} = Vector{Int64}(undef, 0)

        for i in 1:length(positions)
                push!(only_positions, positions[i])
        end

        add_translated_positions(only_positions, all_translated_positions, period, 0.0)
        add_translated_positions(only_positions, all_translated_positions, period, period)
        add_translated_positions(only_positions, all_translated_positions, 0.0, period)
        add_translated_positions(only_positions, all_translated_positions, -period, period)
        add_translated_positions(only_positions, all_translated_positions, -period, 0.0)
        add_translated_positions(only_positions, all_translated_positions, -period, -period)
        add_translated_positions(only_positions, all_translated_positions, 0.0, -period)
        add_translated_positions(only_positions, all_translated_positions, period, -period)

        #print("The length of positions is $(length(positions)), while the length of only_positions is $(length(only_positions))\n")

        no_neighbours::Int64 = Int64(length(positions))
        for i in (no_neighbours+1):length(all_translated_positions)
                index = (i-1)%(no_neighbours)+1
                #print("Index is $index, i is $i\n")
                push!(positions, all_translated_positions[i])
        end
        #print("Length of positions after is $(length(positions))\n")
        return
end
	
