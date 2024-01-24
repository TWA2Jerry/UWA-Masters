function cl(sl::Int32, r::Float64, L::Float64)
	return sl*r/L
end

function pl(rl, cl, alpha)
	return rl-alpha*cl
end

function pl_quick(agent_l::bird, alpha::Float64, model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister})
	l::Int64 = agent_l.id
	rl::Float64 = rl_quick(l, rho, model)
	return rl - alpha*cl(agent_l.collaborator, rho, rect_bound)
end

function wlm(pl::Float64, pm::Float64, beta::Float64 = 1.0)
	return 1/(1+exp((pl-pm)/beta))
end

function change_strat(agent_l::bird, alpha::Float64 = 1.0, beta::Float64 = 1.0, model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister})
	pl::Float64 = pl_quick(agent_l, alpha, model)
	neighbour_pos_vec::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
	for i in 1:no_birds
		if(i != agent_l.id) push!(neighbour_pos_vec, model[i].pos) end
	end	
	neighbour_vec::Vector{Int64} = neighbours_l_r(agent_l.pos, rho, neighbour_pos_vec)
	m::Int64 = rand(neighbour_vec)		
	pm::Float64 = pl_quick(model[m], alpha, model)
	wl::Float64 = wlm(pl, pm, 1.0)
	prob::Float64 = rand(Float64)
	
	return prob < wl ? (1, model[m].collaborator) : (0, model[m].collaborator)
end
