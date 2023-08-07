include("agent_definition.jl")
#Define the data we want

function happiness(agent::bird)
        return abs((agent.A - 1000*sqrt(12))/(pi/2*rho^2-1000*sqrt(12)))
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
