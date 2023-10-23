include("global_vars.jl")
include("move_gradient_file.jl")
include("state_helper.jl")

agent_vals_file = open("agent_vals.txt", "w")

function write_agent_vals(model)
	for i in 1:nagents(model)
		write(agent_vals_file, "$(model.n) ")

		write(agent_vals_file, "$i ")
		
		write(agent_vals_file, "$(model[i].A) ")
		
		kn::Vector{Float64} = [0.0, 0.0, 0.0, 0.0]
        	q::Int64 = 8
        	m::Int64 = 100
        	next_pos_area::Tuple{Tuple{Float64, Float64}, Float64} = move_gradient_alt(model[i], model, kn, q, m, rho, model.target_area) 
		area_next_step::Float64 = next_pos_area[2]
		next_pos::Tuple{Float64, Float64} = next_pos_area[1]
		distance_to_next_pos::Float64 = norm(next_pos .- model[i].pos)
		
		write(agent_vals_file, "$distance_to_next_pos ")
		
		write(agent_vals_file, "$area_next_step\n")
	end
end
