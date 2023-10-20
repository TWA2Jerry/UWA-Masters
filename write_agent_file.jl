include("global_var
include("move_gradient_file.jl")

agent_vals_file = open("agent_vals.txt", "w")

function write_agent_vals(model)
	for i in 1:nagent(model)
		write(agent_vals_file, "$(model.n) ")

		write(agent_vals_file, "$i ")
		
		write(agent_vals_file, "$(model[i].A)\n")

		
	end
end
