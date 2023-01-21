###Function takes in the agent objects. It calculates the center of the group (average of all agent positions). It then calculates the cross product of r_ig and v_i, the relative position of the agent to the group center and its velocity. 
function rot_ord(agents)
	num_agents = length(agents)

	##Calculate the group center
	r_g = (0.0, 0.0)
	for agent in agents
		r_g = r_g .+ agent.pos
	end

	r_g = r_g ./ length(agents)

	##Iterate through the agents again and find the absolute value of the sum of the cross products r_ig \times v_i
	rot_order = 0.0
	for agent in agents
		rot_order += abs(cross((agent.pos .- r_g) ./ norm(agent.pos .- r_g), agent.vel .* agent.speed))
	end
		
	rot_order = rot_order/length(agents)

	return rot_order
end
