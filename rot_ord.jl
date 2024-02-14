###Function takes in the agent objects. It calculates the center of the group (average of all agent positions). It then calculates the cross product of r_ig and v_i, the relative position of the agent to the group center and its velocity. 

function rot_o_alt_generic(r_com::Tuple{Float64, Float64}, velocity::Tuple{Float64, Float64})
	#rcom is the agent position relative to com, velocity is the velocity of the agent (direction they face, not including speed)
	rot_o_raw::Float64 = norm(r_com) > eps ? cross((r_com) ./ norm(r_com), velocity) : 0.0
	return abs(rot_o_raw)
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
	rot_order::Float64 = 0.0
	for agent in agents
		rot_order += cross((agent.pos .- r_g) ./ norm(agent.pos .- r_g), agent.vel .* agent.speed)
	end
		
	rot_order = abs(rot_order)/length(agents)

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
			print("Cross product calculated to be $int_rot, while the radial component was $int_rad. Total was $(int_rot^2 + int_rad^2)\n")
		end
	end

        rot_order = rot_order/num_agents
	#print("Rotational order calculated to be $rot_order\n")

        return rot_order
end


