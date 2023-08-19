###Function that determines the gradient of movement
function move_gradient(agent, model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister},  kn::Vector{Float64}, q::Int64, m::Int64, rho::Float64, target_area::Float64 = 0.0)
	#Calculate the unit vector in the current direction of motion
	dt::Float64 = model.dt
	unit_v::Tuple{Float64,Float64} = agent.vel ./ 1.0
	theta_0::Float64 = atan(unit_v[2], unit_v[1])
	agent_speed::Float64 = 1.0
	vix::Float64 = unit_v[1]
	viy::Float64 = unit_v[2]
	positions::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
	all_agents_iterable = allagents(model)
	for neighbour in all_agents_iterable
		if(neighbour.id == agent.id)
			continue
		end
		pushfirst!(positions, neighbour.pos)	
	end	
	#min_area = inf  #The agent's current DOD area
	min_diff::Float64 = abs(agent.A - target_area)
	min_direction::Tuple{Float64, Float64} = (0.0, 0.0) #This is to set it so that the default direction of move is nowehere (stay in place)
	move_made::Int64 = 0
	pos_area_array::Vector{Tuple{Tuple{Float64,Float64}, Float64}}  = []
	no_angles_considered::Int64 = 0

	#Iterate through all the possible places the agent can move, keeping track of which one minimises area assuming static neighbour positions, though we make sure that if none of the moves optimises the current area, don't move at all
	#print("For agent $(agent.id), its min area is $min_area \n")
	temp_hp::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}} = []

	#For the relic idea, we have a bounding half plane based on the agent's current position and velocity
        relic_x::Float64 = -1.0*(-viy)
        relic_y::Float64 = -vix
        relic_pq::Tuple{Float64, Float64} = (relic_x, relic_y)
        relic_angle::Float64 = atan(relic_y, relic_x)
        relic_is_box::Int64 = -2
        relic_half_plane::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64} = (relic_angle, relic_pq, agent.pos, relic_is_box)
	for i::Int64 in 0:(q-1) #For every direction
		direction_of_move::Tuple{Float64, Float64} = (cos(i*2*pi/q)*vix - sin(i*2*pi/q)*viy, sin(i*2*pi/q)*vix + cos(i*2*pi/q)*viy)
		angle_of_move::Float64 = atan(direction_of_move[2], direction_of_move[1])
		rel_angle::Float64 = ((angle_of_move - theta_0 + pi)+2*pi)%(2*pi) - pi
		angular_conflict::Int64 = 0
		if(abs(rel_angle) > (1)*2*pi/q + eps)
			continue
		end
		no_angles_considered += 1
		for j::Int64 in 1:m #For every position up to m
			if(angular_conflict == 1) 
				break
			end
			conflict::Int64 = 0
			new_agent_pos::Tuple{Float64, Float64} = agent.pos .+ j .* direction_of_move .* agent_speed .* dt
		
			#Check first if there are no other agents in the potential position, note that we don't need to keep updating nearest neighbours since we assume the neighbours of a given agent are static
			for neighbour_position in positions
				if norm(new_agent_pos .- neighbour_position) < 2.0 #If moving in this direction and this m causes a collision, don't consider a move in this direction
					if(j == 1)
						angular_conflict = 1
					end
					conflict = 1
					break
				end			
			end			
			if (conflict == 1 || angular_conflict == 1)		
				continue
			end
			#If there are no other agents in the potential position (no conflicts), go ahead and evaluate the new DOD
                	
			###
			#print("\nThe time to calculate a voronoi cell in move gradient is ")
			agent_voronoi_cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} =  voronoi_cell_bounded(model, new_agent_pos, positions, rho, eps, inf, temp_hp, direction_of_move, relic_half_plane) #Generates the set of vertices which define the voronoi cell
                	new_area::Float64 = voronoi_area(model, new_agent_pos, agent_voronoi_cell, rho) #Finds the area of the agent's voronoi cell
			

			##Some error detection stuff
			if(new_area > pi*rho^2 && abs(new_area-pi*rho^2) > 10^(7))
				print("Conventional area exceeded by agent. For agent position of $new_agent_pos, the cell was $agent_voronoi_cell, with area of $new_area\n")
				print("\n\n\nThe dq for this position was \n")
                        	for i in 1:length(temp_hp)
                                	print("$(temp_hp[i])\n")
                        	end
				AgentsIO.save_checkpoint("simulation_save.jld2", model)
				exit()
			end

			if(rot_ord_check(new_agent_pos, agent_voronoi_cell) != 1)
				print("Rotational order violated for a potential position of $new_agent_pos\n")
				print("\n\n\nThe dq for this position was \n")
                        	for i in 1:length(temp_hp)
                                	print("$(temp_hp[i])\n")
                        	end
				print("The points and angles were\n")
				for i in 1:length(agent_voronoi_cell)
			                vec_to_point = [agent_voronoi_cell[i][1][1] - new_agent_pos[1], agent_voronoi_cell[i][1][2] - new_agent_pos[2]]
                			angle_of_vec = atan(vec_to_point[2], vec_to_point[1])
                        		print("$(agent_voronoi_cell[i]), $(angle_of_vec)")
				end
				AgentsIO.save_checkpoint("simulation_save.jld2", model)
				#exit()
			end

			#=		
			print("\n\n\nThe dq for this position was \n")
			for i in 1:length(temp_hp)
				print("$(temp_hp[i])\n")
			end
			=#

			#print("Potential new area of $new_area\n")
			#=
			#print("The vertices of this convex hull point are\n")
			if(convex_hull_point[agent.id] == 1)
				for i in 1:length(agent_voronoi_cell)
				 vector_to_vertex = agent_voronoi_cell[i][1] .- new_agent_pos
                                 angle_to_vertex = atan(vector_to_vertex[2], vector_to_vertex[1])
                                 print("$angle_to_vertex ")
                                 #print("$(atan(cell[i][1][2], cell[i][1][1])) ")
				end
			end
			print("\n")
			=#

			if (abs(new_area-target_area) < min_diff)
                        	min_diff = abs(new_area-target_area)
				#min_area = new_area
				#print("New min area of $min_area, direction of $direction_of_move\n")
                        	#min_direction = i*2*pi/q < pi ? (i > 1 ? (cos(1*2*pi/q)*vix - sin(1*2*pi/q)*viy, sin(1*2*pi/q)*vix + cos(1*2*pi/q)*viy) : direction_of_move) : (i<q-1 ? (cos(-1*2*pi/q)*vix - sin(-1*2*pi/q)*viy, sin(-1*2*pi/q)*vix + cos(-1*2*pi/q)*viy) : direction_of_move)
				min_direction = direction_of_move
                        	move_made = 1
				#=replace_vector(last_half_planes[Int64(agent.id)], [agent_voronoi_cell, temp_hp, new_agent_pos])
				if(convex_hull_point[agent.id] == 1)
					#print("Min area was lowered for agent $(agent.id), in a potential position of $(new_agent_pos),  here is the temp_hp\n")
				end=#
                	end

		end
		
		#Check area calculation through voronoi package
		#=
		pack_positions = Vector{Point2{Float64}}(undef, nagents(model)) 
		for i in 1:nagents(model)-1
			pack_positions[i+1] = Point2(positions[i])
		end
		pack_positions[1] = Point2(pot_new_pos)
		tess = voronoicells(pack_positions, rect)
		tess_areas = voronoiarea(tess)
		if(abs(new_area-tess_areas[1]) > 0.1)
			print("Area check, our calculated area was $new_area, theirs was $(tess_areas[1])\n")
		end
		=#
		
		#push!(pos_area_array, [angle_of_move, min_area])
	end

	#push!(moves_areas[agent.id], (model.n, agent.A, pos_area_array))
	if(move_made == 1)
		no_move[agent.id] = 0
	else 
		no_move[agent.id] = 1
	end
	
	#Create the noise addition
	epsilon::Vector{Float64} = randn(model.rng, Float64, 2)
	epsilon_prime::Vector{Float64} = randn(model.rng, Float64, 2)
	dW::Vector{Float64} = sqrt(model.dt) .* (epsilon .- epsilon_prime)

	if(move_made==1)
                agent.speed = 1.0
        else 
                #print("No movement made, agent area was $(agent.A)\n")
                turn = rand([-1, 1])
                min_direction = (cos(turn*2*pi/q)*vix - sin(turn*2*pi/q)*viy, sin(turn*2*pi/q)*vix + cos(turn*2*pi/q)*viy)
                agent.speed = 0.0
        end


	#Store the new position for updating in model step
	new_pos[agent.id] = Tuple(min_direction .* agent.speed .* model.dt .+ agent.pos .+ sigma*dW)
	if(new_pos[agent.id][1] > rect_bound || new_pos[agent.id][1] < 0.0 || new_pos[agent.id][2] > rect_bound || new_pos[agent.id][2] < 0.0)
		print("Move gradient file here. Agent $(agent.id) will step overbounds. This is for a rectangle bound of $rect_bound. The position was $(new_pos[agent.id]). This is for time step $(model.n), was the particle part of the convex hull? $(convex_hull_point[agent.id])\n")
		AgentsIO.save_checkpoint("simulation_save.jld2", model)	
		exit()
	end
	
	#This warning might not be as important due to the fact that we may set min_area = inf, and then skip all other possible positions for sampling
	#=if(min_area > pi*rho^2)
		print("Conventional area exceeded by agent $(agent.id)\n")
	end=#

	 #print("The number of angles considered was $no_angles_considered\n")
        #It really doesn't have to be like this, since  at least just for the simple SHH model of Dr.Algar, we can simply return a velocity
        kn[1] = (min_direction .* agent_speed)[1]
        kn[2] = (min_direction .* agent_speed)[2]
	return move_made
end







