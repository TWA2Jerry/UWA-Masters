###Function that determines the gradient of movement
include("record_peripheral_agents.jl")
include("nearest_agents.jl")
include("generate_relic.jl")
include("global_vars.jl")

using StatsBase
using VoronoiCells
function move_gradient(agent::bird, model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister},  kn::Vector{Float64}, q::Int64, m::Int64, rho::Float64, target_area::Float64 = 0.0)
	dt::Float64 = model.dt
	unit_v::Tuple{Float64,Float64} = agent.vel ./ 1.0
	theta_0::Float64 = atan(unit_v[2], unit_v[1])
	agent_speed::Float64 = 1.0
	v::Tuple{Float64,Float64} = agent.vel
	vix::Float64 = unit_v[1]
	viy::Float64 = unit_v[2]
	positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = Vector{Tuple{Tuple{Float64, Float64}, Int64}}(undef, 0)
	all_agents_iterable = allagents(model)
	for neighbour in all_agents_iterable
		if(neighbour.id == agent.id)
			continue
		end
		pushfirst!(positions, (neighbour.pos, neighbour.id))	
	end	
	
	#translate_periodic_quick(positions)	
	#print("Move grad file here. length of positions is $(length(positions))\n")	

	#min_area = inf  #The agent's current DOD area
	min_diff::Float64 = abs(agent.A - target_area)
	min_direction::Tuple{Float64, Float64} = (0.0, 0.0) #This is to set it so that the default direction of move is nowehere (stay in place)
	move_made::Int64 = 0
	pos_area_array::Vector{Tuple{Tuple{Float64,Float64}, Float64}}  = []
	no_angles_considered::Int64 = 0
	num_positions_better::Int32 = 0 #This is to implement and measure what Shannon wants, which is the number of positions the agent considers is better than its current position
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
	
	
	velocity_half_plane = generate_relic_alt(agent.pos, unit_v)
	

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
			for neighbour_position_tup in positions
				neighbour_position::Tuple{Float64, Float64} = neighbour_position_tup[1]
				if norm(new_agent_pos .- neighbour_position) < 2.0 #If moving in this direction and this m causes a collision, don't consider a move in this direction
					if(j == 1)
						angular_conflict = 1
					end
					conflict = 1
					break
				end			
			end			
			
			if (conflict == 1)		
				continue
			end

			#If there are no other agents in the potential position (no conflicts), go ahead and evaluate the new DOD
                	
			### Agent cell calculation
			#print("\nThe time to calculate a voronoi cell in move gradient is ")
			#agent_voronoi_cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} =  voronoi_cell(model, new_agent_pos, positions, rho, eps, inf, temp_hp, direction_of_move, relic_half_plane) #Generates the set of vertices which define the voronoi cell
                	bounded_cell_1 = voronoi_cell_bounded(model, new_agent_pos, positions, rho, eps, inf, temp_hp, direction_of_move, [relic_half_plane])
			bounded_area::Float64 = voronoi_area(model, new_agent_pos, bounded_cell_1, rho) #Finds the area of the agent's voronoi cell
			new_area::Float64 = bounded_area 


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

			if(rot_ord_check(new_agent_pos, bounded_cell_1) != 1)
				#=
				print("Rotational order violated for a potential position of $new_agent_pos\n")
				print("\n\n\nThe dq for this position was \n")
                        	for i in 1:length(temp_hp)
                                	print("$(temp_hp[i])\n")
                        	end
				print("The points and angles were\n")
				for i in 1:length(bounded_cell_1)
			                vec_to_point = [bounded_cell_1[i][1][1] - new_agent_pos[1], bounded_cell_1[i][1][2] - new_agent_pos[2]]
                			angle_of_vec = atan(vec_to_point[2], vec_to_point[1])
                        		print("$(bounded_cell_1[i]), $(angle_of_vec)")
				end
				=#
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
			
			if(abs(new_area-target_area) < abs(agent.A - target_area)) 
				num_positions_better += 1
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
                turn::Int32 = rand([1, -1])
		min_direction = (cos(turn*2*pi/q)*vix - sin(turn*2*pi/q)*viy, sin(turn*2*pi/q)*vix + cos(turn*2*pi/q)*viy)
		agent.speed = 0.0
        end
	#agent.nospots = num_positions_better

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










function move_gradient_alt(agent, model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister},  kn::Vector{Float64}, q::Int64, m::Int64, rho::Float64, target_area::Float64 = 0.0)
	#Calculate the unit vector in the current direction of motion
	dt::Float64 = model.dt
	unit_v::Tuple{Float64,Float64} = agent.vel ./ 1.0
	theta_0::Float64 = atan(unit_v[2], unit_v[1])
	agent_speed::Float64 = 1.0
	vix::Float64 = unit_v[1]
	viy::Float64 = unit_v[2]
	positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = Vector{Tuple{Tuple{Float64, Float64}, Int64}}(undef, 0)
	all_agents_iterable = allagents(model)
	vel_angle::Float64 = atan(viy, vix) #This is for rotating vectors so we know which neighbours are out of view
	for neighbour in all_agents_iterable
		if(neighbour.id == agent.id)
			continue
		end
		
		#=
		rij::Tuple{Float64, Float64} = neighbour.pos .- agent.pos
		rij_angle::Float64 = atan(rij[2], rij[1])
		if(abs(rij_angle - vel_angle) > pi/2) continue end
		=#
		
		pushfirst!(positions, (neighbour.pos, neighbour.id))	
	end	

	#translate_periodic_quick(positions)
        #print("Move grad file here. length of positions is $(length(positions))\n")

	min_area = agent.A  #The agent's current DOD area
	min_diff::Float64 = abs(agent.A - target_area)
	min_direction::Tuple{Float64, Float64} = (0.0, 0.0) #This is to set it so that the default direction of move is nowehere (stay in place)
	min_distance::Float64 = inf
	move_made::Int64 = 0
	pos_area_array::Vector{Tuple{Tuple{Float64,Float64}, Float64}}  = []
	no_angles_considered::Int64 = 0
	num_positions_better::Int32 = 0 #This is to implement and measure what Shannon wants, which is the number of positions the agent considers is better than its current position
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
	best_pos::Tuple{Float64, Float64} = agent.pos
	sampled_positions::Vector{Tuple{Float64, Float64}} = []
	colours = []	
	best_voronoi_cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = []	

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
			#=
			if(angular_conflict == 1) 
				break
			end
			=#
			conflict::Int64 = 0
			new_agent_pos::Tuple{Float64, Float64} = agent.pos .+ j .* direction_of_move .* agent_speed .* dt
			push!(sampled_positions, new_agent_pos)	
			#Check first if there are no other agents in the potential position, note that we don't need to keep updating nearest neighbours since we assume the neighbours of a given agent are static
			for neighbour_position_tup in positions
				neighbour_position::Tuple{Float64, Float64} = neighbour_position_tup[1]
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
			agent_voronoi_cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} =  voronoi_cell_bounded(model, new_agent_pos, positions, rho, eps, inf, temp_hp, direction_of_move, [relic_half_plane]) #Generates the set of vertices which define the voronoi cell
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

			#if (abs(new_area-target_area) < min_diff && conflict != 1)
			lower_area::Float64 = model.lower_area
			upper_area::Float64 = model.upper_area
			if((new_area > lower_area && new_area < upper_area && j < min_distance || (move_made == 0 && abs(new_area-target_area) < min_diff)) && conflict != 1)	
				min_diff = abs(new_area-target_area)
				min_area = new_area
				#print("New min area of $min_area, direction of $direction_of_move\n")
                        	#min_direction = i*2*pi/q < pi ? (i > 1 ? (cos(1*2*pi/q)*vix - sin(1*2*pi/q)*viy, sin(1*2*pi/q)*vix + cos(1*2*pi/q)*viy) : direction_of_move) : (i<q-1 ? (cos(-1*2*pi/q)*vix - sin(-1*2*pi/q)*viy, sin(-1*2*pi/q)*vix + cos(-1*2*pi/q)*viy) : direction_of_move)
				min_direction = direction_of_move
			
				if((new_area > lower_area && new_area < upper_area && j < min_distance))
                                        min_distance = j
                                end	
				
				move_made = 1
				#=replace_vector(last_half_planes[Int64(agent.id)], [agent_voronoi_cell, temp_hp, new_agent_pos])
				if(convex_hull_point[agent.id] == 1)
					#print("Min area was lowered for agent $(agent.id), in a potential position of $(new_agent_pos),  here is the temp_hp\n")
				end=#
				best_pos = new_agent_pos
				best_voronoi_cell = agent_voronoi_cell
                	end
			
			colour = :black
			if(abs(new_area-target_area) < abs(agent.A - target_area)) 
				num_positions_better += 1
				colour = (conflict == 1) ? :orange : :green
				push!(better_positions_vec, new_agent_pos)
			else
				colour = :red
			end
			push!(colours, colour)
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
                turn = rand([1, -1])
                min_direction = (cos(turn*2*pi/q)*vix - sin(turn*2*pi/q)*viy, sin(turn*2*pi/q)*vix + cos(turn*2*pi/q)*viy)
		agent.speed = 0.0
        end
	#agent.nospots = num_positions_better
		
	#Store the new position for updating in model step
	new_pos[agent.id] = Tuple(min_direction .* agent.speed .* model.dt .+ agent.pos .+ sigma*dW)
	correct_x::Float64 = (((new_pos[agent.id])[1])%rect_bound+rect_bound)%rect_bound
        correct_y::Float64 = (((new_pos[agent.id])[2])%rect_bound+rect_bound)%rect_bound
        new_pos[agent.id] = (correct_x, correct_y)
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
	#return Tuple(min_direction .* agent.speed .* model.dt .+ agent.pos .+ sigma*dW)
	#print("Best pos was $best_pos, with a difference of $min_diff, with an area of $min_area\n")
	
	return best_pos, min_area, sampled_positions, colours, move_made, best_voronoi_cell
end




function move_gradient_collab(agent::bird, model, kn::Vector{Float64}, r::Float64 = rho, eta::Float64 = 1.0)
	##Create vector of neighbour positions
	neighbour_positions::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
	for i in 1:no_birds
		push!(neighbour_positions, model[i].pos) 
	end

	##Find neighbours
	neighbour_set::Vector{Int64} = neighbours_l_r(agent.id, r, neighbour_positions)

	##Align velocity with neighbour velocities
	theta_vec::Vector{Float64} = Vector{Tuple{Float64, Float64}}(undef, 0)
	for nid in neighbour_set
		push!(theta_vec, atan(model[nid].vel[2], model[nid].vel[1]))
	end

	theta_tpp::Float64 = mean(theta_vec) + 2* eta .* rand(Float64) .- (eta)

	##Set kn[1,4]
	kn[1] = agent.vel[1]
	kn[2] = agent.vel[2]
	kn[3] = (cos(theta_tpp) .- agent.vel[1])
	kn[4] = (sin(theta_tpp) .- agent.vel[2])

	agent.speed = 1.0
		
	new_pos[agent.id] = agent.pos .+ (agent.vel .+ (kn[3], kn[4])) .* agent.speed		
	correct_x::Float64 = (((new_pos[agent.id])[1])%rect_bound+rect_bound)%rect_bound
        correct_y::Float64 = (((new_pos[agent.id])[2])%rect_bound+rect_bound)%rect_bound
	new_pos[agent.id]  = (correct_x, correct_y)	
	if(new_pos[agent.id][1] > rect_bound || new_pos[agent.id][1] < 0.0 || new_pos[agent.id][2] > rect_bound || new_pos[agent.id][2] < 0.0)
                print("Move gradient file here. Agent $(agent.id) will step overbounds. This is for a rectangle bound of $rect_bound. The position was $(new_pos[agent.id]). This is for time step $(model.n), was the particle part of the convex hull? $(convex_hull_point[agent.id])\n")
                AgentsIO.save_checkpoint("simulation_save.jld2", model)
                exit()
        end

	return 1
end
