include("half_plane_fast.jl")
include("half_plane_bounded.jl")
include("move_gradient_file.jl")
include("draw_circle_part.jl")
include("load_initialise.jl")
function draw_agent_cell(agent_i, model)
	all_agents_iterable =  allagents(model)
	temp_hp::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}} = []
        previous_areas::Vector{Float64} = zeros(nagents(model))
        actual_areas::Vector{Float64} = zeros(nagents(model))

                neighbour_positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = []
                for agent_j in all_agents_iterable
                        if(agent_i.id == agent_j.id)
                                continue
                        end
                        push!(neighbour_positions, (agent_j.pos, agent_j.id))
                end
                ri::Tuple{Float64, Float64} = agent_i.pos
                vix::Float64 = agent_i.vel[1]
                viy::Float64 = agent_i.vel[2]
                relic_x::Float64 = -1.0*(-viy)
                relic_y::Float64 = -vix
                relic_pq::Tuple{Float64, Float64} = (relic_x, relic_y)
                relic_angle::Float64 = atan(relic_y, relic_x)
                relic_is_box::Int64 = 2
                relic_half_plane::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64} = (relic_angle, relic_pq, agent_i.pos, relic_is_box)

                #new_cell_i::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = voronoi_cell(model, ri, neighbour_positions, rho, eps, inf, temp_hp, agent_i.vel, relic_half_plane)
                new_cell_i::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = give_agent_cell(agent_i, model)
		figure = draw_cell(new_cell_i)
		Plots.scatter!(agent_i.pos)
		display(figure)
                new_area::Float64 = voronoi_area(model, ri, new_cell_i, rho)
		 print("Now calculating the voronoi area for agent $(agent_i.id), which was $new_area\n")	
		return new_cell_i 
end

function give_agent_cell(agent_i, model)
        all_agents_iterable =  allagents(model)
        temp_hp::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}} = []
        previous_areas::Vector{Float64} = zeros(nagents(model))
        actual_areas::Vector{Float64} = zeros(nagents(model))

                neighbour_positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = []
                for agent_j in all_agents_iterable
                        if(agent_i.id == agent_j.id)
                                continue
                        end
                        push!(neighbour_positions, (agent_j.pos, agent_j.id))
                end
                ri::Tuple{Float64, Float64} = agent_i.pos
                vix::Float64 = agent_i.vel[1]
                viy::Float64 = agent_i.vel[2]
                relic_x::Float64 = -1.0*(-viy)
                relic_y::Float64 = -vix
                relic_pq::Tuple{Float64, Float64} = (relic_x, relic_y)
                relic_angle::Float64 = atan(relic_y, relic_x)
                relic_is_box::Int64 = 2
                relic_half_plane::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64} = (relic_angle, relic_pq, agent_i.pos, relic_is_box)

                new_cell_i::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = voronoi_cell(model, ri, neighbour_positions, rho, eps, inf, temp_hp, agent_i.vel, relic_half_plane)
        	
		##Procedure for adding the 
		cell_including_circle = []
		#print("Starting\n")
		for i in 1:length(new_cell_i)
			point = new_cell_i[i]
			point_pp = new_cell_i[(i)%length(new_cell_i)+1]
			push!(cell_including_circle, point)
			#print("$(point[3]) $(point_pp[2])\n")
			if(point[3] == 0 && point_pp[2] == 0)
				#print("State helper.jl here. Circle confirmed\n")
				vec_to_point = point[1] .- agent_i.pos
				vec_to_pointpp = point_pp[1] .- agent_i.pos
				theta_1 = atan(vec_to_point[2], vec_to_point[1])
				theta_2 = atan(vec_to_pointpp[2], vec_to_pointpp[1])
				if(theta_2 < theta_1)
					theta_2 += 2*pi
				end
				circle_points = circle_seg(agent_i.pos, rho, theta_1, theta_2)
				for j in 1:length(circle_points[1])
					push!(cell_including_circle, ((circle_points[1][j], circle_points[2][j]), 0, 0))
				end
			end
		end        
		return cell_including_circle
end

function draw_cell(cell)
	print("Draw cell called\n")	
	points::Vector{Tuple{Float64, Float64}} = []
	for i in 1:length(cell)
		push!(points, cell[i][1])
	end		
	push!(points, cell[1][1])
	figure = Plots.plot(points)
	#display(Plots.plot(points))
	return figure
end

function display_model_cell(model)
	figure = draw_model_cell(model)
	display(figure)
end

function draw_model_cell(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister})
	##Scatter the agent positions
	b_positions::Vector{Tuple{Float64, Float64}} = []
	colours::Vector{Float64} = []
        rotations::Vector{Float64} = []
	for i in 1:nagents(model)
		push!(colours, model[i].nospots)
                push!(rotations, atan(model[i].vel[2], model[i].vel[1]))
		push!(b_positions, model[i].pos)
	end
		
	
	figure, ax, colourbarthing = Makie.scatter(b_positions,axis = (; title = "Model state at step $(model.n)", limits = (0, rect_bound, 0, rect_bound)), marker = '→', markersize = 20, rotations = rotations, color = colours, colormap = cgrad(:matter, 300, categorical = true), colorrange = (0, 300))
	
	
	##Draw the cells	
	##For each agent, generate the cells and plot using the normal half plane bounded thingo. 
	for i in 1:nagents(model)
		##Just some colour stuff for the plot
		text!(model[i].pos, text = "$i", align = (:center, :top))
		temp_hp::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}} = []
		positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = Vector{Tuple{Tuple{Float64, Float64}, Int64}}(undef, 0)
		for j in 1:nagents(model)
			if(j == i) continue
			end
			push!(positions, (model[j].pos, model[j].id))	
		end		

		#cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} =  voronoi_cell(model, model[i].pos, positions, rho, eps, inf, temp_hp) 
		cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = give_agent_cell(model[i], model)
		points::Vector{Tuple{Float64, Float64}} = []
		if(length(cell) == 0)
			draw_circle_seg(b_positions[i], rho, 0.0, 2*pi)
			continue
		end		
		for i in 1:length(cell)
        	        push!(points, cell[i][1])
	        end
        	push!(points, cell[1][1])
		Makie.lines!(points, color = :black)
	end 
	Colorbar(figure[1,2], colourbarthing)
	#save("./Cell_Images/shannon_flock_n_=_$(model.n).png", figure)
	#display(figure)
	return figure		
end

function draw_tesselation(positions, model)
	##Scatter the agent positions
        b_positions::Vector{Tuple{Float64, Float64}} = []
        for i in 1:length(positions)
                push!(b_positions, positions[i])
        end


        figure, ax, colourbarthing = Makie.scatter(b_positions,axis = (; limits = (0, rect_bound, 0, rect_bound)), marker = :circle)


        ##Draw the cells
        ##For each agent, generate the cells and plot using the normal half plane bounded thingo.
        for i in 1:length(positions)
                ##Just some colour stuff for the plot
                text!(positions[i], text = "$i", align = (:center, :top))
                temp_hp::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}} = []
                c_positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = Vector{Tuple{Tuple{Float64, Float64}, Int64}}(undef, 0)
                for j in 1:length(positions)
                        if(j == i) continue
                        end
                        push!(c_positions, (positions[j], model[j].id))
                end

                cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} =  voronoi_cell(model, positions[i], c_positions, rho, eps, inf, temp_hp)
                points::Vector{Tuple{Float64, Float64}} = []
		for i in 1:length(cell)
                        push!(points, cell[i][1])
                end
                push!(points, cell[1][1])
                Makie.lines!(points, color = :black)
        end
        #Colorbar(figure[1,2], colourbarthing)
        #save("./Cell_Images/shannon_flock_n_=_$(model.n).png", figure)
        #display(figure)
        return figure

end

function show_move(model, id)
	##First, show the position that the agent with id of id will go to 
	kn::Vector{Float64} = [0.0, 0.0, 0.0, 0.0]
	q::Int64 = 8
	m::Int64 = 100
	pot_pos::Tuple{Float64, Float64} = move_gradient_alt(model[id], model, kn, q, m, rho, model.target_area)[1]   	

	##Next, evaluate and draw the voronoi tesselation of the model given that move of the agent	
	positions::Vector{Tuple{Float64, Float64}} = []
	for i in 1:nagents(model)
		if(i == id) 
			push!(positions, pot_pos)
			continue
		end
		push!(positions, model[i].pos)
	end

	figure = draw_tesselation(positions, model)
	print("Agent $id wanted to move to a new position of $pot_pos from its old position of $(model[id].pos)\n")
	Makie.scatter!(pot_pos, color = :yellow)
	Makie.scatter!(model[id].pos, color = :red)
	display(figure)
end

areas= []
potential_areas = []
function record_moves(model, pos_vels_file, start, agent_ids, no_steps)
###Model is our model, agent_ids should be a vector of the ids of agents you want to track, and no_steps is the number of steps you want to evolve teh model by. We should preset the model to the starting step we want in interactive.			
	empty!(areas)
	empty!(potential_areas)
	for id in agent_ids
		push!(areas, [])
		push!(potential_areas, [])
	end
	
	for i in 0:no_steps
		model = load_initialise(pos_vels_file, start+i)
		#Go through all the thangs, record their 
		for j in 1:length(agent_ids)
			id = agent_ids[j]
			#write(current_areas_file, "($model[agent_ids[id]].A)\n")
			print("The position of agent $(agent_ids[j]) is $(model[id].pos)\n")
			push!(areas[j], model[agent_ids[j]].A) 
		end
	
		#Go through the agents of interest, and use move_gradient_alt to find their best move in current position
		for j in 1:length( agent_ids)
			id = agent_ids[j]
			kn::Vector{Float64} = [0.0, 0.0, 0.0, 0.0]
		        q::Int64 = 8
       			m::Int64 = 100
			best_area =  move_gradient_alt(model[id], model, kn, q, m, rho, model.target_area)[2]
			best_difference = model[id].A - best_area
			#write(best_areas_file, "$best_difference $id\n")
			push!(potential_areas[j], best_area)
		end
		#step!(model, agent_step!, model_step!, 1)
	end
end

function plot_cave_ins(agent_ids)
	Plots.plot([t for t in 1:length(areas[1])], areas, label = transpose(agents_to_track))
end

function give_agent_cell_bounded(agent_i, model)
        all_agents_iterable =  allagents(model)
        temp_hp::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}} = []
        previous_areas::Vector{Float64} = zeros(nagents(model))
        actual_areas::Vector{Float64} = zeros(nagents(model))

                neighbour_positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = []
                for agent_j in all_agents_iterable
                        if(agent_i.id == agent_j.id)
                                continue
                        end
                        push!(neighbour_positions, (agent_j.pos, agent_j.id))
                end
                ri::Tuple{Float64, Float64} = agent_i.pos
                vix::Float64 = agent_i.vel[1]
                viy::Float64 = agent_i.vel[2]
                relic_x::Float64 = -1.0*(-viy)
                relic_y::Float64 = -vix
                relic_pq::Tuple{Float64, Float64} = (relic_x, relic_y)
                relic_angle::Float64 = atan(relic_y, relic_x)
                relic_is_box::Int64 = 2
                relic_half_plane::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64} = (relic_angle, relic_pq, agent_i.pos, relic_is_box)

                new_cell_i::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = voronoi_cell_bounded(model, ri, neighbour_positions, rho, eps, inf, temp_hp, agent_i.vel, relic_half_plane)
        	
		##Procedure for adding the 
		cell_including_circle = []
		print("Starting\n")
		for i in 1:length(new_cell_i)
			point = new_cell_i[i]
			point_pp = new_cell_i[(i)%length(new_cell_i)+1]
			push!(cell_including_circle, point)
			print("$(point[3]) $(point_pp[2])\n")
			if(point[3] == 0 && point_pp[2] == 0)
				print("State helper.jl here. Circle confirmed\n")
				vec_to_point = point[1] .- agent_i.pos
				vec_to_pointpp = point_pp[1] .- agent_i.pos
				theta_1 = atan(vec_to_point[2], vec_to_point[1])
				theta_2 = atan(vec_to_pointpp[2], vec_to_pointpp[1])
				if(theta_2 < theta_1)
					theta_2 += 2*pi
				end
				circle_points = circle_seg(agent_i.pos, rho, theta_1, theta_2)
				for j in 1:length(circle_points[1])
					push!(cell_including_circle, ((circle_points[1][j], circle_points[2][j]), 0, 0))
				end
			end
		end        
		return cell_including_circle
end

function draw_agent_cell_bounded(id, model)
	agent_i = model[id]
	cell = give_agent_cell_bounded(model[id], model)
	 points::Vector{Tuple{Float64, Float64}} = []
                for i in 1:length(cell)
                        push!(points, cell[i][1])
                end
        push!(points, cell[1][1])
	figure = Makie.scatter(agent_i.pos, marker = '→',markersize = 20, rotations = atan(agent_i.vel[2], agent_i.vel[1]))
	Makie.lines!(points, color = :black)
	display(figure)
end

function draw_model_cell_bounded(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister})
	##Scatter the agent positions
	b_positions::Vector{Tuple{Float64, Float64}} = []
	colours::Vector{Float64} = []
        rotations::Vector{Float64} = []
	for i in 1:nagents(model)
		push!(colours, model[i].nospots)
                push!(rotations, atan(model[i].vel[2], model[i].vel[1]))
		push!(b_positions, model[i].pos)
	end
		
	
	figure, ax, colourbarthing = Makie.scatter(b_positions,axis = (; title = "Model state at step $(model.n)", limits = (0, rect_bound, 0, rect_bound)), marker = '→', markersize = 20, rotations = rotations, color = colours, colormap = cgrad(:matter, 300, categorical = true), colorrange = (0, 300))
	
	
	##Draw the cells	
	##For each agent, generate the cells and plot using the normal half plane bounded thingo. 
	for i in 1:nagents(model)
		##Just some colour stuff for the plot
		text!(model[i].pos, text = "$i", align = (:center, :top))
		temp_hp::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}} = []
		positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = Vector{Tuple{Tuple{Float64, Float64}, Int64}}(undef, 0)
		for j in 1:nagents(model)
			if(j == i) continue
			end
			push!(positions, (model[j].pos, model[j].id))	
		end		

		#cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} =  voronoi_cell(model, model[i].pos, positions, rho, eps, inf, temp_hp) 
		cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = give_agent_cell_bounded(model[i], model)
		points::Vector{Tuple{Float64, Float64}} = []
		if(length(cell) == 0)
			draw_circle_seg(b_positions[i], rho, 0.0, 2*pi)
			continue
		end		
		for i in 1:length(cell)
        	        push!(points, cell[i][1])
	        end
        	push!(points, cell[1][1])
		Makie.lines!(points, color = :black)
	end 
	Colorbar(figure[1,2], colourbarthing)
	#save("./Cell_Images/shannon_flock_n_=_$(model.n).png", figure)
	#display(figure)
	return figure		
end

