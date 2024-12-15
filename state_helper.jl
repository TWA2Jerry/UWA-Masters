include("half_plane_fast.jl")
include("half_plane_bounded.jl")
include("move_gradient_file.jl")
include("draw_circle_part.jl")
include("load_initialise.jl")
include("global_vars.jl")
include("give_agent_cell.jl")

###Takes in any cell, of the type returned by voronoi_cell from half_plane bounded (which includes vertex info), and draw lines between vertices of the cell. 
function draw_cell(cell)
	print("Draw cell called\n")	
	points::Vector{Tuple{Float64, Float64}} = []
	for i in 1:length(cell)
		push!(points, cell[i][1])
	end		
	push!(points, cell[1][1])
	figure = Makie.lines(points, axis = (;aspect = 1))
	#display(Plots.plot(points))
	return figure
end

###Same as above, but just draws it on the figure that already exists
function draw_cell!(cell)
        print("state_helper here. Draw cell! called\n")
        points::Vector{Tuple{Float64, Float64}} = []
        for i in 1:length(cell)
                push!(points, cell[i][1])
        end
        push!(points, cell[1][1])
	Makie.lines!(points)
	return
end

function draw_cell_filled!(cell)
	print("state_helper here. Draw cell! called\n")
        points::Vector{Tuple{Float64, Float64}} = []
        for i in 1:length(cell)
                push!(points, cell[i][1])
        end
        Makie.poly!(points)
        return

end

function display_model_cell(model)
	figure = give_model_cell(model)
	display(figure)
end

function give_model(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}; fig_box = ((0.0, 0.0), (rect_bound, rect_bound)), marker_size = 10)
        ##Scatter the agent positions
        b_positions::Vector{Tuple{Float64, Float64}} = []
        colours::Vector{Float64} = []
        rotations::Vector{Float64} = []
        minx = model[1].pos[1]
        maxx = model[1].pos[1]
        miny = model[1].pos[2]
        maxy = model[2].pos[2]
        for i in 1:nagents(model)
                #push!(colours, model[i].nospots)
                push!(rotations, atan(model[i].vel[2], model[i].vel[1]))
                push!(b_positions, model[i].pos)
                minx = min(minx, model[i].pos[1])
                maxx = max(maxx, model[i].pos[1])
                miny = min(miny, model[i].pos[2])
                maxy = max(maxy, model[i].pos[2])
        end


        #figure, ax, colourbarthing = Makie.scatter(b_positions,axis = (; title = "Model state at step $(model.n)", limits = (minx-10, maxx+10, miny-10, maxy+10), aspect = 1), marker = :circle, markersize = 20, rotations = rotations, color = colours, colormap = cgrad(:matter, 300, categorical = true), colorrange = (0, 300))
       #figure, ax, colourbarthing = Makie.scatter(b_positions,axis = (;   limits = (fig_box[1][1], fig_box[2][1], fig_box[1][2], fig_box[2][2]), aspect = 1), marker = '→',  markersize = marker_size, rotations = rotations, color = :black)
        #figure, ax, colourbarthing = Makie.scatter(b_positions,axis = (;title = "Model state at step $(model.n)",  limits = (fig_box[1][1], fig_box[2][1], fig_box[1][2], fig_box[2][2]), aspect = 1, xticklabelsize = 30, yticklabelsize=30), marker = :circle,  rotations = rotations, color = :blue)
	figure, ax, colourbarthing = Makie.scatter(b_positions,axis = (;  limits = (fig_box[1][1], fig_box[2][1], fig_box[1][2], fig_box[2][2]), aspect = 1, xticklabelsize = 30, yticklabelsize=30), marker = :circle,  rotations = rotations, color = :black, markersize =marker_size)	
	#=
	for i in 1:nagents(model)
                text!(model[i].pos .-(0.0, 15.0), text = "$i", align = (:center, :top), fontsize = 40)
 	end
	=#
	#Colorbar(figure[1,2], colourbarthing)
        #save("./Cell_Images/shannon_flock_n_=_$(model.n).png", figure)
        #display(figure)
        return figure, ax
end


###Function that returns a figure of the cells for the entire model. 
function give_model_cell(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}; fig_box = ((0.0, 0.0), (rect_bound, rect_bound)))
	##Scatter the agent positions
	b_positions::Vector{Tuple{Float64, Float64}} = []
	colours::Vector{Float64} = []
        rotations::Vector{Float64} = []
	minx = model[1].pos[1]
	maxx = model[1].pos[1]
	miny = model[1].pos[2]
	maxy = model[2].pos[2]
	for i in 1:nagents(model)
		push!(colours, model[i].nospots)
                push!(rotations, atan(model[i].vel[2], model[i].vel[1]))
		push!(b_positions, model[i].pos)
		minx = min(minx, model[i].pos[1])	
		maxx = max(maxx, model[i].pos[1])
		miny = min(miny, model[i].pos[2])
		maxy = max(maxy, model[i].pos[2])
	end
		
	
	#figure, ax, colourbarthing = Makie.scatter(b_positions,axis = (; title = "Model state at step $(model.n)", limits = (minx-10, maxx+10, miny-10, maxy+10), aspect = 1), marker = :circle, markersize = 20, rotations = rotations, color = colours, colormap = cgrad(:matter, 300, categorical = true), colorrange = (0, 300))
	figure, ax, colourbarthing = Makie.scatter(b_positions,axis = (;  title = "Model state at step $(model.n)", limits = (fig_box[1][1], fig_box[2][1], fig_box[1][2], fig_box[2][2]), aspect = 1), marker = :circle,  rotations = rotations, color = :blue)
	#figure, ax, colourbarthing = Makie.scatter(b_positions,axis = (;title = "Model state at step $(model.n)",  limits = (minx-100, maxx+100, miny-100, maxy+100), aspect = 1), marker = :circle,  rotations = rotations, color = :blue)
	##Draw the cells	
	##For each agent, generate the cells and plot using the normal half plane bounded thingo. 
	for i::Int32 in 1:nagents(model)
		##Just some colour stuff for the plot
		#text!(model[i].pos, text = "$i", align = (:center, :top))
		temp_hp::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}} = []
		positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = Vector{Tuple{Tuple{Float64, Float64}, Int64}}(undef, 0)
		for j::Int32 in 1:nagents(model)
			if(j == i) continue
			end
			push!(positions, (model[j].pos, model[j].id))	
		end		

		#cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} =  voronoi_cell(model, model[i].pos, positions, rho, eps, inf, temp_hp) 
		cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = give_agent_cell(model[i], model)
		points::Vector{Tuple{Float64, Float64}} = []
		for j in 1:length(cell)
        	        push!(points, cell[j][1])
	        end
        	push!(points, cell[1][1])
		Makie.lines!(points, color = :black)
	end 
	#Colorbar(figure[1,2], colourbarthing)
	#save("./Cell_Images/shannon_flock_n_=_$(model.n).png", figure)
	#display(figure)
	return figure		
end

function give_model_cell_circled(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}; fig_box = ((0.0, 0.0), (rect_bound, rect_bound)))
        ##Scatter the agent positions
        b_positions::Vector{Tuple{Float64, Float64}} = []
        colours::Vector{Float64} = []
        rotations::Vector{Float64} = []
        for i in 1:nagents(model)
                #push!(colours, model[i].nospots)
                push!(rotations, atan(model[i].vel[2], model[i].vel[1]))
                push!(b_positions, model[i].pos)
        end


        figure, ax, colourbarthing = Makie.scatter(b_positions,axis = (; limits = (fig_box[1][1], fig_box[2][1], fig_box[1][2], fig_box[2][2]), aspect = 1), marker = '→', markersize = 20, rotations = rotations, color = :blue)


        ##Draw the cells
        ##For each agent, generate the cells and plot using the normal half plane bounded thingo.
        for i in 1:nagents(model)
                ##Just some colour stuff for the plot
                #text!(model[i].pos, text = "$i", align = (:center, :top))
                temp_hp::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}} = []
                positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = Vector{Tuple{Tuple{Float64, Float64}, Int64}}(undef, 0)
                for j in 1:nagents(model)
                        if(j == i) continue
                        end
                        push!(positions, (model[j].pos, model[j].id))
                end

                #cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} =  voronoi_cell(model, model[i].pos, positions, rho, eps, inf, temp_hp)
                #cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = give_agent_cell_circled(model[i], model)
                cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = give_cell_forward_quick(i, model)
		cell = give_cell_circled(cell, model[i].pos)
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


function draw_tesselation(positions, model)
	##Scatter the agent positions
        b_positions::Vector{Tuple{Float64, Float64}} = []
        for i in 1:length(positions)
                push!(b_positions, positions[i])
        end


        figure, ax, colourbarthing = Makie.scatter(b_positions,axis = (; limits = (0, rect_bound, 0, rect_bound), aspect = 1), marker = :circle)


        ##Draw the cells
        ##For each agent, generate the cells and plot using the normal half plane bounded thingo.
        for i in 1:length(positions)
                ##Just some colour stuff for the plot
                text!(positions[i] .+ (5.0, 5.0), text = "$i", align = (:center, :top))
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

function show_move(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}, id::Int64; view_box = ((0.0, 0.0), (rect_bound, rect_bound)), marker_size =10, draw_best_cell_arg = 1)
	##First, show the position that the agent with id of id will go to 
	kn::Vector{Float64} = [0.0, 0.0, 0.0, 0.0]
	q::Int64 = 8
	m::Int64 = 100
	move_tuple = move_gradient_alt(model[id], model, kn, q, m, rho, model.target_area) 
	pot_pos::Tuple{Float64, Float64} = move_tuple[1]
	sampled_positions = move_tuple[3]
	sampled_colours = move_tuple[4]	
	best_area = move_tuple[2]
	best_voronoi_cell = move_tuple[6]
	##Next, evaluate and draw the voronoi tesselation of the model given that move of the agent	
	positions::Vector{Tuple{Float64, Float64}} = []
	for i in 1:nagents(model)
		if(i == id) 
			push!(positions, pot_pos)
			continue
		end
		push!(positions, model[i].pos)
	end

	#figure = give_model_cell_circled(model, fig_box = view_box)
	#figure = give_model_cell(model, fig_box = view_box)
	figure, ax = give_model(model, fig_box = view_box, marker_size= marker_size)
	print("Agent $id wanted to move to a new position of $pot_pos with area of $best_area from its old position of $(model[id].pos) which had an area of $(model[id].A)\n")
        Makie.scatter!(sampled_positions, marker = :utriangle, color = sampled_colours, markersize = marker_size) #The agents sampled positions
	Makie.scatter!(model[id].pos, color = :blue, marker = Circle, markersize = marker_size, rotations = atan(model[id].vel[2], model[id].vel[1])) #The agent of interest's current position
	#Makie.scatter!(pot_pos, markersize = marker_size/2, color = :cyan)
	Makie.scatter!(pot_pos, markersize = marker_size, color = :cyan)
	#Makie.scatter!(pot_pos, markersize = marker_size, color = :blue)
	
	if(draw_best_cell_arg == 1)
		circled_cell = give_cell_circled(best_voronoi_cell, pot_pos)
		draw_agent_cell_bounded!(circled_cell)
	end
	return figure, ax
end

function show_move!(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}, id::Int64; view_box = ((0.0, 0.0), (rect_bound, rect_bound)), marker_size =10, m_arg = 100, m_spacing_arg = 1, qp_arg = 1, draw_best_cell_arg = 1, conflict_dist_arg = 2.0, colorrange_arg = (0, 1), show_calcs = 0)
        #First, show the position that the agent with id of id will go to
        kn::Vector{Float64} = [0.0, 0.0, 0.0, 0.0]
        q::Int64 = 8
        m::Int64 = m_arg
        move_tuple = move_gradient_alt(model[id], model, kn, q, m, rho, model.target_area, m_spacing = m_spacing_arg, qp = qp_arg, conflict_dist_arg = conflict_dist_arg)
        pot_pos::Tuple{Float64, Float64} = move_tuple[1]
        sampled_positions = move_tuple[3]
        sampled_colours = move_tuple[4]
        if(show_calcs == 1)
			for colour in sampled_colours
				print("$colour\n")
			end
		end
		best_area = move_tuple[2]
        best_voronoi_cell = move_tuple[6]
        ##Next, evaluate and draw the voronoi tesselation of the model given that move of the agent
        positions::Vector{Tuple{Float64, Float64}} = []
        for i in 1:nagents(model)
                if(i == id)
                        push!(positions, pot_pos)
                        continue
                end
                push!(positions, model[i].pos)
        end

        print("Agent $id wanted to move to a new position of $pot_pos with area of $best_area from its old position of $(model[id].pos) which had an area of $(model[id].A)\n")
        Makie.scatter!([model[i].pos for i in 1:nagents(model)], marker=:circle, color = :black, markersize = marker_size)
	Makie.scatter!(sampled_positions, marker = :utriangle, color = sampled_colours, markersize = marker_size, colormap = :cool, colorrange = colorrange_arg) #The agents sampled positions
        #Makie.scatter!(model[id].pos, color = :blue, marker = :circle, markersize = marker_size, rotations = atan(model[id].vel[2], model[id].vel[1])) #The agent of interest's current position
        Makie.scatter!(pot_pos, markersize = marker_size/2, color = :cyan)
        if(draw_best_cell_arg == 1)
			circled_cell = give_cell_circled(best_voronoi_cell, pot_pos)
			draw_agent_cell_bounded!(circled_cell)
		end
        #display(figure)
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


function draw_agent_cell_bounded!(bounded_cell)
	cell = bounded_cell
	points::Vector{Tuple{Float64, Float64}} = []
                for i in 1:length(cell)
                        push!(points, cell[i][1])
                end
        push!(points, cell[1][1])
	Makie.lines!(points, color = :black)
	return 
end


function draw_agent_cell_bounded(id, model)
	agent_i = model[id]
	cell = give_agent_cell_circled(model[id], model)
	 points::Vector{Tuple{Float64, Float64}} = []
                for i in 1:length(cell)
                        push!(points, cell[i][1])
                end
        push!(points, cell[1][1])
	figure = Makie.scatter(agent_i.pos, marker = '→', markersize = 20, rotations = atan(agent_i.vel[2], agent_i.vel[1]))
	Makie.lines!(points, color = :black)
	display(figure)
end

	

function record_dod_distributions(pos_vels_file, start, agent_ids, no_steps)
	for i in 0:no_steps
               	model = load_initialise(pos_vels_file, start+i)
		plot_dod_hist(model)		
	end
end



function give_move(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}, id::Int32)
        ##First, show the position that the agent with id of id will go to
        kn::Vector{Float64} = [0.0, 0.0, 0.0, 0.0]
        q::Int64 = 8
        m::Int64 = 100
        move_tuple = move_gradient_alt(model[id], model, kn, q, m, rho, model.target_area)
        pot_pos::Tuple{Float64, Float64} = move_tuple[1]
        sampled_positions = move_tuple[3]
        sampled_colours = move_tuple[4]
        best_area = move_tuple[2]
        ##Next, evaluate and draw the voronoi tesselation of the model given that move of the agent
        positions::Vector{Tuple{Float64, Float64}} = []
        for i in 1:nagents(model)
                if(i == id)
                        push!(positions, pot_pos)
                        continue
                end
                push!(positions, model[i].pos)
        end
	return move_tuple
end

function find_model_limits(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister})
	minx = model[1].pos[1]
        maxx = model[1].pos[1]
        miny = model[1].pos[2]
        maxy = model[1].pos[2]
        for i in 1:nagents(model)
                minx = min(minx, model[i].pos[1])
                maxx = max(maxx, model[i].pos[1])
                miny = min(miny, model[i].pos[2])
                maxy = max(maxy, model[i].pos[2])
        end
	
	return (minx, miny), (maxx, maxy)
end 

better_positions_vec::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)

function draw_cell_context(cell, positions::Vector{Tuple{Float64, Float64}}; fig_box = ((0.0, 0.0), (rect_bound, rect_bound)), filled= 0)
	figure, ax, colourbarthing = Makie.scatter(positions, axis = (; limits = (fig_box[1][1], fig_box[2][1], fig_box[1][2], fig_box[2][2]), aspect = 1), marker = :circle, markersize = 10, color = :black)
	for i in 1:length(positions)
		text!(positions[i], text = "$i", align = (:center, :top))
	end
	
	if(filled == 0)
		draw_cell!(cell)
	else
		draw_cell_filled!(cell)
	end
	return figure, ax
end

function draw_cell_context_quick(id::Int64, model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}; fig_box = ((0.0, 0.0), (rect_bound, rect_bound)), circled = 0, rhop = rho, filled = 0)
	ri::Tuple{Float64, Float64} = model[id].pos
	positions::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
	for i in 1:nagents(model)
		push!(positions, model[i].pos)
	end

	cell_id = give_agent_cell(model[id], model; rhop = rhop)
	if(circled == 1)
		cell_id = give_cell_circled(cell_id, ri; rhop= rhop)
	end
	return draw_cell_context(cell_id, positions; fig_box = fig_box, filled = filled)
end

function draw_cell_forward_context_quick(id::Int64, model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}; fig_box = ((0.0, 0.0), (rect_bound, rect_bound)), circled = 0, rhop = rho, filled = 0)
        ri::Tuple{Float64, Float64} = model[id].pos
        positions::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
        for i in 1:nagents(model)
                push!(positions, model[i].pos)
        end
        cell_id = give_cell_forward_quick(id, model; rhop = rhop)
	if(circled == 1)
		cell_id = give_cell_circled(cell_id, ri; rhop= rhop)
	end
        return draw_cell_context(cell_id, positions; fig_box = fig_box, filled= filled)
end

