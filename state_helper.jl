include("half_plane_fast.jl")
include("move_gradient_file.jl")
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

                new_cell_i::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = voronoi_cell_bounded(model, ri, neighbour_positions, rho, eps, inf, temp_hp, agent_i.vel, relic_half_plane)
                figure = draw_cell(new_cell_i)
		Plots.scatter!(agent_i.pos)
		display(figure)
                new_area::Float64 = voronoi_area(model, ri, new_cell_i, rho)
		 print("Now calculating the voronoi area for agent $(agent_i.id), which was $new_area\n")	
		return new_cell_i 
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

		cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} =  voronoi_cell(model, model[i].pos, positions, rho, eps, inf, temp_hp) 
		points::Vector{Tuple{Float64, Float64}} = []
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
	pot_pos::Tuple{Float64, Float64} = move_gradient_alt(model[id], model, kn, q, m, rho, model.target_area)   	

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
