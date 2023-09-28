include("half_plane_fast.jl")
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

function draw_model_cell(model)
	##Scatter the agent positions
	b_positions::Vector{Tuple{Float64, Float64}} = []
	for i in 1:nagents(model)
		push!(b_positions, model[i].pos)
	end
		
	figure = Plots.scatter(b_positions)

	##Draw the cells	
	##For each agent, generate the cells and plot using the normal half plane bounded thingo. 
	for i in 1:nagents(model)
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
		Plots.plot!(points)
	end 

	display(figure)		
end
