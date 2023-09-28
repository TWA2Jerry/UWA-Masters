include("half_plane_show.jl")

function show_cell_calculation(model, position, neighbour_positions, rho)
	bro = voronoi_cell_show(model, position, neighbour_positions, rho)		
end


function show_cell_agent(model, agent_num, rho)
	neighbour_positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = []
	for i in 1:nagents(model)
		if(i == agent_num)
			continue
		end
		push!(neighbour_positions, (model[i].pos, i))
	end
	show_cell_calculation(model, model[agent_num].pos, neighbour_positions, rho)
end
