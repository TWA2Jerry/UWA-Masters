include("record_peripheral_agents.jl")
include("give_agent_cell.jl")
include("voronoi_area_file.jl")

function plot_dod_hist(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister})
	areas_on_periphery = []
	areas_on_interior = []	
	max_cell_area::Float64 = 0.0

	##Go through and calculate cells for the agents
	for i in 1:nagents(model)
		cell = give_agent_cell(model[i], model)
		status::Int32 = detect_periphery(cell, 1)
		cell_area::Float64 = voronoi_area(model, model[i].pos, cell, rho)
		max_cell_area = max(max_cell_area, cell_area)
		if(status == 1) 
			push!(areas_on_periphery, voronoi_area(model, model[i].pos, cell, rho))
		else 
			push!(areas_on_interior, voronoi_area(model, model[i].pos, cell, rho))
		end							
	end

	b_range = range(0.0, max_cell_area, 100)
	figure = histogram(areas_on_periphery, bins = b_range, label="Areas for agents on the periphery", color = :purple)
	histogram!(areas_on_interior, bins = b_range, label="Areas for agents on the interior", color = :green)	
	return figure
end

function histograms(pos_vels_file)

end
