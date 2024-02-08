using Plots
using InteractiveDynamics
using CairoMakie # choosing a plotting backend
using ColorSchemes
import ColorSchemes.balance

###Function for drawing the plots for model step
function draw_figures(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}, actual_areas::Vector{Float64}, previous_areas::Vector{Float64}, delta_max::Float64, new_pos::Vector{Tuple{Float64, Float64}}, path_points::Vector{Tuple{Float64, Float64}} = [])
        ##Draw the standard figure of the agents with their DODs after the step
        colours::Vector{Float64} = []
        rotations::Vector{Float64} = []
        allagents_iterable = allagents(model)
        target_area::Float64 = model.target_area
        com::Tuple{Float64, Float64} = center_of_mass(model)
        for id in 1:nagents(model)
                #push!(colours, abs(model[id].A-model.target_area)/(delta_max))
                #push!(colours, radial_distance(model[id], com)/200.0)
                #push!(colours, agent_regularity(model[id]))
                push!(rotations, atan(model[id].vel[2], model[id].vel[1]))
        	push!(colours, distance(model[id].pos, best_pos[id])) #This is for helping cave ins. 
	end
        #figure, _ = abmplot(model)
        print("\n\n\nThe number of points in new_pos is $(length(new_pos)), the first element is $(new_pos[1])\n")
        #print("About to do the figure\n")


        #figure, ax, colourbarthing = Makie.scatter([Tuple(point) for point in new_pos], axis = (; title = "Model state at step $(model.n)", limits = (0, rect_bound, 0, rect_bound), aspect = 1), marker = '→', markersize = 20, rotations = rotations, color = colours, colormap = cgrad(:matter, 300, categorical = true))
        #figure, ax, colourbarthing = Makie.scatter([Tuple(point) for point in new_pos], axis = (;  limits = (0, rect_bound, 0, rect_bound), aspect = 1), marker = :circle,  rotations = rotations, color = :black)
	figure, ax, colourbarthing = Makie.scatter([model[i].pos for i in 1:nagents(model)], axis = (;  limits = (0, rect_bound, 0, rect_bound), aspect = 1), marker = '→',  markersize = 20, rotations = rotations, color = colours, colormap = :viridis, colorrange = (0.0, 100.0)) #This is for detecting cave ins better

        for i in 1:nagents(model)
                text!(new_pos[i], text = "$i", align = (:center, :top))
        end

        #print("The number of points in path points is $(length(path_points))\n")
        #draw_path(path_points)
        #title!("Model state at step $(model.n)")
        #text!(model[model.tracked_agent].pos, text = "$(model.tracked_agent)", align = (:center, :top))
        ###tracking the radial distance of each agent from group center
        radial_distances::Vector{Float64}  = []
        #=for i in 1:nagents(model)
                push!(radial_distances, radial_distance(model[i], com))
        end
        for i in 1:nagents(model)
                text!(new_pos[i], text = "$(trunc(radial_distances[i]))", align = (:center, :top))
        end=#
        Colorbar(figure[1,2], colourbarthing)
        save("./Simulation_Images/shannon_flock_n_=_$(model.n).png", figure)
end


function draw_actual_DODs(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}, actual_areas::Vector{Float64}, previous_areas::Vector{Float64}, delta_max::Float64, new_pos::Vector{Tuple{Float64, Float64}}, path_points::Vector{Tuple{Float64, Float64}} = [])
         ##Draw the figure of the agents with their actual DODs
        for id in 1:nagents(model)
                colours[id] = actual_areas[id]/(pi*rho^2)
        end
        figure_actual, ax, colourbarthing = Makie.scatter([Tuple(point) for point in new_pos], axis = (; limits = (0, rect_bound, 0, rect_bound)), marker = '→', markersize = 20, rotations = rotations, color = colours, colormap = :viridis, colorrange = (0.0, 0.250)) #Note that I have no idea what the colorbarthing is for
        #=for i in 1:nagents(model)
                text!(new_pos[i], text = "$i", align = (:center, :top))
        end=#
        Colorbar(figure_actual[1,2], colourbarthing)
        save("./Simulation_Images_Actual_Areas/shannon_flock_n_=_$(model.n).png", figure_actual)

end


function draw_delta_DOD(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}, actual_areas::Vector{Float64}, previous_areas::Vector{Float64}, delta_max::Float64, new_pos::Vector{Tuple{Float64, Float64}}, path_points::Vector{Tuple{Float64, Float64}} = [])
        ##Draw the figure of the agents with their change in DOD
        for id in 1:nagents(model)
                #print("Current A is $(model[id].A), previous areas was $(previous_areas[id])\n")
                colours[id] = (abs(model[id].A - model.target_area)-abs(previous_areas[id]-model.target_area))/(2*delta_max)
        end
        figure_difference, ax, colourbarthing = Makie.scatter([Tuple(point) for point in new_pos], axis = (; limits = (0, rect_bound, 0, rect_bound)), marker = '→', markersize = 20, rotations = rotations, color = colours, colormap = :viridis, colorrange = (-0.1, 0.1)) #Note that I have no idea what the colorbarthing is for
        #=for i in 1:nagents(model)
                text!(new_pos[i], text = "$i", align = (:center, :top))
        end=#
        #Makie.scatter!([Tuple(point) for point in path_points], marker = :circle, color = :black, markersize = 20)
        #draw_path(path_points)
        Colorbar(figure_difference[1,2], colourbarthing)
        save("./Simulation_Images_Difference_Areas/shannon_flock_n_=_$(model.n).png", figure_difference)

end


###Function for drawing future figures and whatnot
function draw_figures_futures(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}, actual_areas::Vector{Float64}, previous_areas::Vector{Float64}, delta_max::Float64, new_pos::Vector{Tuple{Float64, Float64}}, path_points::Vector{Tuple{Float64, Float64}} = [])
        ##Draw the standard figure of the agents with their DODs after the step
        colours::Vector{Float64} = []
        rotations::Vector{Float64} = []
        allagents_iterable = allagents(model)
        target_area::Float64 = model.target_area
        com::Tuple{Float64, Float64} = center_of_mass(model)
        for id in 1:nagents(model)
                #push!(colours, abs(model[id].A-model.target_area)/(delta_max))
                #push!(colours, radial_distance(model[id], com)/200.0)
                push!(colours, distance(model[id].pos, best_pos[id]))
                push!(rotations, atan(model[id].vel[2], model[id].vel[1]))
        end
        #figure, _ = abmplot(model)
        print("\n\n\nThe number of points in new_pos is $(length(new_pos)), the first element is $(new_pos[1])\n")
        #print("About to do the figure\n")


        #figure, ax, colourbarthing = Makie.scatter([Tuple(point) for point in new_pos], axis = (; title = "Model state at step $(model.n)", limits = (0, rect_bound, 0, rect_bound), aspect = 1), marker = '→', markersize = 20, rotations = rotations, color = colours, colormap = cgrad(:matter, 300, categorical = true))
	figure, ax, colourbarthing = Makie.scatter([model[i].pos for i in 1:nagents(model)], axis = (;  limits = (0, rect_bound, 0, rect_bound), aspect = 1), marker = :circle,  rotations = rotations, color = colours, colormap = :viridis, colorrange = (0.0, 100.0))
	Makie.scatter!([Tuple(point) for point in best_pos], marker = :circle,  rotations = rotations, color = :blue)
	
	for i in 1:length(new_pos)
		Makie.lines!([new_pos[i], best_pos[i]], color= :black)
	end

	#=
        for i in 1:nagents(model)
                text!(new_pos[i], text = "$i", align = (:center, :top))
        end =#

        #print("The number of points in path points is $(length(path_points))\n")
        #draw_path(path_points)
        #title!("Model state at step $(model.n)")
        #text!(model[model.tracked_agent].pos, text = "$(model.tracked_agent)", align = (:center, :top))
        ###tracking the radial distance of each agent from group center
        radial_distances::Vector{Float64}  = []
        #=for i in 1:nagents(model)
                push!(radial_distances, radial_distance(model[i], com))
        end
        for i in 1:nagents(model)
                text!(new_pos[i], text = "$(trunc(radial_distances[i]))", align = (:center, :top))
        end=#
        Colorbar(figure[1,2], colourbarthing)
        save("./Better_Positions/shannon_flock_n_=_$(model.n).png", figure)
end

function draw_better_positions(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}, better_positions::Vector{Tuple{Float64, Float64}})
	figure, ax, colourbarthing = Makie.scatter([Tuple(point) for point in better_positions], axis = (;  limits = (0, rect_bound, 0, rect_bound), aspect = 1), marker = :circle, color = (:blue, 0.5))
	Makie.scatter!([model[i].pos for i in 1:nagents(model)], marker = :circle,  color = :black)
	save("./Better_Positions/shannon_flock_n_=_$(model.n).png", figure)
end
