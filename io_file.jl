using Plots
using InteractiveDynamics
using CairoMakie # choosing a plotting backend
using ColorSchemes
import ColorSchemes.balance
#using Agents
using Random

include("agent_definition.jl")

#=
###Animate
model = initialise(1000.0*sqrt(12), 1);
#ac(agent) =  get(balance, abs(agent.A-model.target_area)/(pi*rho^2))
plotkwargs = (; ac = get(balance, 0.7), as  = 10, am = :diamond)

print("Gone past the thang")

abmvideo(
    "Colour_Test.mp4", model, agent_step!, model_step!;
    spf = 1,
        framerate = 12, frames = 24,
    title = "Shannon flock",
        showstep = true,
        ac = get(balance, 0.7), as = 10, am = :diamond
)

print("Finished the vid\n")
=#


function do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)
        close(compac_frac_file)
        close(mean_a_file)
        close(rot_o_file)
        close(rot_o_alt_file)
        close(mean_speed_file)

        compac_frac_file = open("compaction_frac.txt", "r")
        mean_a_file = open("mean_area.txt", "r")
        rot_o_file = open("rot_order.txt", "r")
        rot_o_alt_file = open("rot_order_alt.txt", "r")
        mean_speed_file = open("mean_speed.txt", "r")

cf_array = []
ma_array = []
rot_o_array = []
rot_o_alt_array = []
ms_array = []

for i in 0:no_steps
        push!(cf_array, [])
        push!(ma_array, [])
        push!(rot_o_array, [])
        push!(rot_o_alt_array, [])
        push!(ms_array, [])
end

cf_lines = readlines(compac_frac_file)
ma_lines = readlines(mean_a_file)
rot_o_lines = readlines(rot_o_file)
rot_o_alt_lines = readlines(rot_o_alt_file)
ms_lines = readlines(mean_speed_file)

print("The first thing read from the compac_frac_file was $(cf_lines[1])\n")
for line in cf_lines
        split_line = parse.(Float64, split(line, " "))
        for i in 1:length(split_line)
                #print("The element read was $(split_line[i])\n")
                push!(cf_array[i], split_line[i])
        end
end

for line in ma_lines
        split_line = parse.(Float64, split(line, " "))
        for i in 1:length(split_line)
                push!(ma_array[i], split_line[i])
        end
end

for line in rot_o_lines
        split_line = parse.(Float64, split(line, " "))
        for i in 1:length(split_line)
                push!(rot_o_array[i], split_line[i])
        end
end

for line in rot_o_alt_lines
        split_line = parse.(Float64, split(line, " "))
        for i in 1:length(split_line)
                push!(rot_o_alt_array[i], split_line[i])
                #print("The split line was $(split_line[i])\n")
        end
end

for line in ms_lines
        split_line = parse.(Float64, split(line, " "))
        for i in 1:length(split_line)
                push!(ms_array[i], split_line[i])
                #print("The split line was $(split_line[i])\n")
        end
end

cf_ave_file = open("cf_ave.txt", "w")
ma_ave_file = open("ma_ave.txt", "w")
rot_o_ave_file = open("rot_o_ave.txt", "w")
rot_o_alt_ave_file = open("rot_o_alt_ave.txt", "w")
mean_speed_file = open("mean_speed_ave.txt", "w")

for i in 0:no_steps
        write(cf_ave_file, "$i $(mean(cf_array[i+1]))\n")
        write(ma_ave_file, "$i $(mean(ma_array[i+1]))\n")
        write(rot_o_ave_file, "$i $(mean(rot_o_array[i+1]))\n")
        write(rot_o_alt_ave_file, "$i $(mean(rot_o_alt_array[i+1]))\n")
        write(mean_speed_file, "$i $(mean(ms_array[i+1]))\n")
end
	#=
	close(compac_frac_file)
        close(mean_a_file)
        close(rot_o_file)
        close(rot_o_alt_file)
        close(mean_speed_file)
	=#
	
	close(cf_ave_file)
	close(ma_ave_file)
	close(rot_o_ave_file)
	close(rot_o_alt_ave_file)
	close(mean_speed_file)
end



###Function for drawing the plots for model step
function draw_figures(model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}, actual_areas::Vector{Float64}, previous_areas::Vector{Float64}, delta_max::Float64, new_pos::Vector{Tuple{Float64, Float64}})
	##Draw the standard figure of the agents with their DODs after the step
	colours::Vector{Float64} = []
        rotations::Vector{Float64} = []
        allagents_iterable = allagents(model)
	target_area::Float64 = model.target_area
        for id in 1:nagents(model)
                push!(colours, abs(model[id].A-model.target_area)/(delta_max))
                push!(rotations, atan(model[id].vel[2], model[id].vel[1]))
        end
        #figure, _ = abmplot(model)
        print("\n\n\nThe number of points in new_pos is $(length(new_pos)), the first element is $(new_pos[1])\n")
        #print("About to do the figure\n")
        figure, ax, colourbarthing = Makie.scatter([Tuple(point) for point in new_pos], axis = (; limits = (0, rect_bound, 0, rect_bound)), marker = '→', markersize = 20, rotations = rotations, color = colours, colormap = :viridis, colorrange = (0.0, 1.0)) #Note that I have no idea what the colorbarthing is for
        #=for i in 1:nagents(model)
                text!(new_pos[i], text = "$i", align = (:center, :top))
        end=#
        Colorbar(figure[1,2], colourbarthing)
        save("./Simulation_Images/shannon_flock_n_=_$(model.n).png", figure)


	##Draw the figure of the agents with their actual DODs
	for id in 1:nagents(model)
                colours[id] = actual_areas[id]/(pi*rho^2)
        end
	figure_actual, ax, colourbarthing = Makie.scatter([Tuple(point) for point in new_pos], axis = (; limits = (0, rect_bound, 0, rect_bound)), marker = '→', markersize = 20, rotations = rotations, color = colours, colormap = :viridis, colorrange = (0.0, 1.0)) #Note that I have no idea what the colorbarthing is for
        #=for i in 1:nagents(model)
                text!(new_pos[i], text = "$i", align = (:center, :top))
        end=#
        Colorbar(figure[1,2], colourbarthing)
        save("./Simulation_Images_Actual_Areas/shannon_flock_n_=_$(model.n).png", figure)

	##Draw the figure of the agents with their change in DOD
	for id in 1:nagents(model)
                #print("Current A is $(model[id].A), previous areas was $(previous_areas[id])\n")
		colours[id] = (abs(model[id].A - model.target_area)-abs(previous_areas[id]-model.target_area)+delta_max)/(2*delta_max)
        end
        figure_actual, ax, colourbarthing = Makie.scatter([Tuple(point) for point in new_pos], axis = (; limits = (0, rect_bound, 0, rect_bound)), marker = '→', markersize = 20, rotations = rotations, color = colours, colormap = :viridis, colorrange = (0.0, 1.0)) #Note that I have no idea what the colorbarthing is for
        #=for i in 1:nagents(model)
                text!(new_pos[i], text = "$i", align = (:center, :top))
        end=#
        Colorbar(figure[1,2], colourbarthing)
        save("./Simulation_Images_Difference_Areas/shannon_flock_n_=_$(model.n).png", figure)
end
