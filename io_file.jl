using Plots
using InteractiveDynamics
using CairoMakie # choosing a plotting backend
using ColorSchemes
import ColorSchemes.balance
#using Agents
using Random


include("agent_definition.jl")
include("order_parameters.jl")
include("plot_histograms.jl")
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

#print("The first thing read from the compac_frac_file was $(cf_lines[1])\n")
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
	close(compac_frac_file)
        close(mean_a_file)
        close(rot_o_file)
        close(rot_o_alt_file)
        close(mean_speed_file)
	
	close(cf_ave_file)
	close(ma_ave_file)
	close(rot_o_ave_file)
	close(rot_o_alt_ave_file)
	close(mean_speed_file)
	print("Finished statistics IO\n")
end



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
		push!(colours, 0)
		push!(rotations, atan(model[id].vel[2], model[id].vel[1]))
        end
        #figure, _ = abmplot(model)
        print("\n\n\nThe number of points in new_pos is $(length(new_pos)), the first element is $(new_pos[1])\n")
        #print("About to do the figure\n")
        
	
	#figure, ax, colourbarthing = Makie.scatter([Tuple(point) for point in new_pos], axis = (; title = "Model state at step $(model.n)", limits = (0, rect_bound, 0, rect_bound), aspect = 1), marker = '→', markersize = 20, rotations = rotations, color = colours, colormap = cgrad(:matter, 300, categorical = true), colorrange = (0, 300))
	figure, ax, colourbarthing = Makie.scatter([Tuple(point) for point in new_pos], axis = (;  limits = (0, rect_bound, 0, rect_bound), aspect = 1), marker = :circle,  rotations = rotations, color = :black)


	#=
	for i in 1:nagents(model)
                text!(new_pos[i], text = "$i", align = (:center, :top))
        end
	=#	

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
	#Colorbar(figure[1,2], colourbarthing)
        save("./Simulation_Images/shannon_flock_n_=_$(model.n).png", figure)
end

function draw_path(points)
	Makie.scatter!([Tuple(point) for point in points], marker = :circle, color = :black, markersize = 5)
end


function do_more_io_stuff(adf, mdf)
	rot_o_alt_ave_file = open("ensemble_rot_o_alt.txt", "w")
	for i in 1:no_steps+1
        	write(rot_o_alt_ave_file, "$(i-1) $(mdf[i, 3])\n")
	end

	close(rot_o_alt_ave_file)

	###This section is for ensembles
	mean_happiness_file = open("mean_happiness.txt", "w")
	std_happiness_file = open("std_happiness.txt", "w")

	for step in 1:no_steps+1
        	new_mean::Float64 = 0.0
        	mean_squared::Float64 = 0.0
        	for sim_n in 0:no_simulations-1
                	new_mean += adf[sim_n*(no_steps+1)+step , 3]/no_simulations
                	mean_squared += (adf[sim_n*(no_steps+1)+step, 3])^2/no_simulations
                	#print("New mean was $new_mean, mean squared was $mean_squared\n")
        	end
        	write(mean_happiness_file, "$(step-1) $(new_mean)\n")
        	std_happiness = sqrt(mean_squared - new_mean^2)
        	#std_happiness = sqrt(mean_squared)
        	write(std_happiness_file, "$(step-1) $(std_happiness)\n")
	end

	radial_dist_file = open("circ_radial.txt", "w")
	for i in 1:no_steps+1
        	write(radial_dist_file, "$(i-1) $(mdf[i, 2])\n")
	end
	
	

	rand_happiness_file = open("rand_happpiness.txt", "w")
	for i in 1:no_steps+1
		write(rand_happiness_file, "$(i-1) $(mdf[i,4])\n")
	end

	no_moves_file = open("no_moves.txt", "w")
	for i in 1:no_steps+1
                write(no_moves_file, "$(i-1) $(mdf[i,5])\n")
        end

	polarisation_file = open("polarisation.txt", "w")
	for i in 1:no_steps+1
                write(polarisation_file, "$(i-1) $(mdf[i,6])\n")
        end

	###If I recall correctly, this is basically to write the entire mdf data struct onto the mdf_file. The starting from 2 in j is to account for the fact that the first column is just a column counter which we don't want
	mdf_file = open("mdf_file.txt", "w")
	for i in 1:size(mdf, 1)
		for j in 1:size(mdf, 2)
			if(j < size(mdf,2)) 
				write(mdf_file, "$(mdf[i, j]) ")
			else 
				write(mdf_file, "$(mdf[i, j])") 
			end
		end	
		write(mdf_file, "\n")
	end

	adf_file = open("adf_file.txt", "w")
	for i in 1:size(adf, 1)
                for j in 1:size(adf, 2)
                        if(j < size(adf,2))
                                write(adf_file, "$(adf[i, j]) ")
                        else
                                write(adf_file, "$(adf[i, j])")
                        end
                end
                write(adf_file, "\n")
        end

	close(rot_o_alt_ave_file)
	close(mean_happiness_file)
	close(std_happiness_file)
	close(radial_dist_file)
	close(rand_happiness_file)
	close(no_moves_file)
	close(polarisation_file)
	close(mdf_file)
	close(adf_file)
end

function step_statistics(no_move::Vector{Int64}, model_step::Int64, n::Int64, no_steps::Int64, no_birds::Int64)
	no_no_moves::Int64 = 0
	for i in 1:length(no_move)
		no_no_moves += no_move[i]/no_birds
	end	
	if(n < no_steps)
                write(no_moves_file, "$no_no_moves ")
		#write(compac_frac_file, "$packing_fraction ")
                #write(rot_o_file, "$rot_order ")
                #write(rot_o_alt_file, "$rot_order_alt ")
        else
                write(no_moves_file, "$no_no_moves\n")
		#write(compac_frac_file, "$packing_fraction")
                #write(rot_o_file, "$rot_order")
                #write(rot_o_alt_file, "$rot_order_alt")
        end
	
end

function write_rtave(mdf, tdods)
	###Do some tdod ave recording for radius
        radius_tave_file = open("radius_tave.txt", "w")
        #tdods = Set() #Create a set that will keep track of the tdods that were passed down
        #=for i in 1:size(mdf, 1)
                push!(tdods, mdf[i, size(mdf, 2)-1])
                print("Added to tdods a tdod of $(mdf[i, size(mdf, 2)-1])\n")
        end =#

        no_tdods::Int32 = length(tdods)
        print("The number of tdods registered is $no_tdods\n")
        for i in 1:no_tdods
                average_r::Float64 = 0.0
                r::Float64 = 0.0
                for j in 1:no_simulations
                        for k in trunc(Int, (no_steps+1)/2):no_steps+1
                                #print("The value of i, j k is $i $j $k\n")
                                #r = mdf[((i-1)*no_simulations*(no_steps+1))+(j-1)*(no_steps+1)+(k), 2]
                                r = mdf[(j-1)*no_tdods*(no_steps+1)+(i-1)*(no_steps+1)+(k), 2]
                                #print("The mean radius was registered as $r\n")
                                average_r += r
                        end
                end
                average_r /= (no_simulations*(no_steps+1-(no_steps+1)/2 + 1))
                write(radius_tave_file, "$(tdods[i]) $average_r\n")
        end

end

function write_pos_vel(positions, velocities, pos_vels_file, n)
	for i in 1:length(positions)
		if(i < length(positions))
			write(pos_vels_file, "$(positions[i][1]) $(positions[i][2]) $(velocities[i][1]) $(velocities[i][2]) ")
		else
			write(pos_vels_file, "$(positions[i][1]) $(positions[i][2]) $(velocities[i][1]) $(velocities[i][2]) $n\n")			
		end
	end
end

function parse_model_vals(read_lines)
	parsed_lines = []
	for line in read_lines
		split_line = split(line, " ")
		parsed_line = parse.(Float64, split_line)
		parsed_line[1] = Int64(trunc(parsed_line[1]))
		push!(parsed_lines, parsed_line)
	end
	return parsed_lines
end

function calc_mean_std(values)
	n::Int32 = length(values)
	mean::Float64 = 0.0
	mean_squared::Float64 = 0.0	
	for i in 1:n
		mean += values[i]/n
        	mean_squared += values[i]^2/n
	end
	std_dev::Float64 = sqrt(mean_squared - mean^2)
	return mean, std_dev
end

function average_across_thing(agent_vals_lines, dimension_to_average_across, var_of_interest)
	#Create a set that will hold the 
	s = Set([])
		
	#Iterate through all sample instances 
	for i in 1:length(agent_vals_lines)
		push!(s, agent_vals_lines[i][dimension_to_average_across])
	end

	#Create a map such that associated with each key, the values that the variable along which we average can take, is an array that holds all the values of the variable of interest for that value
	d = Dict([])
	for element in s
		d[element] = [] 
	end	

	#Now go through the lines again and add every lines value
	for line in agent_vals_lines
		xvar = line[dimension_to_average_across]
		push!(d[xvar], line[var_of_interest])
	end
		
	d_ave = Dict([])
	for (key, value) in d
		d_ave[key] = calc_mean_std(value)
	end
	
	return d_ave
end


function average_across_thing_data_frame(data_frame, dimension_to_average_across, var_of_interest)
        #Create a set that will hold the
        s = Set([])

        #Iterate through all sample instances
        for i in 1:size(data_frame, 1)
                push!(s, data_frame[i, dimension_to_average_across])
        end

        #Create a map such that associated with each key, the values that the variable along which we average can take, is an array that holds all the values of the variable of interest for that value
        d = Dict([])
        for element in s
                d[element] = []
        end

        #Now go through the lines again and add every lines value
        for i in 1:size(data_frame, 1)
		xvar = data_frame[i, dimension_to_average_across]
		push!(d[xvar], data_frame[i, var_of_interest])
        end

        d_ave = Dict([])
        for (key, value) in d
                d_ave[key] = calc_mean_std(value)
        end

        return d_ave
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
 
