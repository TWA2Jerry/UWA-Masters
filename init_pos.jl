###Function that assigns an initial position to each agent such that they 
include("global_vars.jl")
include("some_math_functions.jl")
using Random

function assign_positions(cellx, celly, no_agents, spacexspan, spaceyspan, offsetx, offsety, initial_positions, initial_vels)
	init_grid = []

	#Find the number of cells given the minimum exclusion volume for each agent
	no_rows = spaceyspan/celly
	no_cols=  spacexspan/cellx


	#Iterate through the 
	for i in 0:no_rows-1
		#Iterate through the
		for j in 0:no_cols-1
			push!(init_grid, (offsetx+j*cellx + cellx/2, offsety+ i*celly + celly/2))
		end
	end

	#Iterate through all the agent ids
	for id in 1:no_agents
		#Pick a random number i  between 1 and the length of the vector
		rand_num = rand(1:length(init_grid))
		#Assign the i-th coordinate as the initial position of the agent
		push!(initial_positions, init_grid[rand_num])
		#Remove that coordinate from the vector
		deleteat!(init_grid, rand_num)
	end

	for i in 1:no_agents
                rand_vel::Tuple{Float64, Float64} = 2 .* Tuple(rand(Float64, 2)) .- (1.0, 1.0)
                rand_vel = rand_vel ./norm(rand_vel)
                push!(initial_vels, rand_vel)
	end
end


###Function that assigns initial positions based off previously established thingo
function assign_pos_vels(pos_vels_file, initial_positions, initial_vels, step, no_birds)
	print("Assign_pos_vels called\n")
	new_pos_vels_file = open("pos_vels.txt", "r")
	pos_vels_lines = readlines(new_pos_vels_file)
	line_counter::Int32 = 0	
	for line in pos_vels_lines
		print("Line read\n")
		if(line_counter == step)
			split_line = parse.(Float64, split(line, " "))
			for i in 1:no_birds #Yeah this is real dodgy, maybe pass in arg of noagents
				x::Float64 = split_line[(i-1)*4 + 1]
				y::Float64 = split_line[(i-1)*4 + 2]
				vx::Float64 = split_line[(i-1)*4 + 3]
				vy::Float64 = split_line[(i-1)*4 + 4]
				push!(initial_positions, (x, y))
				print("Pushed $((x,y)), initial pos is now $initial_positions\n")
				push!(initial_vels, (vx, vy))
			end
			print("The step number read in was $(split_line[length(split_line)])\n")
			break
		else 
			print("Line thingo incremented\n")
			line_counter += 1
			print("Line counter is $line_counter\n")
			continue
		end
	end
	close(new_pos_vels_file)
end

function init_circle(circle_origin, circle_radius, initial_positions, initial_vels)
	for i in 1:no_birds
                angle_per_bird = 2*pi/(no_birds)
                initial_pos = (circle_radius*cos(angle_per_bird*(i-1)), circle_radius*sin(angle_per_bird*(i-1))) .+ circle_origin
                rand_vel = (-circle_radius*sin(angle_per_bird*(i-1)), circle_radius*cos(angle_per_bird*(i-1))) 
                rand_vel = rand_vel ./norm(rand_vel)
                print("The bird $i's initial position is $initial_pos\n")
                push!(initial_positions, initial_pos)
                push!(initial_vels, rand_vel)
        end	
end
