#Function that spawns agents into cells to prevent married couples. spanx, spany represents the x and y dimension sizes respectively of the spawn space, cellx celly represent the minimum cell dimensions for each agent, offset represents the offset from (0,0) you wanna spawn. 
function assign_positions(cellx, celly, no_agents, spacexspan, spaceyspan, offsetx, offsety, initial_positions)
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
end
