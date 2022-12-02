###Introduction
#In this template, we generate an agent, but now it moves with some velocity and acceleration. However, the acceleration will be randomised at each time step and the velocity updated accordingly
#One thing that we're looking out for in particular is if we can update position, velocity and acceleration, which will be N tuples, with vectors generated from our ODE or movement gradient function. Actually no, in accordance with RK4, or even Euler integration, acceleration isn't an quantity we need to keep track of.

###Preliminaries
using Agents
using Random
using VoronoiCells
using GeometryBasics
using Plots

print("Packages loaded\n")

include("half_plane_alt.jl")
include("convex_hull.jl")

print("Both homemade files included\n")

rho = 100.0
initialised = 0
area_zero = zeros(Int64, 100)
rect = Rectangle(Point2(0,0), Point2(200, 200))
rect_bound = 200.0
moves_areas = [] #This is an array which will allow us to record all the areas and directions considered for each step, for each agent
no_move = ones(Int64, 100) #An array which will allow us to keep track of which agents never move
new_pos = [] #An array that will store the new positions of the agents for movement when we go to the model step
convex_hull_point = zeros(Int64, 100)
last_half_planes = []
D = 9
sigma = 0.0

###Function that takes a vector and calculates the mean of the elements in the vector
function mean(v)
	total = 0.0
	for i in 1:length(v)
		total += v[i]
	end
	
	return total/length(v)
end


###Function that calculates the area of a voronoi cell given the vertices that
#comprise the cell.
function voronoi_area(ri, cell, rho)
       	Area = 0.0
	circle_area = 0.0
	circle_detected = 0
	segment_detected = 0
	balloon_detected = 0
	num_points = length(cell)
	if(num_points == 0)
		Area = pi*rho^2
		return Area
	end
	
	#=
	print(" The vertices for the cell are ")
	for i in 1:num_points
        	vector_to_vertex = cell[i][1] .- ri
        	angle_to_vertex = atan(vector_to_vertex[2], vector_to_vertex[1])
        	print("$angle_to_vertex ")
                #print("$(atan(cell[i][1][2], cell[i][1][1])) ")
        end
                                print("\n")
	=#

	#Iterate through successive pairs of vertices in the cell
	for i in 1:length(cell)
		#Use the shoestring formula to calcualte the area
		j = (i)%num_points+1
		xi = cell[i][1][1]
		yi = cell[i][1][2]
		xj = cell[j][1][1]
		yj = cell[j][1][2]
                Area += 0.5 * (yi + yj)* (xi - xj)
		#If the two vertices are actually intersects with the circle, then in addition to the area calculated from the shoestring formula, you should also add the area of the circle segment 
		if(cell[i][3] == 0 && cell[j][2] == 0) #If the forward line segment for intersect i aand the backward lines segment for intersect j is a circle, then we have a chord
			#print("Circle segments detected\n")
			circle_detected = 1
			chord_length = norm(cell[j][1] .- cell[i][1]) #Calculates the length of the chord between the two vertices lying on the bounding circle
			r = sqrt(rho^2 - (0.5 * chord_length)^2)
			h = rho - r
			circle_segment_area = rho^2*acos((rho-h)/rho) - (rho-h)*sqrt(2*rho*h-h^2) #Calculated according to Wolfram formula 
			
			vec_to_i = cell[i][1] .- ri
			vec_to_ip1 = cell[j][1] .- ri
			angle_to_i = atan(vec_to_i[2], vec_to_i[1])
			angle_to_ip1 = atan(vec_to_ip1[2], vec_to_ip1[1])
			theta = min((angle_to_ip1 - angle_to_i + 2*pi)%(2*pi), (angle_to_i - angle_to_ip1 + 2*pi)%(2*pi))
			alt_circle_segment_area = 0.5* rho^2 * (theta - sin(theta))
			if(abs(alt_circle_segment_area - circle_segment_area)/circle_segment_area > 0.01)
				print("Divergence in area calculated, the angle formula calculated $alt_circle_segment_area, the other $circle_segment_area")
				#exit()
			end

			#Check, if the agent position is inside the chord half plane. 
			chord_vector = cell[j][1] .- cell[i][1]
			chord_point = 0.5 .* chord_vector .+ cell[i][1]
			chord_half_plane = [atan(chord_vector[2], chord_vector[1]), chord_vector, chord_point, 0]
			if(num_points == 2)
                        	return pi*rho^2 - circle_segment_area
                        end
			if(outside(chord_half_plane, ri))
				balloon = pi*rho^2 - circle_segment_area
				balloon_detected = 1
					
				print("Ballon segment detected, balloon area was $balloon.\n")
				#=
				for i in 1:num_points
					vector_to_vertex = cell[i][1] .- ri
					angle_to_vertex = atan(vector_to_vertex[2], vector_to_vertex[1])
					print("$angle_to_vertex ")
				end
				
				print("\n")
				=#
				
				circle_area += balloon
				#exit()
			else 
				segment_detected = 1
				circle_area += circle_segment_area
			end
		end
	end
		#=	
		if(abs(Area)+circle_area > pi*rho^2 && initialised == 0)
			print("Conventional area exceeded for agent, circle detected? $circle_detected. Balloon detected? $balloon_detected. Segment detected? $segment_detected. The circle area was $circle_area and the normal area was $(abs(Area)).\n")
			exit()
		end
		=#
		return  abs(Area)+circle_area
end



###Function that determines the gradient of movement
function move_gradient(agent, model,  kn, q, m, rho)
	#Calculate the unit vector in the current direction of motion
	dt = model.dt
	unit_v = agent.vel ./ 1.0
	theta_0 = atan(unit_v[2], unit_v[1])
	agent_speed = 1.0
	vix = unit_v[1]
	viy = unit_v[2]
	positions = []
	all_agents_iterable = allagents(model)
	for neighbour in all_agents_iterable
		if(neighbour.id == agent.id)
			continue
		end
		pushfirst!(positions, neighbour.pos)	
	end		
	min_area = agent.A #The agent's current DOD area
	min_direction = [0.0, 0.0] #This is to set it so that the default direction of move is nowehere (stay in place)
	move_made = 0
	pos_area_array = []
	no_angles_considered = 0

	#Iterate through all the possible places the agent can move, keeping track of which one minimises area assuming static neighbour positions, though we make sure that if none of the moves optimises the current area, don't move at all
	#print("For agent $(agent.id), its min area is $min_area \n")
	temp_hp = []
	for i in 0:(q-1) #For every direction
		conflict = 0
		direction_of_move = [cos(i*2*pi/q)*vix - sin(i*2*pi/q)*viy, sin(i*2*pi/q)*vix + cos(i*2*pi/q)*viy]
		angle_of_move = atan(direction_of_move[2], direction_of_move[1])
		rel_angle = ((angle_of_move - theta_0 + pi)+2*pi)%(2*pi) - pi
		if(abs(rel_angle) > (1)*2*pi/q + eps)
			continue
		end
		no_angles_considered += 1
		for j in 1:m #For every position up to m
			new_agent_pos = agent.pos .+ j .* direction_of_move .* agent_speed .* dt
		
			#Check first if there are no other agents in the potential position, note that we don't need to keep updating nearest neighbours since we assume the neighbours of a given agent are static
			for neighbour_position in positions
				if norm(new_agent_pos .- neighbour_position) < 2 #If moving in this direction and this m causes a collision, don't consider a move in this direction
					conflict = 1
					break
				end			
			end
			
			if (conflict == 1)		
				continue
			end

			#If there are no other agents in the potential position (no conflicts), go ahead and evaluate the new DOD
                	agent_voronoi_cell = voronoi_cell(new_agent_pos, positions, rho, temp_hp) #Generates the set of vertices which define the voronoi cell
                	new_area = voronoi_area(new_agent_pos, agent_voronoi_cell, rho) #Finds the area of the agent's voronoi cell
			#=		
			print("\n\n\nThe dq for this position was \n")
			for i in 1:length(temp_hp)
				print("$(temp_hp[i])\n")
			end
			=#

			#print("Potential new area of $new_area\n")
			#=
			#print("The vertices of this convex hull point are\n")
			if(convex_hull_point[agent.id] == 1)
				for i in 1:length(agent_voronoi_cell)
				 vector_to_vertex = agent_voronoi_cell[i][1] .- new_agent_pos
                                 angle_to_vertex = atan(vector_to_vertex[2], vector_to_vertex[1])
                                 print("$angle_to_vertex ")
                                 #print("$(atan(cell[i][1][2], cell[i][1][1])) ")
				end
			end
			print("\n")
			=#
			if (new_area < min_area)
                        	min_area = new_area
				print("New min area of $min_area, direction of $direction_of_move\n")
                        	min_direction = direction_of_move
                        	move_made = 1
				replace_vector(last_half_planes[Int64(agent.id)], [agent_voronoi_cell, temp_hp, new_agent_pos])
				if(convex_hull_point[agent.id] == 1)
					print("Min area was lowered for agent $(agent.id), in a potential position of $(new_agent_pos),  here is the temp_hp\n")
					for i in 1:length(temp_hp)
                                		print("$(temp_hp[i])\n")
					end
					print("Also, here's the vertices of the cell\n")
					print("$(agent_voronoi_cell)\n")
				end
                	end

		end
		
		#Check area calculation through voronoi package
		#=
		pack_positions = Vector{Point2{Float64}}(undef, nagents(model)) 
		for i in 1:nagents(model)-1
			pack_positions[i+1] = Point2(positions[i])
		end
		pack_positions[1] = Point2(pot_new_pos)
		tess = voronoicells(pack_positions, rect)
		tess_areas = voronoiarea(tess)
		if(abs(new_area-tess_areas[1]) > 0.1)
			print("Area check, our calculated area was $new_area, theirs was $(tess_areas[1])\n")
		end
		=#
		
		#push!(pos_area_array, [angle_of_move, min_area])
	end

	#push!(moves_areas[agent.id], [model.n, agent.A, pos_area_array])
	if(move_made == 1)
		no_move[agent.id] = 0
	end
	
	print("The number of angles considered was $no_angles_considered\n")
	#It really doesn't have to be like this, since  at least just for the simple SHH model of Dr.Algar, we can simply return a velocity
	kn[1] = (min_direction .* agent_speed)[1]
	kn[2] = (min_direction .* agent_speed)[2]
	
	#Create the noise addition
	epsilon = randn(model.rng, Float64, 2)
	epsilon_prime = randn(model.rng, Float64, 2)
	dW = sqrt(2*D*model.dt) .* (epsilon .- epsilon_prime)

	#Store the new position for updating in model step
	new_pos[agent.id] = Tuple(min_direction .* agent_speed .* model.dt .+ agent.pos .+ sigma*dW)
	if(new_pos[agent.id][1] > rect_bound || new_pos[agent.id][1] < 0.0 || new_pos[agent.id][2] > rect_bound || new_pos[agent.id][2] < 0.0)
		print("Agent $(agent.id) will step overbounds. This is for time step $(model.n), was the particle part of the convex hull? $(convex_hull_point[agent.id])\n")
		exit()
	end

	if(min_area > pi*rho^2)
		print("Conventional area exceeded by agent $(agent.id)\n")
	end
	return move_made
end



###Function for updating the convex hull, returns the points of the convex hull so the area of the convex hull can be calculated
function update_convex_hull(model)
	points = []
	all_agents_iterable = allagents(model)
	for agent in all_agents_iterable
		push!(points, [agent.pos[1], agent.pos[2], agent.id])
		#print("The point added was $([agent.pos[1], agent.pos[2], Int64(agent.id)]), the agent id is $(agent.id)\n")
	end
	
	for i in 1:nagents(model)
		convex_hull_point[i] = 0
	end

	CH = convex_hull(points)
	CH_for_area = []
	
	for point in CH
		convex_hull_point[Int64(point[3])] = 1
		push!(CH_for_area, [[point[1], point[2]], 1, 1]) #Note that the points used for the CH_for_area consists of the actual point itself and 1,1 because the voronoi area calcualtor we're going to use requires a 2 and 3 index, checking them to see if the points are circular. 
	end

	return CH_for_area
end



###Create the agent
mutable struct bird <: AbstractAgent
	id::Int
	pos::NTuple{2, Float64}
	vel::NTuple{2, Float64}
	A::Float64 #The area of the agent's DOD
end
	


print("Agent template created\n")



###Create the initialisation function
using Random #for reproducibility
function initialise(; seed = 123, no_birds = 4)
	#Create the space
	space = ContinuousSpace((200.0, 200.0); periodic = true)
	#Create the properties of the model
	properties = Dict(:t => 0.0, :dt => 1.0, :n => 0, :CHA => 0.0)
	
	#Create the rng
	rng = Random.MersenneTwister(seed)
	
	print("Before model\n")

	#Create the model
	model = ABM(
		bird, space; 
		properties, rng, scheduler = Schedulers.fastest
	)	


	#Generate random initial positions for each bird, then calculate the DoDs
	initial_positions = []
	temp_hp = []
	pack_positions = Vector{Point2{Float64}}(undef, no_birds)
	#=
	for i in 1:no_birds
		rand_position = Tuple(100*rand(Float64, 2)) .+ (50.0, 50.0) 
		push!(initial_positions, rand_position)
		pack_positions[i] = Point2(rand_position)
		push!(moves_areas, [])
		push!(last_half_planes, [])
		push!(new_pos, (0.0, 0.0))
	end
	=#

	rand_position = (50.0, 50.0)
	push!(initial_positions, rand_position)
        pack_positions[1] = Point2(rand_position)
        push!(moves_areas, [])
        push!(last_half_planes, [])
        push!(new_pos, (50.0, 50.0))

	rand_position = (150.0, 50.0)
        push!(initial_positions, rand_position)
        pack_positions[2] = Point2(rand_position)
        push!(moves_areas, [])
        push!(last_half_planes, [])
        push!(new_pos, (150.0, 50.0))

	rand_position = (150.0, 150.0)
        push!(initial_positions, rand_position)
        pack_positions[3] = Point2(rand_position)
        push!(moves_areas, [])
        push!(last_half_planes, [])
        push!(new_pos, (150.0, 150.0))

	rand_position = (50.0, 150.0)
        push!(initial_positions, rand_position)
        pack_positions[4] = Point2(rand_position)
        push!(moves_areas, [])
        push!(last_half_planes, [])
        push!(new_pos, (50.0, 150.0))


	#Calculate the DOD based off the initial positions
	init_tess = voronoicells(pack_positions, rect)
	init_tess_areas = voronoiarea(init_tess)

	#Calculate the DoDs based off the initial positions
	#initial_dods = voronoi_area(initial_positions, rho)
	initial_dods = []
	for i in 1:no_birds
		print("\n\nCalculating initial DOD for agent $i, at position $(initial_positions[i]).")
		ri  = initial_positions[i]
		neighbouring_positions = []
		for j in 1:no_birds
			if(i == j)
				continue 
			end
			push!(neighbouring_positions, initial_positions[j])
		end
		initial_cell = voronoi_cell(ri, neighbouring_positions, rho, temp_hp)
		initial_A = voronoi_area(ri, initial_cell, rho) 
		
		replace_vector(last_half_planes[i], [initial_cell, temp_hp, ri])
			
		print("Initial DOD calculated to be $initial_A\n")
		if(abs(initial_A) > pi*rho^2)
			print("Conventional area exceeded by agent $(i)\n")
			exit()
		elseif initial_A < eps
			print("Effective area of 0. The cell was comprised of vertices $(initial_cell)\n")
			
			area_zero[i] = 1
		end
		if(abs(initial_A-init_tess_areas[i]) > eps)
			print("Difference in area calculated between our code and the voronoi package. Our code calculated $initial_A, theirs $(init_tess_areas[i])\n")
		end
		push!(initial_dods, initial_A)
	
			
		print("The initial half planes for agent $(i) is \n")
		print("$(last_half_planes[i][2])\n")

		print("The initial vertices for agent $i is \n")
		print("$(last_half_planes[i][1])\n")


	end
	#Now make the agents with their respective DoDs and add to the model
	total_area = 0.0
	#=
	for i in 1:no_birds
		agent = bird(i, initial_positions[i], Tuple(rand(Float64, 2)), initial_dods[i])
		agent.vel = agent.vel ./ norm(agent.vel)
		#print("Initial velocity of $(agent.vel) \n")
		add_agent!(agent, initial_positions[i], model)
		total_area += initial_dods[i]
	end	
	=#

	agent = bird(1, initial_positions[1], (-1.0, -1.0), initial_dods[1])
        agent.vel = agent.vel ./ norm(agent.vel)
        #print("Initial velocity of $(agent.vel) \n")
        add_agent!(agent, initial_positions[1], model)
        total_area += initial_dods[1]

	agent = bird(2, initial_positions[2], (1.0, -1.0), initial_dods[2])
        agent.vel = agent.vel ./ norm(agent.vel)
        #print("Initial velocity of $(agent.vel) \n")
        add_agent!(agent, initial_positions[2], model)
        total_area += initial_dods[2]

	agent = bird(3, initial_positions[3], (1.0, 1.0), initial_dods[3])
        agent.vel = agent.vel ./ norm(agent.vel)
        #print("Initial velocity of $(agent.vel) \n")
        add_agent!(agent, initial_positions[3], model)
        total_area += initial_dods[3]


	agent = bird(4, initial_positions[4], (-1.0, 1.0), initial_dods[4])
        agent.vel = agent.vel ./ norm(agent.vel)
        #print("Initial velocity of $(agent.vel) \n")
        add_agent!(agent, initial_positions[4], model)
        total_area += initial_dods[4]


	#Calculate the actual area of the convex hull of the group of birds
	convexhullbro = update_convex_hull(model)
	initial_convex_hull_area = voronoi_area(-1, convexhullbro, rho)
	model.CHA = initial_convex_hull_area
	packing_fraction = nagents(model)*pi*1^2/model.CHA
	print("Packing fraction at n = 0 is $(packing_fraction)\n")
	write(compac_frac_file, "$packing_fraction ")
	average_area = total_area / nagents(model)
        write(mean_a_file, "$average_area ")
	print("Initialisation complete. \n\n\n")
	global initialised = 1
	
	#=
	Plots.scatter(pack_positions, markersize = 6, label = "generators")
annotate!([(pack_positions[n][1] + 0.02, pack_positions[n][2] + 0.03, Plots.text(n)) for n in 1:no_birds])
display(Plots.plot!(init_tess, legend=:topleft))
savefig("voronoi_pack_init_tess.png")
	=#
	#Finally, plot the figure
	#=	
	figure, _ = abmplot(model)	
	save("./Simulation_Images/shannon_flock_n_=_$(0)", figure)
	=#
	figure = Makie.scatter([Tuple(point) for point in initial_positions], axis = (; limits = (0, 200, 0, 200)))
        for i in 1:nagents(model)
                text!(initial_positions[i], text = "$i", align = (:center, :top))
        end
        save("./Simulation_Images/shannon_flock_n_=_$(0).png", figure)


	return model
end  



###Create the agent step function. This will call upon the force or acceleration function. I'm assuming that this function will be applied to each agent in turn
function agent_step!(agent, model)		
	#Update the agent position and velocities, but only if it is a 
	#print("Step!\n", agent.planet)
	dt = model.dt
	k1 = [0.0, 0.0, 0.0, 0.0]

	#Update the slopes, note that technically we should be using the vectorised dot operators, but Julia seems to allow us to be lazy when working with vectors
        #Now, why have we separated the position and velocity as two different vectors unlike PHYS4070? Because the pos is intrinsically a 2D vector for Julia Agents.
        move_made = move_gradient(agent, model, k1, 8, 100, rho)
	
	#Update the agent position and velocity
	new_agent_pos = Tuple(agent.pos .+ dt .* k1[1:2])
        new_agent_vel = Tuple(k1[1:2]) #So note that we're not doing incremental additions to the old velocity anymore, and that's because under Shannon's model, the velocity is just set automatically to whatever is needed to go to a better place. 
	change_in_position = new_agent_pos .- (agent.pos)
	if(move_made==1)
		agent.vel = new_agent_vel
	else 
		#print("No movement made, agent area was $(agent.A)\n")
	end
	#print("New agent pos of $new_agent_pos representing change of $change_in_position\n")
	#print(k1, "\n")
	#print(new_agent_pos, new_agent_vel, "\n")
	#move_agent!(agent, new_agent_pos, model)	
end
	



###Create the model_step function
function model_step!(model)
        all_agents_iterable = allagents(model)
        for agent in all_agents_iterable
                move_agent!(agent, Tuple(new_pos[agent.id]), model)
        end

        #Now recalculate the agent DODs based off their new positions
        total_area = 0.0
	temp_hp = []
	for agent_i in all_agents_iterable
                neighbour_positions = []
                for agent_j in all_agents_iterable
                        if(agent_i.id == agent_j.id)
                                continue
                        end
                        push!(neighbour_positions, agent_j.pos)
                end
                ri = agent_i.pos
                new_cell_i = voronoi_cell(ri, neighbour_positions, rho, temp_hp)
                new_area = voronoi_area(ri, new_cell_i, rho)
                agent_i.A = new_area
		if(agent_i.A > pi*rho^2)
                        print("Conventional area exceeded for agent, circle detected? $circle_detected. Balloon detected? $balloon_detected. Segment detected? $segment_detected. The circle area was $circle_area and the normal area was $(abs(Area)).\n")
                        exit()
                end
		total_area += agent_i.A
        end
        
	#Now update the model's convex hull
	convexhullbro = update_convex_hull(model)
	convex_hull_area = voronoi_area(-1, convexhullbro, rho)
	model.CHA = convex_hull_area
	model.t += model.dt
        model.n += 1

	#Finally, plot the model after the step
	#figure, _ = abmplot(model)
	print("\n\n\nThe number of points in new_pos is $(length(new_pos)), the first element is $(new_pos[1])\n")
	figure = Makie.scatter([Tuple(point) for point in new_pos], axis = (; limits = (0, 200, 0, 200)))
	for i in 1:nagents(model)
		text!(new_pos[i], text = "$i", align = (:center, :top))
	end
	save("./Simulation_Images/shannon_flock_n_=_$(model.n).png", figure)
	packing_fraction = nagents(model)*pi/model.CHA
	print("Packing fraction at n = $(model.n) is $(packing_fraction)\n")
	if(model.n < no_steps)
		write(compac_frac_file, "$packing_fraction ")
	else
		write(compac_frac_file, "$packing_fraction")
	end
	average_area = total_area / nagents(model)
	if(model.n < no_steps)
		write(mean_a_file, "$average_area ")
	else 
		write(mean_a_file, "$average_area")
	end

	last_hp_vert = open("Last_hp_vert.txt", "w")
	for i in 1:nagents(model)
		write(last_hp_vert, "Agent $i, position of $(new_pos[i]), considering position of $(last_half_planes[i][3])\n")
		write(last_hp_vert, "$(last_half_planes[i][1])\n")
		write(last_hp_vert, "$(last_half_planes[i][2])\n")
		write(last_hp_vert, "\n\n")
	end
	close(last_hp_vert)
end

###Test the model has been initialised and works
using InteractiveDynamics
using CairoMakie # choosing a plotting backend


#=
###Initialise the model
model = initialise()
print("Number of agents is $(nagents(model))\n")

#=
###Test the model has been initialised and works
using InteractiveDynamics
using CairoMakie # choosing a plotting backend
=#


figure, _ = abmplot(model)
figure # returning the figure displays it
save("shannon_flock.png", figure)
=#	


###Animate
#model = initialise();

#=
abmvideo(
    "Shannon_flock.mp4", model, agent_step!, model_step!;
    framerate = 4, frames = 120,
    title = "Shannon flock"
)
=#


compac_frac_file = open("compaction_frac.txt", "w")
mean_a_file = open("mean_area.txt", "w")
no_steps = 10
no_simulations = 1
for i in 1:no_simulations
	model = initialise()
	#figure, _ = abmplot(model)
        #save("./Simulation_Images/shannon_flock_n_=_$(0).png", figure)
	step!(model, agent_step!, model_step!, no_steps)
	write(compac_frac_file, "\n")
	write(mean_a_file, "\n")
end

close(compac_frac_file)
close(mean_a_file)

compac_frac_file = open("compaction_frac.txt", "r")
mean_a_file = open("mean_area.txt", "r")


cf_array = []
ma_array = []

for i in 0:no_steps
	push!(cf_array, [])
	push!(ma_array, [])
end

cf_lines = readlines(compac_frac_file)
ma_lines = readlines(mean_a_file)


print("The first thing read from the compac_frac_file was $(cf_lines[1])\n")
for line in cf_lines
	split_line = parse.(Float64, split(line, " "))
        for i in 1:length(split_line)
		print("The element read was $(split_line[i])\n")
		push!(cf_array[i], split_line[i])
        end
end

for line in ma_lines
        split_line = parse.(Float64, split(line, " "))
        for i in 1:length(split_line)
                push!(ma_array[i], split_line[i])
        end
end

cf_ave_file = open("cf_ave.txt", "w")
ma_ave_file = open("ma_ave.txt", "w")

for i in 0:no_steps
	write(cf_ave_file, "$i $(mean(cf_array[i+1]))\n")
	write(ma_ave_file, "$i $(mean(ma_array[i+1]))\n")
end



