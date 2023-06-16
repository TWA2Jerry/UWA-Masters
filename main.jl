###Introduction
#In this template, we generate an agent, but now it moves with some velocity and acceleration. However, the acceleration will be randomised at each time step and the velocity updated accordingly
#One thing that we're looking out for in particular is if we can update position, velocity and acceleration, which will be N tuples, with vectors generated from our ODE or movement gradient function. Actually no, in accordance with RK4, or even Euler integration, acceleration isn't an quantity we need to keep track of.

###Preliminaries
using Agents
using Random
using VoronoiCells
using GeometryBasics
using Plots
using InteractiveDynamics
using CairoMakie # choosing a plotting backend

print("Packages loaded\n")

###Create the agent
mutable struct bird <: AbstractAgent
        id::Int
        pos::NTuple{2, Float64}
        vel::NTuple{2, Float64}
        speed::Float64
        A::Float64 #The area of the agent's DOD, at least in their own eyes
	true_A::Float64 #The true area of the agent's DOD
end

include("half_plane_fast.jl")
include("half_plane_bounded.jl")
include("convex_hull.jl")
include("rot_ord.jl")
include("rot_ord_check.jl")
include("init_pos.jl")
print("All homemade files included\n")

const rho = 100.0
initialised = 0
area_zero = zeros(Int64, 100)
const rect_bound::Float64 = 1000.0
const spawn_dim_x::Float64 = 100.0 #This gives the x dimesnion size of the initial spawning area for the agents
const spawn_dim_y::Float64 = 100.0 #This gives the y dimension size of the initial spawning area for the agents
rect = Rectangle(Point2(0,0), Point2(Int64(rect_bound), Int64(rect_bound)))
moves_areas = [] #This is an array which will allow us to record all the areas and directions considered for each step, for each agent
no_move = ones(Int64, 100) #An array which will allow us to keep track of which agents never move
new_pos = [] #An array that will store the new positions of the agents for movement when we go to the model step
convex_hull_point = zeros(Int64, 100)
last_half_planes = []
const sigma = 0.0

###Function that takes a vector and calculates the mean of the elements in the vector
function mean(v)
	total = 0.0
	for i in 1:length(v)
		total += v[i]
	end
	
	return total/length(v)
end

include("voronoi_area_file.jl")
include("move_gradient_file.jl")

print("Agent template created\n")



###Create the initialisation function
using Random #for reproducibility
function initialise(target_area_arg, simulation_number_arg; seed = 123, no_birds = 100)
	#Create the space
	space = ContinuousSpace((rect_bound, rect_bound); periodic = true)
	#Create the properties of the model
	properties = Dict(:t => 0.0, :dt => 1.0, :n => 0, :CHA => 0.0, :target_area => target_area_arg, :simulation_number => simulation_number_arg)
	
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
	initial_vels = []
	temp_hp::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}}= []
	pack_positions = Vector{Point2{Float64}}(undef, no_birds)
	
	#Initialise the positions based on the spawn-error free function of assign_positions
	assign_positions(2.0, 2.0, no_birds, spawn_dim_x, spawn_dim_y, (rect_bound-spawn_dim_x)/2, (rect_bound-spawn_dim_x)/2, initial_positions)

	for i in 1:no_birds
		#rand_position = Tuple(100*rand(Float64, 2)) .+ (50.0, 50.0) 
		rand_vel = 2 .* Tuple(rand(Float64, 2)) .- (1.0, 1.0)
		rand_vel = rand_vel ./norm(rand_vel)
		#push!(initial_positions, rand_position)
		push!(initial_vels, rand_vel)
		pack_positions[i] = initial_positions[i]
		push!(moves_areas, [])
		push!(last_half_planes, [])
		push!(new_pos, (0.0, 0.0))
	end

	#Calculate the DOD based off the initial positions
	init_tess = voronoicells(pack_positions, rect)
	init_tess_areas = voronoiarea(init_tess)

	#Calculate the DoDs based off the initial positions
	#initial_dods = voronoi_area(model, initial_positions, rho)
	initial_dods = []
	true_initial_dods = []
	for i in 1:no_birds
		print("\n\nCalculatin initial DOD for agent $i, at position $(initial_positions[i]).")
		ri  = Tuple(initial_positions[i])
		neighbouring_positions = Vector{Tuple{Float64, Float64}}(undef, 0)
		for j in 1:no_birds
			if(i == j)
				continue 
			end
			push!(neighbouring_positions, Tuple(initial_positions[j]))
		end
		vix = initial_vels[i][1]
		viy = initial_vels[i][2]
		relic_x = -1.0*(-viy)
        	relic_y = -vix
        	relic_pq = (relic_x, relic_y)
        	relic_angle = atan(relic_y, relic_x)
        	relic_is_box = 2
        	relic_half_plane = (relic_angle, relic_pq, ri, relic_is_box)

		initial_cell = @time voronoi_cell_bounded(model, ri, neighbouring_positions, rho, eps, inf, temp_hp, initial_vels[i], relic_half_plane)
		initial_A = voronoi_area(model, ri, initial_cell, rho) 
	
		true_initial_cell = @time voronoi_cell(model, ri, neighbouring_positions, rho,eps, inf, temp_hp, initial_vels[i])
                true_initial_A = voronoi_area(model, ri, true_initial_cell, rho)


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
		push!(true_initial_dods, true_initial_A)
			
		#print("The initial half planes for agent $(i) is \n")
		#print("$(last_half_planes[i][2])\n")

		#print("The initial vertices for agent $i is \n")
		#print("$(last_half_planes[i][1])\n")


	end
	#Now make the agents with their respective DoDs and add to the model
	total_area = 0.0
	total_speed = 0.0
	for i in 1:no_birds
		agent = bird(i, initial_positions[i], initial_vels[i], 1.0, initial_dods[i], true_initial_dods[i])
		agent.vel = agent.vel ./ norm(agent.vel)
		#print("Initial velocity of $(agent.vel) \n")
		add_agent!(agent, initial_positions[i], model)
		total_area += true_initial_dods[i]/(pi*rho^2)
		total_speed += agent.speed
	end	

	#Calculate the actual area of the convex hull of the group of birds
	convexhullbro = update_convex_hull(model)
	initial_convex_hull_area = voronoi_area(model, -1, convexhullbro, rho)
	model.CHA = initial_convex_hull_area
	packing_fraction = nagents(model)*pi*1^2/model.CHA
	init_rot_ord = rot_ord(allagents(model))
	init_rot_ord_alt = rot_ord_alt(allagents(model))
	print("Packing fraction at n = 0 is $(packing_fraction)\n")
	write(compac_frac_file, "$packing_fraction ")
	average_area = total_area / nagents(model)
        write(mean_a_file, "$average_area ")
	average_speed = total_speed/no_birds
	write(mean_speed_file, "$average_speed ")
	write(rot_o_file, "$init_rot_ord ")
	write(rot_o_alt_file, "$init_rot_ord_alt ")
	print("Initialisation complete. \n\n\n")
	global initialised = 1
	

	###Plotting
	colours::Vector{Float64} = []
	allagents_iterable = allagents(model)
	for agent in allagents_iterable
		push!(colours, abs(agent.A-model.target_area)/(0.5*pi*rho^2))
	end	
	#=
	Plots.scatter(pack_positions, markersize = 6, label = "generators")
annotate!([(pack_positions[n][1] + 0.02, pack_positions[n][2] + 0.03, Plots.text(n)) for n in 1:no_birds])
display(Plots.plot!(init_tess, legend=:topleft))
savefig("voronoi_pack_init_tess.png")
	=#
	#Finally, plot the figure
	figure, ax, colorbarthing = Makie.scatter([Tuple(point) for point in initial_positions], axis = (; limits = (0, rect_bound, 0, rect_bound)), color = colours, colormap = :viridis, colorrange = (0.0, 1.0))
        #=for i in 1:nagents(model) #This is for labelling each dot with the agent number in plot
                text!(initial_positions[i], text = "$i", align = (:center, :top))
        end=#
	Colorbar(figure[1,2], colorbarthing)
        save("./Simulation_Images/shannon_flock_n_=_$(0).png", figure)


	return model
end  



###Create the agent step function. This will call upon the force or acceleration function. I'm assuming that this function will be applied to each agent in turn
function agent_step!(agent, model)		
	#Update the agent position and velocities, but only if it is a 
	#print("Step!\n", agent.planet)
	dt = model.dt
	k1::Vector{Float64} = [0.0, 0.0, 0.0, 0.0]
	target_area::Float64 = model.target_area	

        #Now, why have we separated the position and velocity as two different vectors unlike PHYS4070? Because the pos is intrinsically a 2D vector for Julia Agents.
        move_made = move_gradient(agent, model, k1, 8, 100, rho, target_area)
	
	#Update the agent position and velocity
	new_agent_pos = Tuple(agent.pos .+ dt .* k1[1:2])
        new_agent_vel = Tuple(k1[1:2]) #So note that we're not doing incremental additions to the old velocity anymore, and that's because under Shannon's model, the velocity is just set automatically to whatever is needed to go to a better place. 
	change_in_position = new_agent_pos .- (agent.pos)
	if(move_made==1)
		agent.vel = new_agent_vel
		agent.speed = 1.0
	else 
		#print("No movement made, agent area was $(agent.A)\n")
		agent.vel = new_agent_vel
		agent.speed = 0.0
	end
	#print("New agent pos of $new_agent_pos representing change of $change_in_position\n")
	#print(k1, "\n")
	#print(new_agent_pos, new_agent_vel, "\n")
	#move_agent!(agent, new_agent_pos, model)	
end
	



###Create the model_step function
function model_step!(model)
	#Calculate the rotational order of the agents. After some debate, we've decided that position \times desired_velocity is the way to go. 
        all_agents_iterable = allagents(model)
	rot_order = rot_ord(allagents(model))
        rot_order_alt = rot_ord_alt(allagents(model))
	print("Alternate rotational order returned as $rot_order_alt\n")	
	#Move the agents to their predetermined places 
	for agent in all_agents_iterable
                move_agent!(agent, Tuple(new_pos[agent.id]), model)
		print("Agent position is now $(agent.pos) for a new agent pos of $(new_pos[agent.id])\n")
        end
	
        #Now recalculate the agent DODs based off their new positions
        total_area = 0.0
	total_speed = 0.0
	temp_hp::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}} = []
	for agent_i in all_agents_iterable
                neighbour_positions::Vector{Tuple{Float64, Float64}} = []
                for agent_j in all_agents_iterable
                        if(agent_i.id == agent_j.id)
                                continue
                        end
                        push!(neighbour_positions, agent_j.pos)
                end
                ri = agent_i.pos
		vix = agent_i.vel[1]
		viy = agent_i.vel[2]
		relic_x = -1.0*(-viy)
        	relic_y = -vix
        	relic_pq = (relic_x, relic_y)
        	relic_angle = atan(relic_y, relic_x)
        	relic_is_box = 2
        	relic_half_plane = (relic_angle, relic_pq, agent_i.pos, relic_is_box)

                new_cell_i = voronoi_cell_bounded(model, ri, neighbour_positions, rho, eps, inf, temp_hp, agent_i.vel, relic_half_plane)
                new_area = voronoi_area(model, ri, new_cell_i, rho)
                agent_i.A = new_area
		if(agent_i.A > pi*rho^2)
			print("Conventional area exceeded for agent. Cell was $(new_cell_i), and area was $(new_area)\n")
                        exit()
                end
		#For measuring parameters, we measure the true voronoi cell, which will not use the bounded vision. 
		#print("\n\n\n The time for calulating the voronoi cell in model step is ")
		true_new_cell_i =  voronoi_cell(model, ri, neighbour_positions, rho,eps, inf, temp_hp, agent_i.vel)
                true_new_area = voronoi_area(model, ri, true_new_cell_i, rho)
		#print("The bounded DOD was calculated as $new_area, while the unbounded was calculated as $true_new_area\n")
		agent_i.true_A = true_new_area
		total_area += true_new_area/(pi*rho^2)
		total_speed += agent_i.speed
        end
        
	#Now update the model's convex hull
	convexhullbro = update_convex_hull(model)
	convex_hull_area = voronoi_area(model, -1, convexhullbro, rho)
	model.CHA = convex_hull_area
	model.t += model.dt
        model.n += 1

	#Finally, plot the model after the step
	colours::Vector{Float64} = []
        allagents_iterable = allagents(model)
        for agent in allagents_iterable
                push!(colours, abs(agent.A-model.target_area)/(0.5*pi*rho^2))
        end     
	#figure, _ = abmplot(model)
	#print("\n\n\nThe number of points in new_pos is $(length(new_pos)), the first element is $(new_pos[1])\n")
	figure, ax, colourbarthing = Makie.scatter([Tuple(point) for point in new_pos], axis = (; limits = (0, rect_bound, 0, rect_bound)), color = colours, colormap = :viridis, colorrange = (0.0, 1.0)) #Note that I have no idea what the colorbarthing is for
	#=for i in 1:nagents(model)
		text!(new_pos[i], text = "$i", align = (:center, :top))
	end=#
	Colorbar(figure[1,2], colourbarthing)
	save("./Simulation_Images/shannon_flock_n_=_$(model.n).png", figure)
	
	##Statistics recording
	packing_fraction = nagents(model)*pi/model.CHA
	print("Packing fraction at n = $(model.n) is $(packing_fraction)\n")
	if(model.n < no_steps)
		write(compac_frac_file, "$packing_fraction ")
		write(rot_o_file, "$rot_order ")
		write(rot_o_alt_file, "$rot_order_alt ")
	else
		write(compac_frac_file, "$packing_fraction")
		write(rot_o_file, "$rot_order")
		write(rot_o_alt_file, "$rot_order_alt")
	end
	average_area = total_area / nagents(model)
	average_speed = total_speed/nagents(model)
	if(model.n < no_steps)
		write(mean_a_file, "$average_area ")
		write(mean_speed_file, "$average_speed ")
	else 
		write(mean_a_file, "$average_area")
		write(mean_speed_file, "$average_speed")
	end

	last_hp_vert = open("Last_hp_vert.txt", "w")
	for i in 1:nagents(model)
		write(last_hp_vert, "Agent $i, position of $(new_pos[i]), considering position of $(last_half_planes[i][3])\n")
		write(last_hp_vert, "$(last_half_planes[i][1])\n")
		write(last_hp_vert, "$(last_half_planes[i][2])\n")
		write(last_hp_vert, "\n\n")
	end
	close(last_hp_vert)

	print("Finished step $(model.n)\n\n\n")
end


include("io_file.jl")


no_simulations::Int64 = 1
no_steps::Int64 = 10

function run_ABM()
	global compac_frac_file
        global mean_a_file
        global rot_o_file
        global rot_o_alt_file
	global mean_speed_file
for i in 1:no_simulations
	model = initialise(1000.0*sqrt(12), i)
	#figure, _ = abmplot(model)
        #save("./Simulation_Images/shannon_flock_n_=_$(0).png", figure)
	step!(model, agent_step!, model_step!, no_steps)
	write(compac_frac_file, "\n")
	write(mean_a_file, "\n")
	write(rot_o_file, "\n")
	write(rot_o_alt_file, "\n")
	write(mean_speed_file, "\n")
end

	do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)

end #This should be the end of the function or running the ABM

###This line simulates the model
run_ABM()
