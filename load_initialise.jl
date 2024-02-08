include("half_plane_fast.jl")
include("half_plane_bounded.jl")
include("convex_hull.jl")
include("rot_ord.jl")
include("rot_ord_check.jl")

function load_initialise(pos_vels_file, step; target_area_arg = 1000*sqrt(12), simulation_number_arg = 1, no_bird = 100, seed = 123, tracked_agent_arg = 42, no_moves_arg = 100)
	#Create the space
	space = ContinuousSpace((rect_bound, rect_bound); periodic = true)
	#Create the properties of the model
	properties = Dict(:t => 0.0, :dt => 1.0, :n => step, :CHA => 0.0, :target_area => target_area_arg, :simulation_number => simulation_number_arg, :tracked_agent => tracked_agent_arg, :no_moves => no_moves_arg)
	
	#Create the rng
	rng = Random.MersenneTwister(seed)
	
	print("Before model\n")

	#Create the model
	model = UnremovableABM(
		bird, space; 
		properties, rng, scheduler = Schedulers.fastest
	)	


	#Generate random initial positions for each bird, then calculate the DoDs
	initial_positions::Vector{Tuple{Float64, Float64}} = []
	initial_vels::Vector{Tuple{Float64, Float64}} = []
	temp_hp::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}}= []
	pack_positions = Vector{Point2{Float64}}(undef, no_birds)
	
	

	#Initialise the positions based on the spawn-error free function of assign_positions
	assign_pos_vels(pos_vels_file, initial_positions, initial_vels, step, no_birds) 
	for i in 1:no_bird

		pack_positions[i] = initial_positions[i]
		print("Pack positions i is $(pack_positions[i])\n")
		if(i == tracked_agent)
			push!(tracked_path, initial_positions[i])
		end
	end 

	#Calculate the DOD based off the initial positions
	init_tess = voronoicells(pack_positions, rect)
	init_tess_areas = voronoiarea(init_tess)

	#Calculate the DoDs based off the initial positions
	#initial_dods = voronoi_area(model, initial_positions, rho)
	initial_dods::Vector{Float64} = []
	true_initial_dods::Vector{Float64} = []
	for i::Int32 in 1:no_birds
		print("\n\nCalculatin initial DOD for agent $i, at position $(initial_positions[i]).")
		ri::Tuple{Float64, Float64}  = Tuple(initial_positions[i])
		neighbouring_positions = Vector{Tuple{Tuple{Float64, Float64}, Int64}}(undef, 0)
		for j::Int32 in 1:no_birds
			if(i == j)
				continue 
			end
			push!(neighbouring_positions, (Tuple(initial_positions[j]), j))
		end
		vix::Float64 = initial_vels[i][1]
		viy::Float64 = initial_vels[i][2]
		relic_x::Float64 = -1.0*(-viy)
        	relic_y::Float64 = -vix
        	relic_pq::Tuple{Float64, Float64} = (relic_x, relic_y)
        	relic_angle::Float64 = atan(relic_y, relic_x)
        	relic_is_box::Int64 = -1
        	relic_half_plane::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64} = (relic_angle, relic_pq, ri, relic_is_box)

		initial_cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = @time voronoi_cell_bounded(model, ri, neighbouring_positions, rho, eps, inf, temp_hp, initial_vels[i], [relic_half_plane])
		initial_A::Float64 = voronoi_area(model, ri, initial_cell, rho) 
	
		true_initial_cell::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = @time voronoi_cell(model, ri, neighbouring_positions, rho,eps, inf, temp_hp, initial_vels[i])
                true_initial_A::Float64 = voronoi_area(model, ri, true_initial_cell, rho)

		#=print("The half planes that generated the cell for agent $i were \n")
                        for i in 1:length(temp_hp)
                                print("$(temp_hp[i])\n")
                        end
		=#
		#What is this? last_half_planes[i], (initial_cell, temp_hp, ri)
			
		print("Initial DOD calculated to be $initial_A\n")
		if(abs(initial_A) > pi*rho^2/2)
			print("Main file here. Conventional area exceeded by agent $(i)in position $(initial_positions[i])\n")	
			print("The cell was \n")
			for i in 1:length(initial_cell)
				print("$(initial_cell[i])\n")
			end
			print("The half planes that generated it were \n")
			for i in 1:length(temp_hp)
				print("$(temp_hp[i])\n")
			end
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
	total_area::Float64 = 0.0
	total_speed::Float64 = 0.0
	for i::Int32 in 1:no_birds
		agent = bird(i, initial_positions[i], initial_vels[i], 1.0, initial_dods[i], true_initial_dods[i], target_area_arg, 0, 0.0, 0.0, 0.0, 1, 0.0)
		agent.vel = agent.vel ./ norm(agent.vel)
		#print("Initial velocity of $(agent.vel) \n")
		add_agent!(agent, initial_positions[i], model)
		total_area += true_initial_dods[i]/(pi*rho^2)
		total_speed += agent.speed
	end	
	
	return model
end 
