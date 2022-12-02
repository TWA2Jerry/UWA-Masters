#Preliminaries
using Agents
using Random
using InteractiveDynamics
using CairoMakie

eps = 0.0000000001


#Create the agent 
mutable struct bird <: AbstractAgent
	id::Int
	pos::NTuple{2, Float64}
	vel::NTuple{2, Float64}
	#Create a class characteristic of Int64 type
	group::Int
end

#Initialise the model
function initialise(; seed=123, no_birds = 25)
	#Create the space
        space = ContinuousSpace((15.0, 15.0); periodic = true)
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
	origin_pos = (7.5, 7.5)
	#Generate agent 0. Use 
	agent = bird(0, (7.5, 7.5),(1.0, 0.0), 0)
	add_agent_pos!(agent, model)
	
	#For each q
	vix = agent.vel[1]
        viy = agent.vel[2]
	no_agents = 1
	for q in 0:(8-1)
		#If the q angle is of magnitude greater than so and so, then set the group to 2, otherwise set it to 1 
		direction_of_move = [cos(q*2*pi/8)*vix - sin(q*2*pi/8)*viy, sin(q*2*pi/8)*vix + cos(q*2*pi/8)*viy]
		print("Direction of move is $direction_of_move\n")
		angle_of_move = atan(direction_of_move[2], direction_of_move[1])
		theta_0 = atan(viy, vix)
		rel_angle = ((angle_of_move - theta_0 + pi)+2*pi)%(2*pi) - pi
		agent_class = 1

		#For each m 
		for m in 1:3
			if(q == 2 && m == 3)
				agent_class = 3
			elseif (q== 2 && m == 2)
					agent_class = 2
			end
			new_agent_pos = origin_pos .+ m .* direction_of_move .* 1.0 .* 1.0
			print("Position considered was $new_agent_pos\n")
			new_agent = bird(no_agents, Tuple(new_agent_pos), (1.0, 0.0), agent_class)
			add_agent!(new_agent, Tuple(new_agent_pos), model) 
			no_agents += 1
		end
	end		
	print("Number of agents registered is $no_agents\n")
	return model
end	

#Initialise
model = initialise()
print("The number of agents in the model is $(nagents(model))\n")
#Generate the image
groupcolor(a) = a.group == 0 ? :blue : (a.group == 1 ? :green : (a.group == 2 ? :orange : :black))
groupmarker(a) = :circle
figure, _ = abmplot(model; ac = groupcolor, am = groupmarker, as = 10)
figure
save("./base_ill.png", figure)
			
