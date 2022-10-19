using Agents
using Random

new_position = [(0.0, 0.0), (0.0, 0.0)]

mutable struct agent <: AbstractAgent
       id::Int
       pos::NTuple{2, Float64}
end

function initialise(; seed=123, no_agents=2)
       space=ContinuousSpace((100.0, 100.0); periodic=true)
       rng = Random.MersenneTwister(seed)
       properties = Dict(:t => 0.0)
       model = ABM(agent, space; properties, rng, scheduler = Schedulers.fastest)
       bro = agent(1, (52.0, 0.0))
       bro2 = agent(2, (48.0, 0.0))
       add_agent!(bro, (52.0, 0.0), model)
       add_agent!(bro2, (48.0, 0.0), model)
       return model
end

function agent_step!(agent, model)
       all_agents_iterable = allagents(model)
       vector_to_neighbour = 0.0
       for neighbour in all_agents_iterable
       	if(neighbour.id==agent.id)
       		continue
       	end
	print("Neighbour $(neighbour.id)  position was $(neighbour.pos)\n")
       vector_to_neighbour = neighbour.pos .- agent.pos
       end
       #move_agent!(agent, agent.pos .+ 0.25 .* vector_to_neighbour, model)
       new_position[agent.id] = agent.pos .+ 0.25 .* vector_to_neighbour
end

function model_step!(model)
       print("Model step now\n")
       all_agents_iterable = allagents(model)
       for agent in all_agents_iterable
       move_agent!(agent, new_position[agent.id], model)
       print("Agent $(agent.id) has a new position of $(agent.pos)\n")
       end
end

model = initialise()

step!(model, agent_step!, model_step!, 3)

