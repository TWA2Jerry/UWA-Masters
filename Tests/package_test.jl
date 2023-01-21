using Agents

include("aux.jl")

###Create the agent
mutable struct bird <: AbstractAgent
        id::Int
        pos::NTuple{2, Float64}
        vel::NTuple{2, Float64}
end

using Random #for reproducibility
function initialise(; seed = 123, no_birds = 10)
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
	
	agent = bird(1, Tuple(100*rand(Float64, 2)) .+ (50.0, 50.0), Tuple(rand(Float64, 2)))
	add_agent!(agent, model) 
	return model
end

model = initialise()
print(bro(model))
