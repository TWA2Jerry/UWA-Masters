using Agents

###Create the agent
mutable struct bird <: AbstractAgent
        id::Int
        pos::NTuple{2, Float64}
        vel::NTuple{2, Float64}
        speed::Float64
        A::Float64 #The area of the agent's DOD, at least in their own eyes
        true_A::Float64 #The true area of the agent's DOD
	tdod::Float64
	no_neighbours::Int32
	perimeter_squared::Float64
	rot_o_alt::Float64
	rot_o_alt_corr::Float64
	collaborator::Int32
	rl::Float64
	distance::Float64
	best_A::Float64
	rot_dir::Int32 
end

