using Agents
using Random
using VoronoiCells
using GeometryBasics
using Plots
using InteractiveDynamics
using CairoMakie # choosing a plotting backend

mutable struct bird <: AbstractAgent
        id::Int
        pos::NTuple{2, Float64}
        vel::NTuple{2, Float64}
end


