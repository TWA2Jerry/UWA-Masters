using Agents
using InteractiveDynamics
using CairoMakie
@agent SchellingAgent GridAgent{2} begin
    mood::Bool # whether the agent is happy in its position. (true = happy)
    group::Int # The group of the agent, determines mood as it interacts with neighbors
end

mutable struct SchellingAgent <: AbstractAgent
    id::Int             # The identifier number of the agent
    pos::NTuple{2, Int} # The x, y location of the agent on a 2D grid
    mood::Bool          # ...
    group::Int          # ...
end

model = AgentsIO.load_checkpoint("schelling.jld2"; scheduler = Schedulers.Randomly())
groupcolor(a) = a.group == 1 ? :blue : :red
groupmarker(a) = a.group == 1 ? :circle : :rect
figure, _ = abmplot(model; ac = groupcolor, am = groupmarker, as = 10)
figure # returning the figure displays it
display(figure)

