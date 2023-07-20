Hello future me, or probably presently-confused me. 

This document is just to outline what files are needed to run the project. To sum up, our project is to simulate a group of agents flocking together. 

The main program is in program.jl. This file loads all the necessary files for running the simulations. It also defines the files for input and output.

The file "main.jl" contains all the main functions for the model and agents as well as global variables. 

The type of your model is ::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}, or  Agents.SingleContainerABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, Dict{Int64, bird}, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister}. The former comes from us initialising model as an unremovableABM (model which does not permit removal of agents, but does allow addition as far as I can see (see here: https://juliadynamics.github.io/Agents.jl/stable/examples/schelling/#Using-an-UnremovableABM-1). 

The type of your agent is: 
