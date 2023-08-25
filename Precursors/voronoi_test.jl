include("../agent_definition.jl")
using Random
using VoronoiCells
using GeometryBasics
using Plots
using InteractiveDynamics
using CairoMakie # choosing a plotting backend
using ColorSchemes
import ColorSchemes.balance
print("Packages loaded\n")

const no_birds::Int32 = 100
const rho::Float64 = 100.0
initialised::Int32 = 0
area_zero = zeros(Int32, 100)
const rect_bound::Float64 = 1000.0
const spawn_dim_x::Float64 = 100.0 #This gives the x dimesnion size of the initial spawning area for the agents
const spawn_dim_y::Float64 = 100.0 #This gives the y dimension size of the initial spawning area for the agents
rect = Rectangle(Point2(0,0), Point2(Int64(rect_bound), Int64(rect_bound)))
moves_areas::Vector{Tuple{Int64, Float64, Float64}} = [] #This is an array which will allow us to record all the areas and directions considered for each step, for each agent
no_move = ones(Int32, no_birds) #An array which will allow us to keep track of which agents never move
new_pos::Vector{Tuple{Float64, Float64}} = [(0.0, 0.0) for i in 1:no_birds] #An array that will store the new positions of the agents for movement when we go to the model step
convex_hull_point = zeros(Int32, 100)
last_half_planes::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}} = [(0.0, (0.0, 0.0), (0.0, 0.0), 0) for i in 1:no_birds]
const sigma = 0.0
const tracked_agent::Int64 = rand(1:no_birds)
tracked_path::Vector{Tuple{Float64, Float64}} = []

initial_positions::Vector{Tuple{Float64, Float64}} = []
        initial_vels::Vector{Tuple{Float64, Float64}} = []
        temp_hp::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}}= []
        pack_positions = Vector{Point2{Float64}}(undef, no_birds)


for i in 1:no_birds
                R = 100.0
                #rand_position = Tuple(100*rand(Float64, 2)) .+ (50.0, 50.0)
                angle_per_bird = 2*pi/(no_birds-1)
                initial_pos = i != tracked_agent ? (R*cos(angle_per_bird*(i-1)), R*sin(angle_per_bird*(i-1))) .+ (300.0, 300.0) : (300.0, 300.0)
                rand_vel = (-R*sin(angle_per_bird*(i-1)), R*cos(angle_per_bird*(i-1)))
                #rand_vel = rand_vel ./norm(rand_vel)

                #rand_vel::Tuple{Float64, Float64} = 2 .* Tuple(rand(Float64, 2)) .- (1.0, 1.0)
                #rand_vel = rand_vel ./norm(rand_vel)
                push!(initial_positions, initial_pos)
                push!(initial_vels, rand_vel)
                pack_positions[i] = initial_positions[i]
                #push!(moves_areas, [])
                #push!(last_half_planes, [])
                #=if(model.simulation_number==1)
                        push!(new_pos, (0.0, 0.0))
                end=#
                if(i == tracked_agent)
                        push!(tracked_path, initial_positions[i])
                end
        end
for i in 1:no_birds
                print("the position of agent $i is $(initial_positions[i])\n")
                print("The position of the i-th thang in pack positions is $(pack_positions[i])\n")
        end

        print("Yo\n")
        #Calculate the DOD based off the initial positions
        init_tess = voronoicells(pack_positions, rect)
        init_tess_areas = voronoiarea(init_tess)
        print("Finished the package thang\n")
