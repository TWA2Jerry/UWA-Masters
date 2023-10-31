include("half_plane_fast.jl")
include("half_plane_bounded.jl")
include("move_gradient_file.jl")
include("draw_circle_part.jl")
include("global_vars.jl")

function give_agent_cell(agent_i::bird, model::UnremovableABM{ContinuousSpace{2, true, Float64, typeof(Agents.no_vel_update)}, bird, typeof(Agents.Schedulers.fastest), Dict{Symbol, Real}, MersenneTwister})
        all_agents_iterable =  allagents(model)
        temp_hp::Vector{Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}} = []
        previous_areas::Vector{Float64} = zeros(nagents(model))
        actual_areas::Vector{Float64} = zeros(nagents(model))

                neighbour_positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = []
                for agent_j in all_agents_iterable
                        if(agent_i.id == agent_j.id)
                                continue
                        end
                        push!(neighbour_positions, (agent_j.pos, agent_j.id))
                end
                ri::Tuple{Float64, Float64} = agent_i.pos
                vix::Float64 = agent_i.vel[1]
                viy::Float64 = agent_i.vel[2]
                relic_x::Float64 = -1.0*(-viy)
                relic_y::Float64 = -vix
                relic_pq::Tuple{Float64, Float64} = (relic_x, relic_y)
                relic_angle::Float64 = atan(relic_y, relic_x)
                relic_is_box::Int64 = 2
                relic_half_plane::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64} = (relic_angle, relic_pq, agent_i.pos, relic_is_box)

                new_cell_i::Vector{Tuple{Tuple{Float64, Float64}, Int64, Int64}} = voronoi_cell(model, ri, neighbour_positions, rho, eps, inf, temp_hp, agent_i.vel, relic_half_plane)

                ##Procedure for adding the
                cell_including_circle = []
                #print("Starting\n")
                for i in 1:length(new_cell_i)
                        point = new_cell_i[i]
                        point_pp = new_cell_i[(i)%length(new_cell_i)+1]
                        push!(cell_including_circle, point)
                        #print("$(point[3]) $(point_pp[2])\n")
                        if(point[3] == 0 && point_pp[2] == 0)
                                #print("State helper.jl here. Circle confirmed\n")
                                vec_to_point = point[1] .- agent_i.pos
                                vec_to_pointpp = point_pp[1] .- agent_i.pos
                                theta_1 = atan(vec_to_point[2], vec_to_point[1])
                                theta_2 = atan(vec_to_pointpp[2], vec_to_pointpp[1])
                                if(theta_2 < theta_1)
                                        theta_2 += 2*pi
                                end
                                circle_points = circle_seg(agent_i.pos, rho, theta_1, theta_2)
                                for j in 1:length(circle_points[1])
                                        push!(cell_including_circle, ((circle_points[1][j], circle_points[2][j]), 0, 0))
                                end
                        end
                end
                return cell_including_circle
end

