fig, ax = give_model(model, fig_box = ((200, 200), (600, 600)), marker_size = 30)
hidedecorations!(ax)

##Get the sampled point and plot it
sampled_direction = sampled_direction = rotate_vector(2*pi/8, model[1].vel)
sampled_distance = 50
sampled_point = model[1].pos .+ sampled_direction .* sampled_distance
Makie.scatter!(sampled_point, markersize = 30)

##Getting ready to plot cell for sampled point
relic_half_plane = generate_relic(model[1].pos, model[1].vel)
positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = Vector{Tuple{Tuple{Float64, Float64}, Int64}}(undef, 0)
all_agents_iterable = allagents(model)
for neighbour in all_agents_iterable
               if(neighbour.id == 1)
                   continue
               end
               pushfirst!(positions, (neighbour.pos, neighbour.id))
           end

##Get and plot the sampled cell
sampled_cell = give_cell_forward(sampled_point, positions, model, model[1].vel, relic_half_plane, relic_passed = 1)
#sampled_cell = give_cell(sampled_point, positions, model, model[1].vel, relic_half_plane)
sampled_cell = give_cell_circled(sampled_cell, sampled_point)
draw_cell!(sampled_cell)
Makie.arc!(sampled_point, 100.0, -pi, pi; transparency = true, color = (:red, 0.2))


##Just illustrating the velocity of the current agent
Makie.arrows!([model[1].pos[1]], [model[1].pos[2]], [model[1].vel[1]*20], [model[1].vel[2]*20])

