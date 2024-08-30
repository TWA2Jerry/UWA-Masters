using LaTeXStrings
fig, ax = give_model(model, fig_box = ((200, 200), (600, 600)), marker_size = 40)
hidedecorations!(ax)

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
#sampled_cell = give_cell(model[1].pos, positions, model, model[1].vel, relic_half_plane, relic_passed = 1)
sampled_cell = give_cell(model[1].pos, positions, model, model[1].vel, relic_half_plane)
sampled_cell = give_cell_circled(sampled_cell, model[1].pos)
draw_cell!(sampled_cell)
Makie.arc!(model[1].pos, 100.0, -pi, pi; transparency = true, color = (:red, 0.2))


##Just illustrating the velocity of the current agent
Makie.arrows!([model[1].pos[1]], [model[1].pos[2] + 10.0], [model[1].vel[1]*20], [model[1].vel[2]*20])
text!(model[1].pos .+ model[1].vel .* 20, text = L"\hat{v}_1", align = (:left, :bottom), fontsize = 50)

colsize!(fig.layout, 1, Aspect(1, 1.0))
resize_to_layout!(fig)
