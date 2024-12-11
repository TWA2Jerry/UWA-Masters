using LaTeXStrings
model = AgentsIO.load_checkpoint("./Records/stdod68/stdod68_s50")
fig, ax = give_model(model, fig_box = ((180, 180), (620, 620)), marker_size = 40)
#fig, ax = give_model(model, fig_box = (model[2].pos .- (220, 220), model[2].pos .+ (220, 220)), marker_size = 40)
hidedecorations!(ax)

##Getting ready to plot cell for sampled point
positions::Vector{Tuple{Tuple{Float64, Float64}, Int64}} = Vector{Tuple{Tuple{Float64, Float64}, Int64}}(undef, 0)
all_agents_iterable = allagents(model)
for neighbour in all_agents_iterable
            if(neighbour.id == 1)   
			#if(neighbour.id == 2)
                   continue
               end
               pushfirst!(positions, (neighbour.pos, neighbour.id))
           end

for i in 1:nagents(model)
	text!(model[i].pos .- (0, 10), text = "$i", align = (:center, :top), fontsize = 50)
end

#=
relic_half_plane = generate_relic(model[2].pos, model[2].vel)
perp_vector = rotate_vector(-pi/2, model[2].vel)
perp_angle = atan(perp_vector[2], perp_vector[1])
Makie.arc!(model[2].pos, 200.0, perp_angle, perp_angle + pi; transparency = true, color = (:red, 0.4))
Makie.lines!([model[2].pos .- 200 .* perp_vector, model[2].pos .+ 200 .* perp_vector]; transparency = true, color = (:red, 0.4)) 
Makie.arrows!([model[2].pos[1]], [model[2].pos[2]], [model[2].vel[1]*30], [model[2].vel[2]*30], arrowsize = 40, linewidth = 10)
text!(model[2].pos .+ model[2].vel .* 30 .+ (10, 0), text = L"\hat{v}_2", align = (:left, :bottom), fontsize = 50)
sampled_cell = give_cell_forward(model[2].pos, positions, model, model[2].vel, relic_half_plane, relic_passed = 1, rhop = 200.0)
sampled_cell = give_cell_circled(sampled_cell, model[2].pos, rhop = 200.0)
=#

relic_half_plane = generate_relic(model[1].pos, model[1].vel)
Makie.arc!(model[1].pos, 200.0, 0, pi; transparency = true, color = (:red, 0.4))
Makie.lines!([model[1].pos .- (200, 0), model[1].pos .+ (200.0, 0)]; transparency = true, color = (:red, 0.4)) 
##Just illustrating the velocity of the current agent
Makie.arrows!([model[1].pos[1]], [model[1].pos[2] + 10.0], [model[1].vel[1]*20], [model[1].vel[2]*20], arrowsize = 40, linewidth = 10)
text!(model[1].pos .+ model[1].vel .* 30 .+ (10, 0), text = L"\hat{v}_1", align = (:left, :bottom), fontsize = 50)
sampled_cell = give_cell_forward(model[1].pos, positions, model, model[1].vel, relic_half_plane, relic_passed = 1, rhop = 200.0)
#sampled_cell = give_cell(model[1].pos, positions, model, model[1].vel, relic_half_plane, rhop = 1000.0)
#sampled_cell = give_cell_circled(sampled_cell, model[1].pos)


draw_cell!(sampled_cell)

##Get and plot the sampled cell
colsize!(fig.layout, 1, Aspect(1, 1.0))
resize_to_layout!(fig)

#=
Colorbar(fig[1,2], limits = (0, 1), colormap = :cool)
show_move!(model, 1, view_box = ((200, 200), (600, 600)), m_arg = 2, m_spacing_arg = 100, marker_size = 30, qp_arg = 2)
=#
