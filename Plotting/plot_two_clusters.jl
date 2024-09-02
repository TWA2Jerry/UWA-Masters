###Drawing the tessellations
fig, ax = give_model_cell_circled(model)
hidedecorations!(ax)
colsize!(fig.layout, 1, Aspect(1, 1.0))
resize_to_layout!(fig)

fig2, ax2 = give_model_cell_forward_circled(model)
hidedecorations!(ax2)
colsize!(fig2.layout, 1, Aspect(1, 1.0))
resize_to_layout!(fig2)


###Drawing the adj graphs
positions = Vector{Tuple{Float64, Float64}}(undef, 0)
adj = Array{Vector}(undef, nagents(model))
adj_forward = Array{Vector}(undef, nagents(model))

for i in 1:nagents(model)
	push!(positions, model[i].pos)
end

construct_voronoi_adj_list(model, adj)
construct_forward_voronoi_adj_list(model, adj_forward)

fig3, ax3 = draw_graph(positions, adj)
hidedecorations!(ax3)
colsize!(fig3.layout, 1, Aspect(1, 1.0))
resize_to_layout!(fig3)

fig4, ax4 = draw_graph(positions, adj_forward)
hidedecorations!(ax4)
colsize!(fig4.layout, 1, Aspect(1, 1.0))
resize_to_layout!(fig4)
