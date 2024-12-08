include("../../Plotting/plot_model_template.jl")

model = AgentsIO.load_checkpoint("./Records/stdod67/stdod67_s50")

fig, ax = plot_model(model; limits_arg = (model[101].pos[1]  - 12, model[101].pos[1] + 12, model[101].pos[2]  - 12, model[101].pos[2] + 12), markersize_arg = 40, resize_to_layout_arg = 0)
show_move!(model, 101, view_box = (model[101].pos .- (10, 10), model[101].pos .+ (10, 10)), m_arg = 2, m_spacing_arg = 5, marker_size = 30, draw_best_cell_arg = 0, conflict_dist_arg = 1.0)
Makie.arrows!([model[101].pos[1]], [model[101].pos[2]], [model[101].vel[1]*3], [model[101].vel[2]*3], arrowsize = 40, linewidth = 10)
text!(model[101].pos .+ model[101].vel .* 3, text = L"\hat{v}_1", align = (:left, :top), fontsize = 50)
Colorbar(fig[1,2], limits = (0, 1), colormap = :cool, ticklabelsize = 40, ticks = ([0.0, 0.5, 1.0], [L"0.0", L"A_t", L"A_{\mathrm{max}}"]))

#conf_pos = model[101].pos .+ rotate_vector(2*pi/8, model[101].vel) .* 10 ##This is for illustrating the point that is meant to be conflicted. 
Makie.scatter!(ax, model[102].pos, marker=:circle, markersize = 40)

text!(model[101].pos .+ rotate_vector(2*pi/8, model[101].vel) .* 5 .+ (0.5, 0), text = L"A", align = (:left, :top), fontsize = 30)
text!(model[101].pos .+ model[101].vel .* 5 .+ (0.5, 0), text = L"B", align = (:left, :top), fontsize = 30)
text!(model[101].pos .+ rotate_vector(-2*pi/8, model[101].vel) .* 5 .+ (0.5, 0), text = L"C", align = (:left, :top), fontsize = 30)

hidedecorations!(ax)
