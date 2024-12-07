include("../../Plotting/plot_model_template.jl")

fig, ax = plot_model(model; limits_arg = (model[101].pos[1]  - 12, model[101].pos[1] + 12, model[101].pos[2]  - 12, model[101].pos[2] + 12), markersize_arg = 40)
show_move!(model, 101, view_box = (model[101].pos .- (10, 10), model[101].pos .+ (10, 10)), m_arg = 2, m_spacing_arg = 5, marker_size = 30)
Makie.arrows!([model[101].pos[1]], [model[101].pos[2]], [model[101].vel[1]*3], [model[101].vel[2]*3], arrowsize = 40, linewidth = 10)
