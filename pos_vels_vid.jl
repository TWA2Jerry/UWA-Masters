include("load_randstep.jl")
actual_areas::Vector{Float64} = Vector{Float64}(undef, 0)
previous_areas::Vector{Float64} = Vector{Float64}(undef, 0)
new_pos::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
tracked_path::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)

for i in 1:5000
       model = model = load_initialise(pos_vels_file, i, target_area_arg = 0.0)
       draw_figures(model, actual_areas, previous_areas, 0.5, new_pos,tracked_path)
       end
