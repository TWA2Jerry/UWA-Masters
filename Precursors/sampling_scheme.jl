using CairoMakie
q = 8
m = 3
points::Vector{Tuple{Float64, Float64}} = Vector{Tuple{Float64, Float64}}(undef, 0)
colours = []
push!(points, (0.0, 0.0))
push!(colours, :black)
vix = 1.0
viy = 0.0

for i::Int64 in 0:(q-1)
        direction_of_move::Tuple{Float64, Float64} = (cos(i*2*pi/q)*vix - sin(i*2*pi/q)*viy, sin(i*2*pi/q)*vix + cos(i*2*pi/q)*viy)
        for j::Int64 in 1:m
                push!(points, j.* direction_of_move)
                if(i > 10 && i < 7)
                        push!(colours, :orange)
                else
                        push!(colours, :green)
                end
        end
end

fig, ax, _ = Makie.scatter(points, axis = (; aspect = 1), color=colours, marker=:circle, markersize = 20)
Makie.text!([(0.0,0.0)], text = L"i", fontsize = 40)
Makie.arrows!([0.0], [0.0], [0.85], [0.0], linewidth = 3)
Makie.text!([(0.9,-0.1)], text = L"\vec{v}_i", fontsize = 40, align = (:right, :top)
