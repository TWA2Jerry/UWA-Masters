using CairoMakie

include("polygon_area_file.jl")
include("order_parameters.jl")

#n::Int32 = 3
r::Float64 = 1.0

#vertexes::Vector{Tuple{Float64, Float64}} = Array{Tuple{Float64, Float64}}(undef, 0)


function return_regularities()
	regularities::Vector{Float64} = Array{Float64}(undef, 0)
	
	for n in 3:20
        	vertexes::Vector{Tuple{Float64, Float64}} = Array{Tuple{Float64, Float64}}(undef, 0)
        	for i in 0:n
                	coord::Tuple{Float64, Float64} = (r*cos(2*pi*i/n), r*sin(2*pi*i/n))
                	push!(vertexes, coord)
        	end

		push!(regularities, regularity_metric(vertexes))
	end
	return regularities
end




###Some earlier pretesting
#=
for n in 3:20
        vertexes::Vector{Tuple{Float64, Float64}} = Array{Tuple{Float64, Float64}}(undef, 0)
        for i in 0:n
                coord::Tuple{Float64, Float64} = (r*cos(2*pi*i/n), r*sin(2*pi*i/n))
                push!(vertexes, coord)
        end

        print("$(regularity_metric(vertexes))\n")

end 
=#

#=
for i in 0:n
        coord::Tuple{Float64, Float64} = (r*cos(2*pi*i/n), r*sin(2*pi*i/n))
        push!(vertexes, coord)
end
=#


#figure = Makie.lines(vertexes)
