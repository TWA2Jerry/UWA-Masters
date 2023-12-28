function polygon_area(vertices::Vector{Tuple{Float64, Float64}})
        Area::Float64 = 0.0
        num_points::Int32 = length(vertices)
        for i in 1:length(vertices)
                #Use the shoestring formula to calcualte the area
                j = (i)%num_points+1
                #print("i and j are $i $j\n")
                xi = vertices[i][1]
                yi = vertices[i][2]
                xj = vertices[j][1]
                yj = vertices[j][2]
                Area += 0.5 * (yi + yj)* (xi - xj)
        end

        return Area
end
