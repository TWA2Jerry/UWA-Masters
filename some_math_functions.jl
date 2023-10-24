function norm(v::Tuple{Float64, Float64})
        sum_of_squares::Float64 = 0.0
        for i in 1:length(v)
                sum_of_squares += (v[i])^2
        end

        return sqrt(sum_of_squares)
end

function norm(v::Float64)
        sum_of_squares::Float64 = 0.0
        for i in 1:length(v)
                sum_of_squares += (v[i])^2
        end

        return sqrt(sum_of_squares)
end

function norm(v::Vector{Float64})
        sum_of_squares::Float64 = 0.0
        for i in 1:length(v)
                sum_of_squares += (v[i])^2
        end

        return sqrt(sum_of_squares)
end

###Function that takes a vector and calculates the mean of the elements in the vector
function mean(v)
        total = 0.0
        for i in 1:length(v)
                total += v[i]
        end

        return total/length(v)
end

