const eps::Float64 = 0.0000000001
const inf::Float64 = 1000000000000.0
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

###Function for calculating the intersection between two points
function inter(h1::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}, h2::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64}, eps::Float64, inf::Float64)
        #h1 and h2 represent the half planes we want to calculate the line intersections for
        #print("Calculating the intersection for $(h1[2]) and $(h2[2])\n")
        #m1 = h1[2][2]/h1[2][1]
        #m2 = h2[2][2]/h2[2][1]
        m1::Float64 = 0.0
        m2::Float64 = 0.0
        if(abs(abs(h1[1]) - pi/2) < 0.000001)
                #print("Infinite gradient detected for m1\n")
                m1 = inf
        else
                m1 = h1[2][2]/h1[2][1]
        end
        #print("m1 found to be $m1\n")
        if(abs(abs(h2[1]) - pi/2) < 0.000001)
                m2 = inf
                #print("Infinite gradient detected for m2\n")
        else
                m2 = h2[2][2]/h2[2][1]
        end
        if(abs(m1 - m2) < abs(eps))
                print("Parallel planes yo\n")
                #exit()
                return ((-1.0, -1.0), -1)
        end
        ##print("m1 - m2 found to be $(m1-m2)\n")
        c1::Float64 = h1[3][2] - m1*h1[3][1]
        c2::Float64 = h2[3][2] - m2*h2[3][1]
        xint::Float64 = (c2-c1)/(m1-m2)
        yint::Float64 = -1.0
        if(abs(m1 - inf) < abs(eps))
                yint = m2 * xint + c2
        else
                yint = m1 * xint + c1
        end
        #print("Intersect calculated as $([xint, yint])\n")
        return ((xint, yint), 1)
end

