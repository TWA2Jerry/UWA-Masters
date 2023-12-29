include("some_math_functions.jl")
function generate_artificial_hp(position::Tuple{Float64, Float64}, angle::Float64)
	vel::Tuple{Float64, Float64} = (1.0, 0.0)
	vix::Float64 = vel[1]
	viy::Float64 = vel[2]
        relic_x::Float64 = cos(angle)*vix-sin(angle)*(viy)
        relic_y::Float64 = sin(angle)*vix+cos(angle)*viy
        relic_pq::Tuple{Float64, Float64} = (relic_x, relic_y)
        relic_angle::Float64 = atan(relic_y, relic_x)
        relic_is_box::Int64 = -2
        relic_half_plane::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64} = (relic_angle, relic_pq, position, relic_is_box)
end

###This function generates the half plane for our extended vision. Given the velocity of the current agent, and the angle you want to rotate by, return the appropriate half plane
function generate_artificial_hp_vision(position::Tuple{Float64, Float64}, vel::Tuple{Float64, Float64}, rotate_angle::Float64)  
	rotated_vector::Tuple{Float64, Float64} = rotate_vector(rotate_angle, vel)
	half_plane = generate_artificial_hp(position, atan(rotated_vector[2], rotated_vector[1]))
	return half_plane
end

function generate_relic(position, vel)
	vix::Float64 = vel[1]
        viy::Float64 = vel[2]
	relic_x::Float64 = -1.0*(-viy)
        relic_y::Float64 = -vix
        relic_pq::Tuple{Float64, Float64} = (relic_x, relic_y)
        relic_angle::Float64 = atan(relic_y, relic_x)
        relic_is_box::Int64 = -2
        relic_half_plane::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64} = (relic_angle, relic_pq, position, relic_is_box)
	return relic_half_plane
end

###This function just generates a half plane whose defining points and pq is on the agent position and is the agent vel respectively. Used for calculating the hemispheres in weird_vision. 
function generate_relic_alt(position, vel)
	relic_x::Float64 = vel[1]
        relic_y::Float64 = vel[2]
	relic_pq::Tuple{Float64, Float64} = (relic_x, relic_y)
        relic_angle::Float64 = atan(relic_y, relic_x)
        relic_is_box::Int64 = -2
        relic_half_plane::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64} = (relic_angle, relic_pq, position, relic_is_box)
        return relic_half_plane
end
