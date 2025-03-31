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

###Function that takes in some position, some direction vel, and some angle with which to rotate that direction by, and then returns a relic half plane with a defining point of position, and defining angle of vel rotated by angle. 
function generate_relic_alt(position, vel, angle = 0.0; relic_id = -2)
	rotated_vel = rotate_vector(angle, vel)
		relic_x::Float64 = rotated_vel[1]
        relic_y::Float64 = rotated_vel[2]
	relic_pq::Tuple{Float64, Float64} = (relic_x, relic_y)
        relic_angle::Float64 = atan(relic_y, relic_x)
        relic_is_box::Int64 = -2
        relic_half_plane::Tuple{Float64, Tuple{Float64, Float64}, Tuple{Float64, Float64}, Int64} = (relic_angle, relic_pq, position, relic_id)
        return relic_half_plane
end
