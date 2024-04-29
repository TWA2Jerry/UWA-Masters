include("some_math_functions.jl")
include("rot_ord.jl")

function enforce_rotation(model, desired_rotation)
	n::Int64 = nagents(model)
	com = center_of_mass(model)
	for i in 1:n
		rot_o_raw = rot_o_generic(model[i].pos .- com, model[i].vel)
		if(rot_o_raw * desired_rotation < 0)
			model[i].vel  = rotate_vector((Float64(pi)), model[i].vel)
		end
	end	
	
	return 
end	
