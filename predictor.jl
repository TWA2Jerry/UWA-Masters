function path_collide(ri, rj, veli, velj)
	rix = ri[1]
	riy = ri[2]
	rjx = rj[1]
	rjy = rj[2]
	vix = veli[1]
	viy = veli[2]
	vjx = velj[1]
	vjy = velj[2]
	a = (vix-vjx)^2 + (viy-vjy)^2
	b = 2*(vix - vjx)*(rix-rjx) + 2*(viy-vjy)*(riy-rjy)
	c = (rix-rjx)^2+(riy-rjy)^2-4
	val0 = c
	#print("c was calculated as $c\n")
	val1 = a+b+c
	#print("val1 was calculated as $val1\n")
	tlocal = -b/(2*a)
	#print("tlocal was calculated as $tlocal\n")
	vallocal = a*tlocal^2+b*tlocal + c
	if(tlocal > 1.0 || tlocal < 0.0)
		vallocal = 1.0
	end
	#print("vallocal was calculated as $vallocal\n")
	if(val0 < 0 || val1 < 0 || vallocal < 0)
		return 1
	else
		return 0
	end
end

function collide_predicted(ri,rj,  vi, vj, q, qprime)
	#Define the three possible velocities that agent j can take
	vjx = vj[1]
	vjy = vj[2]
	
	collision_flag = 0
	for i in -div(qprime,2):div(qprime/2)
		velj = [cos(i*2*pi/q)*vjx - sin(i*2*pi/q)*vjy, sin(i*2*pi/q)*vjx + cos(i*2*pi/q)*vjy]
		if(path_collide(ri, rj, vi, vj)==1)
			collision_flag = 1
			break
		end
	end
	#For each velocity, check if the chosen vi will result in a collision with vj
	#If there is a collision with any of the three possible vj, then return collision detected
end
