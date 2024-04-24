function find(a)
	while(a != rep[a])
		a = rep[a]
	end
	return a
end

function merge(a, b)
	rep_a = find(a)
	rep_b = find(b)

	if(rep_a = rep_b)
		return
	end
	
	if(size[rep_a] < size[rep_b])
		size[rep_b] = size[rep_b]+size[rep_a]
		rep[rep_a] = rep_b
	else	
		size[rep_a] = size[rep_a] + size[rep_b]
		rep[rep_b] = rep_a
	end

	return
end	


