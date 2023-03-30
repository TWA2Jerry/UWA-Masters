bro = 1.0

function pants(bro)
	print("$bro\n")
	bro += 5.0
	print("$bro\n")
end

function call_pants()
	pants(bro)
	global bro
	bro = 0.0
end

call_pants()
call_pants()
