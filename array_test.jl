a = zeros(Int64, 100)

function something(i)
	print("Something called\n")
	if(a[i] == 1)
		print("a[i] is currently occupied\n")
	end
end

function main()
	something(1)
	a[1] = 1
	something(1)
end

main()


