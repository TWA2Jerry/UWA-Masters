function bro(a)
	print("yo\n")
	return a + 1
end

a = @time bro(1)
print("$a\n")
