out_file = open("io_test.txt", "w")
a = 5
b = [1,2]
c = (1,2)
print(a ,"\n")
write(out_file, "This is my gun, my only gun", "\n")
write(out_file, "$a \n")
write(out_file, "$b \n")
write(out_file, "$c \n")
write(out_file, "$(c[1]) $(c[2])\n")
write(out_file, "I have no idea\n")
using Random
out_file_3 = open("io_test_2.txt", "w")
for i in 1:20
	write(out_file_3, "$(rand(Float64)) ")
end

for i in 1:3
	bro = open("rewrite_test.txt", "w")
	write(bro, "$i\n")
	close(bro)
end

out_file_4 = open("io_test_3.txt", "a")
write(out_file_4, "Uh, hello\n")
#=
out_file_5 = open("io_test_4.txt", "w+")
write(out_file_5, "Yo\n")
#close(out_file_5)
#thing = open("io_test_4.txt", "r")
#bro = readlines(thing) 
bro = readlines(out_file_5)
#bro = readlines("io_test_4.txt")
print(bro[1])
=#

outfile6 = open("io_test_5.txt", "w")
blah  = 6
include("dude_file.jl")

for i in 1:10	
	global outfile6 = open("io_test_5.txt", "w")
	local_outfile = outfile6
	truncate(outfile6, 0)
	write_dude()
	close(outfile6)
	outfile6 = open("io_test_5.txt", "r")
	bro = readlines(outfile6)
	print(bro[1])
	print(bro[2], "\n")	
end 

