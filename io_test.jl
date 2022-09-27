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
