clear
reset
set key off
set border 3

set yzeroaxis

###This stuff is no longer needed since it was purely for generating bins
#set boxwidth 0.0025 absolute
set style fill solid 1.0 noborder

#bin_width = 0.005
#bin_number(x) = floor(x/bin_width)
#rounded(x) = bin_width * ( bin_number(x) + 0.5 )

#stats 'mh_wall.txt' using 1 name "A"
#x = A_mean
#set arrow from x, graph 0 to x, graph 1 nohead lc rgb "red" 

#stats "../adf_file.txt" using 3 every ::500000::1000100 name "B"
#y = B_mean
#print y

identity(x) = x

plot "../adf_file.txt" using (identity($5)):(1) every ::500000::1000100 smooth frequency with boxes 
set title "Frequency histogram of number of spots available"
set ylabel "Frequency"
set xlabel "Number of spots available"
#set xrange [0.0:1.025]
