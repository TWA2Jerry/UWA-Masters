clear
reset
set key off
set border 3

set yzeroaxis

set boxwidth 10.0 absolute
set style fill solid 1.0 noborder

bin_width = 0.1
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )

#stats 'mh_wall.txt' using 1 name "A"
#x = A_mean
#set arrow from x, graph 0 to x, graph 1 nohead lc rgb "red" 

#stats "../adf_file.txt" using 3 every ::250000::500100 name "B"
#y = B_mean
#print y

plot "../adf_file.txt" using (rounded($4)):(1) every ::250000::500100 smooth frequency with boxes 
set title "Frequency histogram of measured area to tdod"
set ylabel "Frequency"
set xlabel "Ratio of the actual measured area to the tdod"
set xrange [0.0:6000.0]
