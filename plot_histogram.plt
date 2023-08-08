clear
reset
set key off
set border 3

set yzeroaxis

set boxwidth 0.05 absolute
set style fill solid 1.0 noborder

bin_width = 0.005
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )
plot "adf_file.txt" using (rounded($3)):(1) every ::250000::500000 smooth frequency with boxes 
set title "Frequency histogram of (un)happiness"
set ylabel "Frequency"
set xlabel "(un)happiness"

