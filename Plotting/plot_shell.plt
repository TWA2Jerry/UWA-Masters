clear
reset
set key off
set border 3

set yzeroaxis

bin_width = 1.0
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )

set boxwidth 0.5*bin_width absolute
set style fill solid 1.0 noborder

no_lines = floor(system("wc -l ../agent_vals.txt")/1)

plot "../agent_vals.txt" using (rounded($4)):(1) smooth frequency with boxes
set title "Frequency histogram of distances of chosen future positions for agents"
set ylabel "Frequency"
set xlabel "Distance of future position from current position"
#set xrange [0.0:1.025]
