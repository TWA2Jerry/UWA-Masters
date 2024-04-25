#plot "rot_order_alt_tave.txt" using 1:2 title "Rotational order vs target DOD" ps 6 lw 6 with lines
plot "rot_order_alt_tave.txt" using 1:2 ps 6 lw 6 with lines notitle
set logscale x 3
#set title "Rotational order vs Target DOD" font ", 40"
set bmargin 8
set lmargin 15
#set lmargin at screen 0
set rmargin 8
set ylabel "Rotational Order" offset -5, 0 font ", 30"
set xlabel "Target DOD" offset 0,-2 font ", 30"
set ytics font ", 25"
set xtics font ", 22"

