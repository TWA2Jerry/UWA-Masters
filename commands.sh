//Plotting
gnuplot
plot "rot_o_ave.txt" using 1:2 title "Rotational Order vs Time" with lines
set title "Rotational Order vs Time, 3 angles considered"
set xlabel "Time step"
set ylabel "Rotational Order"

