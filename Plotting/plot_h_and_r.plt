set y2tics 0, 10 #This sets the distance between tics, with 0 representing the smallest value, and 10 representing the spaces 
set ytics nomirror
plot "circ_radial.txt" using 1:2 with lines axis x1y2, "ensemble_rot_o_alt.txt" using 1:2 with lines axis x1y1
set title "Mean radial distance and rotational order vs time"

