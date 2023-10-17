clear
reset
set key off
set border 3

set yzeroaxis

bin_width = 100.0
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )

set boxwidth 0.5*bin_width absolute
set style fill solid 1.0 noborder

x = 1000*sqrt(12)
set arrow from x, graph 0 to x, graph 1 nohead lc rgb "red"

no_lines = floor(system("wc -l ../periphery_file.txt")/1)

plot "../periphery_file.txt" using (rounded($1)):($2 ==  1? 1:0) every  ::floor(no_lines/2)::no_lines-1 smooth frequency with boxes, "../periphery_file.txt" using (rounded($1)):($2 ==  0? 1:0) every  ::floor(no_lines/2)::no_lines-1 smooth frequency with boxes
set title "Frequency histogram of DOD for those in the periphery"
set ylabel "Frequency"
set xlabel "DOD"
#set xrange [0.0:1.025]
