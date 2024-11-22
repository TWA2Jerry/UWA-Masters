stats 'random.txt' using 1 name "A"
x = A_mean
print x
set arrow from x, graph 0 to x, graph 1 nohead
plot "random.txt" using 1:2 with lines

