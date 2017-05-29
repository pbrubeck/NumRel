set xlabel "x"
m="./data/data1d.dat"
set terminal x11 0
set nokey
set grid
set title 'u(x)'
plot m using 1:2 with linespoints