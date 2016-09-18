# ===============================================
# plot command for gnuplot
# Exercise 1 Sheet 07 CP 1
# ===============================================

# ===============================================
# FOR TESTING

# set terminal pngcairo dashed
# set output 'test.png'
# test
# set output

# Define Linestyles

# lc: line color
# lt: line type
# dt: dash type
# lw: line width
# pt: point type
# ps: point size
# pi: point interval
# ===============================================

set style func linespoints
set termoption dashed
set style line 1 lc rgb "red" lt 1 lw 2 dt 1 pt 0 pi 2
set style line 2 lc rgb "blue" lt 1 lw 2 dt 1 pt 0 pi 2

# plot commands

set title 'path of electron around proton in 0/0. Initial Condition (x,y) = (0.5,0) (vx,vy) = (0,0.2)'
set xlabel "x"
set ylabel "y"
plot "solution_lf.txt" using 2:4 with lines lc rgb "red" title "path of electron"

