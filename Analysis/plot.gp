# ===============================================
# plot command for gnuplot
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

set title 'test of borisPusher with dt = 0.01 and tEnd = 100'
set xlabel "vx"
set ylabel "vy"
plot "test.txt" using 4:5 with lines lc rgb "red" title "velocity"