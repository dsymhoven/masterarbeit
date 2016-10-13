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
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 5   # blue
set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 7   # red

# plot commands

set xlabel "x"
set ylabel "y"
set xrange [0:30]
set yrange [0:30]
plot "test.txt" using 1:2 with points ls 2, \
	"test.txt" using 3:4 with points ls 2, \
	"test.txt" using 5:6 with points ls 2, \
	"test.txt" using 7:8 with points ls 2