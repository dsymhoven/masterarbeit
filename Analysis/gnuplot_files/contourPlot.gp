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

# line style definitions
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 5   # blue
set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 7   # red

# Set terminal and output
set terminal pngcairo 
set output 'contour.png'
 
# Set various features of the plot
set pm3d
unset surface  # don't need surfaces
set view map
set contour
set key outside
set cntrparam cubicspline  # smooth out the lines
set cntrparam levels 0    # sets the num of contour lines
set pm3d interpolate 20,20 # interpolate the color
 
# Set a nice color palette
set palette rgb 33,13,10
 
# Axes
set xlabel "x"
set ylabel "y"
set xrange [0:30]
set yrange [0:30]
 
# Now plot
splot 'fields.txt' using 1:2:3 notitle with lines lt 1

