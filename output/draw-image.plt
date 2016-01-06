set term png

set xrange [0:*]
set yrange [0:*]
set autoscale fix

set output sprintf("%s.png", filename)
set title filename
plot filename matrix with image