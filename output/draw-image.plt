set term png

set xrange [0:*]
set yrange [0:*]
set output sprintf("%s.png", filename)
set title filename
plot filename matrix with image