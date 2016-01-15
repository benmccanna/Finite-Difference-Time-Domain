set term png enhanced
set pm3d map

set output sprintf("%s.png", filetitle)
set title filetitle
set xlabel "z / c {/Symbol D} t"
set ylabel "time / {/Symbol D} t"

set xrange [0:*]
set yrange [0:*]
set autoscale fix

stats filename matrix using 3 nooutput
cbmax = (abs(STATS_min) > abs(STATS_max) ? abs(STATS_min) : abs(STATS_max))
# set cbrange [-cbmax:cbmax]

splot filename matrix