set term png enhanced
set pm3d map

set output sprintf("%s.png", filetitle)
set title filetitle
set xlabel "z / c {/Symbol D} t"
set ylabel "time / {/Symbol D} t"

if (strstrt(filetitle, "(Ey Field)")) {
  set cblabel "Field strength / V m^{-1}"
} else {
  set cblabel "Field strength / A m^{-1}"
}

set xrange [0:*]
set yrange [0:*]
set autoscale fix

stats filename matrix using 3 nooutput
# cbmax = (abs(STATS_min) > abs(STATS_max) ? abs(STATS_min) : abs(STATS_max))
# set cbrange [-cbmax:cbmax]

# stop cblabel falling off the edge of the image
set rmargin at screen 0.80
set lmargin at screen 0.15

splot filename matrix

# Also output tex files:
set term postscript enhanced color

set title filetitle
unset key
set xlabel "z / c {/Symbol D} t"
set ylabel "time / {/Symbol D} t"

if (strstrt(filetitle, "(Ey Field)")) {
  set cblabel "Field strength / V m^{-1}"
} else {
  set cblabel "Field strength / A m^{-1}"
}

set xrange [0:*]
set yrange [0:*]
set autoscale fix

stats filename matrix using 3 nooutput
# cbmax = (abs(STATS_min) > abs(STATS_max) ? abs(STATS_min) : abs(STATS_max))
# set cbrange [-cbmax:cbmax]

# stop cblabel falling off the edge of the image
set rmargin at screen 0.80
set lmargin at screen 0.15

set output sprintf("%s.eps", filetitle)
splot filename matrix