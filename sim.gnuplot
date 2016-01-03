reset
set terminal gif animate delay 100
set output "animate.gif"
n=24    #25 frames
set xrange [0:200]
set style data lines
do for [i=0:n] {splot  sprintf('sim%d.dat', i*10) using 1 with linespoints; pause 0.5} 
set output