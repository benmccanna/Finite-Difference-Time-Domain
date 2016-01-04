reset
set term png

set output "basic-Ex.png"
plot "basic-Ex.dat" matrix with image

set output "basic-Hy.png"
plot "basic-Hy.dat" matrix with image