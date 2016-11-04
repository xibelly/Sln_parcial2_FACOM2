set terminal png
set output "grafica_fdex_gaussian_noise.png" 
set pointsize 1.5
set xlabe 'x'
set ylabe 'f(x)'
set key top right
plot "fdex_gaussian_noise.dat" w lp title "fdex_gaussian_noise.dat" 
