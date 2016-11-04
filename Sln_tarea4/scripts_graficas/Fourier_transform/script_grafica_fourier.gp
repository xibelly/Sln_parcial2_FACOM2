set terminal png
set output "grafica_fourier0.1.png" 
set pointsize 1.5
set xlabe 'frecuncies'
set ylabe 'Real FT'
set key top right
plot "FT_0.1.dat" u 1 w lp title "FT_0.1.dat" 
