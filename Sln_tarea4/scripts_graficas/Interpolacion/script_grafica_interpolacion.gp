plot "fdex.dat" title  "fdex.dat" 
set terminal png
set output "grafica_fdex_inter_lineal.png" 
set pointsize 1.5
set xlabe 'x'
set ylabe 'f(x)'
set key top right
replot "inter_lineal.dat" title "inter_lineal.dat" 
