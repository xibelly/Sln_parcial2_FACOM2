plot "fdex_random.dat" title  "fdex_random.dat" 
set terminal png
set output "grafica_random_taus.png" 
set pointsize 1.5
set xlabe 'x'
set ylabe 'f(x)'
set key top right
replot "random_taus.dat"w l title "random_taus.dat" 
