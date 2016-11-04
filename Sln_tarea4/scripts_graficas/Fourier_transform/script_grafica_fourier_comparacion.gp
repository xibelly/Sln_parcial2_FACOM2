plot "fdex_fourier0.1.dat" title  "fdex_fourier0.1.dat" 
set terminal png
set output "grafica_fourier0.1_comparacion.png" 
set pointsize 1.5
set xlabe 'x'
set ylabe 'f(x)'
set key top right
replot "IFT_0.1.dat" u 2:1 title "IFT_0.1.dat" 
