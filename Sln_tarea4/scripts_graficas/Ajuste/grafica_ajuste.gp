f(x) = a * x * sin(b*x) #funcion a ajustar

a = 0.91364; b =0.99652 ;    #valores parametros modelo


#fit f(x) 'densidad_halo19.dat' u 1:2 via a,b 

plot "fdex_gaussian_noise.dat" u 1:2 w p lc rgb "red" title 'fdex_gaussian_noise' 

set terminal png
set output 'grafica_ajuste.png'
set pointsize 1.5
set xlabel 'x'
set ylabel 'f(x)'
set title ' Ajuste'
set label 'a = 0.91364' at 5,5
set label 'b = 0.99652' at 5,6



replot  f(x)  lc rgb "blue" title 'Model'
