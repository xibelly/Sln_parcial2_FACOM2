set terminal png
set output "grafica_RealIDFT_surface_density.png" 
set pointsize 1.5
set xlabe ''
set ylabe 'Real IDFT'
set key top left
plot "IDFT_2D.dat" u 1 w lp title  "IDFT_2D.dat" 

 
set terminal png
set output "grafica_ImaIDFT_surface_density.png" 
set pointsize 1.5
set xlabe ''
set ylabe 'Ima IDFT'
set key top left
plot "IDFT_2D.dat" u 2 w lp title  "IDFT_2D.dat"