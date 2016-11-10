set terminal png
set output "grafica_RealDFT_surface_density.png" 
set pointsize 1.5
set xlabe 'Frecuencies'
set ylabe 'Real DFT'
set key top left
plot [0:10000]"DFT_2D.dat" u 1 w lp title  "DFT_2D.dat" 

 
set terminal png
set output "grafica_ImaDFT_surface_density.png" 
set pointsize 1.5
set xlabe 'Frecuencies'
set ylabe 'Ima DFT'
set key top left
plot [0:10000]"DFT_2D.dat" u 2 w lp title  "DFT_2D.dat"