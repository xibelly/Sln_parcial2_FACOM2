set terminal png
set output "grafica_ImaDFT_surface_density.png" 
set pointsize 1.5
set xlabe 'frecuencies'
set ylabe 'Imaginary DFT'
set key top left
plot "DFT_2D.dat" u 2 w d title  "DFT_2D.dat" 
