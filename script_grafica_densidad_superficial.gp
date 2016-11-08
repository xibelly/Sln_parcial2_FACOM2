set terminal png
set output "grafica_surface_density.png" 
set pointsize 1.5
set xlabe 'ID_cell'
set ylabe 'Surface Density'
set key top left
plot "properties_grid.dat" u 1:4 w lp title  "properties_grid.dat" 
