set terminal png
set output "grafica_surface_density_3D.png" 
set pointsize 1.5
set xlabe 'cell_posx'
set ylabe 'cell_posy'
set zlabe 'S(x,y)'
set key top left
splot "datos_a_derivar.dat"  w d title  "properties_grid.dat" 
