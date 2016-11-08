plot "new_positions.dat"  title  "new_positions.dat" 
set terminal png
set output "grafica_grid.png" 
set pointsize 1.5
set xlabe 'x'
set ylabe 'y'
set key top left
replot "grid.dat"  title "grid.dat" 
