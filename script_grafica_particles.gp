plot "positions.dat"  title  "positions.dat" 
set terminal png
set output "grafica_particle_distribution.png" 
set pointsize 1.5
set xlabe 'x'
set ylabe 'y'
set key top left
replot "new_positions.dat"  title "new_positions.dat" 
