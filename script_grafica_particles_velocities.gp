set terminal png
set output "grafica_part_velocities.png" 
set pointsize 1.5
set xlabe 'distance to center'
set ylabe 'Velocities'
set key top left
plot "velocities.dat" u 1:2 title  "velocities.dat" 
