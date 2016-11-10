set terminal png
set output "grafica_part_velocities_x.png" 
set pointsize 1.5
set xlabe 'r'
set ylabe 'Vel_x'
set key top left
plot "velocities.dat" u 1:2 title  "velocities.dat" 

set terminal png
set output "grafica_part_velocities_y.png" 
set pointsize 1.5
set xlabe 'r'
set ylabe 'Vel_y'
set key top left
plot "velocities.dat" u 1:3 title  "velocities.dat" 
