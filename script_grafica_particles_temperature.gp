set terminal png
set output "grafica_part_Temp.png" 
set pointsize 1.5
set xlabe 'r'
set ylabe 'Temperature'
set key top left
plot "velocities.dat" u 1:4 title  "velocities.dat" 

