set terminal png
set output "grafica_gradx.png" 
set pointsize 1.5
set xlabe 'x'
set ylabe 'S'(x,y)'
set key top left
plot "derivates_of_desity.dat"  title "derivates_of_desity.dat" 
