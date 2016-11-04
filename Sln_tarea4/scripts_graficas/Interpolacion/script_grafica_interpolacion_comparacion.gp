plot "fdex.dat" w d  title  "fdex.dat" 
replot "inter_lineal.dat" w l title "inter_lineal.dat" 
replot "inter_poly.dat" w l title "inter_poly.dat" 
replot "inter_cspline.dat" w l title "inter_cspline.dat" 
replot "inter_cspline_per.dat" w l title "inter_cspline_per.dat" 
replot "inter_akima.dat" w l title "inter_akima.dat" 
set terminal png
set output "grafica_interpolacion_comparacion.png" 
set pointsize 1.5
set xlabe 'x'
set ylabe 'f(x)'
set key top right
replot "inter_akima_per.dat" w l title "inter_akima_per.dat" 
