plot "density_interpolacion_5"  title  "xc_5" 
replot "density_interpolacion_10"  title  "xc_10"
replot "density_interpolacion_15"  title  "xc_15"
replot "density_interpolacion_20"  title  "xc_20"
replot "density_interpolacion_25"  title  "xc_25"
replot "density_interpolacion_30"  title  "xc_30"
replot "density_interpolacion_35"  title  "xc_35"
replot "density_interpolacion_40"  title  "xc_40"
replot "density_interpolacion_45"  title  "xc_45"
set terminal png
set output "grafica_inter_density.png" 
set pointsize 1.5
set ylabe 'S(x,y)'
set xlabe 'y'
set key top right
replot "density_interpolacion_50"  title  "xc_50"