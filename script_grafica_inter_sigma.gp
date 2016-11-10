plot "sigmav_interpolacion_5"  title  "xc_5" 
replot "sigmav_interpolacion_10"  title  "xc_10"
replot "sigmav_interpolacion_15"  title  "xc_15"
replot "sigmav_interpolacion_20"  title  "xc_20"
replot "sigmav_interpolacion_25"  title  "xc_25"
replot "sigmav_interpolacion_30"  title  "xc_30"
replot "sigmav_interpolacion_35"  title  "xc_35"
replot "sigmav_interpolacion_40"  title  "xc_40"
replot "sigmav_interpolacion_45"  title  "xc_45"
set terminal png
set output "grafica_inter_sigma.png" 
set pointsize 1.5
set ylabe 'sigma'
set xlabe 'y'
set key top right
replot "sigmav_interpolacion_50"  title  "xc_50"