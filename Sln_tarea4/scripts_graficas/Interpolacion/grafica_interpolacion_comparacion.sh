#!/bin/bash 

echo "Introduce the data1:"
read data1

echo "Introduce the data2:"
read data2

echo "Introduce the data3:"
read data3

echo "Introduce the data4:"
read data4

echo "Introduce the data5:"
read data5

echo "Introduce the data6:"
read data6

echo "Introduce the data7:"
read data7


echo "Introduce the output file:"
read archivo_salida


echo "plot "$data1" w p  title  "$data1" "  > script_grafica_interpolacion_comparacion.gp

echo "replot "$data2" w l title "$data2" "  >> script_grafica_interpolacion_comparacion.gp

echo "replot "$data3" w l title "$data3" "  >> script_grafica_interpolacion_comparacion.gp

echo "replot "$data4" w l title "$data4" "  >> script_grafica_interpolacion_comparacion.gp

echo "replot "$data5" w l title "$data5" "  >> script_grafica_interpolacion_comparacion.gp

echo "replot "$data6" w l title "$data6" "  >> script_grafica_interpolacion_comparacion.gp


echo "set terminal png" >> script_grafica_interpolacion_comparacion.gp
echo "set output "$archivo_salida" " >> script_grafica_interpolacion_comparacion.gp
echo "set pointsize 1.5" >> script_grafica_interpolacion_comparacion.gp
echo "set xlabe 'x'"  >> script_grafica_interpolacion_comparacion.gp
echo "set ylabe 'f(x)'"  >> script_grafica_interpolacion_comparacion.gp
echo "set key top right" >> script_grafica_interpolacion_comparacion.gp


echo "replot "$data7" w l title "$data7" "  >> script_grafica_interpolacion_comparacion.gp

gnuplot script_grafica_interpolacion_comparacion.gp
