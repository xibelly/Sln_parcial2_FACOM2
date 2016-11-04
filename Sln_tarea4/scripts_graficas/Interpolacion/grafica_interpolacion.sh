#!/bin/bash 

echo "Introduce the data1:"
read data1

echo "Introduce the data2:"
read data2

echo "Introduce the output file:"
read archivo_salida


echo "plot "$data1" title  "$data1" "  > script_grafica_interpolacion.gp

echo "set terminal png" >> script_grafica_interpolacion.gp
echo "set output "$archivo_salida" " >> script_grafica_interpolacion.gp
echo "set pointsize 1.5" >> script_grafica_interpolacion.gp
echo "set xlabe 'x'"  >> script_grafica_interpolacion.gp
echo "set ylabe 'f(x)'"  >> script_grafica_interpolacion.gp
echo "set key top right" >> script_grafica_interpolacion.gp

echo "replot "$data2" title "$data2" "  >> script_grafica_interpolacion.gp

gnuplot script_grafica_interpolacion.gp
