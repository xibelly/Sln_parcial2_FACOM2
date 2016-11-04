#!/bin/bash 

echo "Introduce the data1:"
read data1

echo "Introduce the data2:"
read data2

echo "Introduce the output file:"
read archivo_salida


echo "plot "$data1" title  "$data1" "  > script_grafica_fourier_comparacion.gp

echo "set terminal png" >> script_grafica_fourier_comparacion.gp
echo "set output "$archivo_salida" " >> script_grafica_fourier_comparacion.gp
echo "set pointsize 1.5" >> script_grafica_fourier_comparacion.gp
echo "set xlabe 'x'"  >> script_grafica_fourier_comparacion.gp
echo "set ylabe 'f(x)'"  >> script_grafica_fourier_comparacion.gp
echo "set key top right" >> script_grafica_fourier_comparacion.gp

echo "replot "$data2" u 2:1 title "$data2" "  >> script_grafica_fourier_comparacion.gp

gnuplot script_grafica_fourier_comparacion.gp
