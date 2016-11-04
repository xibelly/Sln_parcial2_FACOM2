#!/bin/bash 

echo "Introduce the data1:"
read data1

echo "Introduce the output file:"
read archivo_salida


echo "set terminal png" > script_grafica_fourier.gp
echo "set output "$archivo_salida" " >> script_grafica_fourier.gp
echo "set pointsize 1.5" >> script_grafica_fourier.gp
echo "set xlabe 'frecuncies'"  >> script_grafica_fourier.gp
echo "set ylabe 'Real FT'"  >> script_grafica_fourier.gp
echo "set key top right" >> script_grafica_fourier.gp

echo "plot "$data1" u 1 w lp title "$data1" "  >> script_grafica_fourier.gp

gnuplot script_grafica_fourier.gp
