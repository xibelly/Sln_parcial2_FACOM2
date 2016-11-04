#!/bin/bash 

echo "Introduce the data1:"
read data1

echo "Introduce the output file:"
read archivo_salida


echo "set terminal png" > script_grafica_gaussian_noise.gp
echo "set output "$archivo_salida" " >> script_grafica_gaussian_noise.gp
echo "set pointsize 1.5" >> script_grafica_gaussian_noise.gp
echo "set xlabe 'x'"  >> script_grafica_gaussian_noise.gp
echo "set ylabe 'f(x)'"  >> script_grafica_gaussian_noise.gp
echo "set key top right" >> script_grafica_gaussian_noise.gp

echo "plot "$data1" w lp title "$data1" "  >> script_grafica_gaussian_noise.gp

gnuplot script_grafica_gaussian_noise.gp
