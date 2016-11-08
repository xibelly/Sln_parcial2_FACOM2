#!/bin/bash 

echo "Introduce the data1:"
read data1


echo "Introduce the output file:"
read archivo_salida



echo "set terminal png" > script_grafica_densidad_superficial.gp
echo "set output "$archivo_salida" " >> script_grafica_densidad_superficial.gp
echo "set pointsize 1.5" >> script_grafica_densidad_superficial.gp
echo "set xlabe 'ID_cell'"  >> script_grafica_densidad_superficial.gp
echo "set ylabe 'Surface Density'"  >> script_grafica_densidad_superficial.gp
echo "set key top left" >> script_grafica_densidad_superficial.gp

echo "plot "$data1" u 1:4 w lp title  "$data1" "  >> script_grafica_densidad_superficial.gp

gnuplot script_grafica_densidad_superficial.gp
