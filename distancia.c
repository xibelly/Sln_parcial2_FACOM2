/*Xibelly Eliseth Mosquera Escobar
 * 
 * programa: distancia.c
 * 
 * 
 *Este programa se usa para calcular la distancia entre particulas
 * 
 * 
 * La funcion "distancia" recibe como entradas
 * las coordenadas (xi,yi) y (xj,yj).
 * Ambas entradas se le pasan directamente a la funcion 
 * la cual retorna el valor de la distancia ecludea entre dos puntos dados
 *  
 */

#include<math.h>

double distancia(double xi, double xj, double yi, double yj)
{

  double dist, dx, dy;
  
  dx = (xi-xj)*(xi-xj);
  dy = (yi-yj)*(yi-yj);
  
  
  dist = sqrt(dx + dy); 

  return dist;

}
