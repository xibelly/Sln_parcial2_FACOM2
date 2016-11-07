/*Xibelly Eliseth Mosquera Escobar
 * 
 * programa: derivadas.c
 * 
 * 
 *
 * La funcion "deriva" interpola los datos del campo densidad para luego
 * derivarlos
 *
 * NOTA: interpolamos con akima_periodico 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_double.h>


double deriva(double z, void * params)
{
  int i, nread;
  int N=1024;
  int Nlines;
  double xo, function;
  double *x,*y, *xx, *yy;

  gsl_interp_accel *acelerador;
  gsl_spline *interpolador;
  
  size_t p[N];
  
  int k = N;

    
  y = (double *) malloc(N *sizeof(double));
  x = (double *) malloc(N *sizeof(double));
  
  yy = (double *) malloc(N *sizeof(double));
  xx = (double *) malloc(N *sizeof(double));
  
        
  for(i=0; i< N; i++)
    {
      xx[i] = celda[i].r;
      yy[i] = celda[i].den_super;
      
      //function += celda[i].den_super;
      

    } 
   
  gsl_sort_smallest_index (p, k, yy, 1, N);//Para asegurar la interpolacion ordenamos en forma creciente
  
  Nlines = 0;

  for(i=0; i< N; i++)
    {         
      if( xx[p[i]] == xx[p[i+1]] )continue; //se impone para que los datos sean estrictamente crecientes

      x[i] = xx[p[i]];
      y[i] = yy[p[i]];

      printf("%lf %lf\n", x[i], y[i]);
      Nlines = Nlines +1 ; 
      
   }

  printf("%d\n", Nlines);

  acelerador   = gsl_interp_accel_alloc();
  interpolador = gsl_spline_alloc(gsl_interp_akima_periodic, (size_t) Nlines);
  gsl_spline_init(interpolador, x, y, (size_t) Nlines);

  
  function = gsl_spline_eval(interpolador, 0, acelerador);//se interpolan los datos del campo de densidad
  

  gsl_spline_free(interpolador);
  gsl_interp_accel_free(acelerador);
  
  
  return function;
  
}


 
 
  
  
  
  
