/*Xibelly Eliseth Mosquera Escobar
 * 
 * programa: derivadas.c
 * 
 * 
 *
 * La funcion "deriva" interpola los datos del campo densidad para luego
 * derivarlos -> campo densidad superficial como funcion de la masa
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


double deriva(double z, void * data)
{
  int i, j, l, nread;
  size_t n = ((struct data *)data)->n;

  
  int Nlines;
  double xo, function;
  double *x,*y, *xx, *yy;
  double a, b;
  
  gsl_interp_accel *acelerador;
  gsl_spline *interpolador;
  
  FILE *write = NULL;  
  FILE *read = NULL;  
 
  yy = (double *) malloc(n *sizeof(double));
  xx = (double *) malloc(n *sizeof(double));
  
        
  for(i=0; i< n; i++)
    {
      xx[i] = celda[i].masa;  //Campo densidad superficial en funcion de la masa
      yy[i] = celda[i].den_super;  
      

    } 
   
 //Para asegurar la interpolacion ordenamos en forma creciente
  gsl_sort(xx, 1, n);

  gsl_sort(yy, 1, n);

  Nlines = 0;

  write = fopen("datos_a_derivar.dat","w");

  for(i=0; i<n; i++)
    {         
      if((xx[i] < xx[i+1]) && (yy[i] < yy[i+1]) )//condicion de estrictamente crecientes
	{
	  fprintf(write,"%lf %lf\n", xx[i],yy[i]); 
	  Nlines = Nlines +1 ; 
	}
    }
  fclose(write);
  
  
  y = (double *) malloc(Nlines *sizeof(double));
  x = (double *) malloc(Nlines *sizeof(double));
  
  
  read = fopen("datos_a_derivar.dat","r");
  
  for(j=0; j<Nlines; j++)
    { 
      nread = fscanf(read,"%lf %lf",&x[j], &y[j]);//datos que seran interpolados
            
    }
  
  fclose(read);

  
  acelerador   = gsl_interp_accel_alloc();
  interpolador = gsl_spline_alloc(gsl_interp_akima_periodic, (size_t) Nlines);
  gsl_spline_init(interpolador, x, y, (size_t) Nlines);

  
  function = gsl_spline_eval(interpolador, 0, acelerador);//se interpolan los datos del campo de densidad
  

  gsl_spline_free(interpolador);
  gsl_interp_accel_free(acelerador);
  
  
  return function; //funcion que sera derivada en cada punto de la malla
  
}


 
 
  
  
  
  
