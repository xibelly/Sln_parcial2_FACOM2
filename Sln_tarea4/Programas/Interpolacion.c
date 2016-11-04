/*
Xibelly Eliseth Mosquera Escobar

Solucion punto 1 Tarea 4 COMPUTACION CIENTIFICA AVANZADA

1) Generar una tabla de datos con 20 puntos entre x=0 y x=2*M_PI

2) Interpolar -> vamos a generar una funcion para cada tipo de interpolador
   evluar la funcion en puntos igualmente espaciados, con frecuencia de 0.01

   vamos a imprimir en disco cada tipo de interpolacion

3) Graficar los resultados de todos los tipos de interpolador

4) Evaluar la diferencia cuadratica media entre los diferentes metodos de interpolacion usados. Para ello vamos a crear una funcion que recibe las interpolaciones que van a ser comparadas y vamos imprimir en pantalla cada resultado talque el usuario decida cual es el metodo que presenta menos error.

 */


/*
  NOTA: ESTE PROGRAMA SE ENCARGA DE LOS PUNTOS 1), 2)

 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort_float.h>


//Parametros globales//

int N;
double *x, *y;

FILE *in = NULL;
FILE *out1 = NULL;
FILE *out2 = NULL;
FILE *out3 = NULL;
FILE *out4 = NULL;
FILE *out5 = NULL;
FILE *out6 = NULL;
FILE *out7 = NULL;

double value;
  
gsl_interp_accel *acelerador;
gsl_spline *interpolador;
  

//Funciones///

double interpolador_lineal(double xo)
{

  
  acelerador   = gsl_interp_accel_alloc();
  interpolador = gsl_spline_alloc(gsl_interp_linear, (size_t) N);
  gsl_spline_init(interpolador, x, y, (size_t) N);
  
  value = gsl_spline_eval(interpolador, xo, acelerador);
  
  gsl_spline_free(interpolador);
  gsl_interp_accel_free(acelerador);
  
  return value;
  
}


double interpolador_poli(double xo)
{
    
  acelerador   = gsl_interp_accel_alloc();
  interpolador = gsl_spline_alloc(gsl_interp_polynomial, (size_t) N);
  gsl_spline_init(interpolador, x, y, (size_t) N);
  
  value = gsl_spline_eval(interpolador, xo, acelerador);
  
  gsl_spline_free(interpolador);
  gsl_interp_accel_free(acelerador);
  
  return value;
  
}

double interpolador_cspline(double xo)
{
  
  acelerador   = gsl_interp_accel_alloc();
  interpolador = gsl_spline_alloc(gsl_interp_cspline, (size_t) N);
  gsl_spline_init(interpolador, x, y, (size_t) N);
  
  value = gsl_spline_eval(interpolador, xo, acelerador);
  
  gsl_spline_free(interpolador);
  gsl_interp_accel_free(acelerador);
  
  return value;
  
}

double interpolador_cspline_per(double xo)
{
    
  acelerador   = gsl_interp_accel_alloc();
  interpolador = gsl_spline_alloc(gsl_interp_cspline_periodic, (size_t) N);
  gsl_spline_init(interpolador, x, y, (size_t) N);
  
  value = gsl_spline_eval(interpolador, xo, acelerador);
  
  gsl_spline_free(interpolador);
  gsl_interp_accel_free(acelerador);
  
  return value;
  
}

double interpolador_akima(double xo)
{
    
  acelerador   = gsl_interp_accel_alloc();
  interpolador = gsl_spline_alloc(gsl_interp_akima, (size_t) N);
  gsl_spline_init(interpolador, x, y, (size_t) N);
  
  value = gsl_spline_eval(interpolador, xo, acelerador);
  
  gsl_spline_free(interpolador);
  gsl_interp_accel_free(acelerador);
  
  return value;
  
}


double interpolador_akima_per(double xo)
{
  
  
  acelerador   = gsl_interp_accel_alloc();
  interpolador = gsl_spline_alloc(gsl_interp_akima_periodic, (size_t) N);
  gsl_spline_init(interpolador, x, y, (size_t) N);
  
  value = gsl_spline_eval(interpolador, xo, acelerador);
  
  gsl_spline_free(interpolador);
  gsl_interp_accel_free(acelerador);
  
  return value;
  
}


int main (int argc, char *argv[])
{
  
  int i;
  double xo, yo;
  double dx;
  
  //Carga de parametros//
  
  N = 20;

  y = (double *) malloc(N *sizeof(double));

  x = (double *) malloc(N *sizeof(double));

  //Se crea la funcion f(x) = x*sin(x)//

  dx = 2 * M_PI/ N;

  in=fopen("fdex.dat","w");
  
  for (i = 0; i < N; i++)
    {
      x[i] = dx * i;
      y[i] = (x[i]) *sin (x[i]);
      fprintf(in,"%g %g\n", x[i], y[i]);
    }
  fclose(in);
  
  //Se evaluan las interpolaciones//

  //INTERPOLACION LINEAL
  
  out1=fopen("inter_lineal.dat","w");
  
  for (xo = x[0]; xo < x[N-1]; xo += 0.01)
    {
      yo = interpolador_lineal(xo);
      fprintf (out1,"%g %g\n", xo, yo);
    }

  fclose(out1);
  
   //INTERPOLACION POLINIMICA
  
  out2=fopen("inter_poly.dat","w");
  
  for (xo = x[0]; xo < x[N-1]; xo += 0.01)
    {
      yo = interpolador_poli(xo);
      fprintf (out2,"%g %g\n", xo, yo);
    }

  fclose(out2);

  //INTERPOLACION CSPLINE
  
  out3=fopen("inter_cspline.dat","w");
  
  for (xo = x[0]; xo < x[N-1]; xo += 0.01)
    {
      yo = interpolador_cspline(xo);
      fprintf (out3,"%g %g\n", xo, yo);
    }

  fclose(out3);
  
  //INTERPOLACION CSPLINE PERIODICO
  
  out4=fopen("inter_cspline_per.dat","w");
  
  for (xo = x[0]; xo < x[N-1]; xo += 0.01)
    {
      yo = interpolador_cspline_per(xo);
      fprintf (out4,"%g %g\n", xo, yo);
    }

  fclose(out4);

  //INTERPOLACION AKIMA
  
  out5=fopen("inter_akima.dat","w");
  
  for (xo = x[0]; xo < x[N-1]; xo += 0.01)
    {
      yo = interpolador_akima(xo);
      fprintf (out5,"%g %g\n", xo, yo);
    }

  fclose(out5);

  //INTERPOLACION AKIMA PERIODICA
  
  out6=fopen("inter_akima_per.dat","w");
  
  for (xo = x[0]; xo < x[N-1]; xo += 0.01)
    {
      yo = interpolador_akima_per(xo);
      fprintf (out6,"%g %g\n", xo, yo);
    }
  
 
  fclose(out6);

  //Liberamos memoria//

  free(x);
  free(y);
  
  return 0;
}

