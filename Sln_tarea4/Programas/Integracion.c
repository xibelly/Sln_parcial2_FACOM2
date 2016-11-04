/*
Xibelly Eliseth Mosquera Escobar

Solucion punto 3 Tarea 4 COMPUTACION CIENTIFICA AVANZADA

Usar la misma funcion interpolante para calcular la integral de f(x) usando
cuadraturas adaptativas de Gauss-konrod. Cree una funcion cuyo integrando se la funcion interpolante y usela para hace la integracion.

Pase los demas parametros del proceso de interpolacion como parametros a la funcion. Compare el resultado de la integracion  con lo que se obtine de la rutina de interpolacion de gsl que evalua la integral del conjunto de puntos a interpolar y con el valor analito exacto.

1) Crear una funcion que lea los datos interpolados que presentan la menor diferencia cuadratica media 

2) Integrar dichos datos de 0 a 2*M_PI usando cuadraturas adaptativas. 
 */



#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>


//Variables globales//

int N;
char *file;

//Estructuras//

struct data   
{
  double x, f;
  
}; 
struct data *datos; 
 
//Llamados a subrutinas//

#include "input.c"

//Funciones//

double f ()//function x*sin(x) 
{
  int i;
  double f;

  for(i=0; i<N; i++)
    {
      f = datos[i].f;

    }
     
  return f;
}


//Programa Principal//

int main (int argc, char *argv[])
{
  
  int i;
  double a, b;
  double result, error, expected;
  size_t limit;
  gsl_function F;
  gsl_integration_workspace * w;
  
  printf("%d\n",argc);
  
  if(argc != 3)
    {
      printf("ERROR--> use as:\n");
      printf("%s Nline file\n",argv[0]);
      exit(0);  
    }

  //Carga de parametros//

  N = atoi(argv[1]);

  file = argv[2];

  printf("%s %d %s\n",argv[0], N, file);

  a = 0.0;

  b = 2.0 * M_PI;

  expected = -6.2832;

  limit = 10000;
  
  F.function = &f; //f(x)=x*sin(x)

  w = gsl_integration_workspace_alloc(limit);//store memory for all things needed to integrate

  datos  = (struct data *) malloc(N *sizeof(struct data));
  
  //Lectura archivos -carga de datos-//
  
  read_file(file,N);

  //Integracion//

  gsl_integration_qags(&F, a,b, 0, 1e-7, 10000,w, &result, &error);//function that integarte the function f(x) using Qags. 

  printf ("result          = % .18f\n", result);
  printf ("exact result    = % .18f\n", expected);
  printf ("estimated error = % .18f\n", error);
  printf ("actual error    = % .18f\n", result - expected);
  printf ("intervals       = %zu\n", w->size);

  
  gsl_integration_workspace_free (w);//free the memory used

  return 0;
}

