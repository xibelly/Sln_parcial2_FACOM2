/*
Xibelly Eliseth Mosquera Escobar

Solucion punto 4 Tarea 4 COMPUTACION CIENTIFICA AVANZADA

Tome la misma función anterior y sumele un ruido distribuido gausianamente con una varianza de
0.5
Use las rutinas de ajuste de gsl (no lineal) para ajustar la función g(x)=a*x*sin(b*x) a los datos.
Obtenga los valores de los parámetros a y b que mejor describen el modelo.

1) leer los datos interpolados y sumarle el ruidio gaussiano

2) usar las librerias de GSL para hacer un ajuste -no lineal- a dichos datos a la fucion g(x) = a * x * sin(b * x) y determinar los valores de los parametros a y b, que mejor describen el modelo.

 */

/*
 Este programa se encarga de resolver el punto 1)
 */

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


FILE *write = NULL;

//Variables Globales

int N;
char *file;

struct data   
{
  double x, f;
  
}; 
struct data *datos; 
 
//Llamados a subrutinas//

#include "input.c"


//PROGRAMA PRINCIPAL//

int main (int argc, char *argv[])
{
  int i;
  double *y, *x;
  const gsl_rng_type * T;
  gsl_rng * r;
  unsigned int k;

  double mu = 0.5;

   
  //Carga de Parametros//
  
  printf("%d\n",argc);
  
  if(argc != 3)
    {
      printf("ERROR--> use as:\n");
      printf("%s Nline file\n",argv[0]);
      exit(0);  
    }

  N =  atoi(argv[1]);

  file = argv[2];

  printf("%s %d %s\n",argv[0], N, file);

  datos  = (struct data *) malloc(N *sizeof(struct data));
  
  y = (double *) malloc(N *sizeof(double));

  x = (double *) malloc(N *sizeof(double));
  
  /* Lectra de los datos a ser ajustados */
  
  read_file(file,N);

  /* Se genera el ruido gaussiano a la funcion interpolada */

  // create a generator chosen by the environment variable GSL_RNG_TYPE 

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  //print N random variates chosen from the poisson distribution with mean parameter mu

  write = fopen("fdex_gaussian_noise.dat","w");

  for (i = 0; i < N; i++)
    {
          
      k = gsl_ran_poisson (r, mu); //Gaussian distribution

      y[i] = datos[i].f + k;  //Funcion interpolada con ruido gaussiano
                
      x[i] = datos[i].x;

      fprintf (write,"%g %g \n", x[i] , y[i]);
    }

  printf("STATE OF RANDOM NOISE GENERATION IS: SUCESS\n");
  
  gsl_rng_free (r);
  
  return 0;
}


