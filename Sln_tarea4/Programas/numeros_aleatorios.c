#include <stdio.h>
/*
Xibelly Eliseth Mosquera Escobar

Solucion punto 2 Tarea 4 COMPUTACION CIENTIFICA AVANZADA

1) Vamos a generar 1000 numeros aleatorios entre 0 y 2*M_PI

2) Con dichos puntos genramos la funcion f(x) = x*sin(x) -> se escribe en disco

3)Interpolamos la funcion anterior con la interpolacion que tiene la menor diferencia cuadratica media -> akima_periodico. -> se imprime en disco.

Nota: Para interpolar debemos ordenar de forma creciente los datos, por lo que usaremos gsl_sort.

NOta: Se hacen grafiacas de la funcion original con la interpolacion obtenida a traves de cada uno de los metodos de generacion de numeros aleatorios utilizados.

 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include<gsl/gsl_rng.h>
#include <gsl/gsl_sort_int.h>
#include <gsl/gsl_sort_double.h>


#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#define N 1000

//Variables globales//

FILE *out = NULL;
FILE *out1 = NULL;
FILE *out2 = NULL;
FILE *out3 = NULL;
FILE *out4 = NULL;
FILE *out5 = NULL;
FILE *out6 = NULL;
FILE *out7 = NULL;
FILE *out8 = NULL;
FILE *out9 = NULL;
FILE *out10 = NULL;

double *x, *y;
double value;
  
gsl_interp_accel *acelerador;
gsl_spline *interpolador;


//Estructuras//



struct data_result   
{
  
  double *value;
  
}; 
struct data_result RESULT; 

 
//Funciones//

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


int random_gfsr4()
{
  
  int i;
  long seed;
  gsl_rng *rng;  //generador de # aleatorios

   
  rng = gsl_rng_alloc (gsl_rng_gfsr4);     
  seed = time(NULL)*getpid();    //semilla que cambia con el tiempo
  gsl_rng_set (rng, seed);       //se establce la semilla
  
  for(i=0; i<N; i++)
    {
      RESULT.value[i] = gsl_rng_uniform (rng);

    }
 
  gsl_rng_free (rng);

  printf("STATE OF RANDOM NUMBER GENERATION BY 'GFSR4': SUCCESS\n");
 
  return 0;
  
}

int random_taus()
{
  int i;
  long seed;
  gsl_rng *rng;  //generador de # aleatorios

  
  rng = gsl_rng_alloc (gsl_rng_taus);     
  seed = time(NULL)*getpid();    //semilla que cambia con el tiempo
  gsl_rng_set (rng, seed);       //se establce la semilla
  
  for(i=0; i<N; i++)
    {
      RESULT.value[i] =  gsl_rng_uniform (rng);

    }
 
  gsl_rng_free (rng);

  printf("STATE OF RANDOM NUMBER GENERATION BY 'TAUS': SUCCESS\n");
 
  return 0;
  
}

int random_rand48()
{
  int i;
  long seed;
  gsl_rng *rng;  //generador de # aleatorios

  
  rng = gsl_rng_alloc (gsl_rng_rand48);     
  seed = time(NULL)*getpid();    //semilla que cambia con el tiempo
  gsl_rng_set (rng, seed);       //se establce la semilla
  
  for(i=0; i<N; i++)
    { 
      RESULT.value[i] =gsl_rng_uniform (rng);

    }
 
  gsl_rng_free (rng);                           

  printf("STATE OF RANDOM NUMBER GENERATION BY 'RAND48': SUCCESS\n");
  
  return 0;
  
}

int random_mt19937()
{
  int i;
  long seed;
  gsl_rng *rng;  //generador de # aleatorios

  
  rng = gsl_rng_alloc (gsl_rng_mt19937);     
  seed = time(NULL)*getpid();    //semilla que cambia con el tiempo
  gsl_rng_set (rng, seed);       //se establce la semilla
  
  for(i=0; i<N; i++)
    {
      RESULT.value[i] = gsl_rng_uniform (rng);

    }
 
  gsl_rng_free (rng);                           

  printf("STATE OF RANDOM NUMBER GENERATION BY 'MT19937': SUCCESS\n");
  
  return 0;
  
}

int random_ranlxs0()
{
  int i;
  long seed;
  gsl_rng *rng;  //generador de # aleatorios

  
  rng = gsl_rng_alloc (gsl_rng_ranlxs0);     
  seed = time(NULL)*getpid();    //semilla que cambia con el tiempo
  gsl_rng_set (rng, seed);       //se establce la semilla
  
  for(i=0; i<N; i++)
    {
      RESULT.value[i] =gsl_rng_uniform (rng);

    }
 
  gsl_rng_free (rng);                           

  printf("STATE OF RANDOM NUMBER GENERATION BY 'RANLXS0': SUCCESS\n");
  
  return 0;
  
}

int random_ranlxs1()
{
  int i;
  long seed;
  gsl_rng *rng;  //generador de # aleatorios

  
  rng = gsl_rng_alloc (gsl_rng_ranlxs1);     
  seed = time(NULL)*getpid();    //semilla que cambia con el tiempo
  gsl_rng_set (rng, seed);       //se establce la semilla
  
  for(i=0; i<N; i++)
    {
      RESULT.value[i] = gsl_rng_uniform (rng);

    }
 
  gsl_rng_free (rng);

  printf("STATE OF RANDOM NUMBER GENERATION BY 'RANLXS1': SUCCESS\n");
 
  return 0;
  
}

int random_mrg()
{
  int i;
  long seed;
  gsl_rng *rng;  //generador de # aleatorios

  
  rng = gsl_rng_alloc (gsl_rng_mrg);     
  seed = time(NULL)*getpid();    //semilla que cambia con el tiempo
  gsl_rng_set (rng, seed);       //se establce la semilla
  
  for(i=0; i<N; i++)
    {
      RESULT.value[i] = gsl_rng_uniform (rng);

    }
 
  gsl_rng_free (rng);                           

  printf("STATE OF RANDOM NUMBER GENERATION BY 'MRG': SUCCESS\n");
  
  return 0;
  
}

int random_ranlux()
{
  int i;
  long seed;
  gsl_rng *rng;  //generador de # aleatorios

  
  rng = gsl_rng_alloc (gsl_rng_ranlux);     
  seed = time(NULL)*getpid();    //semilla que cambia con el tiempo
  gsl_rng_set (rng, seed);       //se establce la semilla
  
  for(i=0; i<N; i++)
    {
      RESULT.value[i] = gsl_rng_uniform (rng);

    }
 
  gsl_rng_free (rng);                           

  printf("STATE OF RANDOM NUMBER GENERATION BY 'RANLUX': SUCCESS\n");
  
  return 0;
  
}

int random_ranlxd1()
{
  int i;
  long seed;
  gsl_rng *rng;  //generador de # aleatorios

  
  rng = gsl_rng_alloc (gsl_rng_ranlxd1);     
  seed = time(NULL)*getpid();    //semilla que cambia con el tiempo
  gsl_rng_set (rng, seed);       //se establce la semilla
  
  for(i=0; i<N; i++)
    {
      RESULT.value[i] = gsl_rng_uniform (rng);

    }
 
  gsl_rng_free (rng);                           
 
  return 0;
  
}

int random_cmrg()
{
  int i;
  long seed;
  gsl_rng *rng;  //generador de # aleatorios

  
  rng = gsl_rng_alloc (gsl_rng_cmrg);     
  seed = time(NULL)*getpid();    //semilla que cambia con el tiempo
  gsl_rng_set (rng, seed);       //se establce la semilla
  
  for(i=0; i<N; i++)
    {
      RESULT.value[i] = gsl_rng_uniform (rng);

    }
 
  gsl_rng_free (rng);                           

  printf("STATE OF RANDOM NUMBER GENERATION BY 'CMRG': SUCCESS\n");
  
  return 0;
  
  }

//Programa Principal//

int main (int argc, char *argv[]) 
{
  int i;
  double xo, yo;
  double *xx, *ff;
  double dx;
  double *X, *F;
  
  //Carga de parametros//

  
  RESULT.value  = (double *) malloc(N *sizeof(double));

  X  = (double *) malloc(N *sizeof(double));

  F  = (double *) malloc(N *sizeof(double));
  
  
  y = (double *) malloc(N *sizeof(double));

  x = (double *) malloc(N *sizeof(double));

  
  xx  = (double *) malloc(N *sizeof(double));
  
  ff  = (double *) malloc(N *sizeof(double));

  

  size_t p[N];

  int k = N;

  dx = 2 * M_PI / N;

  out = fopen("fdex_random.dat","w");

  for(i=0; i<N; i++)
    {

      xx[i] =(dx) * i;
      
      ff[i] = (xx[i])* sin (xx[i]);

      fprintf (out,"%g %g\n", xx[i], ff[i]);
      
    }

  fclose(out);
  
  
  //Generacion de numeros aleatorios//

  //GSFR4

  random_gfsr4();

  out1 = fopen("random_gfsr4.dat","w");

  for(i=0; i<N; i++)
    {

      X[i] = 2.0 * M_PI * RESULT.value[i];
      
      F[i] = X[i] * sin (X[i]);
      
    }

  
  gsl_sort_smallest_index (p, k, X, 1, N);//Se ordenan los datos de forma creciente

  //INTERPOLACION AKIMA PERIODICA

  for(i=0; i<N; i++)
    {
      x[i] = X[p[i]];

      y[i] = F[p[i]];
    }
    
    
  for (xo = x[0]; xo < x[N-1]; xo = xo + 0.01)
    {
      yo = interpolador_akima_per(xo);
      fprintf (out1,"%g %g\n", xo, yo);
    }
  
  
  fclose(out1);

  
  
  //TAUS

  random_taus();

  out2 = fopen("random_taus.dat","w");

  for(i=0; i<N; i++)
    {

      X[i] = 2.0 * M_PI * RESULT.value[i];
      
      F[i] = X[i] * sin (X[i]);
      
    }

  
  gsl_sort_smallest_index (p, k, X, 1, N);//Se ordenan los datos de forma creciente

  //INTERPOLACION AKIMA PERIODICA

  for(i=0; i<N; i++)
    {
      x[i] = X[p[i]];

      y[i] = F[p[i]];
    }
    
    
  for (xo = x[0]; xo < x[N-1]; xo = xo + 0.01)
    {
      yo = interpolador_akima_per(xo);
      fprintf (out2,"%g %g\n", xo, yo);
    }
  
  
  fclose(out2);

  
  
  //RAND48

  random_rand48();

  out3 = fopen("random_rand48.dat","w");

  for(i=0; i<N; i++)
    {

      X[i] = 2.0 * M_PI * RESULT.value[i];
      
      F[i] = X[i] * sin (X[i]);
      
    }

  
  gsl_sort_smallest_index (p, k, X, 1, N);//Se ordenan los datos de forma creciente

  //INTERPOLACION AKIMA PERIODICA

  for(i=0; i<N; i++)
    {
      x[i] = X[p[i]];

      y[i] = F[p[i]];
    }
    
    
  for (xo = x[0]; xo < x[N-1]; xo = xo + 0.01)
    {
      yo = interpolador_akima_per(xo);
      fprintf (out3,"%g %g\n", xo, yo);
    }
  
  
  fclose(out3);

  //MT19937

  random_mt19937();
  
  out4 = fopen("random_mt19937.dat","w");

  for(i=0; i<N; i++)
    {

      X[i] = 2.0 * M_PI * RESULT.value[i];
      
      F[i] = X[i] * sin (X[i]);
      
    }

  
  gsl_sort_smallest_index (p, k, X, 1, N);//Se ordenan los datos de forma creciente

  //INTERPOLACION AKIMA PERIODICA

  for(i=0; i<N; i++)
    {
      x[i] = X[p[i]];

      y[i] = F[p[i]];
    }
    
    
  for (xo = x[0]; xo < x[N-1]; xo = xo + 0.01)
    {
      yo = interpolador_akima_per(xo);
      fprintf (out4,"%g %g\n", xo, yo);
    }
  
  
  fclose(out4);

  
  
  //RANLXS0

  random_ranlxs0();
  
  out5 = fopen("random_ranlxs0.dat","w");

  for(i=0; i<N; i++)
    {

      X[i] = 2.0 * M_PI * RESULT.value[i];
      
      F[i] = X[i] * sin (X[i]);
      
    }

  
  gsl_sort_smallest_index (p, k, X, 1, N);//Se ordenan los datos de forma creciente

  //INTERPOLACION AKIMA PERIODICA

  for(i=0; i<N; i++)
    {
      x[i] = X[p[i]];

      y[i] = F[p[i]];
    }
    
    
  for (xo = x[0]; xo < x[N-1]; xo = xo + 0.01)
    {
      yo = interpolador_akima_per(xo);
      fprintf (out5,"%g %g\n", xo, yo);
    }
  
  
  fclose(out5);


  //RANLXS1
  
  random_ranlxs1();
  
  out6 = fopen("random_ranlxs1.dat","w");

  for(i=0; i<N; i++)
    {

      X[i] = 2.0 * M_PI * RESULT.value[i];
      
      F[i] = X[i] * sin (X[i]);
      
    }

  
  gsl_sort_smallest_index (p, k, X, 1, N);//Se ordenan los datos de forma creciente

  //INTERPOLACION AKIMA PERIODICA

  for(i=0; i<N; i++)
    {
      x[i] = X[p[i]];

      y[i] = F[p[i]];
    }
    
    
  for (xo = x[0]; xo < x[N-1]; xo = xo + 0.01)
    {
      yo = interpolador_akima_per(xo);
      fprintf (out6,"%g %g\n", xo, yo);
    }
  

  fclose(out6);

  //MRG

  random_mrg();
  
  out7 = fopen("random_mrg.dat","w");

  for(i=0; i<N; i++)
    {

      X[i] = 2.0 * M_PI * RESULT.value[i];
      
      F[i] = X[i] * sin (X[i]);
      
    }

  
  gsl_sort_smallest_index (p, k, X, 1, N);//Se ordenan los datos de forma creciente

  //INTERPOLACION AKIMA PERIODICA

  for(i=0; i<N; i++)
    {
      x[i] = X[p[i]];

      y[i] = F[p[i]];
    }
    
    
  for (xo = x[0]; xo < x[N-1]; xo = xo + 0.01)
    {
      yo = interpolador_akima_per(xo);
      fprintf (out7,"%g %g\n", xo, yo);
    }
  
  
  fclose(out7);

  //RANLUX
  
  random_ranlux();
  
  out8 = fopen("random_ranlux.dat","w");

 for(i=0; i<N; i++)
    {

      X[i] = 2.0 * M_PI * RESULT.value[i];
      
      F[i] = X[i] * sin (X[i]);
      
    }

  
  gsl_sort_smallest_index (p, k, X, 1, N);//Se ordenan los datos de forma creciente

  //INTERPOLACION AKIMA PERIODICA

  for(i=0; i<N; i++)
    {
      x[i] = X[p[i]];

      y[i] = F[p[i]];
    }
    
    
  for (xo = x[0]; xo < x[N-1]; xo = xo + 0.01)
    {
      yo = interpolador_akima_per(xo);
      fprintf (out8,"%g %g\n", xo, yo);
    }
  
  
  fclose(out8);


  //RANLXD1

  random_ranlxd1();
  
  out9 = fopen("random_ranlxd1.dat","w");

  for(i=0; i<N; i++)
    {

      X[i] = 2.0 * M_PI * RESULT.value[i];
      
      F[i] = X[i] * sin (X[i]);
      
    }

  
  gsl_sort_smallest_index (p, k, X, 1, N);//Se ordenan los datos de forma creciente

  //INTERPOLACION AKIMA PERIODICA

  for(i=0; i<N; i++)
    {
      x[i] = X[p[i]];

      y[i] = F[p[i]];
    }
    
    
  for (xo = x[0]; xo < x[N-1]; xo = xo + 0.01)
    {
      yo = interpolador_akima_per(xo);
      fprintf (out9,"%g %g\n", xo, yo);
    }
  
  
  fclose(out9);


  //CMRG

  random_cmrg();
  
  out10 = fopen("random_cmrg.dat","w");
  
  for(i=0; i<N; i++)
    {
      
      X[i] = 2.0 * M_PI * RESULT.value[i];
      
      F[i] = X[i] * sin (X[i]);
      
    }
  
  
  gsl_sort_smallest_index (p, k, X, 1, N);//Se ordenan los datos de forma creciente
  
  //INTERPOLACION AKIMA PERIODICA
  
  for(i=0; i<N; i++)
    {
      x[i] = X[p[i]];
      
      y[i] = F[p[i]];
    }
    
    
  for (xo = x[0]; xo < x[N-1]; xo = xo + 0.01)
    {
      yo = interpolador_akima_per(xo);
      fprintf (out10,"%g %g\n", xo, yo);
    }
  
  
  fclose(out10);

 
  
  return 0;
}
