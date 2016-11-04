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


#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>


FILE *lee = NULL;
FILE *lee1 = NULL;
FILE *lee2 = NULL;


//Variables Globales



unsigned int N;
char *file;

//Estructura de datos
struct data {
  size_t n;
  double * y;
  double * sigma;
};


//Funciones
int expb_f (const gsl_vector * x, void *data, gsl_vector * f) //Funcion que retorna un puntero el cual guarda los valores de f(t,params)
{
  int nread;
  
  size_t n = ((struct data *)data)->n;
  double *y = ((struct data *)data)->y;

  double *sigma = ((struct data *) data)->sigma;
  
  double A = gsl_vector_get (x, 0);
  
  double b = gsl_vector_get (x, 1);

  double X, F_x;  
  
  size_t i;
          

  lee1 = fopen(file,"r");
  
  for (i =0; i < n; i++)
    {
      nread = fscanf(lee1,"%lf %lf",&X, &F_x);

      //printf("%lf %lf\n",X, F_x);
      
      /* Model : Yi =  A * t * sin(b * t)*/

      double t = X;
     
      double Yi = A  * t * sin(b * t) ;
      
      gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
     
    }
    
  fclose(lee1);

  printf("THE MODEL IS CREATED CORRECTLY\n");
  
  return GSL_SUCCESS;

  
}

int expb_df (const gsl_vector * x, void *data, gsl_matrix * J)//Funcion que almacena la matriz resultante J
{

  int nread;
  size_t n = ((struct data *)data)->n;
  double *sigma = ((struct data *) data)->sigma;
  
  double A = gsl_vector_get (x, 0);
  double b = gsl_vector_get (x, 1);

  double X, F_x;  
  size_t i;
   
  
  
  lee2 = fopen(file,"r");
  
  for (i = 0; i <n; i++)
    {
      nread = fscanf(lee2,"%lf %lf",&X,&F_x);
      
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /*       Yi = A * t * sin(B * t)   */
      /* and the xj are the parameters (A,b) */

      //double pos =i;
      double t = X;
      double s = sigma[i];
      double e = t * sin(b * t); //Derivada respecto al parametro A
      double e2 = t * t * cos(b * t); //Derivada respecto al parametro b
      
      gsl_matrix_set (J, i, 0, e/s); 
      

      gsl_matrix_set (J, i, 1, A * e2 /s);  
      
      
      
      
    }

  printf("THE MATRICES ARE STORED CORRECTLY\n");
  fclose(lee2);
  return GSL_SUCCESS;
}

int expb_fdf (const gsl_vector * x, void *data,gsl_vector * f, gsl_matrix * J)//Se crea la funcion que sera minimizada
{
  expb_f (x, data, f);
  expb_df (x, data, J);

   
  return GSL_SUCCESS;
}

//The main part of the program sets up a Levenberg-Marquardt solver and some simulated data. with standard deviation = 0.1 over a range of all dataset timesteps. The initial guess for the parameters is chosen as (2.0, 1.0).


void print_state (size_t iter, gsl_multifit_fdfsolver * s);


//PROGRAMA PRINCIPAL//

int main (int argc, char *argv[])
{

  int nread;
 
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int i, iter = 0;
   
  const size_t p = 2;
  
  gsl_matrix *covar = gsl_matrix_alloc (p, p); //Se establece la matriz de covarianza(Da el error estadistico sobre los mejores parametros de ajuste)

  double pos;
  double *y, *sigma, X, F_x, *xi;  
 
  gsl_multifit_function_fdf f;
  double x_init[2] = { 2.0, 1.0};
  gsl_vector_view x = gsl_vector_view_array (x_init, p);
  const gsl_rng_type * type;
  gsl_rng * r;

  //Carga de Parametros//
  
  printf("%d\n",argc);
  
  if(argc != 3)
    {
      printf("ERROR--> use as:\n");
      printf("%s Nline file\n",argv[0]);
      exit(0);  
    }

  N =  atoi(argv[1]);

  const size_t n = N;
  
  file = argv[2];

  printf("%s %d %s\n",argv[0], N, file);

  
  y = (double *) malloc(N *sizeof(double));

  xi = (double *) malloc(N *sizeof(double));

  sigma = (double *) malloc(N *sizeof(double));

  struct data d = {n, y, sigma};
   
  gsl_rng_env_setup();
  
  type = gsl_rng_default;
  r = gsl_rng_alloc (type);
  
  f.f = &expb_f;
  f.df = &expb_df;
  f.fdf = &expb_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

    
  /* Lectra de los datos a ser ajustados */
  
  lee = fopen(file,"r");
 
  for (i = 0; i < n; i++)
    {
     
      nread = fscanf(lee,"%lf %lf",&X, &F_x);

      xi[i] = X;
      
      y[i] = F_x;

      sigma[i] = 0.1;
            
      
    }

 
  
  T = gsl_multifit_fdfsolver_lmsder; //Se aloca memoria para una instacia de tipo Levenberg-Marquardt solver para los  puntos y los 2 par
  s = gsl_multifit_fdfsolver_alloc (T, n, p);//Se define el solver
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);//Se inicializa el solver
  
  print_state (iter, s);

  printf("hasta aqui bien\n");
    
  do //Se analixa el estado de las iteraciones para continuar con el proceso o detenerlo segun el caso
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);//funcion que realiza una sola iteracion del solver, conduce la iteracion de cada algoritmo
      
      
      printf ("status = %s\n", gsl_strerror (status)); //Selecciona el tipo de error de la libreria strerror
      
      print_state (iter, s);
      
      if (status)
        break;

      status = gsl_multifit_test_delta (s->dx, s->x,1e-4, 1e-4);//Funcion que estima de manera adecuada los parametros
    }
  while (status == GSL_CONTINUE && iter < 1000);

  
  
  gsl_multifit_covar (s->J, 0.0, covar);//Funcion que retorna la matriz de covarianza

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  { 
    double chi = gsl_blas_dnrm2(s->f); //Funcion que pormedio de BLAS,realizando operciones aritmeticas entre matrices calcula la distribuccion chi squared
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 

    printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);

    printf ("A      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
    printf ("b      = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
    
  }

  printf ("status = %s\n", gsl_strerror (status));

  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
  gsl_rng_free (r);
   fclose(lee);
  return 0;
}

void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3u x = % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
	  gsl_blas_dnrm2 (s->f));
}
