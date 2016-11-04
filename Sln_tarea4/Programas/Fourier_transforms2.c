/*
Xibelly Eliseth Mosquera Escobar

Solucion punto 5 Tarea 4 COMPUTACION CIENTIFICA AVANZADA

1) Muestrear la funcion f(x) = x*sin(x) entre (0,2*M_PI)
con una frecuencia de

a) 1

b) 0.5

c) 0.1

2) Calcular la DFT 1D en cada caso -> usamos librerias de FFTW3
   e imprimos en disco 

3) Comparar resultado numerico y analitico
 

  A modo de prueba se calcula tanto la DFT con la IDFT y se comparan los resultados obtenidos, ver graficas.

 */


/*
  NOTA: Este programa resuelve el CASO b)
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>

//Variables globales//

int N;

FILE *input = NULL;
FILE *output = NULL;
FILE *output2 = NULL;
fftw_complex *in1=NULL, *in2=NULL, *in3=NULL;

//Estructuras//

struct data   
{
  double *x, *f;
  
}; 
struct data datos; 


//Funciones//

void fourier(int N)
{
  int i;
  fftw_complex *in=NULL, *out=NULL, *out2;
  fftw_plan my_plan, my_plan2;
  
  in  = (fftw_complex *) malloc((size_t) N*sizeof(fftw_complex));
  out  = (fftw_complex *) malloc((size_t) N*sizeof(fftw_complex));
  out2  = (fftw_complex *) malloc((size_t) N*sizeof(fftw_complex));
  
  for(i=0; i<N; i++)
    {
      in[i][0] = datos.f[i];       //Parte real
      in[i][1] = 0.0;              //Parte compleja
      
      //printf("%lf %lf\n", in[i][0], in[i][1]);            
    }
  
  
  my_plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE); //DFT:Discrete Fourier Transform

  fftw_execute(my_plan);

  output = fopen("FT_0.5.dat","w");

  for(i=0; i<N; i++)
    {
      fprintf(output,"%lf %lf\n", out[i][0] , out[i][1]);
      //printf("%lf %lf\n", out[i][0], out[i][1]);
    }

  fclose(output);
  
  printf("THE STATE OF DFT IS: SUCESS\n");
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  my_plan2 = fftw_plan_dft_1d(N, out, out2, FFTW_BACKWARD, FFTW_ESTIMATE);//IDFT:Inverse Discrete Fourier Transform

  fftw_execute(my_plan2);

  output2 = fopen("IFT_0.5.dat","w");

  for(i=0; i<N; i++)
    {
      fprintf(output2,"%lf %lf %lf\n", out2[i][0]/N, datos.x[i],out2[i][1]/N);
      //printf("%lf %lf\n", out2[i][0]/N, out2[i][1]/N);
    }

  fclose(output2);
  
  fftw_destroy_plan(my_plan);
  
  fftw_free(out);
  
   
  printf("THE STATE OF IDFT IS: SUCESS\n");
  
}



//Programa Principal


int main (int argc, char *argv[])
{
  
  int i;
  double xo, yo;
  double dx;
    

  
  //Carga de parametros//
  
  N = 20;

  datos.x = (double *) malloc(N *sizeof(double));

  datos.f = (double *) malloc(N *sizeof(double));

    
  dx = 2.0 * M_PI/ (1.0 * N);
  
  
  //SE MUESTREA LA FUNCION  f(x) = x*sin(x)//

  //CASO b)

 

  input=fopen("fdex_fourier0.5.dat","w");
  
  for (i = 0; i < N; i++)
    {
      datos.x[i] = dx * (i+0.5);

      datos.f[i] =  datos.x[i] *sin (datos.x[i]);

      fprintf(input,"%g %g\n", datos.x[i], datos.f[i]);
      
            
    }
  fclose(input);

  printf("THE STATE OF SAMPLE DATA FOR FRECUENCY 0.5 IS: SUCESS\n");

  //Se calculan la respectiva transformada de Fourier//

  fourier(N);

  //Se libera memoria usada
  

  free(datos.x);
  free(datos.f);
  fftw_free(in1);
  fftw_free(in2);
  fftw_free(in3);

  return 0;
}

