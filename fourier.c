/*Xibelly Eliseth Mosquera Escobar
 * 
 * programa: forier.c
 * 
 *
 * La funcion "fourier" calcula la 2D DFT y la 2D IDFT del campo
 * densidad supeficial de masa.
 *
 * NOTA: En principio el campo densidad superficial
 * no es periodico, por lo que se deben imponer PBC's adecuadas.
 * O hacer un padding -> rellenar con ceros regiones del arreglo 
 * que contiene los datos.
 *
 * En este caso se decidio hacer un padding.
 */


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>

struct PAD
{
  double *density;
};
struct PAD pad;


void fourier(int n, double L)
{
  FILE *input = NULL;
  FILE *output = NULL;
  FILE *output2 = NULL;
  FILE *salida= NULL;

  int i, N;
  fftw_complex *in=NULL, *dft=NULL, *idft=NULL;
  fftw_plan my_plan, my_plan2;
  
  
  N = n*n;

  pad.density = (double*) malloc((size_t)N* sizeof(double));
  
 
  in  = (fftw_complex *) malloc((size_t)N* sizeof(fftw_complex));
  dft  = (fftw_complex *) malloc((size_t)N* sizeof(fftw_complex));
  idft  = (fftw_complex *) malloc((size_t)N*sizeof(fftw_complex));
  
  printf("THE PADDING FOR SURFACE DENSITY FIELD IS APPLIED\n");

  int top1 = floor ((N*0.5) - n);

  int top2 = floor (N*0.5);

  for(i=0; i<top1; i++) 
    {
         
      pad.density[i] = 0.0;
     
    }

  for(i=top1; i<top2; i++) //Padding
    {
         
      pad.density[i] = celda[i-top1].den_super;
      
    }

  for(i=top2; i<N; i++)
    {
      pad.density[i] = 0.0;
    }

  salida = fopen("pbc_data.dat","w");

  for(i=0; i<N; i++) 
    {
      in[i][0] = pad.density[i];
      in[i][1] = 0.0;

      fprintf(salida,"%lf\n", in[i][0]);
    }

    
  my_plan = fftw_plan_dft_2d(n, n, in, dft, FFTW_FORWARD, FFTW_ESTIMATE); //2D DFT:Discrete Fourier Transform

  fftw_execute(my_plan);

  output = fopen("DFT_2D.dat","w");

  for(i=0; i<N; i++)
    {
      fprintf(output,"%lf %lf\n", dft[i][0] , dft[i][1]);
      
    }

  fclose(output);
  
  printf("THE STATE OF DFT IS: SUCESS\n");
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  my_plan2 = fftw_plan_dft_2d(n, n, dft, idft, FFTW_BACKWARD, FFTW_ESTIMATE);//2D IDFT:Inverse Discrete Fourier Transform

  fftw_execute(my_plan2);

  output2 = fopen("IDFT_2D.dat","w");

  for(i=0; i<N; i++)
    {
      fprintf(output2,"%lf %lf\n", idft[i][0]/N, idft[i][1]/N);
    }

  fclose(output2);
  
  fftw_destroy_plan(my_plan2);
    
  printf("THE STATE OF IDFT IS: SUCESS\n");
  
  fftw_free(in);
  fftw_free(dft);
  fftw_free(idft);
  
}
