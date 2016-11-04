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
#include <math.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_int.h>
#include <gsl/gsl_sort_double.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <distancia.h>

#define m 1.0
#define t 1.0
#define C 1.0

//Variables Globales
int N_part, Ncell;
double mu_x, mu_y;

FILE *out1=NULL;

//Estructuras
struct particles   
{
  double *x;
  double *y;
  
  double *vx;
  double *vy;
}; 
struct particles part; 

//Funciones
int random_number()
{
  int i;
  const gsl_rng_type * T;
  gsl_rng * r;

  FILE *write=NULL;
  
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  
  write = fopen("positions.dat","w");

  for (i = 0; i < N_part; i++)
    {
          
      part.x[i] = gsl_ran_gaussian (r, mu_x); //Gaussian distribution

      part.y[i] = gsl_ran_gaussian (r, mu_y); //Gaussian distribution
      
      fprintf (write,"%lf %lf \n", part.x[i] , part.y[i]);
    }

 
  printf("STATE OF PARTICLES GENERATION IS: SUCESS\n");
  
  gsl_rng_free (r);

  fclose(write);
  
  return 0;
}



int main(int argc, char *argv[])
{

  int i, j, NtotalCells, cell_ID;
  double xmin, xmax, ymin, ymax, zmin, zmax, BoxSize;
  double CellSize, CellVol;
  double x, y;
  double xc, yc;
  double M;
  double *r;
  double Temp;
  double Lbox;
 
  
  FILE *file = NULL, *file1=NULL;
  FILE *datospos = NULL, *sigma= NULL;

  double Xmax;
  double Ymax;
  double Zmax;
  double Xmin;
  double Ymin;
  double Zmin;

  //Carga de Parametros//

  N_part = atoi(argv[1]); 
             
  Ncell = atoi(argv[2]);

  mu_x = atof(argv[3]);

  mu_y = atof(argv[4]);

  if(argc != 5)
    {
      printf("ERROR--> use as:\n");
      printf(" USE:  <Nparticles> <Ncells>  <mu_x>  <mu_y>\n");
      exit(0);
    }

  printf("%d %d %lf %lf\n", N_part, Ncell, mu_x, mu_y );

  part.x = (double *) malloc(N_part *sizeof(double));  /*Particles*/

  part.y = (double *) malloc(N_part *sizeof(double));

  part.vx = (double *) malloc(N_part *sizeof(double)); /*Velocities*/

  part.vy = (double *) malloc(N_part *sizeof(double));

  r = (double *) malloc(N_part *sizeof(double));

  size_t p[N_part];

  int k = N_part;


  //-------------------------------------------------------------------------GENERACION DE PARTICULAS

  random_number();


  //--------------------------------------------------------------------------CALCULO CENTRO DISTRIBUCION

  /*coordenadas centro de masa*/

  M = N_part * m ;  // -> la masa total 

  for (i = 0; i < N_part; i++)
    {
      
     xc += ( part.x[i] * m ) / M;   //Coordenada X 
     
     yc += ( part.y[i] * m ) / M;   //Coordenada Y
            
    }

  printf("THE CENTER OF MASS COORDINATES ARE:\n");
  printf("xc :%lf yc: %lf\n", xc, yc);


  /*distancia al centro de la distribucion - Temperatura - Velocidad-*/

  out1= fopen("velocities.dat","w");
  
  for (i = 0; i < N_part; i++)
    {
          
      r[i] = distancia(xc, part.x[i], yc, part.y[i]);

      Temp = C *r[i];                        //Temperatura como una funcion de la distancia

      part.vx[i] = (r[i] * part.x[i]) /t;    //velocidad como una funcion de la distancia

      part.vy[i] = (r[i] * part.y[i]) /t;

      fprintf (out1,"%lf %lf %lf\n", r[i], part.vx[i], part.vy[i]);
      
    }
      
  fclose(out1);

  /*Se calcula la distancia maxima -Lbox-*/
  
  gsl_sort_largest_index (p, k, r, 1,N_part);//Se ordenan las particulas en forma decreciente

  Lbox = 1.01 * r[p[0]];

  printf("THE BOX SIZE IS:\n");
  printf("%lf\n", Lbox);
  
  //-------------------------------------------------------------------------CONSTRUCCION DEL GRID                                                                      
  /*

  CellSize    = Header.BoxSize/Ncell;
  CellVol     = pow(CellSize,3);
  NtotalCells = Ncell*Ncell*Ncell;



  celda = malloc((size_t) NtotalCells*sizeof(struct cells));
  //asignar tamañao a variable part tipo struct particles                                                                                                               

  file=fopen("centros.dat","w");

  for(k=0; k<Ncell; k++) // z                                                                                                                                           
    {

      z=  CellSize*0.5 +k*CellSize ;

      for(j=0; j<Ncell; j++) //y                                                                                                                                        
        {

          y=  CellSize*0.5 +j*CellSize ;

          for(i=0; i<Ncell; i++) // x                                                                                                                                   
            {

              //coordenada del centro en el eje x                                                                                                                       
              x =   CellSize *0.5 +i*  CellSize ;
	      //identificador de celda                                                                                                                                  
              cell_ID=(k*Ncell+j)*Ncell+i;
              celda[cell_ID].cell_ID = cell_ID;

              //coordenadas del centro de cada celda                                                                                                                    
              celda[cell_ID].x = x;
              celda[cell_ID].y = y;
              celda[cell_ID].z = z;

              celda[cell_ID].Np = 0;

              fprintf(file,"%f\t %f\t %f\n",x,y,z);

            }
        }
    }

  fclose(file);

  //---------------------------------------------------------------ASIGNACION DE PARTICULAS A CADA CELDITA                                                              

  for(i=0; i<NtotPart; i++)
    {
      part[i].v = sqrt(part[i].vel[0]*part[i].vel[0] + part[i].vel[1]*part[i].vel[1] + part[i].vel[2]*part[i].vel[2]);
    }

  Nguess = 100*NtotPart/NtotalCells;
  printf(" *** Using %lfGb of ram for vel in grid\n",Nguess*sizeof(double)/(1024*1024*1024.0));

  sigma=fopen("sigma_pos.dat","w");

  for(k=0; k<NtotalCells; k++)
    {

      velocities = (double *) malloc((size_t) Nguess*sizeof(double));
      if(velocities == NULL)
        {
          printf("do not use such a large amount of memory! (alloc velocity)\n");
          exit(0);
        }

      Xmax = celda[k].x + CellSize*0.5;
      Xmin = celda[k].x - CellSize*0.5;

      Ymax = celda[k].y + CellSize*0.5;
      Ymin = celda[k].y - CellSize*0.5;

      Zmax = celda[k].z + CellSize*0.5;
      Zmin = celda[k].z - CellSize*0.5;
      for(i=0; i<NtotPart; i++)
        {

          if( (part[i].pos[0]>=Xmin) && (part[i].pos[0]<Xmax) )
            if( (part[i].pos[1]>=Ymin) && (part[i].pos[1]<Ymax) )
              if( (part[i].pos[2]>=Zmin) && (part[i].pos[2]<Zmax) )
                {
                  velocities[celda[k].Np] = 1.0*part[i].v;
                  celda[k].Np = celda[k].Np+1;

                  if(celda[k].Np >= Nguess)
                    {
                      Nguess = 2*Nguess;
                      velocities = (double *) realloc(velocities, (size_t) Nguess*sizeof(double));
                      printf("WARNING: Increasing the value of Nguess...\n");
                      //exit(0);                                                                                                                                        
                    }
		}

        }

      velocities = (double *) realloc(velocities, (size_t) celda[k].Np*sizeof(double));

      //gsl_sort(velocities, 1, celda[k].Np);                                                                                                                           
      //median = gsl_stats_mean_from_sorted_data(v, 1, celda[cell_ID].Np); //va a calcular la dispersion de velocidades en cada celda                                   

      celda[k].sigmaV = gsl_stats_sd(velocities, 1, celda[k].Np);
      if(celda[k].Np <= 1)
        celda[k].sigmaV = 0;

      if((k%100) == 0)
        printf("%10d %10d %16.8e %12.8f %12.8f %12.8f\n",k, celda[k].Np, celda[k].sigmaV, celda[k].x, celda[k].y, celda[k].z);

      fprintf(sigma,"%10d %10d %16.8e %12.8f %12.8f %12.8f\n",k, celda[k].Np, celda[k].sigmaV, celda[k].x, celda[k].y, celda[k].z);

      free(velocities);

    }

  */

 
  return 0;

}

