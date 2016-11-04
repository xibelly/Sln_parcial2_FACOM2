/*
Xibelly Eliseth Mosquera Escobar

Solucion parcial 2  COMPUTACION CIENTIFICA AVANZADA

1) Generar una distribucion de particulas, por medio de una dsitribucion gaussiana. A su vez se asigna la velocidad y la temperatura como una funcion de la distancia de las particulas al centro de la distribucion.

2) Se crea un grid rectangular con un tamaño Lbox=101 la distancia maxima de las particulas al centro. Se carga el grid, identificando las particulas que pertenecen a una celda dada y se calculan el siguiente conjunto de propiedades para la malla:

a)la masa total contenida en cada celda
b)la densidad superficial de masa
c)el vector velocidad
d)posicion del centro de masa 
e)valor de la dispersion de velocidades y cuartiles

3)

4) Se calcula el campo de densidad superficial y se le calcula tanto la primera como la segunda derivada -gradiente-.

5) A dicho campo de densidad se le calcula la 2D DFT -Discrete Fourier Transform- y la 2D IDFT -Inverse Discrete Fourier Transform-
 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_double.h>

#include <gsl/gsl_statistics.h>
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
FILE *out2=NULL;
FILE *out3=NULL;
FILE *out4=NULL;

//Estructuras
struct particles   
{
  double *x;           /*coordenadas particulas*/
  double *y;
  
  double *vr;          /*vector velocidad*/
  

  double *new_posx;    /*coordenadas particulas re-escaladas*/
  double *new_posy;
  double *new_r;

  double new_xc;         /*Coordenadas del nuevo centro*/
  double new_yc;
}; 
struct particles part; 

struct cells    
{

  int cell_ID;           /*ID de la celda*/
  int Np;                /*# particulas en la celda*/
  
  double x;             /*coordenadas de la celda*/
  double y;

  double xc;            /*coordenadas del centro de masa*/
  double yc;
  
  double vc;            /* vector velocidad del centro de masa*/
  double vr;            /* vector velocidad radial*/

  double masa;          /*masa contenida en la celda*/
  double den_super;     /*densidad superficial de la celda*/

  double r;             /*distancia al centro de masa de la celda*/
  double sigma;         /*dispersion de velocidades*/

  double cuartil;       /* cuartil de la celda*/
  
}; 
struct cells *celda; 



//Funciones
int random_number()
{
  int i;
  const gsl_rng_type * T;
  gsl_rng * r;

  FILE *write=NULL;
  
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  printf("WRITING FLIE: positions.dat\n");
  
  write = fopen("positions.dat","w");

  if(write == NULL)
    printf("THE FILE: positions.dat CAN NOT BE OPENED\n");
  

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
  double *vc_celda, *sigma_cell;
  double CellSize;
  double x, y;
  double xc, yc;
  double M;
  double *r;
  double Temp;
  double Lbox;
 
  
  FILE *file = NULL, *file1=NULL;
  FILE *datospos = NULL, *sigma= NULL;

  double *dist, distmin, distmax;
  double xi, yi, xj, yj;
  double X_centro, Y_centro;
  
  double Xmax;
  double Ymax;
  double Zmax;
  double Xmin;
  double Ymin;
  double Zmin;

  double median, upperq, lowerq;
  

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

  part.vr = (double *) malloc(N_part *sizeof(double)); /*Velocities*/
  
  
  part.new_posx = (double *) malloc(N_part *sizeof(double));  /*New sistem of Particles*/

  part.new_posy = (double *) malloc(N_part *sizeof(double));  

  part.new_r = (double *) malloc(N_part *sizeof(double));

    
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

  
  /*Con las coordenadas del centro de masa - Se calcula la distancia maxima */

  /*Re-escalamos el sistema y calculamos Temperatura y Velocidad como una funcion de la distancia al nuevo centro-*/
  
  distmin = 10000*N_part;

  distmax = 0.0;  
 

  out1= fopen("velocities.dat","w");

  if(out1 == NULL)
    printf("THE FILE: velocities.dat  CAN NOT BE OPENED\n");


  out2= fopen("new_positions.dat","w");

  if(out2 == NULL)
    printf("THE FILE: new_positions.dat  CAN NOT BE OPENED\n");
  
  
  for (i = 0; i < N_part; i++)
    {

      xj = part.x[i] ;
      yj = part.y[i] ;

      
      r[i] = distancia(xc, xj, yc, yj);  //distancia de las particulas al centro
      
      
      if(distmin>r[i])//se calcula la distancia min y max
	{
	  distmin=r[i];
	}
      
      
      if(r[i]>distmax)
	{
	  distmax=r[i];
	}
      

      //re-escalamos la distribucion
      
      part.new_posx[i] = xj + distmax;
      part.new_posy[i] = yj + distmax;


      //calculamos el nuevo centro de la distribucion
      
      part.new_xc += ( part.new_posx[i] * m ) / M;   //Coordenada X del nuevo centro 
     
      part.new_yc += ( part.new_posy[i] * m ) / M;   //Coordenada Y del nuevo centro 
     

      part.new_r[i] = distancia(part.new_xc, part.new_posx[i], part.new_yc,  part.new_posy[i]);  //distancia de las particulas al nuevo centro
      

      Temp = C * part.new_r[i] ;           //Temperatura como una funcion de la distancia

      part.vr[i] = (part.new_r[i]) /t;    //velocidad como una funcion de la distancia
     

      fprintf (out1,"%lf %lf %lf\n", r[i], part.vr[i], Temp);
      fprintf (out2, "%lf %lf\n", part.new_posx[i] , part.new_posy[i] ); 
    }
   
  fclose(out1);
  fclose(out2);

  printf("THE CENTER OF MASS COORDINATES ARE:\n");
  printf("xc :%lf yc: %lf\n",part.new_xc, part.new_yc );


  /*-Lbox -> 101% dist_maxima -*/

  printf("Valor min de dist %.9lf\n",distmin);
  printf("Valor max de dist: %.9lf\n",distmax);
 
    
  Lbox = 1.01 * 2 * distmax;

  printf("THE BOX SIZE IS:\n");
  printf("%lf\n", Lbox);
  
  //-------------------------------------------------------------------------CONSTRUCCION DEL GRID                                                                      
  

  CellSize    = Lbox/(1.0 * Ncell);

  NtotalCells = Ncell*Ncell;

  
  celda = malloc((size_t) NtotalCells*sizeof(struct cells));
                                                                                                 vc_celda = (double *) malloc((size_t) NtotalCells *sizeof(double));
												 sigma_cell = (double *) malloc((size_t) NtotalCells *sizeof(double));
												 

  out3=fopen("malla.dat","w");

  if(out3 == NULL)
    printf("THE FILE: malla.dat  CAN NOT BE OPENED\n");
  

  for(j=0; j<Ncell; j++)
    
    {    
      
      for(i=0; i<Ncell; i++)

	{   
	  cell_ID=Ncell * j + i;       //identificador de celda
	  celda[cell_ID].cell_ID = cell_ID;
	  
	  //coordenadas del centro de cada celda                                                                                                                    
	  celda[cell_ID].x = CellSize * 0.5 +i * CellSize;  //coordenada en el eje x
	  celda[cell_ID].y = CellSize * 0.5 +j * CellSize;   //coordenada en el eje y

	  //numero de particulas en la celda
	  celda[cell_ID].Np = 0;
	  
	  fprintf(out3,"%lf %lf\n",celda[cell_ID].x, celda[cell_ID].y);
	  
	}
        
    }
  
  fclose(out3);

  //---------------------------------------------------------------ASIGNACION DE PARTICULAS A CADA CELDA                                                              

  /*Calculamos las propiedades para la malla -> masa total, 
    densidad superficial, posicion centro de masa, vector velocidad,
    dispersion de velocidadesy cuartil.*/


  out4=fopen("properties_grid.dat","w");

  if(out4 == NULL)
    printf("THE FILE: malla.dat  CAN NOT BE OPENED\n");
  
  
  for(i=0; i<NtotalCells; i++)
    {

      Xmax = celda[i].x + CellSize*0.5;  //Coordenadas extremos de la celda
      Xmin = celda[i].x - CellSize*0.5;

      Ymax = celda[i].y + CellSize*0.5;
      Ymin = celda[i].y - CellSize*0.5;

      
      for(j=0; j<N_part; j++)
        {

          if( (part.new_posx[j]>=Xmin) && (part.new_posx[j]<Xmax) )
            if( (part.new_posy[j]>=Ymin) && (part.new_posy[j]<Ymax) ) 
              {
	
		celda[i].Np = celda[i].Np + 1;

		celda[i].masa = celda[i].masa + 1;

		celda[i].xc += (  celda[i].x * celda[i].masa ) ;

		celda[i].yc += (  celda[i].y * celda[i].masa ) ;

		celda[i].r = distancia(celda[i].xc, celda[i].x, celda[i].yc, celda[i].y );

		celda[i].vr = celda[i].r / t;
		
		celda[i].vc +=  (  celda[i].vr * celda[i].masa ) ;

		vc_celda[i] = celda[i].vc;
	
		celda[i].sigma = gsl_stats_sd(vc_celda, 1, NtotalCells);

		sigma_cell[i] = celda[i].sigma;
	      } 
	}                                

      fprintf(out4,"%d %d %lf %lf %lf %lf %lf %lf %lf\n", celda[i].cell_ID, celda[i].Np, celda[i].masa, celda[i].xc, celda[i].yc, celda[i].r, celda[i].vr, celda[i].vc, celda[i].sigma);
      
    }
  fclose(out4);
  
  /*Calculo de los cuantiles -sigma-*/

  /*
  
  gsl_sort(sigma_cell, 1, NtotalCells);

  upperq = gsl_stats_quantile_from_sorted_data(sigma_cell, 1, NtotalCells, 0.5);

  lowerq = gsl_stats_quantile_from_sorted_data(sigma, 1, NtotalCells, 0.5);

  printf("The upper quartile is %g\n", upperq);
  
  printf("The lower quartile is %g\n", lowerq);

  */
  
  return 0;

}

