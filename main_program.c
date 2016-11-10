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

3) Se interpolan las cantidades anteriores considerando x-> xc y y-> variando  
entre (0,Lbox) -se generan # aleatorios en este rango-.

4) Se calcula el campo de densidad superficial y se le calcula tanto la primera como la segunda derivada -gradiente-.

5) A dicho campo de densidad se le calcula la 2D DFT -Discrete Fourier Transform- y la 2D IDFT -Inverse Discrete Fourier Transform-
 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_double.h>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <distancia.h>

#define m 1.0
#define t 1.0
#define C 1.0

//Variables Globales
int N_part, Ncell, NtotalCells;
double mu_x, mu_y;
double L;
double value;
double *xx, *yy;
double CellSize;  

FILE *out1=NULL;
FILE *out2=NULL;
FILE *out3=NULL;
FILE *out4=NULL;
FILE *out5=NULL;
FILE *out6=NULL;
FILE *out7=NULL;
FILE *out8=NULL;
FILE *out9=NULL;
FILE *extra=NULL;
FILE *extra2=NULL;



//Estructuras
struct particles   
{
  double *x;           /*coordenadas particulas*/
  double *y;
  
  double *vrx;          /*vector velocidad*/
  double *vry;

  double *new_posx;    /*coordenadas particulas re-escaladas*/
  double *new_posy;
  double *new_r;

  double new_xc;         /*Coordenadas del nuevo centro*/
  double new_yc;
}; 
struct particles part; 

struct cells    
{
  int N_cells;
  
  int cell_ID;           /*ID de la celda*/
  int Np;                /*# particulas en la celda*/

  double *part_posx;     /*Coordenas particulas en la celda*/
  double *part_posy;
  
  double *part_vx;      /*Velocidades particulas en la celda*/
  double *part_vy;  
  
  double *mag_v;        /*Magnitud de la velocidades de las particulas*/

  double x;             /*coordenadas de la celda*/
  double y, new_y;

  double xc;            /*coordenadas del centro de masa*/
  double yc;
  
  double vcx;          /* vector velocidad del centro de masa*/
  double vcy;          /* vector velocidad radial*/
     

  double masa;          /*masa contenida en la celda*/
  double den_super;     /*densidad superficial de la celda*/

  double r, new_r;                 /*distancia al centro de masa de la celda*/
  double sigma, new_sigma;         /*dispersion de velocidades*/

  double cuartil;                  /* cuartil de la celda*/
  
}; 
struct cells *celda; 

struct DATA
{
  double *result;      /*Numeros aleatorios generados entre 0 y Lbox*/
  size_t n;
  double *F;
  
};
struct DATA datos;

struct data
{
  size_t n;
  double Lbox;
};


//Funciones

#include "distribution.c"

#include "random.c"

#include "interpolador.c"

#include "derivadas.c"

#include "fourier.c"



/*----------------------------------------------------------------PROGRAMA PRINCIPAL*/
 
int main(int argc, char *argv[])
{

  int i, j, cell_ID;
  double xmin, xmax, ymin, ymax, zmin, zmax, BoxSize;
  double *vc_celda, *sigma_cell;
  
  double x, y;
  double xc, yc;
  double M;
  double *r;
  double Temp;
  double Lbox;

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
  double skewness;
  double kurtosis;
  double *X,*Y,*F,*Xo;
  double xo, yo;

    

  //Carga de Parametros//

  L = atof(argv[1]);        //Longitud de la caja

  N_part = atoi(argv[2]);   //# de particlas
             
  Ncell = atoi(argv[3]);    //# de celdas

  mu_x = atof(argv[4]);     //varianza para el eje x

  mu_y = atof(argv[5]);     //varaianza para el eje y

  if(argc != 6)
    {
      printf("ERROR--> use as:\n");
      printf(" USE:<Lbox>  <Nparticles> <Ncells>  <mu_x>  <mu_y>\n");
      exit(0);
    }

  printf("%lf %d %d %lf %lf\n", L, N_part, Ncell, mu_x, mu_y );

  NtotalCells = Ncell*Ncell;

   
  part.x = (double *) malloc(N_part *sizeof(double));  /*Particles*/

  part.y = (double *) malloc(N_part *sizeof(double));

  part.vrx = (double *) malloc(N_part *sizeof(double)); /*Velocities*/
  
  part.vry = (double *) malloc(N_part *sizeof(double));

  
  part.new_posx = (double *) malloc(N_part *sizeof(double));  /*New sistem of Particles*/
  part.new_posy = (double *) malloc(N_part *sizeof(double));  

  part.new_r = (double *) malloc(N_part *sizeof(double));

    
  r = (double *) malloc(N_part *sizeof(double));


  celda = malloc((size_t) NtotalCells*sizeof(struct cells));         /*Cells*/
  vc_celda = (double *) malloc((size_t) NtotalCells *sizeof(double));
  sigma_cell = (double *) malloc((size_t) NtotalCells *sizeof(double));

  
  for(i=0; i<NtotalCells; i++)
    {
      celda[i].part_posx = (double *) malloc(N_part *sizeof(double));  /*Particles coordinates into the cell*/
      celda[i].part_posy = (double *) malloc(N_part *sizeof(double));
      
      celda[i].part_vx = (double *) malloc(N_part *sizeof(double));  /*Particles velocities into the cell*/
      celda[i].part_vy = (double *) malloc(N_part *sizeof(double));

      celda[i].mag_v = (double *) malloc(N_part *sizeof(double));

    }

   double * mag_v = (double *) malloc(N_part *sizeof(double));

  datos.result = (double *) malloc((size_t) NtotalCells *sizeof(double)); /*Random number generation*/

  datos.F = (double *) malloc((size_t) NtotalCells *sizeof(double)); 


  X  = (double *) malloc(NtotalCells *sizeof(double)); /*Arrays Interpolation*/
  Xo  = (double *) malloc(NtotalCells *sizeof(double));  
  Y  = (double *) malloc(NtotalCells *sizeof(double));
  F  = (double *) malloc(NtotalCells *sizeof(double));
  
  yy = (double *) malloc(NtotalCells *sizeof(double));
  xx = (double *) malloc(NtotalCells *sizeof(double));

  

  


  //-------------------------------------------------------------------------GENERACION DE PARTICULAS

  random_distribution(N_part, mu_x, mu_y); //-> programa : distribution.c 

  printf("THE PARTICLE DISTRIBUTION IS CREATED\n"); 

  //--------------------------------------------------------------------------CALCULO CENTRO DISTRIBUCION

  /*coordenadas centro de masa*/

  M = N_part * m ;  // -> la masa total 

  for (i = 0; i < N_part; i++)
    {
      
     xc += ( part.x[i] * m ) / M;   //Coordenada X 
     
     yc += ( part.y[i] * m ) / M;   //Coordenada Y
            
    }

  
  /*Con las coordenadas del centro de masa - Se calcula la distancia maxima */

  /*Re-escalamos el sistema y calculamos Temperatura y el vector Velocidad como una funcion
    de la distancia al nuevo centro-*/
  
  distmin = 100000*N_part;

  distmax = 0.0;  
 

  out1= fopen("velocities.dat","w");

  if(out1 == NULL)
    printf("THE FILE: velocities.dat  CAN NOT BE OPENED\n");


  out2= fopen("new_positions.dat","w");

  if(out2 == NULL)
    printf("THE FILE: new_positions.dat  CAN NOT BE OPENED\n");


  
  printf("WRITING FLIE: new_positions.dat\n");

  printf("WRITING FLIE: velocities.dat\n");


  
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
      
      
      
      part.new_posx[i] = part.x[i] + distmax ;
      part.new_posy[i] = part.y[i] + distmax ; 
      
      
      //calculamos el nuevo centro de la distribucion
      
      part.new_xc += ( part.new_posx[i] * m ) / M;   //Coordenada X del nuevo centro 
      
      part.new_yc += ( part.new_posy[i] * m ) / M;   //Coordenada Y del nuevo centro 
      
      
      part.new_r[i] = distancia(part.new_xc, part.new_posx[i], part.new_yc,  part.new_posy[i]);  //distancia de las particulas al nuevo centro
      
      
      Temp = C * part.new_r[i] ;           //Temperatura como una funcion de la distancia
      
      part.vrx[i] = (part.new_r[i] * part.new_posx[i]) /t;    //velocidad como una funcion de la distancia
      part.vry[i] = (part.new_r[i] * part.new_posy[i]) /t;
      

      fprintf (out1,"%lf %lf %lf %lf\n", r[i], part.vrx[i], part.vry[i], Temp);
      fprintf (out2, "%lf %lf\n", part.new_posx[i] , part.new_posy[i] ); 
    }
   
  fclose(out1);
  fclose(out2);

  printf("THE CENTER OF MASS COORDINATES ARE:\n");
  printf("xc :%lf yc: %lf\n",part.new_xc, part.new_yc );


  /*-Lbox -> 101% dist_maxima -*/

  printf("MINIMUM VALUE OF DISTANCE: %.9lf\n",distmin);
  printf("MAXIMUM VALUE OF DISTANCE: %.9lf\n",distmax);
 
    
  Lbox = 1.01 * 2.0 * distmax;

  printf("THE BOX SIZE IS:\n");
  printf("%lf\n", Lbox);
  
  //-------------------------------------------------------------------------CONSTRUCCION DEL GRID                                                                      
  

  CellSize    = Lbox/(1.0 * Ncell);

    
    
  out3=fopen("grid.dat","w");

  if(out3 == NULL)
    printf("THE FILE: grid.dat  CAN NOT BE OPENED\n");
  
  printf("WRITING FLIE: grid.dat\n");

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

  printf("THE GRID IS CREATED\n"); 

  //---------------------------------------------------------------ASIGNACION DE PARTICULAS A CADA CELDA                                                              

  /*Calculamos las propiedades para la malla -> masa total, 
    densidad superficial, posicion centro de masa, vector velocidad,
    dispersion de velocidades y cuartil, de cada una de las celdas que la
    conforman.*/

  double sumaxc, sumayc, sumavcx, sumavcy;

  out4=fopen("properties_grid.dat","w");

  if(out4 == NULL)
    printf("THE FILE: properties_grid.dat  CAN NOT BE OPENED\n");
  
  
  printf("WRITING FLIE: properties_grid.dat\n");


  extra=fopen("particles_in_cell.dat","w");

  if(extra == NULL)
    printf("THE FILE: particles_in_cell.dat  CAN NOT BE OPENED\n");
  
  
  printf("WRITING FLIE: particles_in_cell.dat\n");


  extra2=fopen("cuartiles_curtosis_skewness.dat","w");

  if(extra2 == NULL)
    printf("THE FILE: cuartiles_curtosis_skewness.dat  CAN NOT BE OPENED\n");
  
  
  printf("WRITING FLIE: cuartiles_curtosis_skewness.dat\n");

  int suma = 0;

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

		/*Particulas dentro de cada celda*/
		
		celda[i].part_posx[j] = part.new_posx[j];

		celda[i].part_posy[j] = part.new_posy[j];

		celda[i].part_vx[j] = part.vrx[j];

		celda[i].part_vy[j] = part.vry[j];

		/*Calculo del # total de particulas, masa y densidad superficial*/

		celda[i].Np = celda[i].Np + 1;

		celda[i].masa = celda[i].masa + 1;

		celda[i].den_super = celda[i].masa / (CellSize * CellSize);

		
		 /*Calculo de las coordenadas centro de masa y distancia de este*/

		sumaxc +=   celda[i].part_posx[j] ;

		sumayc +=   celda[i].part_posy[j]  ;

		celda[i].r = distancia(celda[i].xc, part.new_posx[j], celda[i].yc, part.new_posy[j] );

		 /*Calculo del vector velocidad del centro de masa*/

			
		sumavcx +=   celda[i].part_vx[j]   ;

		sumavcy +=   celda[i].part_vy[j]  ;

		/*Magnitud de la velocidad de las particulas*/


		celda[i].mag_v[j] = sqrt( (celda[i].part_vx[j]*celda[i].part_vx[j]) + (celda[i].part_vy[j]*celda[i].part_vx[j])   );
		
		
		//mag_v = (double *) realloc(mag_v, (size_t) suma*sizeof(double));		
		mag_v[j] = celda[i].mag_v[j];
			
		/*Calculo de la dispersion de velocidades */
      
		celda[i].sigma = gsl_stats_sd(celda[i].mag_v, 1, N_part);

		 /*Calculo de curtosis-skewness */

		skewness = gsl_stats_skew(celda[i].mag_v, 1, N_part);
		
		kurtosis = gsl_stats_kurtosis(celda[i].mag_v, 1, N_part);
      

		 
		fprintf(extra,"%d %lf %lf %lf %lf \n",celda[i].cell_ID, celda[i].part_posx[j], celda[i].part_posy[j], celda[i].part_vx[j], celda[i].part_vy[j]);

		
		
	      } 
	  suma = celda[i].Np + 1;
	}                  

      celda[i].xc = sumaxc /celda[i].masa;
	  
      celda[i].yc = sumayc /celda[i].masa;
      
      celda[i].vcx = sumavcx /celda[i].masa; 
      
      celda[i].vcy = sumavcy /celda[i].masa; 

      
      if(celda[i].masa == 0.0)
	{
      
	  celda[i].xc = 0.0;
	  
	  celda[i].yc = 0.0;
	  
	  celda[i].vcx = 0.0; 
	  
	  celda[i].vcy = 0.0; 

	}
	  
      /*Calculo de cuartiles */

      gsl_sort(mag_v, 1, suma);  //se ordena en forma creciente
      
      median = gsl_stats_median_from_sorted_data(mag_v, 1, suma);
      
      upperq = gsl_stats_quantile_from_sorted_data(mag_v, 1, suma, 0.75);
      
      lowerq = gsl_stats_quantile_from_sorted_data(mag_v, 1, suma, 0.25);
     
      
      fprintf(out4,"%d %d %lf %lf %lf %lf %lf %lf %lf\n",
	      celda[i].cell_ID, celda[i].Np, celda[i].masa, celda[i].den_super, 
	      celda[i].xc, celda[i].yc, 
	      celda[i].vcx, celda[i].vcy, celda[i].sigma);
      

      fprintf(extra2,"%d %d %lf %lf %lf %lf %lf\n",celda[i].cell_ID, celda[i].Np,median, upperq, lowerq, skewness, kurtosis);

    }
  fclose(extra);
  fclose(out4);
  fclose(extra2);

  
  
  printf("THE GRID IS LOADED\n"); 
  
  //--------------------------------------------------------------------------CALCULO INTERPOLACION CANTIDADES ANTERIORES

  /*Interpolamos del conjunto anterior de propiedades
    de cada celda.

    Para ello calculamos de nuevo dichas propiedades
    dejando fijo    x = xcentro_celda
                    y = variando de (0, Lbox) -> generamos # aleatorios en
			       este rango.
   
			       
  */

  random_gfsr4(Ncell); //Genracion de # alatorios

  for(j=0; j<Ncell; j++)
    {    
      Xo[cell_ID] = Lbox * datos.result[cell_ID];
    }


   //////////////////////////////////Interpolacion Masa///////////////////////////////////
  
  
  out5=fopen("inter_mass.dat","a");

  if(out5 == NULL)
    printf("THE FILE: inter_mass.dat  CAN NOT BE OPENED\n");
  
  
  printf("WRITING FLIE: inter_mass.dat\n");

  for(j=0; j<Ncell; j++)
    {    
      
      for(i=0; i<Ncell; i++)
	{   
	  cell_ID=Ncell * j + i;	  
	  
	  xx[cell_ID] =  celda[cell_ID].yc; //eje_y celda
	  
	  yy[cell_ID] =  celda[cell_ID].masa;  //funcion a interpolar
       	  
	    
	}
      gsl_sort(xx, 1, NtotalCells); //ordena en forma creciente yc_celda para cada columna de celdas
      
      for (xo = Xo[0]; xo < Xo[NtotalCells-1]; xo = xo + 10)
	{
	  
	  yo = interpolador_akima_per(xo, Ncell*Ncell, xx, yy);
	  fprintf (out5,"%g %g\n", xo, yo);
	}
      
    }
  fclose(out5);
  
  printf("STATE OF MASS INTERPOLATION IS: SUCESS\n"); 

  
  //////////////////////////////////Interpolacion Densidad superficial///////////////////////////////////

  out6=fopen("inter_density.dat","w");
  
  if(out6 == NULL)
    printf("THE FILE: inter_density.dat  CAN NOT BE OPENED\n");
  
  
  printf("WRITING FLIE: inter_density.dat\n");

   
 
  for(j=0; j<Ncell; j++)
    {    
      
      for(i=0; i<Ncell; i++)
	{   
	  cell_ID=Ncell * j + i;	  
	  
	  xx[cell_ID] =  celda[cell_ID].yc; //eje_y celda
	  
	  yy[cell_ID] =  celda[cell_ID].den_super;  //funcion a interpolar
	  
	  
	}
      gsl_sort(xx, 1, NtotalCells); //ordena en forma creciente yc_celda para cada columna de celdas
      
      for (xo = Xo[0]; xo < Xo[Ncell-1]; xo = xo + 10)
	{
	  
	  yo = interpolador_akima_per(xo, Ncell, xx, yy);
	  fprintf (out6,"%g %g\n", xo, yo);
	  
	}
      
    }
  fclose(out6);
  
  printf("STATE OF SURFACE-DENSITY INTERPOLATION IS: SUCESS\n");



  
  
  //--------------------------------------------------------------------------CALCULO DERIVADAS CAMPO DE DENSIDAD

  /*Derivadas del campo densidad superficial. Para obtener el gradiente de
    dicho campo, interpolamos talque tengamos una funcion que podamos derivar en los
    puntos (x,y) de cada celda.*/
  
  

  deriva(Ncell);
  

  //--------------------------------------------------------------------------CALCULO FOURIER TRANSFORMS

  /*Calculamos la 2D DFT -Discrete Fourier Transform- y la 2D IDFT -Inverse Discrete Fourier Transform-
    
   */

  
  fourier(NtotalCells, Lbox);

  return 0;

}


