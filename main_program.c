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

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
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
int N_part, Ncell, NtotalCells;
double mu_x, mu_y;
double L;
double value;
double *xx, *yy;
  

FILE *out1=NULL;
FILE *out2=NULL;
FILE *out3=NULL;
FILE *out4=NULL;
FILE *out5=NULL;
FILE *out6=NULL;
FILE *out7=NULL;
FILE *out8=NULL;
FILE *out9=NULL;

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
  int N_cells;
  
  int cell_ID;           /*ID de la celda*/
  int Np;                /*# particulas en la celda*/
  
  double x;             /*coordenadas de la celda*/
  double y, new_y;

  double xc;            /*coordenadas del centro de masa*/
  double yc;
  
  double vc, new_vc;    /* vector velocidad del centro de masa*/
  double vr, new_vr;    /* vector velocidad radial*/

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
};
struct DATA datos;

struct data
{
  size_t n;
};


//Funciones

#include "distribution.c"

#include "random.c"

#include "interpolador.c"

#include "derivadas.c"


/*----------------------------------------------------------------PROGRAMA PRINCIPAL*/
 
int main(int argc, char *argv[])
{

  int i, j, cell_ID;
  double xmin, xmax, ymin, ymax, zmin, zmax, BoxSize;
  double *vc_celda, *sigma_cell;
  double CellSize;
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
  double *X,*Y,*F;
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

  part.vr = (double *) malloc(N_part *sizeof(double)); /*Velocities*/
  
  
  part.new_posx = (double *) malloc(N_part *sizeof(double));  /*New sistem of Particles*/
  part.new_posy = (double *) malloc(N_part *sizeof(double));  

  part.new_r = (double *) malloc(N_part *sizeof(double));

    
  r = (double *) malloc(N_part *sizeof(double));


  celda = malloc((size_t) NtotalCells*sizeof(struct cells));         /*Cells*/
  vc_celda = (double *) malloc((size_t) NtotalCells *sizeof(double));
  sigma_cell = (double *) malloc((size_t) NtotalCells *sizeof(double));
  

  datos.result = (double *) malloc((size_t) NtotalCells *sizeof(double)); /*Random number generation*/


  X  = (double *) malloc(NtotalCells *sizeof(double)); /*Arrays Interpolation*/
  Y  = (double *) malloc(NtotalCells *sizeof(double));
  F  = (double *) malloc(NtotalCells *sizeof(double));
  
  yy = (double *) malloc(NtotalCells *sizeof(double));
  xx = (double *) malloc(NtotalCells *sizeof(double));

  
  size_t p[NtotalCells]; /*Sort of cell porperties*/

  int k = NtotalCells;



  //-------------------------------------------------------------------------GENERACION DE PARTICULAS

  random_distribution(N_part, mu_x, mu_y);

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

  /*Re-escalamos el sistema y calculamos Temperatura y Velocidad como una funcion de la distancia al nuevo centro-*/
  
  distmin = 10000*N_part;

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

  out4=fopen("properties_grid.dat","w");

  if(out4 == NULL)
    printf("THE FILE: properties_grid.dat  CAN NOT BE OPENED\n");
  
  
  printf("WRITING FLIE: properties_grid.dat\n");

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

		 /*Calculo del # total de particulas, masa y densidad superficial*/

		celda[i].Np = celda[i].Np + 1;

		celda[i].masa = celda[i].masa + 1;

		celda[i].den_super = celda[i].masa / (CellSize * CellSize);

		
		 /*Calculo de las coordenadas centro de masa y distancia de este*/

		celda[i].xc += (  celda[i].x * celda[i].masa ) ;

		celda[i].yc += (  celda[i].y * celda[i].masa ) ;

		celda[i].r = distancia(celda[i].xc, part.new_posx[j], celda[i].yc, part.new_posy[j] );

		 /*Calculo del vector velocidad y velocidad del centro de masa*/

		celda[i].vr = celda[i].r / t;
		
		celda[i].vc +=  (  celda[i].vr * celda[i].masa ) ;

		vc_celda[i] = celda[i].vc;


		/*Calculo de la dispersion de velocidades -velocidad del centro de masa-*/

		celda[i].sigma = gsl_stats_sd(vc_celda, 1, NtotalCells);

		sigma_cell[i] = celda[i].sigma;


		 /*Calculo de los cuantiles -sigma-*/

		gsl_sort(sigma_cell, 1, NtotalCells);


		median = gsl_stats_median_from_sorted_data(sigma_cell, 1, NtotalCells);
		
		upperq = gsl_stats_quantile_from_sorted_data(sigma_cell, 1, NtotalCells, 0.75);
  
		lowerq = gsl_stats_quantile_from_sorted_data(sigma_cell, 1, NtotalCells, 0.25);

	      } 
	}                        
      
     
      fprintf(out4,"%d %d %lf %lf %lf %lf %lf %lf %g %g %g\n",
	      celda[i].cell_ID, celda[i].Np, celda[i].masa, celda[i].den_super, 
	      celda[i].xc, celda[i].yc, celda[i].vc, 
	      celda[i].sigma, median, upperq, lowerq);
      
    }
  fclose(out4);

  printf("THE GRID IS LOADED\n"); 
  
  //--------------------------------------------------------------------------CALCULO INTERPOLACION CANTIDADES ANTERIORES

  /*Interpolamos del conjunto anterior de propiedades
    de cada celda.

    Para ello calculamos de nuevo dichas propiedades
    dejando fijo    x = xcentro_celda
                    y = variando de (0, Lbox) -> generamos # aleatorios en
			       este rango.
   
			       
  */



  random_gfsr4(NtotalCells); //Genracion de # alatorios

  
  for(i=0; i<NtotalCells; i++)  //f(xc,y)
    {
      X[i] = celda[i].xc;
      
      Y[i] = Lbox * datos.result[i];
      
    }

  //////////////////////////////////Interpolacion Masa///////////////////////////////////
  
  
  out5=fopen("inter_mass.dat","w");

  if(out5 == NULL)
    printf("THE FILE: inter_mass.dat  CAN NOT BE OPENED\n");
  
  
  printf("WRITING FLIE: inter_mass.dat\n");

  
   for(i=0; i<NtotalCells; i++)
    {
      F[i] = celda[i].masa;  //f(xc,y) -> mass
      
    }
 

  gsl_sort_smallest_index (p, k, Y, 1, NtotalCells);//Para asegurar la interpolacion se ordenan los datos de forma creciente

   //INTERPOLACION AKIMA PERIODICA

  for(i=0; i<NtotalCells; i++)
    {
      xx[i] = Y[p[i]];

      yy[i] = F[p[i]];
    }
    
  for(i=0; i<100; i=i+10)
    {
      for (xo = xx[0]; xo < xx[NtotalCells-1]; xo = xo + 10)
	{
	  yo = interpolador_akima_per(xo, NtotalCells, xx, yy);
	  fprintf (out5,"%d %g\n", i, yo);
	}
    }
  
  
  fclose(out5);

  printf("STATE OF MASS INTERPOLATION IS: SUCESS\n"); 

  //////////////////////////////////Interpolacion Densidad superficial///////////////////////////////////
  
  

  out6=fopen("inter_density.dat","w");

  if(out6 == NULL)
    printf("THE FILE: inter_density.dat  CAN NOT BE OPENED\n");
  
  
  printf("WRITING FLIE: inter_density.dat\n");

  
   for(i=0; i<NtotalCells; i++)
    {
      F[i] = celda[i].den_super; //f(xc,y) -> surface density
      
    }
 

  gsl_sort_smallest_index (p, k, Y, 1, NtotalCells);//Para asegurar la interpolacion se ordenan los datos de forma creciente

   //INTERPOLACION AKIMA PERIODICA

  for(i=0; i<NtotalCells; i++)
    {
      xx[i] = Y[p[i]];

      yy[i] = F[p[i]];
    }
    
  
  for(i=0; i<100; i=i+10)
    {

      for (xo = xx[0]; xo < xx[NtotalCells-1]; xo = xo + 10)
	{
	  yo = interpolador_akima_per(xo, NtotalCells, xx, yy);
	  fprintf (out6,"%d %g\n", i, yo);
	}
    }
      
  fclose(out6);

  printf("STATE OF SURFACE-DENSITY INTERPOLATION IS: SUCESS\n");


  //////////////////////////////////Interpolacion Vector velocidad del centro de masa///////////////////////////////////
   
  

  out7=fopen("inter_velocity_center.dat","w");

  if(out7 == NULL)
    printf("THE FILE: inter_velocity_center.dat  CAN NOT BE OPENED\n");
  
  
  printf("WRITING FLIE: inter_velocity_center.dat\n");

  
   for(i=0; i<NtotalCells; i++)
    {
      celda[i].new_y = Y[i];
      
      celda[i].new_r = distancia(celda[i].xc, celda[i].xc, celda[i].yc, celda[i].new_y);


      /*Calculo del vector velocidad y velocidad del centro de masa*/
      
      celda[i].new_vr = celda[i].new_r / t;
      
      celda[i].new_vc +=  (  celda[i].new_vr * celda[i].masa ) ;

      
      F[i] = celda[i].new_vc; //f(xc,y) -> velocity center vector
      
    }



  gsl_sort_smallest_index (p, k, Y, 1, NtotalCells);//Para asegurar la interpolacion se ordenan los datos de forma creciente

   //INTERPOLACION AKIMA PERIODICA

  for(i=0; i<NtotalCells; i++)
    {
      xx[i] = Y[p[i]];

      yy[i] = F[p[i]];
    }
    
  for(i=0; i<100; i=i+10)
    {
      
      for (xo = xx[0]; xo < xx[NtotalCells-1]; xo = xo + 10)
	{
	  yo = interpolador_akima_per(xo, NtotalCells,xx ,yy);
	  fprintf (out7,"%d %g\n", i, yo);
	}
    }
  
  fclose(out7);

  printf("STATE OF VELOCITY-VECTOR INTERPOLATION IS: SUCESS\n");

  
  //////////////////////////////////Interpolacion dispersion de velocidades///////////////////////////////////
   
  

  out8=fopen("inter_velocity_dispersion.dat","w");

  if(out8 == NULL)
    printf("THE FILE: inter_velocity_dispersion.dat  CAN NOT BE OPENED\n");
  
  
  printf("WRITING FLIE: inter_velocity_dispersion.dat\n");

  
   for(i=0; i<NtotalCells; i++)
    {
      vc_celda[i] = celda[i].new_vc;

      celda[i].new_sigma = gsl_stats_sd(vc_celda, 1, NtotalCells);

      F[i] = celda[i].new_sigma; //f(xc,y) -> velocity dispersion
      
    }
 

  gsl_sort_smallest_index (p, k, Y, 1, NtotalCells);//Para asegurar la interpolacion se ordenan los datos de forma creciente

   //INTERPOLACION AKIMA PERIODICA

  for(i=0; i<NtotalCells; i++)
    {
      xx[i] = Y[p[i]];

      yy[i] = F[p[i]];
    }
    
  for(i=0; i<100; i=i+10)
    {
      
      for (xo = xx[0]; xo < xx[NtotalCells-1]; xo = xo + 10)
	{
	  yo = interpolador_akima_per(xo, NtotalCells, xx, yy);
	  fprintf (out8,"%d %g\n", i, yo);
	}
      
    }
  fclose(out8);
  
  printf("STATE OF VELOCITY DISPERSION INTERPOLATION IS: SUCESS\n");

  //--------------------------------------------------------------------------CALCULO DERIVADAS CAMPO DE DENSIDAD

  /*Derivadas del campo densidad superficial. Para obtener el gradiente de
    dicho campo, interpolamos talque tengamos una funcion que podamos derivar en los
    puntos (x,y) de cada celda.*/
  
  const size_t n = NtotalCells;
  gsl_function H;
  double result, result1, result2, abserr1, abserr2;

  struct data d = {n};

  H.function = &deriva; //funcion interpolacion del campo densidad 
  H.params = &d;
 
  
 out9=fopen("derivates_of_desity.dat","w");

  if(out9 == NULL)
    printf("THE FILE: derivates_of_desity.dat  CAN NOT BE OPENED\n");
  
  
  printf("WRITING FLIE: derivates_of_desity.dat\n");
 
  for(i=0; i<NtotalCells; i++)
    {
      gsl_deriv_central(&H, celda[i].x, 1e-8, &result1, &abserr1); //Componente x del gradiente
      
      gsl_deriv_central(&H, celda[i].y, 1e-8, &result2, &abserr2); //Componente y del gradiente
      
      
      fprintf(out9,"%lf %lf %lf %lf\n", celda[i].x, result1, celda[i].y, result2);
      
    }
  
    fclose(out9);

  printf("STATE OF DERIVATES OF SURFACE DENSITY IS: SUCESS\n");
 

  return 0;

}


