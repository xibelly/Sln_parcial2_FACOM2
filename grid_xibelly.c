/*
Xibelly Eliseth Mosquera Escobar

Solucion Parcial I COMPUTACION CIENTIFICA AVANZADA
 */


/*
 Parametros ideales:

Lbox = 3000

Nbox = 25

 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <gsl/gsl_sort_int.h>
#include <gsl/gsl_sort_double.h>

#include <distancia.h>

//////VARIABLE GLOBALES//////

int Nparticles, Nboxes;
double Lbox;
double *dist;

///////intento

FILE *pf=NULL;
FILE *pf2=NULL;

struct particle   
{
  double pos[2];
  
  
}; 
struct particle *part; 

struct malla   
{
  double *X, *Y;  
  double *pos_x;
  double *pos_y;
  double *centro_x;
  double *centro_y;
  int *part_cells;
  
}; 
struct malla grid;

struct esquina
{
  double x;
  double y;

};
struct esquina *Esquina;


////////////////////////////////////

///////////LLAMADOS A fUNCIONES///////////

#include"input.c"

#include"output.c"

//////////////////////////////////

int main(int argc, char **argv) 
{
  int n, i, j, l, m;
  int index;
  int Nesquinas;
  int suma;
  
  char *filename;
  double box_delta;
  double delta;

  double *dist_grid;
  double distmin, distmax;
  double X_centro, Y_centro;
  double xi, yi, xj, yj;
  

  printf("%d\n",argc);

  if(argc != 5)
    {
      printf("ERROR--> use as:\n");
      printf("%s Npartilces Lbox Nboxes datafile\n",argv[0]);
      exit(0);  
    }

  /*Carga de parametros*/

  Nparticles = atoi(argv[1]);
  
  Lbox       = atof(argv[2]);   
    
  Nboxes     = atoi(argv[3]);   

  filename   = argv[4];

  printf("%s %d %lf %d %s\n",argv[0], Nparticles, Lbox, Nboxes, filename);

  box_delta = Lbox/(1.0*Nboxes); 

  delta = 10* box_delta; //tamaño de las celda
  
  n = Nboxes*Nboxes;

  Nesquinas = (4*n);

  grid.X = (double*) malloc(n *sizeof(double));  //malla

  grid.Y = (double *) malloc(n *sizeof(double));

  grid.pos_x = (double *) malloc(n *sizeof(double));

  grid.pos_y = (double *) malloc(n *sizeof(double));

  grid.centro_x = (double *) malloc(n *sizeof(double));

  grid.centro_y = (double *) malloc(n *sizeof(double));

  dist_grid = (double *) malloc(Nparticles *sizeof(double));

  grid.part_cells = (int *) malloc(Nparticles *sizeof(int));

  Esquina = (struct esquina *) malloc((size_t) Nesquinas* sizeof(struct esquina));
  


  dist = (double *) malloc(Nparticles *sizeof(double));  //particulas

  part = (struct particle *) malloc(Nparticles *sizeof(struct particle));

  

  size_t p[Nparticles];

  int k = Nparticles;

  
  /*Carga de datos*/
  
  read_file(filename, Nparticles);

  pf2 = fopen("coordenadas2.dat","w");

  for(i=0; i<Nparticles; i++)//lectura
    {
      fprintf(pf2,"%lf %lf\n",part[i].pos[0], part[i].pos[1]);
    }

  
  /* Construir el grid*/

  pf = fopen("malla.dat","w");

  for (l=0; l<Nboxes;l++)
    {
      for (j=0; j<Nboxes; j++)
	{
	  
	  index = Nboxes * l + j;
	  //printf("%d\n",index);
	  grid.X[index]= l * box_delta;
	  grid.Y[index] = j * box_delta;
	  
	  fprintf(pf,"%lf %lf\n",grid.X[index], grid.Y[index]);
	}
	
    }

  fclose(pf);

  
  /* Llenamos el grid*/
  suma = 0;
  
  for (l=0; l<Nboxes;l++)
    {
      for (j=0; j<Nboxes; j++)
	{
	  index = Nboxes * l + j;
	  
	  //printf("%d\n",index); //centro de las celdas

	  //grid.centro_y[index] = distancia(grid.X[index], grid.X[index+1], grid.Y[index], grid.Y[index+1]);

	  //printf("%lf\n", grid.centro_y[index]);
	  //grid.centro_x[index] = distancia(grid.X[index], grid.X[index+1], grid.Y[index], grid.Y[index+1]);

	  
	  /*for(m=Nboxes; m<Nesquinas; m++)
	    {
	      
	      Esquina[m].x = grid.X[index];
	      
	      Esquina[m].y = grid.Y[index];
	      
	      Esquina[m+1].x = grid.X[index+1];
	      
	      Esquina[m+1].y = grid.Y[index+1];
	      
	      grid.centro_x[m] = distancia(Esquina[m].x, Esquina[m+1].x, Esquina[m].y, Esquina[m+1].y);
	      
	      grid.centro_y[m] = distancia(Esquina[m].x, Esquina[m+1].x, Esquina[m].y, Esquina[m+1].y);

	      printf("%lf %lf\n", grid.centro_x[m], grid.centro_y[m]);
	  
	      }*/
	  for(i=0; i<Nparticles; i++)
	    {
	      if((part[i].pos[0] >= grid.X[index]) && (part[i].pos[0] <= grid.X[index+Nboxes]))
		{
		  
		  if((part[i].pos[1] >= grid.Y[index]) && (part[i].pos[1] <= grid.Y[index+1]))
		    {
		      suma = suma + 1;  
		    }
		  
		}

	      grid.part_cells[index] = suma;
	      //printf("particles/cell[%d]: %d\n", index, grid.part_cells[index]);
	    }
	  
	}
    }

  printf("%lf %lf %d\n", grid.X[624], grid.X[624],  grid.part_cells[624]);
  
  /*Calculamos los vecinos -los 2 mas cercanos- */


  /*Distancia  al centro del sistema */
  
  for(i=0; i<Nparticles; i++)
    {
      dist[i] = distancia(1500, part[i].pos[0], 1500, part[i].pos[1]);
    }
  
  gsl_sort_smallest_index (p, k, dist, 1,Nparticles);//Se ordenan las particulas en forma creciente
  
  free(dist);
  
  /* Coordenadas del centro geometrico*/

  X_centro =  part[p[0]].pos[0];

  Y_centro =  part[p[0]].pos[1];

  printf("coordenadas del centro\n");

  printf("Xc: %lf Yc: %lf\n", X_centro, Y_centro);
  
  
  /*Trasladamos el sistema al centro (0,0)*/

  dist = (double *) malloc(Nparticles *sizeof(double));

  distmin = 100*Nparticles;

  distmax = 0.0;
  
  for(i=0; i<Nparticles; i++)
    {
      xi = X_centro - X_centro;
      yi = Y_centro - Y_centro;

      xj =  part[i].pos[0] - X_centro;
      yj =  part[i].pos[1] - Y_centro;
      
      dist[i] = distancia(xi, xj, yi, yj); //distancia de las particulas al nuevo centro
      
      if(distmin>dist[i])//se calcula la distancia min y max
	{
	  distmin=dist[i];
	}
      
      
      if(dist[i]>distmax)
	{
	  distmax=dist[i];
	}

      //write_file("distancia_centro.dat", Nparticles, part, dist);
    }
  printf("Valor min de dist %.9lf\n",distmin);
  printf("Valor max de dist: %.9lf\n",distmax);

  
  
  
  
  /* Imprimir resultados y salir */

  return 0;

}
