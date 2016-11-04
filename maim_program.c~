#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
//#include "allvars.c"
//#include "read_binary.c"

int main(int argc, char *argv[])
{

  int Ncell, i, j, k, NtotalCells, cell_ID, p;
  float xmin, xmax, ymin, ymax, zmin, zmax, BoxSize;
  float CellSize, CellVol;
  float x, y, z;
  double median,Nguess, *velocities;
  char *input_file=NULL;
  FILE *file = NULL, *file1=NULL;
  FILE *datospos = NULL, *sigma= NULL;

  float Xmax;
  float Ymax;
  float Zmax;
  float Xmin;
  float Ymin;
  float Zmin;

  struct cells *celda=NULL;

  input_file = argv[1];       //recibe el archivo con las posiciones                                                                                                    
  Ncell = atoi(argv[2]); //recibe el numero de particiones que haremos sobre cada eje.   

  if(argc != 3)
    {
      printf(" USE: ./exec <1DNcells> <Nparticles> <input_file>\n");
      exit(0);
    }

  printf(" running %d %s\n",Ncell, input_file);

  readGadgetBinary(input_file);

  //-------------------------------------------------------------------------CONSTRUCCION DEL GRID                                                                      


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
  return 0;

}

