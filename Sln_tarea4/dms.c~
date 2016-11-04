/*
Xibelly Eliseth Mosquera Escobar

Solucion punto 1 Tarea 4 COMPUTACION CIENTIFICA AVANZADA

1) Generar una tabla de datos con 20 puntos entre x=0 y x=2*M_PI

2) Interpolar -> vamos a generar una funcion para cada tipo de interpolador
   evluar la funcion en puntos igualmente espaciados, con frecuencia de 0.01

   vamos a imprimir en disco cada tipo de interpolacion

3) Graficar los resultados de todos los tipos de interpolador

4) Evaluar la diferencia cuadratica media entre los diferentes metodos de interpolacion usados. Para ello vamos a crear una funcion que recibe las interpolaciones que van a ser comparadas y vamos imprimir en pantalla cada resultado talque el usuario decida cual es el metodo que presenta menos error.

 */

/*
  NOTA: ESTE PROGRAMA SE ENCARGA DEL PUMTO 4)

 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>


//Variables globales//

int N;
char *file;
FILE *out=NULL;

//Estructuras//

struct data   
{
  double x, f;
  
}; 
struct data *datos; 
 
//Llamados a subrutinas//

#include "input.c"

//Funciones//

double dms(double f1, double f2)
{
  double  DMS;

  DMS = (f1 - f2) * (f1 - f2);

  return DMS;

}

//Programa Principal//

int main (int argc, char *argv[])
{
  
  int i;
  double xo, yo;
  double dx, suma;
  double f1, f2;
  double R, x;


  printf("%d\n",argc);
  
  if(argc != 3)
    {
      printf("ERROR--> use as:\n");
      printf("%s Nline file\n",argv[0]);
      exit(0);  
    }
  //Carga de parametros//

  N = atoi(argv[1]);

  file = argv[2];

  printf("%s %d %s\n",argv[0], N, file);

  datos  = (struct data *) malloc(N *sizeof(struct data));

  //Lectura archivos -carga de datos-//
  
  read_file(file,N);

  //Se calcula la diferencia cuadratica media//

  suma = 0.0;

  dx = 2 * M_PI/ N;

  out = fopen("fdex_sample.dat","w");
  
  for(i=0; i<N; i++)
    {
      x = dx * i;

      f1 = (x) *sin (x);
     
      f2 = datos[i].f;

      R = dms(f1, f2);

      suma = suma + R ;

      fprintf(out,"%lf %lf\n",x, f1);
      
    }
 printf("THE DIFERENCE MEAN SQUARE OF %s IS:\n", file);
 printf("%lf\n", suma/N);
  

  return 0;
}

