/*Xibelly Eliseth Mosquera Escobar
 * 
 * programa: random.c
 * 
 * 
 *Se encarga  generar numeros aleatorios uniformemente distribuidos usando
 *el generador gsl_rng_gfsr4 
 * 
 * La funcion "read_file" recibe como entrada
 * el # de datos, la asignasion se hace a traves de
 * la variable ARGS en el makefile incluido con
 * el programa.
 *
 * 
 
 */



int random_gfsr4(int N)
{
  
  int i;
  long seed;
  gsl_rng *rng;  //generador de # aleatorios

   
  rng = gsl_rng_alloc (gsl_rng_gfsr4);     
  seed = time(NULL)*getpid();    //semilla que cambia con el tiempo
  gsl_rng_set (rng, seed);       //se establce la semilla
  
  for(i=0; i<N; i++)
    {
      datos.result[i] = gsl_rng_uniform (rng);
      

    }
 
  gsl_rng_free (rng);

  printf("STATE OF RANDOM NUMBER GENERATION BY 'GFSR4': SUCCESS\n");
 
  return 0;
  
}


