/*Xibelly Eliseth Mosquera Escobar
 * 
 * programa: distribution.c
 * 
 * 
 *Se encarga  generar la distribucion de particulas en 2D
 * 
 * La funcion "random_distribution" recibe los siguinetes parametros:
 * el # de datos que se quieren generar, la vaianza para el eje x y 
 * la varianza para el eje y. Esta rutina usa gsl_ran_gaussian
 * para generar una distribucion gaussiana.
 * 
 
 */


int random_distribution(int N_part,double mu_x, double mu_y)  //Genera la distribucion de particulas
{
  int i;
  double x, y;
  const gsl_rng_type * T;
  gsl_rng * r;

  FILE *write=NULL;

  part.x = (double *) malloc(N_part *sizeof(double));  /*Particles*/

  part.y = (double *) malloc(N_part *sizeof(double));

  
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  printf("WRITING FLIE: positions.dat\n");
  
  write = fopen("positions.dat","w");

  if(write == NULL)
    printf("THE FILE: positions.dat CAN NOT BE OPENED\n");
  

  for (i = 0; i < N_part; i++)
    {
          
      x = gsl_ran_gaussian (r, mu_x); //Gaussian distribution for axis x

      part.x[i] = x * L;
      
      y = gsl_ran_gaussian (r, mu_y); //Gaussian distribution for axis y
      
      part.y[i] = y * L;

      fprintf (write,"%lf %lf \n", part.x[i] , part.y[i]);
    }

 
  printf("STATE OF PARTICLES GENERATION IS: SUCESS\n");
  
  gsl_rng_free (r);

  fclose(write);
  
  return 0;
}
