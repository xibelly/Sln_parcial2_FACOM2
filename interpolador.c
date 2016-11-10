/*Xibelly Eliseth Mosquera Escobar
 * 
 * programa: interpolador.c
 * 
 * 
 *Se encarga  generar de interpolar un conjunto de datos
 * 
 * La funcion "interpolador_akima_per" recibe los siguinetes parametros:
 * el punto de partida para la interpolacion, el # de datos a interpolar,
 * el conjunto x, y a ser interpolados.
 *
 */

int interpolador(char *namefunc, int k,int Npoints,double Lbox, double *y, double *f)

{

  int i;
  double yi,fi;
 
  //random generator stuff 

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
 
  //gsl interpolation stuff

  gsl_interp_accel *acc;
  acc = gsl_interp_accel_alloc();

  gsl_interp *interpolator;
  interpolator = gsl_interp_alloc(gsl_interp_cspline, Npoints); 
  
  gsl_interp_init(interpolator,y,f,Npoints);
  
  //write in a file

  char buff[200];
  sprintf(buff,"%s_%s_%d",namefunc,"interpolacion",k);
 
  FILE *pf=NULL;
  pf = fopen(buff,"w");
  
  for(i=0;i<100;i++)
    {
      
      yi = gsl_rng_uniform(r) * (Lbox);  //eje y -> distribucion de # aletorios 
      fi = gsl_interp_eval(interpolator,y,f,yi,acc);
      fprintf (pf,"%e %e\n", yi, fi);

    }
  
  gsl_interp_free (interpolator);
  gsl_interp_accel_free (acc);  
  gsl_rng_free (r);                               

  fclose(pf);
  
  return 0;
}

int curvas(int Ndiv,double Ntotcells, double Lbox, double *ycenters)
{

  int i,j,k;
  double *Mass = NULL, *Sigma = NULL, *Dens = NULL;

  Mass = (double *)malloc((size_t)Ndiv*sizeof(double)); 
  Sigma = (double *)malloc((size_t)Ndiv*sizeof(double)); 
  Dens = (double *)malloc((size_t)Ndiv*sizeof(double)); 

  for(k=0; k<Ndiv; k++) //cuenta sobre los centros en x
    {

      j=0;
      
      for(i=k; i<Ntotcells; i+=Ndiv)
	{

	  Mass[j] = celda[i].masa;
	  Sigma[j] = celda[i].sigma;
	  Dens[j] = celda[i].den_super;
	 

	  j++;
	}

      //Interpolation

      interpolador("mass",k,Ndiv,Lbox, ycenters, Mass);
      interpolador("sigmav",k,Ndiv,Lbox, ycenters, Sigma);
      interpolador("density",k,Ndiv,Lbox, ycenters, Dens);
     
    }


  free(Mass);
  free(Sigma);
  free(Dens);

  return 0;

}
