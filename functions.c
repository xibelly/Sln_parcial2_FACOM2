#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "allvars.c"

int pos_generator(double Lbox, double sigmax, double sigmay)

{

  int i;

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

do
  {
    part[i].pos[0] = gsl_ran_gaussian(r,sigmax);
    part[i].pos[1] = gsl_ran_gaussian(r,sigmay);
    
    if(part[i].pos[0]>= -Lbox && part[i].pos[0] <= Lbox && part[i].pos[1] >= -Lbox && part[i].pos[1] <= Lbox)
      
      {
	i++;
      }
    
  }while(i<Ntotpart);

  gsl_rng_free(r);

  return 0;
}

int interp_func(char *namefunc, int k,int Npoints,double Lbox, double *y, double *f,const gsl_interp_type *a)

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
  interpolator = gsl_interp_alloc(a, Npoints); 
  
  gsl_interp_init(interpolator,y,f,Npoints);
  
  //write in a file

  char buff[200];
  sprintf(buff,"%s_%s_%s_%d",namefunc,"interp",gsl_interp_name(interpolator),k); 
  FILE *pf=NULL;
  pf = fopen(buff,"w");
  
  for(i=0;i<100;i++)
    {
      
      yi = gsl_rng_uniform(r) * (Lbox);
      fi = gsl_interp_eval(interpolator,y,f,yi,acc);
      fprintf (pf,"%e %e\n", yi, fi);

    }
  
  gsl_interp_free (interpolator);
  gsl_interp_accel_free (acc);  
  gsl_rng_free (r);                               

  fclose(pf);
  
  return 0;
}



int xfix_curves(int Ndiv,double Ntotcells, double Lbox, const gsl_interp_type *a, double *ycenters)

{

  int i,j,k;
  double *shareMass = NULL, *shareSigmaV = NULL, *shareDens = NULL;

  shareMass = (double *)malloc((size_t)Ndiv*sizeof(double)); 
  shareSigmaV = (double *)malloc((size_t)Ndiv*sizeof(double)); 
  shareDens = (double *)malloc((size_t)Ndiv*sizeof(double)); 

  for(k=0; k<Ndiv; k++) //cuenta sobre los centros en x

    {

      j=0;
      //  printf("************\t k = %d\n", k);
      //printf("==array of masses===\n");
      for(i=k; i<Ntotcells; i+=Ndiv)
	{

	  shareMass[j] = celda[i].Mass;
	  shareSigmaV[j] = celda[i].sigmaV;
	  shareDens[j] = celda[i].density;
	  //  printf("\t%lf %lf\n",ycenters[j],shareMass[j]);

	  j++;
	}

      //printf("=========\n");
      
      //printf("Interpolation=====\n");

      //ycenters[0] = 0.0;
      // ycenters[Ndiv-1] = Lbox;

      //interpolating
      interp_func("mass",k,Ndiv,Lbox, ycenters, shareMass,a);
      interp_func("sigmav",k,Ndiv,Lbox, ycenters, shareSigmaV,a);
      interp_func("density",k,Ndiv,Lbox, ycenters, shareDens,a);
     
    }


  free(shareMass);
  free(shareSigmaV);
  free(shareDens);

  return 0;

}


int gsl_derivatives(int Ndiv, double *xcenters, double *ycenters, double *DensField)
{

  size_t i, j;

  double xi;
  double yj;

  FILE *pf = NULL;
  pf = fopen("derivatives_density_field.dat","w");

  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  const size_t N = 10;             /* number of points to interpolate */
  
  gsl_interp2d *interpolator = gsl_interp2d_alloc(T, (size_t)Ndiv, (size_t)Ndiv);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();

  gsl_interp2d_init(interpolator, xcenters, ycenters, DensField, Ndiv, Ndiv);

  /*
  // interpolate N values in x and y and print out grid for plotting
  for (i = 0; i < N; ++i)
    {
      double xi = i / (N - 1.0);

      for (j = 0; j < N; ++j)
        {
          double yj = j / (N - 1.0);
          double densij = gsl_interp2d_eval(interpolator,xcenters,ycenters,DensField, xi, yj, xacc, yacc);

	  double dfdx = gsl_interp2d_eval_deriv_x(interpolator,xcenters, ycenters, DensField, xi,  yj, xacc, yacc);
	  double d2fdx = gsl_interp2d_eval_deriv_xx(interpolator,xcenters, ycenters, DensField, xi,  yj, xacc, yacc);
	  double dfdy = gsl_interp2d_eval_deriv_y(interpolator,xcenters, ycenters, DensField, xi,  yj, xacc, yacc);
	  double d2fdy = gsl_interp2d_eval_deriv_yy(interpolator,xcenters, ycenters, DensField, xi,  yj, xacc, yacc);

          printf("%lf %lf %lf %lf\n", xi, yj, densij,dfdx);
        }
      //printf("\n");
    }
 */ 

  for (i=0; i < Ndiv; i++)
    {

      xi = xcenters[i];

      for (j = 0; j < Ndiv; ++j)
        {

	  yj = ycenters[j];
     	  
	  double densij = gsl_interp2d_eval(interpolator,xcenters,ycenters,DensField, xi, yj, xacc, yacc);

	  double dfdx = gsl_interp2d_eval_deriv_x(interpolator,xcenters, ycenters, DensField, xi,  yj, xacc, yacc);
	  double d2fdx = gsl_interp2d_eval_deriv_xx(interpolator,xcenters, ycenters, DensField, xi,  yj, xacc, yacc);
	  double dfdy = gsl_interp2d_eval_deriv_y(interpolator,xcenters, ycenters, DensField, xi,  yj, xacc, yacc);
	  double d2fdy = gsl_interp2d_eval_deriv_yy(interpolator,xcenters, ycenters, DensField, xi,  yj, xacc, yacc);

          fprintf(pf,"%lf %lf %lf %lf %lf %lf %lf\n", xi, yj, densij,dfdx,dfdy,d2fdx,d2fdy);
        }

    }

  gsl_interp2d_free(interpolator);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
 
  fclose(pf);

  return 0;
}


