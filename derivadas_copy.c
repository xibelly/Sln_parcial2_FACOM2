/*Xibelly Eliseth Mosquera Escobar
 * 
 * programa: derivadas.c
 * 
 * 
 *
 * La funcion "deriva" interpola los datos del campo densidad para luego
 * derivarlos -> campo densidad superficial como funcion de la masa
 *
 * NOTA: interpolamos con akima_periodico 
 */
//const size_t n = NtotalCells;
  //gsl_function H;
  //double result, result1, result2, abserr1, abserr2;

  //struct data d = {n, Lbox};

  //H.function = &deriva; //funcion interpolacion del campo densidad 
  //H.params = &d;
 
  
 
  //for(i=0; i<NtotalCells; i++)
  //{
      //gsl_deriv_central(&H, celda[i].xc, 1e-8, &result1, &abserr1); //Componente x del gradiente
      
      //gsl_deriv_central(&H, celda[i].yc, 1e-8, &result2, &abserr2); //Componente y del gradiente
      
      
      //fprintf(out9,"%lf %lf %lf %lf\n", celda[i].xc, result1, celda[i].yc, result2);

      
      
  //}  
  //fclose(out9);


double deriva(double ff, void * data)
{
  int i, j, l, nread;
  size_t n = ((struct data *)data)->n;
  size_t Lbox = ((struct data *)data)->Lbox;
  
  FILE *write = NULL;  
  FILE *read = NULL;  
   
  int Nlines;
  double xo, function;
  double *x,*y, *z, *xx, *yy, *zz;
  double xi, yi;
  
  yy = (double *) malloc(n *sizeof(double));
  xx = (double *) malloc(n *sizeof(double));
  zz = (double *) malloc(n *sizeof(double));

    
  for(i=0; i< n; i++)
    {
      xx[i] = celda[i].xc;  
      yy[i] = celda[i].yc;

      zz[i] = celda[i].den_super;  //Campo densidad superficial en funcion de la masa
   
    } 
  
  gsl_sort(xx, 1, n);
  gsl_sort(yy, 1, n);
  
  write = fopen("datos_a_derivar.dat","w");

  Nlines = 0;

  for(i=0; i<n; i++)
    {         
      if((xx[i] < xx[i+1]) && (yy[i] < yy[i+1]) )//condicion de estrictamente crecientes
	{
	  fprintf(write,"%lf %lf %lf\n", xx[i],yy[i],zz[i]); 
	  Nlines = Nlines +1 ; 
	}
    }
  fclose(write);
  
  
  y = (double *) malloc(Nlines *sizeof(double));
  x = (double *) malloc(Nlines *sizeof(double));
  z = (double *) malloc(Nlines *sizeof(double));

  const gsl_interp2d_type *T = gsl_interp2d_bilinear;

  const size_t nx = Nlines;//sizeof(xx) / sizeof(double); /* x grid points */
  const size_t ny = Nlines;//sizeof(yy) / sizeof(double); /* y grid points */
  
  gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();
  

  read = fopen("datos_a_derivar.dat","r");
  
  for(j=0; j<Nlines; j++)
    { 
      nread = fscanf(read,"%lf %lf %lf",&x[j], &y[j], &z[j]);//datos que seran interpolados
            
    }
  
  fclose(read);

  

  /* initialize interpolation */
  gsl_spline2d_init(spline, x, y, z, nx, ny); 
 
  for(i=0; i< nx; i++)
    {
           
      xi = i / (n - 1.0);
      
      for(j=0; j< ny; j++)
	{
	  
	  yi = j / (n - 1.0);
	    
	  function = gsl_spline2d_eval(spline, xi, yi, xacc, yacc);
	  
	
	} 
    }
   
   
  gsl_spline2d_free(spline);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);

 
  return function;
  

  /*
 //Para asegurar la interpolacion ordenamos en forma creciente

 //gsl_interp_accel *acelerador;
//gsl_spline *interpolador;
  gsl_sort(xx, 1, n);

  gsl_sort(yy, 1, n);

  Nlines = 0;

  write = fopen("datos_a_derivar.dat","w");

  for(i=0; i<n; i++)
    {         
      if((xx[i] < xx[i+1]) && (yy[i] < yy[i+1]) )//condicion de estrictamente crecientes
	{
	  fprintf(write,"%lf %lf\n", xx[i],yy[i]); 
	  Nlines = Nlines +1 ; 
	}
    }
  fclose(write);
  
  
  y = (double *) malloc(Nlines *sizeof(double));
  x = (double *) malloc(Nlines *sizeof(double));
  
  
  read = fopen("datos_a_derivar.dat","r");
  
  for(j=0; j<Nlines; j++)
    { 
      nread = fscanf(read,"%lf %lf",&x[j], &y[j]);//datos que seran interpolados
            
    }
  
  fclose(read);

  
  acelerador   = gsl_interp_accel_alloc();
  interpolador = gsl_spline_alloc(gsl_interp_akima_periodic, (size_t) Nlines);
  gsl_spline_init(interpolador, x, y, (size_t) Nlines);

  
  function = gsl_spline_eval(interpolador, 0, acelerador);//se interpolan los datos del campo de densidad
  

  gsl_spline_free(interpolador);
  gsl_interp_accel_free(acelerador);
  

  */

 
}


 
 
  
  
  
  
