/*Xibelly Eliseth Mosquera Escobar
 * 
 * programa: derivadas.c
 * 
 * 
 *
 * La funcion "deriva" interpola los datos del campo densidad para luego
 * derivarlos -> campo densidad superficial respecto a las coordenadas de las
 * celdas
 *
 * NOTA: interpolamos con gsl_interp2d_bilinear
 */


void deriva(int N)
{
  int i, j, nread;
  int cell_ID;
 
  FILE *write = NULL;  
  FILE *read = NULL;  
  
  int n = N*N;
  
  int Nlines;
  double dx, dy, ddx, ddy;
  double *x,*y, *z, *xx, *yy, *zz;
  double xi, yi;
  
  yy = (double *) malloc(n *sizeof(double));
  xx = (double *) malloc(n *sizeof(double));
  zz = (double *) malloc(n *sizeof(double));

  out9=fopen("derivates_of_desity.dat","w");
  
  if(out9 == NULL)
    printf("THE FILE: derivates_of_desity.dat  CAN NOT BE OPENED\n");
  
  
  printf("WRITING FLIE: derivates_of_desity.dat\n");
  
  write = fopen("datos_a_derivar.dat","w");

  
  for(j=0; j<N; j++)
    
    {    
      
      for(i=0; i<N; i++)

	{   
	  cell_ID=Ncell * j + i;       //identificador de celda
	  celda[cell_ID].cell_ID = cell_ID;
	  
	  xx[cell_ID] = celda[cell_ID].x;  
	  yy[cell_ID] = celda[cell_ID].y;
	  
	  zz[cell_ID] = celda[cell_ID].den_super;  //Campo densidad superficial en funcion de la masa
	  fprintf(write,"%lf %lf %lf\n", xx[cell_ID],yy[cell_ID],zz[cell_ID]);

	 	 
	} 
  
    }
  fclose(write);
  	  
  y = (double *) malloc(N *sizeof(double));
  x = (double *) malloc(N *sizeof(double));
  

  for(i=0; i<N; i++)
    {
      x[i] = CellSize * 0.5 + i * CellSize; //primera fila del grid  
      y[i] = CellSize * 0.5 + i * CellSize; //primera columna del grid     
    }

  
  
      
  
  const gsl_interp2d_type *T = gsl_interp2d_bilinear;

  const size_t nx = N; 
  const size_t ny = N; 
  
  gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();
  
  
  /* initialize interpolation */
  gsl_spline2d_init(spline, x, y, zz, nx, ny); 
 

  
  /*interpolate and derivate*/
  for(i=0; i< nx; i++)
    {
           
      xi = 0.5 * CellSize+i / (N - 1.0);
      
      for(j=0; j< ny; j++)
	{
	 
	  yi = 0.5 * CellSize +j / (N - 1.0);
	    
	  dx = gsl_spline2d_eval(spline, xi, yi, xacc, yacc);
	  dy = gsl_spline2d_eval_deriv_y(spline, xi, yi, xacc, yacc);

	  ddx = gsl_spline2d_eval_deriv_xx(spline, xi, yi, xacc, yacc);
	  ddy = gsl_spline2d_eval_deriv_yy(spline, xi, yi, xacc, yacc);
	  
	  fprintf(out9,"%e %e %e %e %e %e\n",j*CellSize,j*CellSize ,dx, dy, ddx, ddy);
	} 
    }
   
   
  gsl_spline2d_free(spline);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);

 
  printf("STATE OF DERIVATES OF SURFACE DENSITY IS: SUCESS\n");
}


 
 
  
  
  
  
