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

double interpolador_akima_per(double xo, int N,  double *xx, double *yy) //Interpola usando akima periodico
{  
  gsl_interp_accel *acelerador;
  gsl_spline *interpolador;

 
  acelerador   = gsl_interp_accel_alloc();
  interpolador = gsl_spline_alloc(gsl_interp_akima_periodic, (size_t) N);
  gsl_spline_init(interpolador, xx, yy, (size_t) N);
  
  value = gsl_spline_eval(interpolador, xo, acelerador);
  
  gsl_spline_free(interpolador);
  gsl_interp_accel_free(acelerador);

  
  
  return value;
  
}
