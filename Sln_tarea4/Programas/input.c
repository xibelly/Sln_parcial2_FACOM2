/*Xibelly Eliseth Mosquera Escobar
 * 
 * programa: input.c
 * 
 * 
 * Se encarga de leer los datos de
 * un archivo dado le longitud arbitraria
 * usando malloc.
 * 
 * 
 * La funcion "read_file" recibe como entradas
 * la longitud del archivo y el nombre del archivo
 * a leer. Ambas entradas se hacen usando 
 * la variable ARGS en el makefile incluido con
 * el programa.
 *
 * 
 
 */



int read_file(char *filename, int N)
{
  int i, nread;
  double X, Y;

  FILE *pf = NULL;

  printf("READING FILE: %s\n", filename);
  
  pf = fopen(filename,"r");

  if(pf == NULL)
    printf("THE %s CAN NOT BE OPENED\n",filename);
  
  datos = (struct data*) malloc(N *sizeof(struct data));

  printf(" ->In read_file:\n");
    
  for(i=0; i<N; i++)
    {
      nread = fscanf(pf,"%lf %lf",&X, &Y);

      datos[i].x = X ; 

      datos[i].f = Y ;
      
    }

  printf("READ FILE STATE: SUCESSFUL\n");
  
  return 0;
}
