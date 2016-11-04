CC = gcc

CFLAGS = -c -O3  -I/home/$(USER)/local/include/ -I/usr/include/ -I.

CFLAGSDEBUG = -g -Wall -c -I/home/$(USER)/local/include/ -I/usr/include/ 

LFLAGS = -lm -L/home/$(USER)/local/lib -Wl,-R /home/$(USER)/local/lib 


Interpolacion:
	echo Estoy compilando $@.c
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) -lm -lgsl -lgslcblas -o  $@.x		

dms:
	echo Estoy compilando $@.c
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) -lm  -o  $@.x		

numeros_aleatorios:
	echo Estoy compilando $@.c
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) -lm -lgsl -lgslcblas -o  $@.x		

Fourier_transforms:
	echo Estoy compilando $@.c
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) -lm -lfftw3 -o  $@.x

Fourier_transforms2:
	echo Estoy compilando $@.c
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) -lm -lfftw3 -o  $@.x

Fourier_transforms3:
	echo Estoy compilando $@.c
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) -lm -lfftw3 -o  $@.x


Integracion:
	echo Estoy compilando $@.c
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) -lm -lgsl -lgslcblas -o  $@.x		

ajuste:
	echo Estoy compilando $@.c
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) -lm -lgsl -lgslcblas -o  $@.x		

random_noise:
	echo Estoy compilando $@.c
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) -lm -lgsl -lgslcblas -o  $@.x		



debug: 
	echo Estoy compilando para debugear $@.c
	$(CC) $(CFLAGSDEBUG) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS)  -o $@.x


clean:
	rm -rf *.~
	rm -rf *.out
	rm -rf *#
	rm -rf *.o
	rm -rf *.a
	rm -rf *.so
	rm *.x	
