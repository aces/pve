# Makefile for pve and pve3

CC      = cc
INCLUDE = -I/usr/local/unstable/include -I/usr/local/unstable/mni/include 
LDFLAGS = -L/usr/local/unstable/lib -L/usr/local/unstable/mni/lib -lvolume_io -lminc -lnetcdf -lm -lmalloc
CFLAGS  = -O2 -xansi $(INCLUDE)
OBJ     =  pve.o pve_aux.o least_trimmed_squares.o minvarellipsoid.o 
OBJ3    =  pve3.o pve3_aux.o matrix3.o minvarellipsoid3.o
PROG    = ../bin/pve
PROG3   = ../bin/pve3


all  : pve pve3

pve  : $(OBJ)
	$(CC) $(CFLAGS)  $(OBJ) -o $(PROG)  $(LDFLAGS)

pve3 : $(OBJ3)
	$(CC) $(CFLAGS) $(OBJ3) -o $(PROG3)  $(LDFLAGS)

pve.o                  : pve.c pve_aux.h 
pve_aux.o              : pve_aux.c pve_aux.h
minvarellipsoid.o      : minvarellipsoid.c minvarellipsoid.h
least_trimme_squares.o : least_trimmed_squares.c least_trimmed_squares.h
pve3.o                 : pve3.c pve3_aux.h matrix3.h
pve3_aux.o             : pve3_aux.c pve3_aux.h matrix3.h
matrix3.o              : matrix3.c matrix3.h
minvarellipsoid3.o     : minvarellipsoid3.c minvarellipsoid3.h matrix3.h 

clean:
	\rm -f $(OBJ) $(OBJ3) $(PROG) $(PROG3)
