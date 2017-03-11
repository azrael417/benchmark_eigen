FC = ftn

FFLAG = -qopenmp

.SUFFIXES : .o .f90
.f90.o :
	$(FC) $(FFLAG) -c $<

MAIN = eigen_main.o

OBJCT = \
        parameters.o \
        matrix_module.o \
        pdsyevd_module.o \
        elpa_module.o

all: ./eigen.x

./eigen.x: $(OBJCT) $(MAIN)
	$(FC) $(FFLAG) -o ./eigen.x $(MAIN) $(OBJCT)

clean:
	rm -f *.o *.mod

cleanall:
	rm -f *.o *.mod eigen.x
