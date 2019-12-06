COMPILER= gfortran

FLAGS = -Wall

EXEC = nn_fit

SRC = $(wildcard *.f90) 

OBJ = $(SRC:.f90=.o)

$(EXEC): $(OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^ $(LAPACK)

precisions.o: precisions.f90
	$(COMPILER) $(FLAGS) -c $<

av18.o: av18.f90 precisions.o
	$(COMPILER) $(FLAGS) -c $<

main.o: main.f90 precisions.o av18.o
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod

mrproper: clean
	rm -rf $(EXEC)
