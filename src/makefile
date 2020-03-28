COMPILER= gfortran

FLAGS = -Wall

EXEC = nn_fit

SRC = $(wildcard *.f90)

OBJ = $(SRC:.f90=.o)

$(EXEC): $(OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^ $(LAPACK)

precisions.o: precisions.f90
	$(COMPILER) $(FLAGS) -c $<

coul90.o: coul90.f90
	$(COMPILER) $(FLAGS) -c $<

num_recipes.o: num_recipes.f90 precisions.o constants.o
	$(COMPILER) $(FLAGS) -c $<

test_derivatives.o: test_derivatives.f90 precisions.o num_recipes.o av18.o
	$(COMPILER) $(FLAGS) -c $<

constants.o: constants.f90 precisions.o
	$(COMPILER) $(FLAGS) -c $<

av18.o: av18.f90 precisions.o constants.o num_recipes.o
	$(COMPILER) $(FLAGS) -c $<

nn_phaseshifts.o: nn_phaseshifts.f90 precisions.o num_recipes.o constants.o coul90.o av18.o
	$(COMPILER) $(FLAGS) -c $<

amplitudes.o: amplitudes.f90 precisions.o nn_phaseshifts.o
	$(COMPILER) $(FLAGS) -c $<

observables.o: observables.f90 constants.o amplitudes.o nn_phaseshifts.o precisions.o
	$(COMPILER) $(FLAGS) -c $<

read_write.o: read_write.f90 amplitudes.o precisions.o constants.o nn_phaseshifts.o
	$(COMPILER) $(FLAGS) -c $<

main.o: main.f90 precisions.o
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod

mrproper: clean
	rm -rf $(EXEC)
