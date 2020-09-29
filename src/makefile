COMPILER = gfortran

FLAGS = -fcheck=all -Wall -O3

MAIN_EXEC = nn_fit

MC_PHASES_EXEC = write_mc_phases

SRC = $(wildcard *.f90)

OBJ = $(SRC:.f90=.o)

MAIN_OBJ = $(filter-out mc_phases.o, $(OBJ))

MC_PHASES_OBJ = $(filter-out main.o, $(OBJ))

all: $(MC_PHASES_EXEC) $(MAIN_EXEC)

$(MC_PHASES_EXEC): $(MC_PHASES_OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^ $(LAPACK)

$(MAIN_EXEC): $(MAIN_OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^ $(LAPACK)

precisions.o: precisions.f90
	$(COMPILER) $(FLAGS) -c $<

coul90.o: coul90.f90
	$(COMPILER) $(FLAGS) -c $<

string_functions.o: string_functions.f90
	$(COMPILER) $(FLAGS) -c $<

utilities.o: utilities.f90 precisions.o
	$(COMPILER) $(FLAGS) -c $<

num_recipes.o: num_recipes.f90 precisions.o constants.o
	$(COMPILER) $(FLAGS) -c $<

test_derivatives.o: test_derivatives.f90 precisions.o num_recipes.o
	$(COMPILER) $(FLAGS) -c $<

constants.o: constants.f90 precisions.o
	$(COMPILER) $(FLAGS) -c $<

st_basis_2_partial_waves.o: st_basis_2_partial_waves.f90 precisions.o
	$(COMPILER) $(FLAGS) -c $<

em_nn_potential.o: em_nn_potential.f90 precisions.o constants.o st_basis_2_partial_waves.o
	$(COMPILER) $(FLAGS) -c $<

av18.o: av18.f90 precisions.o constants.o em_nn_potential.o st_basis_2_partial_waves.o
	$(COMPILER) $(FLAGS) -c $<

pion_exchange.o : pion_exchange.f90 precisions.o em_nn_potential.o
	$(COMPILER) $(FLAGS) -c $<

delta_shell.o: delta_shell.f90 precisions.o constants.o utilities.o pion_exchange.o
	$(COMPILER) $(FLAGS) -c $<

nn_phaseshifts.o: nn_phaseshifts.f90 precisions.o num_recipes.o constants.o coul90.o delta_shell.o
	$(COMPILER) $(FLAGS) -c $<

amplitudes.o: amplitudes.f90 precisions.o nn_phaseshifts.o constants.o num_recipes.o utilities.o
	$(COMPILER) $(FLAGS) -c $<

observables.o: observables.f90 constants.o amplitudes.o nn_phaseshifts.o precisions.o deuteron.o delta_shell.o
	$(COMPILER) $(FLAGS) -c $<

read_write.o: read_write.f90 precisions.o amplitudes.o observables.o nn_phaseshifts.o constants.o  delta_shell.o utilities.o
	$(COMPILER) $(FLAGS) -c $<

exp_data.o: exp_data.f90 precisions.o string_functions.o utilities.o
	$(COMPILER) $(FLAGS) -c $<

deuteron.o: deuteron.f90 precisions.o delta_shell.o constants.o
	$(COMPILER) $(FLAGS) -c $<

derivatives.o: derivatives.f90 av18.o num_recipes.o nn_phaseshifts.o amplitudes.o observables.o deuteron.o delta_shell.o pion_exchange.o
	$(COMPILER) $(FLAGS) -c $<

chi_square.o: chi_square.f90 precisions.o exp_data.o observables.o delta_shell.o amplitudes.o
	$(COMPILER) $(FLAGS) -c $<

main.o: main.f90 exp_data.o
	$(COMPILER) $(FLAGS) -c $<

mc_phases.o: mc_phases.f90 precisions.o delta_shell.o read_write.o
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod

mrproper: clean
	rm -rf $(MAIN_EXEC) $(MC_PHASES_EXEC)
