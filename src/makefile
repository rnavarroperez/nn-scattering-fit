COMPILER = gfortran

FLAGS = -O3 -fopenmp -fcheck=all -Wall -pedantic -fimplicit-none -std=f2008 -fwhole-file -march=native

#FLAGS = -g -O0 -fcheck=all -Wall -Wextra -Wconversion -pedantic -fimplicit-none -fbacktrace -std=f2008 -fwhole-file

LAPACK = -llapack -lblas

MAIN_EXEC = nn_fit

MC_PHASES_EXEC = write_mc_phases

PLOT_POT_EXEC = plot_potential

SRC = $(wildcard *.f90)

OBJ = $(SRC:.f90=.o)

MAIN_OBJ = $(filter-out mc_phases.o plot_potential.o, $(OBJ))

MC_PHASES_OBJ = $(filter-out main.o plot_potential.o, $(OBJ))

PLOT_POT_OBJ = $(filter-out main.o mc_phases.o, $(OBJ))

all: $(MC_PHASES_EXEC) $(MAIN_EXEC) $(PLOT_POT_EXEC)

$(MAIN_EXEC): $(MAIN_OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^ $(LAPACK)

$(MC_PHASES_EXEC): $(MC_PHASES_OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^ $(LAPACK)

$(PLOT_POT_EXEC): $(PLOT_POT_OBJ)
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

long_range_chiral_potentials.o: long_range_chiral_potentials.f90 precisions.o constants.o special_functions.o quadrature.o
	$(COMPILER) $(FLAGS) -c $<

short_range_chiral_potentials.o: short_range_chiral_potentials.f90 precisions.o constants.o
	$(COMPILER) $(FLAGS) -c $<

special_functions.o: special_functions.f90 precisions.o
	$(COMPILER) $(FLAGS) -c $<
	
basis_change.o: basis_change.f90 precisions.o
	$(COMPILER) $(FLAGS) -c $<

em_nn_potential.o: em_nn_potential.f90 precisions.o constants.o basis_change.o quadrature.o
	$(COMPILER) $(FLAGS) -c $<

av18.o: av18.f90 precisions.o constants.o em_nn_potential.o basis_change.o string_functions.o
	$(COMPILER) $(FLAGS) -c $<

av18_compatibility.o: av18_compatibility.f90 precisions.o
	$(COMPILER) $(FLAGS) -c $<

pion_exchange.o : pion_exchange.f90 precisions.o constants.o basis_change.o
	$(COMPILER) $(FLAGS) -c $<

delta_shell.o: delta_shell.f90 precisions.o constants.o utilities.o em_nn_potential.o pion_exchange.o string_functions.o
	$(COMPILER) $(FLAGS) -c $<

nn_phaseshifts.o: nn_phaseshifts.f90 precisions.o num_recipes.o constants.o coul90.o delta_shell.o em_nn_potential.o quadrature.o
	$(COMPILER) $(FLAGS) -c $<

amplitudes.o: amplitudes.f90 precisions.o nn_phaseshifts.o constants.o num_recipes.o utilities.o
	$(COMPILER) $(FLAGS) -c $<

observables.o: observables.f90 constants.o amplitudes.o nn_phaseshifts.o precisions.o deuteron.o delta_shell.o
	$(COMPILER) $(FLAGS) -c $<

read_write.o: read_write.f90 precisions.o amplitudes.o observables.o nn_phaseshifts.o constants.o  delta_shell.o utilities.o av18.o av18_compatibility.o long_range_chiral_potentials.o short_range_chiral_potentials.o
	$(COMPILER) $(FLAGS) -c $<

exp_data.o: exp_data.f90 precisions.o string_functions.o utilities.o amplitudes.o
	$(COMPILER) $(FLAGS) -c $<

deuteron.o: deuteron.f90 precisions.o delta_shell.o constants.o ode_solver.o quadrature.o
	$(COMPILER) $(FLAGS) -c $<

derivatives.o: derivatives.f90 precisions.o num_recipes.o av18.o delta_shell.o observables.o
	$(COMPILER) $(FLAGS) -c $<

chi_square.o: chi_square.f90 precisions.o exp_data.o observables.o delta_shell.o
	$(COMPILER) $(FLAGS) -c $<

optimization.o: optimization.f90 precisions.o exp_data.o delta_shell.o chi_square.o av18.o read_write.o
	$(COMPILER) $(FLAGS) -c $<

quadrature.o: quadrature.f90 precisions.o
	$(COMPILER) $(FLAGS) -c $<

ode_solver.o: ode_solver.f90 precisions.o delta_shell.o
	$(COMPILER) $(FLAGS) -c $<

main.o: main.f90 precisions.o av18.o av18_compatibility.o delta_shell.o exp_data.o optimization.o
	$(COMPILER) $(FLAGS) -c $<

mc_phases.o: mc_phases.f90 precisions.o delta_shell.o read_write.o
	$(COMPILER) $(FLAGS) -c $<

random_num.o: random_num.f90 constants.o precisions.o
	$(COMPILER) $(FLAGS) -c $<

randomize_exp.o: randomize_exp.f90 exp_data.o precisions.o random_num.o
	$(COMPILER) $(FLAGS) -c $<

plot_potential.o: plot_potential.f90
	$(COMPILER) $(FLAGS) -c $<	

chiral_potential.o: chiral_potential.f90 precisions.o delta_shell.o long_range_chiral_potentials.o short_range_chiral_potentials.o
	$(COMPILER) $(FLAGS) -c $<	

clean:
	rm -rf *.o *.mod

mrproper: clean
	rm -rf $(MAIN_EXEC) $(MC_PHASES_EXEC) $(PLOT_POT_EXEC)
