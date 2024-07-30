# nn-scattering-fit

Fit a NN potential to reproduce experimental scattering data

4 different potentials have been implemented. The original AV18, the Delta-Shell (DS) potential (with OPE starting at 3.0 fm), the modified version of the DS potential that also adjust the pion-nucleon coupling constants allowing for charge symmetry breaking, and the fully local version of the Norfolk potentials at N3LO.

## Compilation

The code can be quickly compiled using the `make` command inside the `src/` directory. This will generate three executables `nn_fit`, `write_mc_phases` and `plot_potentials`

## Running the main code.

The `nn_fit` binary is the main executable of the code. If executed without inputs, it takes the original AV18 potential with its original fitting parameters and readjusts them to the Granada data base using the Levenberg-Marquardt algorithm

### Input namelist file

A namelist file can be given as an input when executing on command line (e.g. `./nn_fit av18.namelist`) to specify what kind of potential to adjust along with other options for the optimization. Several namelist files are included as examples. The namelist groups that can be included in the file are

1. `data_base`
2. `nn_potential`
3. `local_integration`
4. `delta_shell_integration`
5. `deuteron`
6. `potential_parameters`
7. `adjust_parameter`
8. `output`

Not every group has to be included in the namelist file. Not including a group can mean that either that group is not relevant to that type of interaction (e.g. `delta_shell_integration` for a fully local potential) or that default values for the variables in that namelist group will be used. Similarly if a variable within a namelist group is not included a default value will be used for that variable. Each namelist group is explained below.

#### `data_base`

To which database will the interaction be adjusted. The one variable in this group is `database_file` and should be the file containing the data base to be fitted. There are currently two options; `'database/granada_database.dat'` (default value) and `'database/phases_database.dat'`. The first one corresponds to the self-consistent database by the Granada group; the second one fits to a collection of phase-shifts instead of observables. Fitting to phase-shifts is useful when fitting a new interaction to obtain parameters that can be used as starting values for a fit to observables.

#### `nn_potenial`

Contains two variables to specify which type of potential will be fitted and up to what laboratory energy include in the database.

The first variable is `name` and can have four possible values

1. `AV18`. The original AV18 potential and its related variations
2. `N3LO`. The fully local version of the Norfolk potential at N3LO
3. `ds_ope30`. The Delta-Shell (DS) potential (with OPE starting at 3.0 fm)
4. `ds_ope30_fff`. The modified version of the DS potential that also adjusts the pion-nucleon coupling constants allowing for charge symmetry breaking

The second variable is `t_lab_limit` and corresponds to the maximum laboratory energy in MeV that will be included in the data to fit. The default value is `350.0` MeV

#### `local_integration`

Contains two variables to specify how a local interaction will be integrated when calculating the corresponding phase-shifts.

The first variable is `r_max` and indicates the upper integration limit in fm. Here one usually wants to use a value where the potential is so small that values above it do not contribute to the phase-shift. Recommended values are between `12.5` fm and `13.0` fm. Using much larger values will simply add numerical noise to the calculation of the phase-shifts.

The second variable is `delta_r` and indicates the integration step in fm. A recommended value is `0.0078125` fm. Bigger values may not be precise enough and smaller values can unnecessarily increase computation time.

#### `delta_shell_integration`

Contains four variables to specify how a Delta-Shell interaction will be integrated when calculating the corresponding phase-shifts.

The first variable is `r_max` and indicates the upper integration limit in fm. Here one usually wants to use a value where the potential is so small that values above it do not contribute to the phase-shift. Recommended values are between `12.5` fm and `13.0` fm. Using much larger values will simply add numerical noise to the calculation of the phase-shifts.

The second variable is `n_lambdas` and indicates the number of delta-shells used between 0.0 and the start of the OPE potential. For the standard DS potential the correct value is `5`.

The third variable is `dr_core` and indicates the distance in fm between the delta-shells. The first delta-shell is located at the value of `dr_core`. For the standard DS potential the correct value is `0.6`.

The fourth variable is `dr_tail` and indicates the integration step in fm once the OPE potential starts at `n_lambdas*dr_core`. For the standard DS potential the correct value is `0.5` fm. This may seem not small enough, but due to the soft changing nature of OPE after 3.0 fm the integration is precise and smaller values can unnecessarily increase computation time.

#### `deuteron`

Contains two variables to specify whether and how the deuteron binding energy will be included among the observables to fit. The experimental error bar of the deuteron binding energy is so small that it can make it hard for a potential to leave a local minimum during the optimization procedure. This is especially true when new adjustable parameters are added to an already optimized interaction. A not so elegant workaround is to fit the new interaction without including the deuteron binding energy, and then repeat the optimization including the deuteron binding energy using the optimized parameters as a starting point.

The first variable is `fit_deuteron` and its default value is `.true.`.

The second variable is `relativistic` and indicates whether or not to use relativistic kinematics when calculating the deuteron binding energy (see the `binding_energy` subroutine in `deuteron.f90` for details). The AV18 and Norfolk potentials use non-relativistic kinematics. The DS potentials use relativistic kinematics.

#### `potential_paramters`

Contains a single variable `parameters` which is an array of reals with the adjustable parameters. The number of parameters in the array depend on the type of potential that will be adjusted.

#### `adjust_parameters`

Contains a single variable `mask` which is an array of logical variables. It has the same size as `parameters`. Parameters marked as `.true.` will be adjusted during the optimization. Parameters marked as `.false.` will be left unchanged during the optimization. 

#### `output`

Contains two variables to specify whether and how to write results into output files.

The first variable is `save_results` and its default value is `.true.`.

The second variable is `output_name` and its value will be used as a prefix in the name of all of the output files described below. Its default value is `'results'`

 
### Output

If the `save_results` variable within the `output` namelist is set to `.true.` several output files are created after the optimization is completed. The value of the `output_name` variable is used as a prefix when naming all of the output files. The output files created by the program are

#### `'[output_name]_parameters.txt'`

It lists the potential setup and initial parameters as specified by the given namelist file or default values. Then lists the final adjusted parameters along with their error bars. Parameters that were kept fixed during the optimization are marked with a `*`. It also lists the total chi square, the number of data included in the fit and the chi square per number of data.

The number of data can have small fluctuations due to the fact that some normalization may or may not be included depending in their value.

#### `'[output_name]_plots.dat'`

If the potential is of the AV18 type, the nuclear part of the potential is saved (in MeV) as a function of radius (in fm) in the operators basis.

#### `'[output_name].in'`

If the potential is of the AV18 type, the adjusted parameters are written in the format used as input files in the codes used by Maria Piarulli and Bob Wiring. See `write_marias_format` subroutine in `av18_compatibility.f90` for details.

#### `'[output_name]_pp_v_partial_wave.dat'` and `'[output_name]_np_v_partial_wave.dat'`

If the potential is of the local type, the potential is saved (in MeV) as a function of radius (in fm) in different partial waves.


#### `'[output_name]_phases.txt'`

A file with 3 tables of phase-shifts with error bars. The phase-shifts are calculated at the 11 'canonical' energies 

#### `'[output_name]_pp_phases_t_lab.dat'` and `'[output_name]_np_phases_t_lab.dat'`

The potential phase-shifts are saved (in degrees) as a function of laboratory energy (in MeV) in different partial waves.

#### `'[output_name]_pp_phases_k_cm.dat'` and `'[output_name]_np_phases_k_cm.dat'`

The potential phase-shifts are saved (in degrees) as a function of center of mass momentum (in MeV) in different partial waves.

## Other Executables

### `plot_potentials`

This executable can be used when you simply want to plot the inner part (from `0.0` to `2.0` fm) of an already optimized local potential in operator basis with error bars. The already optimized parameters are given with a namelist file, the code will calculate the total chi square once to obtain the covariance matrix and then plot the potential with error bars.

### `write_mc_phases`

Reads a set of monte-carlo generated set of parameters for the DS potential (generated with a different older code) and writes the corresponding phase-shifts.

## Using the original constant values for the AV18 potential.

Some fundamental constants like hbar times c, masses, magnetic moments, fine structure, etc currently have values that are different from the values used in the original AV18 potential. In order to check for consistency with previous results the original values can be used by uncommenting the corresponding lines in the `constants.f90` file and recompiling the code.