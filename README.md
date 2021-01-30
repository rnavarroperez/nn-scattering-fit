# nn-scattering-fit
Fit a NN potential to reproduce experimental scattering data

3 different potential have been implemented. The original AV18, the Delta-Shell (DS) potential (with OPE starting at 3.0 fm) and the modified version of the DS potential that also adjust the pion-nucleon coupling constants allowing for charge symmetry breaking.

## Compilation

The code can be quickly compiled using the `make` command insde the `src/` directory. This will generate two executables `nn_fit` and `write_mc_phases`

## Running the main code.

The `nn_fit` binary is the main executable of the code. In it's default form it takes the original AV18 potential with it's original fitting parameters and readjusts them to the Granada data base using the Levenberg-Marquardt algorithm
