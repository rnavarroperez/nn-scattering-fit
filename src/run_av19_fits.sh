#!/bin/bash

./nn_fit av19.namelist > av19_fit.out
./nn_fit av19_cutoff11.namelist > av19_cutoff11_fit.out
./nn_fit av19_cutoff13.namelist > av19_cutoff13_fit.out
./nn_fit av19_cutoff15.namelist > av19_cutoff15_fit.out
./nn_fit av19_cutoff17.namelist > av19_cutoff17_fit.out
./nn_fit av19_cutoff19.namelist > av19_cutoff19_fit.out
./nn_fit av19_cutoff23.namelist > av19_cutoff23_fit.out

