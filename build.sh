#!/bin/bash
# FRLDM + Shell + pairing (main.exe). Needs FFTW3 and LAPACK.

rm -f *.mod
rm -f *.exe
rm -f *.o

export OMP_NUM_THREADS=4
FC="gfortran -O2 -fopenmp"

$FC -c constant_mod.f90
$FC -c nucleus_mod.f90
$FC -c shape_mod.f90
$FC -c grid_mod.f90
$FC -I/usr/include -c calculate_potential_hockney.f90
$FC -c frldm_mod.f90
$FC -c micro_constant_mod.f90
$FC -c ho_basis_mod.f90
$FC -c potential_gh_mod.f90
$FC -c shell_bcs_mod.f90
$FC -c main.f90

$FC \
    constant_mod.o \
    nucleus_mod.o \
    shape_mod.o \
    grid_mod.o \
    calculate_potential_hockney.o \
    frldm_mod.o \
    micro_constant_mod.o \
    ho_basis_mod.o \
    potential_gh_mod.o \
    shell_bcs_mod.o \
    main.o \
    -lfftw3 -llapack \
    -o main.exe
