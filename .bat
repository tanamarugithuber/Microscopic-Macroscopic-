
gfortran -c constant_mod.f90
gfortran -c nucleus_mod.f90
gfortran -c grid_mod.f90
gfortran -c CG_method_mod.f90
gfortran -c frldm_mod.f90
gfortran -c main.f90
gfortran constant_mod.o nucleus_mod.o grid_mod.o CG_method_mod.o frldm_mod.o main.o -o main.exe
