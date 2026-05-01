del /Q *.mod 2>nul
del /Q *.exe 2>nul
del /Q *.o   2>nul
gfortran -fopenmp -c constant_mod.f90
gfortran -fopenmp -c nucleus_mod.f90
gfortran -fopenmp -c grid_mod.f90
gfortran -fopenmp -c CG_method_mod.f90
gfortran -fopenmp -c frldm_mod.f90
gfortran -fopenmp -c main.f90

gfortran -fopenmp constant_mod.o nucleus_mod.o grid_mod.o CG_method_mod.o frldm_mod.o main.o -o main.exe
