# 1D-SE-Shooting-algorithm-
In this repository a 1D schreudinger equation shooting algorithm is developed for the course "Scientific Programming".

A picture of the approximate workflow can be found in the main repository. Three seperate subroutines are called, which perform the tasks of gathering (and ordening) the input data, running the actual calculation and deleting the allocated memmory for the input data. This way these three seperate tasks can be tested easely in unit test. 

compiling order:
gfortran blas.f lapack.f diagonalization.F90 trapez.f90 WriterModule.f90 PotentialModule.f90 ThreePointSolverModule.f90 ShootingAlgorithmModule.f90 ProgramModule.f90 TestModule.f90 TestProgram.f90 
