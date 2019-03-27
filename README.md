# 1D-SE-Shooting-algorithm-
In this repository a 1D schreudinger equation shooting algorithm is developed for the course "Scientific Programming".

A picture of the approximate workflow can be found in the main repository. Three seperate subroutines are called, which perform the tasks of gathering (and ordening) the input data, running the actual calculation and deleting the allocated memmory for the input data. This way these three seperate tasks can be tested easely in unit test. 


Note after Update on 27 march: Due to changing of some file names, a few non-updated/ used files are still present in the main folder. I will try to remove them, but cant figure it out yet. for now just use the following compiling code: 
gfortran lapack.f blas.f diagonalization.F90 trapez.f90 PotentialModule.f90 ThreePointSolverModule.f90 ShootingAlgorithmModule.f90 ProgramModule.f90 TestModule.f90 TestProgram.f90 
