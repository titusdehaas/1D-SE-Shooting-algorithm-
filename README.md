# 1D-SE-Shooting-algorithm-
In this repository a 1D schreudinger equation shooting algorithm is developed for the course "Scientific Programming".

compiling order:
gfortran blas.f lapack.f diagonalization.F90 trapez.f90 PotentialModule.f90 WriterModule.f90 ThreePointSolverModule.f90 ShootingAlgorithmModule.f90 ProgramModule.f90 TestModule.f90 ShootingTwoEnds.f90

after compiling the program execurted by running the executable file (a.out if no name is specified). Output can be found in the
files: eigenvectors_shooting.out, eigenvalues_shooting.out, eigenvectors_threepoint.out, eigenvalues_threepoint.out. 
The eigenvectors_shooting.out and eigenvectors_threepoint.out contain the eigenvector array values in columns. These can be plotted 
using, for example, Gnuplot. 
Also a log file is genrated which contains more information about each cycle in the shooting algorithm.

the main program can be found in the ShootingTwoEnds.f90 file. There two subroutines are called, which are defined in the
ProgramModule.f90 file. One initialises the calculation, the other one runs the calculation. 
blas.f lapack.f diagonalization.F90 trapez.f90: where supplied by the course coordinator and contains subroutines for diagonalising
matrices and integrating.
ThreePointSolverModule.f90 contains the subroutine for performing the three-point scheme calculation 
ShootingAlgorithmModule.f90 contains the subroutine for performing the shooting method calculation
TestModule.f90 contains two unit tests for the three-point scheme subroutine and the shooting method subroutine
WriterModule.f90 contains the subroutines that are used for writing the output files and the logfile

###################################################

IMPORTANT NOTE: 
Unfortunately the program still contains a bug. The shooting method does not (yet) produce the right eigenvalues. Somehow, the even
values of n do not converge to an eigenvalue. Also, the produced eigenvectors are clearly not right. For example, the eigenvector 
for n = 1 has 10 nodes, and not contiues. Also its starts aimed in the wrong direction. I think I've made and error somewhere in
the formulas used in the shooting method. I've actually spend a lot of time trying to find the bug and im a bit sorrow that I was
not able to find it. Therefore I was wondering if you maybe could have one final look at it before I officially hand it in? If not,
it's also fine. In that case this is the final version. 

Kind Regards,
Titus 

##################################################

