program ShootingTwoEnds 
   use TestModule
   use ProgramModule
   use PotentialModule
   use ProgramModule
   implicit none 
   save 
   type (Calculation) :: calc
   
   ! Running the calculations
   call NewCalc(calc) !reads parameters and initialises log file 
   call RunCalc(calc) !performs the calculations 
  
   !============================================================
   ! Following can be uncommented for testing 
   !call TestShootingAlgorithm(calc) 
   !call TestThreePointSolver(calc)  
   !===========================================================
end program 


