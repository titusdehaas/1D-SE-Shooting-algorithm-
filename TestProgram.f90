program TestProgram 
   use TestModule
   use PotentialModule
   implicit none 
   save 
   type (Calculation) :: test 
   call TestNewCalc(test)
   call TestPotential(test)  
   call TestShootingAlgorithm(test) 
   call TestThreePointSolver(test)  
end program 


