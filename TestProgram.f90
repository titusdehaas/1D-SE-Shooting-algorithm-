program TestProgram 
   use ShootingTestModule
   use ShootingModule  
   implicit none 
   save 
   type (Calculation) :: test 
   call TestNewCalc(test) 
end program 

