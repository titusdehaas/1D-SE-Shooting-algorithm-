module ShootingAlgorithmModule
   implicit none  
   save
   public 
   private 
   
contains 
   subroutine ShootingAlgorithm(N, h, convergence, potential, boundaryconditions, eigenvalues, eigenvectors) 
      real(8), intent(in) :: h, convergence, potential(:), bondaryconditions 
      integer, intent(in) :: N 
      real(8), allocatable, intent(out):: eigenvalues, eigenvectors 

