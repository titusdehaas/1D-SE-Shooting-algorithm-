module ProgramModule 
   
   use PotentialModule
   use ThreePointSolverModule   
   implicit none 
   save 
   public NewCalc   
   private 

   ! perhaps better to store the interface and convergence trashold as optional arguments in the "grid type", and alter the name to
   ! calc type. 

contains

   subroutine NewCalc(self)
      type (Calculation), intent (inout) :: self
      real(8):: h, h2, convergence 
      integer:: v, N, N2 
      character(32) :: firstline, secondline 

      open (7, file = "input.txt") 
      read(7,*) firstline
      read(7,*) secondline  
      read(7,*) N, h, N2, h2, v,  convergence
      close(7)

      self%ThreePointGrid%N = N
      self%ThreePointGrid%h = h
      self%ShootingGrid%N = N2
      self%ShootingGrid%h = h2
      self%potential = v 
      self%convergence = convergence 

      ! this can be writen in a quicker way once done 
   end subroutine 

   subroutine RunCalc(self) 
      type (Calculation), intent (inout) :: self 
      real(8)::  h
      real(8), allocatable :: potential(:), eigenvalues(:), eigenvectors(:,:) 
      integer:: N 

      N = self%ThreePointGrid%N 
      h = self%ThreePointGrid%h

      Call NewGrid(N, h, Potential) 
      call ThreePointSolver(N, h, potential, eigenvalues, eigenvectors) 


!      call ShooterSolver(self, eigenvaluesshooter, eigenvctorsshooter) 
!      
!      call writer(eigenvalues, eigenvectors) 
!      call writer(eigenvaluesshooter, eigenvectorsshooter)  
   end subroutine 
   



end module 


