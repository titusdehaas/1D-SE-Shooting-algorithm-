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
      real(8):: init_mesh, mesh, convergence, boundary_conditions(4), alpha, D, x_e  
      integer:: init_gridsize, gridsize, number_solutions  
      character(32) :: dummyline, potential_type 

      open (7, file = "input.txt") 
      read(7,*) dummyline
      read(7,*) dummyline  
      read(7,*) init_gridsize, init_mesh, gridsize, mesh, convergence, potential_type, D, alpha, x_e  
      read(7,*) dummyline
      read(7,*) boundary_conditions(1), boundary_conditions(2), boundary_conditions(3), boundary_conditions(4), number_solutions
      close(7)

      self%init_grid%N = init_gridsize 
      self%init_grid%h = init_mesh
      self%calculation_grid%N = gridsize 
      self%calculation_grid%h = mesh 
      self%potential_type  = potential_type  
      self%convergence = convergence 
      self%boundary_conditions = boundary_conditions 
      self%number_solutions = number_solutions 
      self%alpha = alpha 
      self%D = D
      self%x_e = x_e  

      ! this can be writen in a quicker way once done 
   end subroutine 

   subroutine RunCalc(self) 
      type (Calculation), intent (inout) :: self 
      real(8)::  h
      real(8), allocatable :: potential(:), eigenvalues(:), eigenvectors(:,:) 
      integer:: N 

      N = self%init_grid%N 
      h = self%init_grid%h

      Call NewGrid(N, h, Potential) 
      call ThreePointSolver(N, h, potential, eigenvalues, eigenvectors) 


!      call ShooterSolver(self, eigenvaluesshooter, eigenvctorsshooter) 
!      
!      call writer(eigenvalues, eigenvectors) 
!      call writer(eigenvaluesshooter, eigenvectorsshooter)  
   end subroutine 
   



end module 


