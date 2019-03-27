module TestModule 
   use PotentialModule 
   use ProgramModule
   use ThreePointSolverModule
   use ShootingAlgorithmModule 
   implicit none 
   save 
   public TestNewCalc, TestThreePointSolver, TestPotential, TestShootingAlgorithm, AnalyticalSolutionPIAB
   private 

contains 
   subroutine TestNewCalc(self) 
      type (calculation), intent(inout) :: self
      call NewCalc(self) 

      write(*,*) "------ Testing the New Calculation set up and parameters ------"
      Write(*,*) "you set three point mesh size to", self%init_grid%N 
      write(*,*) "you set three point mesh size to:", self%init_grid%h
      write(*,*) "you set Shooting grid size to:", self%calculation_grid%N
      write(*,*) "you set Shooting grid size to:", self%calculation_grid%h
      write(*,*) "you set convergence accuracy to:", self%convergence
      write(*,*) "initialising for", self%potential_type, "type of potential"
   end subroutine

   subroutine TestPotential(self) 
      type (Calculation), intent(inout) :: self
      real(8), allocatable :: potential_array(:) 
      integer :: i 
      character(32) :: potential_type 

      write(*,*) "------ Testing potential ------"
      potential_type = "particle in a box"

      call NewPotential(self%init_grid%N, self%init_grid%h, potential_type, potential_array) 

      
      write(*,*) "first ten potential values"
      do i = 1, 10 
         write(*,*) potential_array(i)
      end do 
   end subroutine  

   subroutine TestThreePointSolver(self) 
      type (Calculation), intent(inout) :: self
      real(8), allocatable :: potential_array(:), eigenvalues(:), eigenvectors(:,:), analytical_solutions(:) 
      character(32) :: potential_type
      integer :: i 
      real(8), parameter :: treshold = 0.05
    
      potential_type = "particle in a box"

      write(*,*) "------ Testing the 3-point scheme algorithm ------"
      call NewPotential(self%calculation_grid%N, self%calculation_grid%h, potential_type, potential_array) 

      print*, "the first ten potential values:"
      do i = 1, 10 
         write(*,*) potential_array(i)
      end do 

      call ThreePointSolver(self%calculation_grid%N, self%calculation_grid%h, potential_array, eigenvalues, eigenvectors)
      write(*,*) "these are the first ten eigenvalues:"
      do i = 1, 20 
         write(*,*) eigenvalues(i) 
      end do

      call AnalyticalSolutionPIAB(self, analytical_solutions)
      do i = 1, 10
         if (analytical_solutions(i) - eigenvalues(i) > treshold) then
            print*, "Error detected in three point scheme calculation"
         else
            write(*,*) "three point scheme calculation succesfully finished"
         end if
      end do


   end subroutine 
   
   subroutine TestShootingAlgorithm(self) 
      type(calculation), intent(inout) :: self 
      real(8), allocatable :: potential_array(:), eigenvalues(:), eigenvectors(:,:)  
      character(32) :: potential_type 
      integer :: i 
      potential_type = "particle in a box"

      write(*,*) "------ Testing shooting algorithm ------"
      call NewPotential(self%calculation_grid%N, self%calculation_grid%h, potential_type, potential_array) 
      
      call ShootingAlgorithm(self%calculation_grid%N, self%calculation_grid%h, self%init_grid%N, self%init_grid%h, & 
         &self%convergence, potential_array, self%boundary_conditions, self%number_solutions, eigenvalues, eigenvectors) 
     
      write(*,*) "Shooting algorithm first ten eigenvalues:"
      do i = 1, 10
         write(*,*) eigenvalues(i)
      end do 
      
   end subroutine

   subroutine AnalyticalSolutionPIAB(self, energy) 
      type(calculation), intent(inout) :: self 
      real(8), parameter :: pi = 3.14159265359d0, h = 6.626070150e-34, m = 1.0d0  
      real(8), allocatable, intent(out)  :: energy(:) 
      real(8) :: l, mesh
      integer :: gridsize, number_solutions, i 
     
      gridsize = self%calculation_grid%N 
      mesh =  self%calculation_grid%h 
      number_solutions = self%number_solutions 
      l = mesh* gridsize 
      allocate(energy(number_solutions)) 

      write(*,*) "first ", number_solutions, "analytical solutions for a partice in a box:"
      do i = 1, number_solutions 
         energy(i) = (i**2)*(pi**2)*(1.0d0/(1.0d0*m*(l**2)))
      end do 
      
      do i = 1, number_solutions 
         write(*,*) energy(i)
      end do 
   end subroutine 

end module 
