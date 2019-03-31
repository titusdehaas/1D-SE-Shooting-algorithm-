module TestModule 
   use PotentialModule 
   use ProgramModule
   use ThreePointSolverModule
   use ShootingAlgorithmModule
   use WriterModule 
   implicit none 
   save 
   public TestThreePointSolver, TestShootingAlgorithm, AnalyticalSolutionPIAB
   private 

contains

   subroutine TestThreePointSolver(self)
      ! subroutine that tests the three point scheme by comparing the eigenvalues obtained with the analytical solution  
      type (Calculation), intent(inout) :: self
      type (OutputFile) :: output 
      real(8), allocatable :: potential_array(:), eigenvalues(:), eigenvectors(:,:), analytical_solutions(:) 
      character(32) :: potential_type
      integer :: i, gridsize 
      real(8) :: mesh, alpha, D, x_e
      real(8), parameter :: treshold = 0.05
      logical :: error = .false.
      
      gridsize = self%calculation_grid%gridsize
      mesh = self%calculation_grid%mesh 
      potential_type = "Particle in a box"
      
      !testing the three point scheme eigenvalues by comparing them to anaylitcal solutions 
      write(*,*) "------ Testing the 3-point scheme algorithm ------"
      call NewPotential(gridsize, mesh, potential_type, potential_array)
      call ThreePointSolver(self%calculation_grid%gridsize, self%calculation_grid%mesh, potential_array, eigenvalues)
      call AnalyticalSolutionPIAB(self, analytical_solutions)
      do i = 1, self%number_solutions 
         if (abs(analytical_solutions(i)) - abs(eigenvalues(i)) > treshold) then
            print*, "Error detected in three point scheme calculation: solution", i
            error = .true. 
         end if  
      end do

      if (error) then
         write(*,*) "------An error occured during three-point scheme calc.-----"
         write(*,*) "-----------------------------------------------------------"
      else
         write(*,*) "------Three point scheme calculation was succesfull------"
         write(*,*) "---------------------------------------------------------"
      endif  
      
      !writing output to file       
   end subroutine 
   
   subroutine TestShootingAlgorithm(self) 
      !subroutine that test the shooting algorithm by comaring the obtained eigenvalues for a particle in a box with the analytical
      !solution 
      type(calculation), intent(inout) :: self 
      type(OutputFile) :: output 
      real(8), allocatable :: potential_array(:), eigenvalues(:), eigenvectors(:,:), analytical_solutions(:)   
      character(32) :: potential_type 
      integer :: i, gridsize, number_solutions 
      real(8) :: convergence, mesh
      real(8), parameter :: treshold = 0.05 
      logical :: error = .false. 

      !reading input 
      gridsize = self%calculation_grid%gridsize 
      mesh = self%calculation_grid%mesh 
      convergence = self%convergence
      number_solutions = self%number_solutions
      potential_type = "Particle in a box"

      ! calling the shooting algorithm and comparing optained eigenvalues with analytical solution
      write(*,*) "------ Testing shooting algorithm ------"
      call NewPotential(gridsize, mesh, potential_type, potential_array) 
      
      call ShootingAlgorithm(gridsize, mesh, gridsize, mesh, convergence,potential_array, &
              &number_solutions, eigenvalues, eigenvectors) 
      call AnalyticalSolutionPIAB(self, analytical_solutions)
      do i = 1, self%number_solutions
         if (abs(analytical_solutions(i)) - abs(eigenvalues(i)) > treshold) then
            print*, "Error detected in three point scheme calculation: solution", i
            error = .true.
         end if  
      end do

      if (error) then
         write(*,*) "------An error occured during shootin method calcculation--"
         write(*,*) "-----------------------------------------------------------"
      else
         write(*,*) "------shooting method calculation was succesfull---------"
         write(*,*) "---------------------------------------------------------"
      endif


   end subroutine

   subroutine AnalyticalSolutionPIAB(self, energy) 
      !subroutine that calculates the analytical solutions for a particle in a box (in atomic units) 
      type(calculation), intent(inout) :: self 
      real(8), parameter :: pi = 3.14159265359d0
      real(8), allocatable, intent(out)  :: energy(:) 
      real(8) :: l, mesh
      integer :: gridsize, number_solutions, i 
      
      !reading input 
      gridsize = self%calculation_grid%gridsize 
      mesh =  self%calculation_grid%mesh 
      number_solutions = self%number_solutions 
      l = mesh* gridsize 
      allocate(energy(number_solutions)) 

      write(*,*) "analytical partice in a box first:", self%number_solutions, "solutions"
      do i = 1, self%number_solutions
         energy(i) = (i**2)*(pi**2)*(1.0d0/(2.0d0*(l**2))) !analytical solution 
      end do 
      
      do i = 1, self%number_solutions
         write(*,*) energy(i)
      end do 
   end subroutine 

end module 
