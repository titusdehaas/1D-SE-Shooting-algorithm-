module TestModule 
   use PotentialModule 
   use ProgramModule
   use ThreePointSolverModule
   use ShootingAlgorithmModule
   use WriterModule 
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
      write(*,*) "you set boundary conditions to:", self%boundary_conditions 
      write(*,*) "you set the number of solutions to:", self%number_solutions 
   end subroutine

   subroutine TestPotential(self) 
      type (Calculation), intent(inout) :: self
      real(8), allocatable :: potential_array(:) 
      integer :: i 
      character(32) :: potential_type 

      write(*,*) "------ Testing potential ------"
      potential_type = self%potential_type

      call NewPotential(self%init_grid%N, self%init_grid%h, potential_type, potential_array, self%alpha, self%D, self%x_e) 

      
      write(*,*) "first ten potential values"
      do i = 1, 10
         write(*,*) potential_array(i)
      end do 
   end subroutine  

   subroutine TestThreePointSolver(self) 
      type (Calculation), intent(inout) :: self
      type (OutputFile) :: output 
      real(8), allocatable :: potential_array(:), eigenvalues(:), eigenvectors(:,:), analytical_solutions(:) 
      character(32) :: potential_type
      integer :: i, gridsize 
      real(8) :: mesh, alpha, D, x_e
      real(8), parameter :: treshold = 0.05
      logical :: error = .false.
      
      gridsize = self%calculation_grid%N
      mesh = self%calculation_grid%h
      alpha = self%alpha
      D = self%D
      x_e = self%x_e
      potential_type = self%potential_type
      
      write(*,*) "------ Testing the 3-point scheme algorithm ------"
      call NewPotential(gridsize, mesh, potential_type, potential_array, alpha, D, x_e)
      call ThreePointSolver(self%calculation_grid%N, self%calculation_grid%h, potential_array, eigenvalues, eigenvectors)
      write(*,*) "three point scheme first:", self%number_solutions, "solutions"

      do i = 1, self%number_solutions  
         write(*,*) eigenvalues(i) 
      end do

      call AnalyticalSolutionPIAB(self, analytical_solutions)
      do i = 1, self%number_solutions 
         if (abs(analytical_solutions(i)) - abs(eigenvalues(i)) > treshold) then
            print*, "Error detected in three point scheme calculation: solution", i
            error = .true. 
         end if  
      end do

            !writing output to files
      call NewOutputFile(7, "EIGENVECTORS_3PS.out", output)
      call Write(output, self%calculation_grid%h, self%number_solutions, eigenvectors)
      call NewOutputFile(8, "EIGENVALUES_3PS.out", output)
      call Write(output, eigenvalues)

      if (error) then
         write(*,*) "------Three point scheme calculation was has an error------"
         write(*,*) "-----------------------------------------------------------"
      else
         write(*,*) "------Three point scheme calculation was succesfull------"
         write(*,*) "---------------------------------------------------------"
      endif  
      
      !writing output to files
      call NewOutputFile(7, "EIGENVECTORS_3PS.out", output)
      call Write(output, self%calculation_grid%h, self%number_solutions, eigenvectors)  
      call NewOutputFile(8, "EIGENVALUES_3PS.out", output) 
      call Write(output, eigenvalues)
       
   end subroutine 
   
   subroutine TestShootingAlgorithm(self) 
      type(calculation), intent(inout) :: self 
      type(OutputFile) :: output 
      real(8), allocatable :: potential_array(:), eigenvalues(:), eigenvectors(:,:)  
      character(32) :: potential_type 
      integer :: i, gridsize 
      real(8) :: mesh, alpha, D, x_e 
      gridsize = self%calculation_grid%N 
      mesh = self%calculation_grid%h 
      alpha = self%alpha 
      D = self%D
      x_e = self%x_e 
      potential_type = self%potential_type 


      write(*,*) "------ Testing shooting algorithm ------"
      call NewPotential(gridsize, mesh, potential_type, potential_array, alpha, D, x_e) 
      
      call ShootingAlgorithm(self%calculation_grid%N, self%calculation_grid%h, self%init_grid%N, self%init_grid%h, & 
         &self%convergence, potential_array, self%boundary_conditions, self%number_solutions, eigenvalues, eigenvectors) 
     
      write(*,*) "Shooting algorithm first:", self%number_solutions," eigenvalues:"
      do i = 1, self%number_solutions
         write(*,*) eigenvalues(i)
      end do

            !writing output to files
      call NewOutputFile(9, "eigenvectors_shooting.out", output)
      call Write(output, self%calculation_grid%h, self%number_solutions, eigenvectors)  
      call NewOutputFile(10, "eigenvalues_shooting.out", output)      
      call Write(output, eigenvalues) 
      write(*,*) "------Shooting algorithm calculation was completed------"
 
      
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

      write(*,*) "analytical partice in a box first:", self%number_solutions, "solutions"
      do i = 1, self%number_solutions
         energy(i) = (i**2)*(pi**2)*(1.0d0/(2.0d0*m*(l**2)))
      end do 
      
      do i = 1, self%number_solutions
         write(*,*) energy(i)
      end do 
   end subroutine 

end module 
