module ProgramModule 
   
   use PotentialModule
   use ShootingAlgorithmModule
   use ThreePointSolverModule   
   use WriterModule
   implicit none 
   save 
   public NewCalc, RunCalc, Calculation     
   private 

   ! a type that stores all data needed for a calculation
   type Calculation
      type (grid) :: init_grid
      type (grid) :: calculation_grid
      real (8) :: convergence                   !delta E at which convegergence is accepted 
      real (8) :: boundary_conditions(4)        !specifies the boundary conditions 
      integer :: number_solutions               !specifies the number of solutions
      character(32)  :: potential_type          !speciefies potential type   
      real(8) :: alpha, D, x_e                  !morse parameters: 
 
   end type
 

contains

   subroutine NewCalc(self)
      !subroutine that reads all input paramters from the input file and stores it in a calculation type
      type (Calculation), intent (inout) :: self
      real(8):: init_mesh, mesh, convergence, boundary_conditions(4), alpha, D, x_e  
      integer:: init_gridsize, gridsize, number_solutions  
      character(32) :: dummyline, potential_type
      type(OutputFile):: log_type 
      
      !preparing logfile program output
      call NewOutputFile(10, "logfile.txt", log_type) ! this empties the logfile for a new calculation 
      call EmptyOutputfile(log_type) 
      write(*,*) "========================================================"


      ! reads the input parameters
      open (7, file = "input.txt") 
      read(7,*) dummyline
      read(7,*) dummyline  
      read(7,*) init_gridsize, init_mesh, gridsize, mesh, convergence, potential_type, D, alpha, x_e  
      read(7,*) dummyline
      read(7,*) number_solutions
      close(7)

      ! writes the input parameters
      self%init_grid%gridsize = init_gridsize 
      self%init_grid%mesh = init_mesh
      self%calculation_grid%gridsize = gridsize 
      self%calculation_grid%mesh = mesh 
      self%potential_type  = potential_type  
      self%convergence = convergence 
      self%number_solutions = number_solutions 
      self%alpha = alpha 
      self%D = D
      self%x_e = x_e       
       
      ! logging 
      call Logger(log_type, "=======================================================")
      call Logger(log_type,"coarse mesh size was set to ", self%init_grid%mesh)
      call Logger(log_type, "coarse grid size was set to:", self%init_grid%gridsize)
      Call Logger(log_type,"Shooting mesh size was set to:", self%calculation_grid%mesh)
      call Logger(log_type,"Shooting grid size was set to ", self%calculation_grid%gridsize)
      call Logger(log_type,"Convergence accuracy was set to:", self%convergence)
      call Logger(log_type,"Number of solutions was set to:", self%number_solutions)
      call Logger(log_type,"Initialising for:", self%potential_type)
      call Logger(log_type, "=======================================================")


   end subroutine 

   subroutine RunCalc(self)
      ! subroutine that performs the calculation. The input parameters are passed in the calculation data type. 
      ! the subroutine calls both the 3 point scheme and shooting method to calculate the eigenvalues and eigenvectors for the
      ! specified input data. 

      type (Calculation), intent(inout) :: self
      type (OutputFile) :: output  
      real(8), allocatable :: potential_array(:), eigenvalues(:), eigenvectors(:,:)
      character(32) :: potential_type
      integer :: i, gridsize, init_gridsize, number_solutions 
      real(8) :: mesh, init_mesh, convergence, alpha, D, x_e
      

      !reading input parameters 
      gridsize = self%calculation_grid%gridsize
      mesh = self%calculation_grid%mesh
      init_gridsize = self%init_grid%gridsize
      init_mesh = self%init_grid%mesh
      convergence = self%convergence
      number_solutions = self%number_solutions 
      alpha = self%alpha
      D = self%D
      x_e = self%x_e
      potential_type = self%potential_type
      
      !generating potential grid 
      Write(*,*) "------ Setting up new potential-------------------------" 
      call NewPotential(gridsize, mesh, potential_type, potential_array, alpha, D, x_e)
      
      !performing shooting algorithm 
      write(*,*) "------ Performing shooting method-----------------------" 
      call ShootingAlgorithm(gridsize, mesh, init_gridsize, init_mesh, convergence, potential_array, &
         &number_solutions, eigenvalues, eigenvectors)
      
      !writing output to files
      call NewOutputFile(11, "eigenvectors_shooting.out", output)
      call Write(output, mesh, number_solutions, eigenvectors)
      call NewOutputFile(12, "eigenvalues_shooting.out", output)
      call Write(output, eigenvalues)
      write(*,*) "------Shooting algorithm calculation was completed------"

      !performing three-point scheme calculation 
      write(*,*) "------Performing Three-Point Scheme Calculation---------"
      call ThreePointSolver(gridsize, mesh, potential_array, eigenvalues, eigenvectors)
      !writing output to files
      call NewOutputFile(7, "eigenvectors_threepoint.out", output)
      call Write(output, mesh, number_solutions, eigenvectors)
      call NewOutputFile(8, "eigenvalues_threepoint.out", output)
      call Write(output, eigenvalues)
      write(*,*) "------Three-point scheme calculation completed----------"
      write(*,*) "========================================================"
   end subroutine 
   



end module 


