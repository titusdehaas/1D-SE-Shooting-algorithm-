module TestModule 
   use PotentialModule 
   use ProgramModule
   use ThreePointSolverModule 
   implicit none 
   save 
   public TestNewCalc, TestThreePointSolver, TestPotential
   private 

contains 
   subroutine TestNewCalc(self) 
      type (calculation), intent(inout) :: self
      call NewCalc(self) 

      write(*,*) "------ Testing the New Calculation set up and parameters ------"
      Write(*,*) "you set three point mesh size to", self%ThreePointGrid%N 
      write(*,*) "you set three point mesh size to:", self%ThreePointGrid%h
      write(*,*) "you set Shooting grid size to:", self%ShootingGrid%N
      write(*,*) "you set Shooting grid size to:", self%ShootingGrid%h
      write(*,*) "you set convergence accuracy to:", self%convergence
      
   end subroutine

   subroutine TestPotential(self) 
      type (Calculation), intent(inout) :: self
      real(8), allocatable :: PotentialArray(:) 
      integer :: i  
      write(*,*) "------ Testing potential ------"
  
      call NewPotential(self%ThreePointGrid%N, self%ThreePointGrid%h, ParticleInABox, PotentialArray) 

      
      write(*,*) "first ten potential values"
      do i = 1, 10 
         write(*,*) PotentialArray(i)
      end do 
   end subroutine  

   subroutine TestThreePointSolver(self) 
      type (Calculation), intent(inout) :: self
      real(8), allocatable :: potential(:), eigenvalues(:), eigenvectors(:,:) 
      integer :: i 

      write(*,*) "------ Testing the 3-point scheme algorithm ------"
      call NewPotential(self%ThreePointGrid%N, self%ThreePointGrid%h, ParticleInABox, Potential) 

      print*, "the first ten potential values:"
      do i = 1, 10 
         write(*,*) potential(i)
      end do 

      call ThreePointSolver(self%ThreePointGrid%N, self%ThreePointGrid%h,&
      &Potential, eigenvalues, eigenvectors)
      write(*,*) "these are the first ten eigenvalues:"
      do i = 1, 10 
         write(*,*) eigenvalues(i) 
      end do

   end subroutine 


end module 
