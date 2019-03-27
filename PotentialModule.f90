module PotentialModule

   implicit none 
   save 
   public NewGrid, NewPotential, Grid, Calculation, ParticleInABox, Gaussian  
   private  
   
   type Grid
      real(8) :: h
      integer :: N
   end type

   type Calculation
      type (grid) :: init_grid 
      type (grid) :: calculation_grid 
      real (8) :: convergence ! delta E at which convegergence is accepted 
      real (8) :: boundary_conditions(4) 
      integer :: number_solutions
      character(32)  :: potential_type ! speciefies potential type   

   end type

 

contains 
   subroutine NewGrid(N, h, grid) 
      real(8), intent(in) :: h 
      real(8), allocatable, intent(out):: grid(:)  
      integer :: i, N 

      allocate(grid(N))
      grid = 0 
      do i = 1, N 
         grid(i) = grid(i) + i*h 
      end do 
   end subroutine

   subroutine NewPotential(N,h, potential_type, potential_array, v0, alpha) 
      !subroutine that genrates an array containing the values of a speciefied potential for a given gridsize and meshsize 
      real(8), intent(in):: h ! meshsize
      integer, intent(in):: N ! Gridsize 
      character(32), intent(in) :: potential_type !pass type of potential as a string
      real(8), intent(in), optional:: v0, alpha
     
      real(8), allocatable, intent(out) :: potential_array(:)    
      if (potential_type == "particle in a box") then 
         call ParticleInAbox(N, h, potential_array) 
      elseif (potential_type == "gaussian") then 
         call Gaussian(N, h, V0, alpha, potential_array) 
      else 
         write(*,*)"No valid potential was entered"
      endif 

   end subroutine
   
   subroutine  ParticleInABox(N, h, box_array)
      real(8), intent(in) :: h
      integer, intent(in) :: N
      real(8), allocatable, intent(out) :: box_array(:) 
      call NewGrid(N, h, box_array) 
      box_array = 0.0d0       
   end subroutine   
   
   subroutine  Gaussian(N, h, v0, alpha, gaussian_array) 
      real(8), intent(in) :: h, v0, alpha 
      integer, intent(in) :: N 
      real(8), allocatable, intent(out) :: gaussian_array(:)
      integer :: i 
      call NewGrid(N, h, gaussian_array) 
      do i = 1, N 
         gaussian_array(i) = -v0*exp(-alpha*(gaussian_array(i))**2)
      end do  
   end subroutine  

   ! instead of using interface with functions, use if (string == ...) then 
   ! call ...() 
           

end module 
