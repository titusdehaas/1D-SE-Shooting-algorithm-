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
      type (grid) :: ThreePointGrid
      type (grid) :: ShootingGrid
      real (8) :: convergence
      integer :: potential !1=partilce in a box, 2 = gaussian etc.  
   end type

 

contains 
   subroutine NewGrid(N, h, Grid) 
      real(8), intent(in) :: h 
      real(8), allocatable, intent(out):: Grid(:)  
      integer :: i, N 

      allocate(Grid(N))
      Grid = 0 
      do i = 1, N 
         Grid(i) = Grid(i) + i*h 
      end do 
   end subroutine

   subroutine NewPotential(N,h, PotentialFunc, PotentialArray) 
      real(8), intent(in):: h 
      integer, intent(in):: N 
      real(8), allocatable, intent(out) :: potentialArray(:) 
      real(8), allocatable :: Grid(:) 
      integer :: i 
      interface 
         real(8) function PotentialFunc(x) 
            real(8), intent(in) :: x 
         end function 
      end interface 

      allocate(PotentialArray(N)) 
      call NewGrid(N, h, Grid) 
      do i = 1, N 
         potentialArray(i) = PotentialFunc(Grid(i)) 
      end do 
   end subroutine
   
   real(8) function ParticleInABox(x)
      real(8), intent(in) :: x 
      ParticleInABox = x 
   end function  
   
   real(8) function Gaussian(v0, alpha, x) 
      real(8), intent(in) :: v0, alpha, x 
      Gaussian = -v0*exp(-alpha*x**2) 
   end function 


           

end module 
