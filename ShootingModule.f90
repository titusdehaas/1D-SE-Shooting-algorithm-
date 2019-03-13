module ShootingModule 
   
   use BuildGridModule  
   implicit none 
   save 
   public NewCalc, Calculation, Grid  
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

end module 


