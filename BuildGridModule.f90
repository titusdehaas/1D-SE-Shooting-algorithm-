module BuildGridModule

   implicit none 
   save 
   public NewGrid 
   private  
 

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

end module 
