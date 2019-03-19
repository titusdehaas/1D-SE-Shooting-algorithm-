module ThreePointSolverModule
   use Diagonalization 
   implicit none 
   save
   
contains 
   subroutine ThreePointSolver(N, h, potential, eigenvalues, eigenvectors)
      !declaration of variables
      real(8), intent(in) :: potential(:), h 
      integer, intent(in) :: N 
      real(8), allocatable, intent(out):: eigenvalues(:), eigenvectors(:,:)
      real(8), allocatable :: s(:,:), v(:,:), l(:,:)  
      integer :: i 
      
      ! Matrix sizes are allocated 
      allocate(s(N,N), v(N,N), l(N,N), eigenvalues(N), eigenvectors(N,N)) 
      s = 0
      v = 0
      
      ! setup the values not reached by the do loop; perhaps this can be done more elegantly
      s(1,1) = -2.0d0 
      s(1,2) = 1.0d0
      s(2,1) = 1.0d0 
      s(N,N) = -2.0d0 
      s(1,N) = 1.0d0 
      s(N,1) = 1.0d0 
      s(N, N-1) = 1.0d0 
      
      ! setting up the tripple diagonal matrix 
      do i = 2, N-1 
         s(i,i) = -2.0d0  
         s(i,i-1) = 1.0d0  
         s(i,i+1) = 1.0d0
      end do

      do i = 1, N 
         v(i, i) = potential(i) 
      end do

      l = (1/h**2)*(s+v)
      call diagonalize(l, eigenvectors, eigenvalues) 
   end subroutine  
end module 

      
      
       
      
         

