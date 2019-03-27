module ThreePointSolverModule
   use Diagonalization 
   implicit none 
   save
   
contains 
   subroutine ThreePointSolver(gridsize, mesh, potential, eigenvalues, eigenvectors)
      !declaration of variables
      real(8), intent(in) :: potential(:), mesh 
      integer, intent(in) :: gridsize
      real(8), allocatable, intent(out):: eigenvalues(:)
      real(8), allocatable, optional, intent (out) :: eigenvectors(:,:)
      real(8), allocatable :: s(:,:), v(:,:), l(:,:), optional_eigenvectors(:,:)   
      integer :: i 
      
      
      ! Matrix sizes are allocated 
      allocate(s(gridsize, gridsize), v(gridsize, gridsize), l(gridsize, gridsize), eigenvalues(gridsize),&
         &optional_eigenvectors(gridsize, gridsize))
      if(present(eigenvectors)) then 
         allocate(eigenvectors(gridsize, gridsize)) 
      endif 
       
      s = 0.0d0 
      v = 0.0d0 



      do i = 1, gridsize
         s(i,i) = -2.0d0 
      end do 

      do i = 2, gridsize
         s(i, i-1) = 1.0d0 
      end do 

      do i = 1, gridsize-1
         s(i, i+1) = 1.0d0 
      end do  
      
      ! setting up the tripple diagonal matrix 

      do i = 1, gridsize 
         v(i, i) = potential(i) 
      end do

      l = (-1.0d0/mesh**2)*s+v
       
       

      call diagonalize(l, optional_eigenvectors, eigenvalues) 
      if(present(eigenvectors)) eigenvectors = optional_eigenvectors 
   end subroutine  
end module 

      
      
       
      
         

