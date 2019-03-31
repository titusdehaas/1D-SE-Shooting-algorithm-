module ThreePointSolverModule
   use Diagonalization 
   implicit none 
   save
   
contains 
   subroutine ThreePointSolver(gridsize, mesh, potential, eigenvalues, eigenvectors)
      !subroutine that uses three point scheme to solve differential equation and given eigenvectors and eigenvalues as output. 
      !requires specified gridsize, meshspace and an array of the same dimensions with nummerical values for a given potential.  
      real(8), intent(in) :: potential(:), mesh 
      integer, intent(in) :: gridsize
      real(8), allocatable, intent(out):: eigenvalues(:)
      real(8), allocatable, optional, intent (out) :: eigenvectors(:,:) !eigenvectors is optional argument 
      real(8), allocatable :: s(:,:), v(:,:), l(:,:), optional_eigenvectors(:,:) ! l = -1/2s+v, the matrix to be diagonalised.   
      integer :: i 
      
      
      ! Matrix sizes are allocated 
      allocate(s(gridsize, gridsize), v(gridsize, gridsize), l(gridsize, gridsize), eigenvalues(gridsize),&
         &optional_eigenvectors(gridsize, gridsize))
      if(present(eigenvectors)) then 
         allocate(eigenvectors(gridsize, gridsize)) 
      endif 
       
      s = 0.0d0 
      v = 0.0d0 


      !genrating matrix triple diagonal matrix 
      do i = 1, gridsize
         s(i,i) = -2.0d0 
      end do 

      do i = 2, gridsize
         s(i, i-1) = 1.0d0 
      end do 

      do i = 1, gridsize-1
         s(i, i+1) = 1.0d0 
      end do  
      
      !generating diagonal matrix with potential values 
      do i = 1, gridsize 
         v(i, i) = potential(i) 
      end do
      
      ! addition of the two matrices
      l = (-0.5d0/(mesh**2))*s+v
       
      !diagonalisation of the matrix, obtaining the eigenvalues and eigenvectors in a one and two dimensional array respectively
      call diagonalize(l, optional_eigenvectors, eigenvalues) 
      if(present(eigenvectors)) eigenvectors = optional_eigenvectors
      
   end subroutine  
end module 

      
      
       
      
         

