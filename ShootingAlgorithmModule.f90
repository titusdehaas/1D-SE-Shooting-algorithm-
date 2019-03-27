module ShootingAlgorithmModule
   use integration_module
   use ThreePointSolverModule 
   implicit none  
   save
   public ShootingAlgorithm 
   private 
   
contains 
   subroutine ShootingAlgorithm(gridsize, mesh, init_gridsize, init_mesh, convergence, potential, boundaryconditions, &
         &number_solutions, eigenvalues, eigenvectors) 

      ! Subroutine that calculates eigenvalues of a passed potential, using the shooting algorithm.
      ! It uses an  
      real(8), intent(in) :: mesh, init_mesh, convergence, potential(:), boundaryconditions(4) 
      integer, intent(in) :: init_gridsize, number_solutions   
      integer, intent(inout) :: gridsize 
      real(8), allocatable, intent(out):: eigenvalues(:), eigenvectors(:,:)
      real(8), allocatable :: lambda(:), inwards_array(:), outwards_array(:)  
      integer :: i, j, x_m, cyclecount, maxcyclecount   
      real(8) :: oldlambda, newlambda, deltalambda, normfactor
      real(8) :: derivative_inward, derivative_outward, integral_inward, integral_outward  
      real(8) :: difference_fractions, weighed_integral_in, weighed_integral_out 

      if (mod(gridsize, 2) == 1) then
         print*, "gridsize has the right dimension (odd) "
         gridsize = gridsize 
      else 
         gridsize = gridsize + 1
         print*, "Gridsize was increased to an odd number to minimise integration error"
      end if 

      maxcyclecount = 4 
      cyclecount = 1
      j = 1
      x_m = (gridsize+1)/2


      call ThreePointSolver(init_gridsize, init_mesh, potential, lambda)
      allocate(inwards_array(x_m+1), outwards_array(x_m+1)) ! allocation of the inwards and outwards solutions 
      allocate(eigenvalues(number_solutions), eigenvectors(number_solutions, gridsize)) ! alocation of the final output 

      newlambda = lambda(j)

      inwards_array = 0 
      outwards_array = 0 

      ! initialising boundary conditions
      inwards_array(1) = boundaryconditions(1)
      inwards_array(2) = boundaryconditions(2)
      outwards_array(gridsize-1) = boundaryconditions(3)
      outwards_array(gridsize) = boundaryconditions(4)

      do  ! set up convergence loop 
         
         oldlambda = newlambda

         ! calculation of the inwards and outwards solution at a common eigenvalue j. 
         do i = 2, x_m+1  
            inwards_array(i+1) = inwards_array(i-1)+(mesh**2)*(oldlambda - potential(i)+(2.0d0/(mesh**2.0d0)))*inwards_array(i) 
         end do 

         do i = gridsize-1, x_m-1, -1 
            outwards_array(i-1) = outwards_array(i+1)+(mesh**2)*(oldlambda - potential(i)+(2.0d0/(mesh**2.0d0)))*outwards_array(i)
         end do

         !testing 
         write(*,*) "mesh", mesh 
         write(*,*) "boundary conditions", boundaryconditions(:)
         write(*,*) "inwards_array(2)", inwards_array(2)
         write(*,*) " integeral_inward", integral_inward
         write(*,*) "oldlambda", oldlambda
         write(*,*) "deltalambda", deltalambda
         print*, inwards_array(:) 
         print*, outwards_array(:) 
         ! normalization of both solutions
         call newton_cotes(inwards_array, mesh, 1, gridsize, normfactor) 
         inwards_array = inwards_array/normfactor
         call newton_cotes(inwards_array, mesh, 1, gridsize, normfactor)  
         outwards_array = outwards_array/normfactor 

         ! determine integralvalue over  inwards and outwards functions
         call newton_cotes(inwards_array**2, mesh, 1, x_m, integral_inward) 
         call newton_cotes(outwards_array**2, mesh, x_m, gridsize, integral_outward) 

         

         ! determine derivatives of inwards and outwards solution 
         
         derivative_inward  = (inwards_array(x_m-1) + inwards_array(x_m + 1))/(2.0d0*mesh)
         derivative_outward  = (outwards_array(x_m - 1) + outwards_array(x_m + 1))/(2.0d0*mesh) 
         
         ! grouping terms together
         difference_fractions = (derivative_inward/inwards_array(x_m))-(derivative_outward/outwards_array(x_m))
         weighed_integral_in = (1.0d0/(inwards_array(x_m)**2.0d0))*integral_inward
         weighed_integral_out = (1.0d0/(outwards_array(x_m)**2.0d0))*integral_outward  
         
         ! calculation of the correction to the eigenvalue
         deltalambda = 0.5d0*(difference_fractions)/(weighed_integral_in*weighed_integral_out)
         
         ! calculating new eigenvalue and checking convergence 
         
         
         newlambda  = oldlambda - deltalambda 
         cyclecount = cyclecount + 1
         if (newlambda - oldlambda >= convergence) then 
            write(*,*) "shooting algorithm has converged in:", cyclecount, "cycles"
         end if  
         if (cyclecount >= maxcyclecount) then 
            write(*,*) "shooting algorithm does not converge within", cyclecount,"cycles"
            exit 
         end if 
         
      end do 
       
      eigenvalues(j) = newlambda
      eigenvectors(j, 1:x_m) = inwards_array ! eigenvalue array is constructed from inwards solutions as the first slice
      eigenvectors(j, x_m: gridsize) = outwards_array ! and the outwards solution as the second slice

       

   end subroutine 
end module  
