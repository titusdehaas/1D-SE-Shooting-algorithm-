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
      real(8) :: oldlambda, newlambda, deltalambda, normfactor, normfactor2
      real(8) :: derivative_inward, derivative_outward, integral_inward, integral_outward  
      real(8) :: difference_fractions, weighed_integral_in, weighed_integral_out 

      if (mod(gridsize, 2) == 1) then
         print*, "gridsize has the right dimension (odd) "
         gridsize = gridsize 
      else 
         gridsize = gridsize + 1
         print*, "Gridsize was increased to an odd number to minimise integration error"
      end if 

      maxcyclecount = 3
      cyclecount = 0
      j = 0
      x_m = (gridsize+1)/2


      call ThreePointSolver(init_gridsize, init_mesh, potential, lambda)
      allocate(inwards_array(x_m+1), outwards_array(x_m+1)) ! allocation of the inwards and outwards solutions 
      allocate(eigenvalues(number_solutions), eigenvectors(number_solutions, gridsize)) ! alocation of the final output
      inwards_array = 0 
      outwards_array = 0 
      eigenvalues = 0 
      eigenvectors = 0 

      ! initialising boundary conditions
      inwards_array(1) = boundaryconditions(1)
      inwards_array(2) = boundaryconditions(2)
      outwards_array(x_m) = boundaryconditions(3)
      outwards_array(x_m+1) = boundaryconditions(4)
      
      ! find the jth eigenvalue, for 1 to the number of solutions requested
      do j = 1, number_solutions

      newlambda = lambda(j) 

      do  ! set up convergence loop 
         
         oldlambda = newlambda
          
         ! calculation of the inwards and outwards solution at a common eigenvalue j. 
         do i = 2, x_m  
            inwards_array(i+1) = -inwards_array(i-1)+(mesh**2)*(oldlambda - potential(i)+(2.0d0/(mesh**2.0d0)))*inwards_array(i) 
         end do 

         do i = x_m, 2, -1 
            outwards_array(i-1) = -outwards_array(i+1)+(mesh**2)*(oldlambda - potential(i)+(2.0d0/(mesh**2.0d0)))*outwards_array(i)
         end do
          
         ! normalization of both solutions
         call newton_cotes(inwards_array, mesh, 1, x_m+1, normfactor) 
         inwards_array = inwards_array/normfactor
         call newton_cotes(outwards_array, mesh, 1, x_m+1, normfactor2)  
         outwards_array = outwards_array/normfactor2 

         ! determine integralvalue over  inwards and outwards functions
         call newton_cotes(inwards_array**2, mesh, 1, x_m+1, integral_inward) 
         call newton_cotes(outwards_array**2, mesh, 1, x_m+1, integral_outward) 

         

         ! determine derivatives of inwards and outwards solution 
         
         derivative_inward  = (inwards_array(x_m-1) + inwards_array(x_m+1))/(2.0d0*mesh)
         derivative_outward  = (outwards_array(1) + outwards_array(3))/(2.0d0*mesh) 
         
         ! grouping terms together
         difference_fractions = (derivative_inward/inwards_array(x_m))-(derivative_outward/outwards_array(2))
         weighed_integral_in = (1.0d0/(inwards_array(x_m)**2.0d0))*integral_inward
         weighed_integral_out = (1.0d0/(outwards_array(2)**2.0d0))*integral_outward  
         
         ! calculation of the correction to the eigenvalue
         deltalambda = 0.5d0*(difference_fractions)/(weighed_integral_in*weighed_integral_out)
         
         ! calculating new eigenvalue and checking convergence
         newlambda  = oldlambda - deltalambda

         !testing 

         write(*,*) "-------------- solution:", j, "cycle:",cyclecount,"-----------------------------"
         write(*,*) "mesh:", mesh
         write(*,*) "boundary conditions:", boundaryconditions(:)
         write(*,*) "normfactor1:", normfactor
         write(*,*) "normfactor2:", normfactor2 
         write(*,*) "integeral_inward:", integral_inward
         write(*,*) "integeral_outward:", integral_outward
         write(*,*) "derivative inwards:", derivative_inward
         write(*,*) "derivative outwards:", derivative_outward 
         write(*,*) "difference fractions:", difference_fractions
         write(*,*) "weighed intgeral inwards:", weighed_integral_in 
         write(*,*) "weighed intgeral outwards:", weighed_integral_out 
         write(*,*) "deltalambda:", deltalambda
         write(*,*) "oldlambda:", oldlambda 
         write(*,*) "newlambda:", newlambda 
        ! write(*,*) "----inwards array:-----"
        ! do i = 1, x_m+10
        !    write(*,*) i, ":", inwards_array(i) 
        ! end do 
        ! write(*,*)"---- outwards array----"
        ! do i = X_m + 10 , 1 , -1 
        !    write(*,*) i, ":",outwards_array(i) 
        ! end do 

         write(*,*) "-------------------------------------------"

         
          
         cyclecount = cyclecount + 1
         if (newlambda - oldlambda <= convergence) then 
            write(*,*) "shooting algorithm has converged in:", cyclecount, "cycles"
            exit
         end if  
         if (cyclecount >= maxcyclecount) then 
            write(*,*) "shooting algorithm does not converge within", cyclecount,"cycles" 
         end if 
         
      end do  
      eigenvalues(j) = newlambda
      eigenvectors(j, 1:x_m) = inwards_array ! eigenvalue array is constructed from inwards solutions as the first slice
      eigenvectors(j, x_m: gridsize) = outwards_array ! and the outwards solution as the second slice

      write(*,*)"Eigenvalue:", j, "=", eigenvalues(j) 
      cyclecount = 0 
   end do 
   end subroutine 
end module  
