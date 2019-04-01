module ShootingAlgorithmModule
   use integration_module
   use ThreePointSolverModule 
   use WriterModule
   implicit none  
   save
   public ShootingAlgorithm 
   private 
   
contains 
   subroutine ShootingAlgorithm(gridsize, mesh, init_gridsize, init_mesh, convergence, potential, &
         &number_solutions, eigenvalues, eigenvectors) 

      ! Subroutine that calculates eigenvalues of a passed potential, using the shooting algorithm.
      ! User must suply gridsize and mesh for both a course extimation calculation and for shooting method. 
      ! also a convergence criterion, a potential array, boundary conditions and the number of solutions must be passed. 

      real(8), intent(in) :: mesh, init_mesh, convergence, potential(:) 
      integer, intent(in) :: init_gridsize, number_solutions   
      integer, intent(inout) :: gridsize 
      real(8), allocatable, intent(out):: eigenvalues(:), eigenvectors(:,:)
      real(8), allocatable :: lambda(:), inwards_array(:), outwards_array(:), boundaryconditions(:,:)  
      integer :: i, j, x_m, cyclecount, maxcyclecount   
      real(8) :: oldlambda, newlambda, deltalambda, normfactor, phisqrt,  phisqrt1, phisqrt2 !lambda values and normalisation 
      real(8) :: derivative_inward, derivative_outward, integral_inward, integral_outward  ! derivatives and integral values
      real(8) :: difference_fractions, weighed_integral_in, weighed_integral_out, step_inward, step_outward !groupin together terms  
      type(OutputFile) :: log_type


      ! setting up the logfile
      Call NewOutputFile(10, "logfile.txt", log_type)
      call Logger(log_type, "------------------------------------------------")
      Call Logger(log_type, "-----------Shooting Algorithm Logfile-----------")
      call Logger(log_type, "------------------------------------------------")
      

      ! setting the gridsize to an odd number, to make sure there is a median(x_m) 
      if (mod(gridsize, 2) == 1) then
         print*, "gridsize has the right dimension (odd) "
         gridsize = gridsize 
      else 
         gridsize = gridsize + 1
         print*, "Gridsize was increased to an odd number to minimise integration error"
      end if 
      
      !initialising calculation 
      maxcyclecount = 40 ! sets the max cycles before algoritm is exited. 
      cyclecount = 0 
      j = 0
      x_m = (gridsize+1)/2 !setting median 
     
      ! using the threepoint scheme to obtain an initial estimation for the eigenvalues 
      call ThreePointSolver(init_gridsize, init_mesh, potential, lambda, boundaryconditions)
      allocate(inwards_array(x_m+1), outwards_array(x_m+1)) ! allocation of the inwards and outwards solutions 
      allocate(eigenvalues(number_solutions), eigenvectors(number_solutions, gridsize)) ! alocation of the final output
      inwards_array = 0 
      outwards_array = 0 
      eigenvalues = 0 
      eigenvectors = 0 
 
      ! find the jth eigenvalue, for 1 to the number of solutions requested
      do j = 1, number_solutions

      call logger(log_type, "calculation of the:", j, "solution") 
      newlambda = lambda(j) 

      do  ! set up convergence loop 
         
         oldlambda = newlambda 
         inwards_array = 0
         outwards_array = 0

         ! initialising boundary conditions
         outwards_array(1) = boundaryconditions(1,j)
         outwards_array(2) = boundaryconditions(2,j)
         inwards_array(x_m) = boundaryconditions(init_gridsize-1,j)
         inwards_array(x_m+1) = boundaryconditions(init_gridsize,j)

          
         ! calculation of the inwards and outwards solution at a common eigenvalue j. 
         do i = 2, x_m 
            step_outward = (mesh**2.0d0)*(-2.0d0*oldlambda + 2.0d0*potential(i))*outwards_array(i) 
            outwards_array(i+1)=-outwards_array(i-1) + 2*outwards_array(i) + step_outward
         end do 
         do i = x_m, 2, -1 
            step_inward = (mesh**2.0d0)*(-2.0d0*oldlambda + 2.0d0*potential(i))*inwards_array(i)
            inwards_array(i-1)=-inwards_array(i+1) + 2*inwards_array(i) + step_inward 
         end do 
          
         ! normalization of both solutions     
         call newton_cotes(inwards_array**2, mesh, 1, x_m+1, phisqrt1) 
         call newton_cotes(outwards_array**2, mesh, 1, x_m+1, phisqrt2)  
         normfactor = 1/sqrt(phisqrt1+phisqrt2) 
         inwards_array = inwards_array*normfactor
         outwards_array = outwards_array*normfactor  

         ! determine integralvalue over  inwards and outwards functions
         call newton_cotes((inwards_array**2), mesh, 1, x_m+1, integral_inward) 
         call newton_cotes((outwards_array**2), mesh, 1, x_m+1, integral_outward)        

         ! determine derivatives of inwards and outwards solution 
         derivative_inward  = (inwards_array(3) - inwards_array(1))/(2.0d0*mesh)
         derivative_outward  = (outwards_array(x_m+1) - outwards_array(x_m-1))/(2.0d0*mesh) 
         
         ! grouping terms together
         difference_fractions = (derivative_inward/inwards_array(2))-(derivative_outward/outwards_array(x_m))
         weighed_integral_in = integral_inward/(inwards_array(2)**2.0d0)
         weighed_integral_out = integral_outward/(outwards_array(x_m)**2.0d0)

         ! calculation of the correction to the eigenvalue
         deltalambda = 0.5d0*(difference_fractions)/(weighed_integral_out + weighed_integral_in)
         
         ! calculating new eigenvalue and checking convergence
         newlambda  = oldlambda - deltalambda

         !logging/ debugging information 
         call Logger(log_type, "------------------------------------------------") 
         call Logger(log_type, "---cylce:", cyclecount) 
         call Logger(log_type,"normfactor", normfactor) 
         call Logger(log_type,"newlambda", newlambda) 
         call Logger(log_type,"oldlambda", oldlambda) 
         call Logger(log_type," deltalambda", deltalambda)  
         

         !============================================================================================
         ! following can be uncommented for dubugging 

         !write(*,*) "-------------- solution:", j, "cycle:",cyclecount,"-----------------------------"
         !write(*,*) "gridsize, x_m:", gridsize, x_m
         !write(*,*) "size(outwards_array):", size(outwards_array(:))
         !write(*,*) "mesh:", mesh
         !write(*,*) "boundary conditions:", boundaryconditions(1,j), boundaryconditions(2,j),&
         !   &boundaryconditions(init_gridsize-1,j), boundaryconditions(init_gridsize,j) 
         !write(*,*) "normfactor:", normfactor
         !write(*,*) "phisqrt 1 & 2", phisqrt1, phisqrt2 
         !write(*,*) "integeral_inward:", integral_inward
         !write(*,*) "integeral_outward:", integral_outward
         !write(*,*) "derivative inwards:", derivative_inward
         !write(*,*) "derivative outwards:", derivative_outward 
         !write(*,*) "inwards_array(2), outwards_array(x_m):", inwards_array(2), outwards_array(x_m) 
         !write(*,*) "derivative_inward/inwards_array(2):", derivative_inward/inwards_array(2)
         !write(*,*) "derivative_outward/outwards_array(x_m):", derivative_outward/outwards_array(x_m)
         !write(*,*) "difference fractions:", difference_fractions
         !write(*,*) "weighed intgeral inwards:", weighed_integral_in 
         !write(*,*) "weighed intgeral outwards:", weighed_integral_out 
         !write(*,*) "deltalambda:", deltalambda
         !write(*,*) "oldlambda:", oldlambda 
         !write(*,*) "newlambda:", newlambda 
         !write(*,*) "----inwards array:-----"
              
         !do i = 1, x_m+1
         !   write(*,*) i, ":", inwards_array(i) 
         !end do 
         !write(*,*)"---- outwards array----"
         !do i = x_m+1 , 1 , -1 
         !   write(*,*) i, ":",outwards_array(i) 
         !end do 

         !write(*,*) "-------------------------------------------"
         !========================================================================================== 
         
         ! checking the convergence criteria 
         cyclecount = cyclecount + 1
         if (abs(deltalambda) <= convergence) then
            write(*,*) "--------------------------------------------------------" 
            write(*,*) "shooting algorithm has converged in:", cyclecount, "cycles"
            call Logger(log_type, "------------------------------------------------")
            call Logger(log_type, "shooting algorithm has converged in:", cyclecount, "cycles") !writing number of cycles to logfile 
            exit
         end if  

         ! exiting the program if convergence is not reached within a given amount of cycles
         if (cyclecount >= maxcyclecount) then 
            write(*,*) "--------------------------------------------------------"
            write(*,*) "shooting algorithm does not converge within", cyclecount,"cycles"
            call Logger(log_type, "shooting algorithm does not converge within:", cyclecount, "cycles")
            exit  
         end if 
      end do  

      !pasting the obtained eigenvalues and eigenfunctions in the therefore specified arrays 
      eigenvalues(j) = newlambda
      eigenvectors(j, 1:x_m) = outwards_array(1:x_m) ! eigenvalue array is constructed from inwards solutions as the first slice
      eigenvectors(j, x_m+1: gridsize) = inwards_array(2:x_m) ! and the outwards solution as the second slice
      
      !normalisation
      print*, size(eigenvectors(j,:)) 
      call newton_cotes((eigenvectors(j,:)**2), mesh, 1, gridsize, phisqrt)
      normfactor = 1/phisqrt 
      eigenvectors(j,:)  = eigenvectors(j,:)*normfactor 



      !logging the results of the shooting algorithm calculation
      write(*,*)"3 point eigenvalue:", j, "=", lambda(j) 
      write(*,*)"Shooting algorithm eigenvalue:", j, "=", eigenvalues(j) 
      
      call Logger(log_type, "3 point eigenvalue:", lambda(j)) 
      call Logger(log_type, "Shooting algorithm eigenvalue:", eigenvalues(j)) 
      call Logger(log_type, "------------------------------------------------")
      call Logger(log_type, "------------------------------------------------")
      
      
      cyclecount = 0 !resetting cycle count
   end do
   deallocate(inwards_array, outwards_array) ! deallocation of the inwards and outwards solutions 
   end subroutine 
end module  
