module PotentialModule

   implicit none 
   save 
   public NewGrid, NewPotential, Grid, Calculation, ParticleInABox, MorsePotential
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
      real(8) :: alpha 
      real(8) :: D 
      real(8) :: x_e
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

   subroutine NewPotential(N,h, potential_type, potential_array, D, alpha, x_e) 
      !subroutine that genrates an array containing the values of a speciefied potential for a given gridsize and meshsize 
      real(8), intent(in):: h ! meshsize
      integer, intent(in):: N ! Gridsize 
      character(32), intent(in) :: potential_type !pass type of potential as a string
      real(8), intent(in), optional :: D, alpha, x_e
      real(8), allocatable, intent(out) :: potential_array(:)

      print*, alpha     
      if (potential_type == "Particle in a box") then 
         call ParticleInAbox(N, h, potential_array) 
      elseif (potential_type == "Morse potential" .and. present(D) .and. present(alpha) .and. present(x_e)) then 
         call MorsePotential(N, h, D, alpha, x_e, potential_array) 
      else 
              write(*,*)"No valid potential or not enought parameters where added"
      endif 

   end subroutine
   
   subroutine  ParticleInABox(N, h, box_array)
      real(8), intent(in) :: h
      integer, intent(in) :: N
      real(8), allocatable, intent(out) :: box_array(:) 
      call NewGrid(N, h, box_array) 
      box_array = 0.0d0       
   end subroutine   
   
   subroutine  MorsePotential(N, h, D, alpha, x_e, morse_array) 
      real(8), intent(in) :: h, D, alpha, x_e 
      integer, intent(in) :: N 
      real(8), allocatable, intent(out) :: morse_array(:)
      integer :: i 
      call NewGrid(N, h, morse_array)
      
      do i = 1, N 
         morse_array(i) = D*(1-exp(-alpha*((morse_array(i)-x_e)**2)))
      end do  

   end subroutine  

   ! instead of using interface with functions, use if (string == ...) then 
   ! call ...() 
           

end module 
