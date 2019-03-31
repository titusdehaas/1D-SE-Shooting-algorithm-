module PotentialModule

   implicit none 
   save 
   public NewGrid, NewPotential, Grid, ParticleInABox, MorsePotential
   private  
   
   ! a type that stores the data needed to generate a grid
   type Grid
      real(8) :: mesh
      integer :: gridsize
   end type
  
contains 
   subroutine NewGrid(gridsize, mesh, grid)
      ! subroutine that generates a grid for a given gridsize and mesh 
      real(8), intent(in) :: mesh
      real(8), allocatable, intent(out):: grid(:)  
      integer :: i, gridsize 

      allocate(grid(gridsize))
      grid = 0 
      do i = 1, gridsize 
         grid(i) = grid(i) + i*mesh 
      end do 
   end subroutine

   subroutine NewPotential(gridsize, mesh, potential_type, potential_array, D, alpha, x_e) 
      !subroutine that genrates an array containing the values of a specified potential for a given gridsize and meshsize 
      real(8), intent(in):: mesh ! meshsize
      integer, intent(in):: gridsize ! Gridsize 
      character(32), intent(in) :: potential_type !pass type of potential as a string
      real(8), intent(in), optional :: D, alpha, x_e !parameters for a Morse potential
      real(8), allocatable, intent(out) :: potential_array(:) !output array 
      if (potential_type == "Particle in a box") then 
         call ParticleInAbox(gridsize, mesh, potential_array) 
      elseif (potential_type == "Morse potential" .and. present(D) .and. present(alpha) .and. present(x_e)) then 
         call MorsePotential(gridsize, mesh, D, alpha, x_e, potential_array) 
      else 
              write(*,*)"No valid potential or not enought parameters where added"
      endif 
   end subroutine
   
   subroutine  ParticleInABox(gridsize, mesh, box_array)
      ! subroutine that generates potential values for a particle in a box 
      real(8), intent(in) :: mesh
      integer, intent(in) :: gridsize
      real(8), allocatable, intent(out) :: box_array(:) 
      call NewGrid(gridsize, mesh, box_array) 
      box_array = 0.0d0       
   end subroutine   
   
   subroutine  MorsePotential(gridsize, mesh, D, alpha, x_e, morse_array) 
      ! subroutine that calculates potential values for a Morse potential 
      real(8), intent(in) :: mesh, D, alpha, x_e 
      integer, intent(in) :: gridsize 
      real(8), allocatable, intent(out) :: morse_array(:)
      integer :: i 
      call NewGrid(gridsize, mesh, morse_array)
      !calculate the potential array values 
      do i = 1, gridsize 
         morse_array(i) = D*(1-exp(-alpha*((morse_array(i)-x_e)**2)))
      end do  

   end subroutine            

end module 
