module ShootingTestModule 
   use shootingModule 
   implicit none 
   save 
   public TestNewCalc
   private 

contains 
   subroutine TestNewCalc(self)
      type (calculation), intent(inout) :: self
      call NewCalc(self) 
      if (self%ThreePointGrid%N == 10 ) then 
          write(*,*) "you set three point grid size to:", self%ThreePointGrid%N
          write(*,*) "correct data type"
          ! I want to write a kind test, how do I do that? 
      end if 
      write(*,*) "you set three point mesh size to:", self%ThreePointGrid%h
      write(*,*) "you set Shooting grid size to:", self%ShootingGrid%N
      write(*,*) "you set Shooting grid size to:", self%ShootingGrid%h
      write(*,*) "you set convergence accuracy to:", self%convergence

    end subroutine 
end module 
