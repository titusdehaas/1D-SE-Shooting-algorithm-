module WriterModule
   use PotentialModule 
   implicit none 
   save 
   public Write, NewOutputFile, OutputFile   
   private

   type OutputFile 
      character(32) :: file_name 
      integer :: file_number
   end type 

   interface Write 
      module procedure WriteEigenvectors, WriteArray
   end interface 

contains 
   subroutine WriteEigenvectors(self, mesh,  number_solutions, eigenvectors)
      type(OutputFile), intent(in) :: self
      integer, intent (in) :: number_solutions  
      real(8), intent(in) :: eigenvectors(:,:), mesh 
      integer:: i 
      real(8), allocatable :: x_axis(:)      
      call NewGrid(size(eigenvectors(1,:)), mesh, x_axis) 
      open(self%file_number, file=self%file_name)
      do i = 1, size(eigenvectors(1,:)) 
         write(self%file_number, '(f4.2, 3f10.6)') x_axis(i), eigenvectors(i,1:number_solutions)
      end do

      close (self%file_number) 
   end subroutine 

   subroutine WriteArray(self, values) 
      type(OutputFile), intent(in) :: self 
      real(8), intent(in) :: values(:) 
      integer :: i 
      
      open(self%file_number, file=self%file_name) 
      write(self%file_number,*) "eigenvalues obtained with three point scheme"
      do i = 1, size(values) 
         write(self%file_number,*) i, values(i)
      end do
      close(self%file_number)       

   end subroutine 

   subroutine NewOutputFile(file_number, file_name, self) 
      integer, intent(in) :: file_number 
      character(len=*), intent(in) :: file_name
      type(OutputFile), intent(out) :: self 

      self%file_number = file_number 
      self%file_name = file_name 
   end subroutine 
end module 
