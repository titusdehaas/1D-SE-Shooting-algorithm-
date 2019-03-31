module WriterModule
   use PotentialModule 
   implicit none 
   save 
   public Write, Logger, NewOutputFile, EmptyOutputFile, OutputFile  
   private
   
   type OutputFile !data type containing information for writing of outputfile 
      character(32) :: file_name 
      integer :: file_number
   end type 


   ! overloading the eigenvectors/eigenvalue and logging procedures 
   interface Write 
      module procedure WriteEigenvectors, WriteArray
   end interface 
 
   interface Logger
      module procedure LogReal, LogInt, LogSentence 
   end interface 

contains 
   subroutine WriteEigenvectors(self, mesh,  number_solutions, eigenvectors)
      !subroutine for writing the obtained eigenvectors and a x-axis in columns on a specified file 
      type(OutputFile), intent(in) :: self !data to specify output file 
      integer, intent (in) :: number_solutions  
      real(8), intent(in) :: eigenvectors(:,:), mesh !mesh is needed to print out a column of x-values 
      integer:: i 
      real(8), allocatable :: x_axis(:)      
      call NewGrid(size(eigenvectors(1,:)), mesh, x_axis) !to generate a x-axis 
      open(self%file_number, file=self%file_name)

      !writing the columns
      do i = 1, size(eigenvectors(1,:)) 
         write(self%file_number, '(f6.3, 100f14.10)') x_axis(i), eigenvectors(i,1:number_solutions)
      end do                 
      !note: to write more then 100 eigenvector columns, update data specification
      close (self%file_number) 
   end subroutine 

   subroutine WriteArray(self, values) 
      !subroutine for writing an array of values in a specified file
      type(OutputFile), intent(in) :: self 
      real(8), intent(in) :: values(:) 
      integer :: i 
      
      open(self%file_number, file=self%file_name) 
      write(self%file_number,*) "# the obtained eigenvalues"
      do i = 1, size(values) 
         write(self%file_number,'(i4, f18.10)') i, values(i)
      end do
      close(self%file_number)       

   end subroutine 

   subroutine LogReal(self, statement, logvalue, unit) 
      !subroutine to add a real variable name + value to a (log) file 
      type(Outputfile) :: self
      character(*), intent(in):: statement
      character(*), intent(in), optional:: unit 
      real(8), intent(in) :: logvalue
      open(self%file_number, file= self%file_name, position = 'append')   
      write(self%file_number,*) statement, logvalue 
      close(self%file_number)
   end subroutine 

   subroutine LogInt(self, statement, logvalue, unit)
      !subroutine to add an integer variable name + value to a (log) file 
      type(Outputfile) :: self
      character(*), intent(in):: statement
      character(*), intent(in), optional:: unit
      integer, intent(in) :: logvalue
      open(self%file_number, file= self%file_name, position = 'append')
      if (present(unit)) then 
         write(self%file_number,*) statement, logvalue, unit
      else 
         write(self%file_number,*) statement, logvalue
      end if 
      close(self%file_number)
   end subroutine
   
   subroutine LogSentence(self, statement, statement2, statement3)
      !subroutine to add a sentence, containing max 3 strings, to a (log) file 
      type(Outputfile) :: self
      character(*), intent(in):: statement
      character(*), intent(in), optional:: statement2, statement3 
      open(self%file_number, file= self%file_name, position = 'append')
      if (present(statement2) .and. present(statement3)) then 
         write(self%file_number,*) statement, statement2, statement3 
       else if (present(statement2)) then 
         write(self%file_number,*) statement, statement2
      else  
         write(self%file_number,*) statement
      end if 
      close(self%file_number)
   end subroutine


   subroutine NewOutputFile(file_number, file_name, self)
      !subroutine to store output specifcations(name and number) into a data type: OutputFile  
      integer, intent(in) :: file_number 
      character(len=*), intent(in) :: file_name
      type(OutputFile), intent(out) :: self 
      self%file_number = file_number 
      self%file_name = file_name 
   end subroutine 
   
   subroutine EmptyOutputFile(self)
      !subroutine to empty an output file 
      type(outputFile) :: self
      open(self%file_number, file= self%file_name)
      write(self%file_number,*) 
      close(self%file_number)
   end subroutine
 





end module 
