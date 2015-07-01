!--------------------------------------------------------
! Functions and subroutines for delivering error messages
!--------------------------------------------------------
module error_messages

contains

subroutine error_allocate(istat)
!---------------------------
! ERROR in allocating memory
!---------------------------
  implicit none
  integer :: istat

  if (istat /= 0) then
    write(*,'(a)') 'ERROR: in allocating memory'
    stop
  end if

end subroutine error_allocate

subroutine error_open_file(filename)
!------------------------
! ERROR in opening a file
!------------------------
  use global_constants
  use string_handling
  implicit none
  character(len=chr80) :: filename

  write(*,'(2(a))') 'ERROR: in opening file ', trim(quote(filename))
  stop

end subroutine error_open_file

subroutine error_open_scratch()
!------------------------
! ERROR in opening a file
!------------------------
  use global_constants
  use string_handling
  implicit none

  write(*,'(a)') 'ERROR: in opening a temporary file'
  stop

end subroutine error_open_scratch

subroutine error_read_file(filename,line_number)
!-----------------------------
! ERROR in reading from a file
!-----------------------------
  use global_constants
  use string_handling
  implicit none
  integer, optional :: line_number
  character(len=chr80) :: filename

  if (.not. present(line_number)) then
    write(*,'(2(a))') 'ERROR: in reading file ', trim(quote(filename))
    stop
  else
    write(*,'(4(a))') 'ERROR: in reading file ', trim(quote(filename)) , &
      ', line ', trim(int2str(line_number))
    stop
  end if

end subroutine error_read_file

subroutine error_write_file(filename)
!---------------------------
! ERROR in writing to a file
!---------------------------
  use global_constants
  use string_handling
  implicit none
  character(len=chr80) :: filename

  if (len_trim(filename) /= 0) then
    write(*,'(2(a))') 'ERROR: in writing file ', trim(quote(filename))
    stop
  else
    write(*,'(a)') 'ERROR: in writing to a temporary file'
    stop
  end if

end subroutine error_write_file

end module error_messages
