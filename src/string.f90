!--------------------------------------------------------
! Functions and subroutines for character string handling
!--------------------------------------------------------
module string_handling

contains

function int2str(i) result(s) 
!-------------------------------------------
! Convert an integer into a character string
! - Add exception handling
!-------------------------------------------
  use global_constants
  implicit none

  integer :: i
  character(len=chr80) :: s

  write(s,'(i80)') i
  s=adjustl(trim(s))

end function int2str

function str2int(s) result(i) 
!-------------------------------------------
! Convert a character string into an integer
! - Add exception handling
!-------------------------------------------
  use global_constants
  implicit none

  integer :: i,ios
  character(len=chr80) :: s

  read(s,*,iostat=ios) i 
  if (ios /= 0) then
    write(*,'(3(a))') 'ERROR: in converting character string ', &
      trim(quote(s)), ' to integer'
    stop
  end if

end function str2int

function str2float(s) result(x) 
!-------------------------------------------
! Convert a character string into an integer
! - Add exception handling
!-------------------------------------------
  use global_constants
  implicit none

  integer :: ios
  real(kind=rk) :: x
  character(len=chr80) :: s

  read(s,*,iostat=ios) x 
  if (ios /= 0) then
    write(*,'(3(a))') 'ERROR: in converting character string ', &
      trim(quote(s)), ' to a floating point number'
    stop
  end if

end function str2float

function quote(strin) result(strout)
!----------------------------------------
! Add single quotes to a character string
! - Add exception handling
!----------------------------------------
  use global_constants
  implicit none

  character(len=chr80) :: strin
  character(len=chr82) :: strout

  strout=char(39) // trim(adjustl(strin)) // char(39)

end function quote

function double_quote(strin) result(strout)
!----------------------------------------
! Add double quotes to a character string
! - Add exception handling
!----------------------------------------
  use global_constants
  implicit none

  character(len=chr80) :: strin
  character(len=chr82) :: strout

  strout=char(34) // trim(adjustl(strin)) // char(34)

end function double_quote

function uppercase(str1) result(str2)
!----------------------------------------------
! Change lowercase letters to uppercase letters
!----------------------------------------------
  use global_constants

  integer :: i
  character(len=chr80) :: str1,str2

  str2=''
  do i=1,len_trim(str1)
    if (iachar(str1(i:i)) >= 97 .and. iachar(str1(i:i)) <= 122) then
      str2(i:i)=char(iachar(str1(i:i))-32)
    else
      str2(i:i)=str1(i:i)
    end if
  end do

end function uppercase

function lowercase(str1) result(str2)
!----------------------------------------------
! Change uppercase letters to lowercase letters
!----------------------------------------------
  use global_constants

  integer :: i
  character(len=chr80) :: str1,str2

  str2=''
  do i=1,len_trim(str1)
    if (iachar(str1(i:i)) >= 65 .and. iachar(str1(i:i)) <= 90) then
      str2(i:i)=char(iachar(str1(i:i))+32)
    else
      str2(i:i)=str1(i:i)
    end if
  end do

end function lowercase

function lowercase_30(str1) result(str2)
!----------------------------------------------
! Change uppercase letters to lowercase letters
!----------------------------------------------
  use global_constants

  integer :: i
  character(len=chr30) :: str1,str2

  str2=''
  do i=1,len_trim(str1)
    if (iachar(str1(i:i)) >= 65 .and. iachar(str1(i:i)) <= 90) then
      str2(i:i)=char(iachar(str1(i:i))+32)
    else
      str2(i:i)=str1(i:i)
    end if
  end do

end function lowercase_30

function format_real_sdf(ncols) result(fstring)
!---------------------------------------------------
! Generate a format string for space delimited reals
!---------------------------------------------------
  use global_constants
  implicit none
  integer :: ncols
  character(len=30) :: ctmp
  character(len=80) :: fstring

  select case(ncols)
  case (1) 
    fstring='(es15.7e3)'
  case (2)
    fstring='(es15.7e3,1x,es15.7e3)'
  case (3:)
    write(ctmp,'(i30)') ncols-1
    fstring='('  // trim(adjustl(ctmp)) // '(es15.7e3,1x),es15.7e3)'
  case default
    write(*,'(a)') 'ERROR: in generating format string'
    stop
  end select

end function format_real_sdf

function format_char_sdf(ncols) result(fstring)
!---------------------------------------------------------------
! Generate a format string for comma separated character strings
!---------------------------------------------------------------
  use global_constants
  implicit none
  integer :: ncols
  character(len=30) :: ctmp
  character(len=80) :: fstring

  select case(ncols)
  case (1) 
    fstring='(a)'
  case (2)
    fstring='(a,1x,a)'
  case (3:)
    write(ctmp,'(i30)') ncols-1
    fstring='('  // trim(adjustl(ctmp)) // '(a,1x),a)'
  case default
    write(*,'(a)') 'ERROR: in generating format string'
    stop
  end select

end function format_char_sdf

function format_real_csv(ncols) result(fstring)
!---------------------------------------------------
! Generate a format string for comma separated reals
!---------------------------------------------------
  use global_constants
  implicit none
  integer :: ncols
  character(len=30) :: ctmp
  character(len=80) :: fstring

  select case(ncols)
  case (1) 
    fstring='(es15.7e3)'
  case (2)
    fstring='(es15.7e3,",",es15.7e3)'
  case (3:)
    write(ctmp,'(i30)') ncols-1
    fstring='('  // trim(adjustl(ctmp)) // '(es15.7e3,","),es15.7e3)'
  case default
    write(*,'(a)') 'ERROR: in generating format string'
    stop
  end select

end function format_real_csv

function format_char_csv(ncols) result(fstring)
!---------------------------------------------------------------
! Generate a format string for space delimited character strings
!---------------------------------------------------------------
  use global_constants
  implicit none
  integer :: ncols
  character(len=30) :: ctmp
  character(len=80) :: fstring

  select case(ncols)
  case (1) 
    fstring='(a)'
  case (2)
    fstring='(a,",",a)'
  case (3:)
    write(ctmp,'(i30)') ncols-1
    fstring='('  // trim(adjustl(ctmp)) // '(a,","),a)'
  case default
    write(*,'(a)') 'ERROR: in generating format string'
    stop
  end select

end function format_char_csv

function format_char_cssv(ncols) result(fstring)
!---------------------------------------------------------------
! Generate a format string for space delimited character strings
! Revision 2. Insert a space after each comma
!---------------------------------------------------------------
  use global_constants
  implicit none
  integer :: ncols
  character(len=30) :: ctmp
  character(len=80) :: fstring

  select case(ncols)
  case (1) 
    fstring='(a)'
  case (2)
    fstring='(a,", ",a)'
  case (3:)
    write(ctmp,'(i30)') ncols-1
    fstring='('  // trim(adjustl(ctmp)) // '(a,", "),a)'
  case default
    write(*,'(a)') 'ERROR: in generating format string'
    stop
  end select

end function format_char_cssv

function count_sdf_columns(input_line) result(ncols)
!--------------------------------------------------------
! Count the number of columns in an SDF-format input line
!--------------------------------------------------------
  use global_constants
  implicit none

  integer :: i,ncols
  logical :: preceding_space
  character(len=input_line_length) :: input_line,input_line_tmp

  ! Remove space characters from the beginning of the input line
  input_line_tmp=adjustl(input_line)

  if (input_line_tmp(1:1) == char(32)) then
    ! No columns
    ncols=0
  else
    ! One or more columns
    ncols=1; preceding_space=.false.
    do i=1,len_trim(input_line_tmp)
      if (input_line_tmp(i:i) == char(32)) then
        preceding_space=.true.
      else
        if (preceding_space) then
          ncols=ncols+1
        end if
        preceding_space=.false.
      end if
    end do
  end if

end function count_sdf_columns

function numerical(stmp_in) result(answer)
!-------------------------------------------------------------------
! Determine whether a character string contains only integers or not
!-------------------------------------------------------------------
  use global_constants
  implicit none

  integer :: i
  logical :: answer
  character(len=chr80) :: stmp
  character(len=chr80) :: stmp_in

  stmp=adjustl(stmp_in)

  if (len_trim(stmp) == 0) then
    answer=.false.
  else
    answer=.true.
    loop: do i=1,len_trim(stmp)
      if ((iachar(stmp(i:i)) < 48) .or. (iachar(stmp(i:i)) > 57)) then
        if (i == 1) then
          answer=.false.
        else
          if (iachar(stmp(i:i)) /= 32 .and. &
              iachar(stmp(i:i)) /= 0) then
            answer=.false.
          end if
        end if
      end if
    end do loop
  end if

end function numerical

function numerical_30(stmp_in) result(answer)
!-------------------------------------------------------------------
! Determine whether a character string contains only integers or not
!-------------------------------------------------------------------
  use global_constants
  implicit none

  integer :: i
  logical ::answer
  character(len=chr30) :: stmp
  character(len=chr30) :: stmp_in

  stmp=adjustl(stmp_in)

  answer=.true.
  loop: do i=1,len_trim(stmp)
    if ((iachar(stmp(i:i)) < 48) .or. (iachar(stmp(i:i)) > 57)) then
      if (i == 1) then
        answer=.false.
      else
        if (iachar(stmp(i:i)) /= 32 .and. &
            iachar(stmp(i:i)) /= 0) then
          answer=.false.
        end if
      end if
    end if
  end do loop

end function numerical_30

function numerical_20(stmp_in) result(answer)
!-------------------------------------------------------------------
! Determine whether a character string contains only integers or not
!-------------------------------------------------------------------
  use global_constants
  implicit none

  integer :: i
  logical ::answer
  character(len=chr20) :: stmp
  character(len=chr20) :: stmp_in

  stmp=adjustl(stmp_in)

  answer=.true.
  loop: do i=1,len_trim(stmp)
    if ((iachar(stmp(i:i)) < 48) .or. (iachar(stmp(i:i)) > 57)) then
      if (i == 1) then
        answer=.false.
      else
        if (iachar(stmp(i:i)) /= 32 .and. &
            iachar(stmp(i:i)) /= 0) then
          answer=.false.
        end if
      end if
    end if
  end do loop

end function numerical_20

function basename(filename) result(base)
!----------------------
! Extract filename base
!----------------------
  use global_constants
  implicit none

  integer :: i,j
  character(len=chr80) :: filename,base

  base=''; j=1
  loop: do i=1,len_trim(filename)
    if (filename(i:i) == '.') then
      exit loop
    else
      base(j:j)=filename(i:i)
      j=j+1
    end if
  end do loop

end function basename

end module string_handling
