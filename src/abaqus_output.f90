!--------------------------------------------
! Functions and subroutines for ABAQUS output
!--------------------------------------------
module abaqus_output 

contains

subroutine abaqus_output_module()
!---------
! Main sub
!---------
  use global_constants
  use global_variables
  implicit none

  write(*,'(a)') ''
  write(*,'(a)') 'ABAQUS output module'

  !----------------------------------------------
  ! Create external files for boundary conditions
  ! and add include-keywords to ABAQUS input file
  !----------------------------------------------
  if (fem_data_available) then

    select case(trim(transfer_quantity))
    case ('wall_temperature')
      
      call dump_abaqus_amplitudes()
      call dump_abaqus_boundaries()

    case ('net_heat_flux')
      
      write(*,'(t3,a)') 'Warning: net heat flux dump not supported yet'

    case ('adiabatic_surface_temperature')
      
      call dump_abaqus_physical()
      call dump_abaqus_amplitudes()
      call dump_abaqus_cradiate()
      call dump_abaqus_cfilm()

    end select

  else
    write(*,'(t3,a)') 'Nothing to be done'

  end if

end subroutine abaqus_output_module

!***************
! Auxiliary subs
!***************

subroutine dump_abaqus_physical()
!-------------------------------------
! Add *Amplitude lines to ABAQUS input
!-------------------------------------
  use abaqus_reader
  use error_messages
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,ios,nlines
  logical :: omit_lines,replace
  character(len=input_line_length) :: keyword,input_line
  character(len=input_line_length), dimension(:), allocatable :: abaqus_input

  open(unit=iochannel(1),file=trim(fem_input_file),status='old',iostat=ios)
  if (ios /= 0) call error_open_file(fem_input_file)

  !------------------------------------------
  ! Count the number of lines in ABAQUS input
  !------------------------------------------
  nlines=0; replace=.false.
  count_lines: do
    read(iochannel(1),'(a)',iostat=ios) input_line
    select case(ios)
    case (:-1)
      exit count_lines
    case(1:)
      call error_read_file(fem_input_file)
    case default
      if (input_line(1:26) == '** fds2fem-physical-marker') then
        replace=.true.
      end if
      nlines=nlines+1
    end select
  end do count_lines

  rewind(iochannel(1))
  allocate(abaqus_input(nlines),stat=ios); call error_allocate(ios)

  !---------------------
  ! Read in ABAQUS input
  !---------------------
  read_lines: do i=1,nlines
    read(iochannel(1),'(a)',iostat=ios) abaqus_input(i)
    select case(ios)
    case (:-1)
      exit read_lines
    case(1:)
      call error_read_file(fem_input_file)
    end select
  end do read_lines

  rewind(iochannel(1))

  if (.not. replace) then
    write(*,'(t3,2(a))') 'Adding entry for physical constants in file ', &
      trim(quote(fem_input_file))
  else
    write(*,'(t3,2(a))') 'Replacing old entry for physical constants in file ', &
      trim(quote(fem_input_file))
  end if

  !----------------------------
  ! Add physical constants line
  !----------------------------
  omit_lines=.false.
  add_lines: do i=1,nlines

    input_line=adjustl(abaqus_input(i)) 

    if (input_line(1:26) == '** fds2fem-physical-marker') then
      omit_lines=.true.
      cycle add_lines
    end if

    if (.not. omit_lines) then
      write(iochannel(1),'(a)',iostat=ios) trim(abaqus_input(i))
    end if

    if (input_line(1:1) == '*' .and. &
        input_line(2:2) /= '*') then
        
        keyword=get_keyword(input_line)
        if (lowercase(keyword) == 'end assembly') then
              
          write(iochannel(1),'(a)') '** fds2fem-physical-marker (do not remove)'
          write(iochannel(1),'(2(a))') '*Physical Constants, absolute zero=-273.15, stefan boltzmann=5.669E-8' 
          ! FDS SIGMA=5.670373E-8_EB
          ! https://physics.nist.gov/cuu/Constants/, Source: 2014 CODATA 5.670367 E-8
          write(iochannel(1),'(2(a))') '*Physical Constants, absolute zero=-273.15, stefan boltzmann=5.670E-8' 

        end if

    end if
    
    if (omit_lines) omit_lines=.false.

  end do add_lines

  close(unit=iochannel(1))

  deallocate(abaqus_input,stat=ios); call error_allocate(ios) 
 
end subroutine dump_abaqus_physical

subroutine dump_abaqus_amplitudes()
!-------------------------------------
! Add *Amplitude lines to ABAQUS input
!-------------------------------------
  use abaqus_reader
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,k,ios,nlines
  logical :: omit_lines,replace
  character(len=chr80) :: amplitude_file
  character(len=input_line_length) :: keyword,input_line
  character(len=input_line_length), dimension(:), allocatable :: abaqus_input

  amplitude_file=trim(fem_input_file) // '.amp'
  open(unit=iochannel(1),file=trim(amplitude_file),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(amplitude_file)

  do j=1,nnodes_fem
    write(iochannel(1),'(2(a))') '*Amplitude, name=AMP-', trim(int2str(j)) 
    do k=1,ntimes_fem-1
      write(iochannel(1),'(es15.7e3,a,es15.7e3,a)') fem_time(k), ', ', fem_data(k,j), ','
    end do
    write(iochannel(1),'(es15.7e3,a,es15.7e3)') fem_time(ntimes_fem), ', ', fem_data(ntimes_fem,j)
  end do

  if (read_hcoeff) then
    ! Provide amplitudes for heat transfer coefficients
    do j=1,nnodes_fem
      write(iochannel(1),'(2(a))') '*Amplitude, name=AMP-H-', trim(int2str(j)) 
      do k=1,ntimes_fem-1
        write(iochannel(1),'(es15.7e3,a,es15.7e3,a)') fem_time(k), ', ', fem_hcoeff(k,j), ','
      end do
      write(iochannel(1),'(es15.7e3,a,es15.7e3)') fem_time(ntimes_fem), ', ', fem_hcoeff(ntimes_fem,j)
    end do
  end if

  close(unit=iochannel(1))

  write(*,'(t3,5(a))') trim(int2str(ntimes_fem)), &
    ' time steps of amplitude data written in file ', &
    trim(quote(amplitude_file))

  open(unit=iochannel(1),file=trim(fem_input_file),status='old',iostat=ios)
  if (ios /= 0) call error_open_file(fem_input_file)

  !------------------------------------------
  ! Count the number of lines in ABAQUS input
  !------------------------------------------
  nlines=0; replace=.false.
  count_lines: do
    read(iochannel(1),'(a)',iostat=ios) input_line
    select case(ios)
    case (:-1)
      exit count_lines
    case(1:)
      call error_read_file(fem_input_file)
    case default
      if (input_line(1:27) == '** fds2fem-amplitude-marker') then
        replace=.true.
      end if
      nlines=nlines+1
    end select
  end do count_lines

  rewind(iochannel(1))
  allocate(abaqus_input(nlines),stat=ios); call error_allocate(ios)

  !---------------------
  ! Read in ABAQUS input
  !---------------------
  read_lines: do i=1,nlines
    read(iochannel(1),'(a)',iostat=ios) abaqus_input(i)
    select case(ios)
    case (:-1)
      exit read_lines
    case(1:)
      call error_read_file(fem_input_file)
    end select
  end do read_lines

  rewind(iochannel(1))

  if (.not. replace) then
    write(*,'(t3,2(a))') 'Adding include entry for amplitudes in file ', &
      trim(quote(fem_input_file))
  else
    write(*,'(t3,2(a))') 'Replacing old include entry for amplitudes in file ', &
      trim(quote(fem_input_file))
  end if

  !--------------------
  ! Add amplitude lines
  !--------------------
  omit_lines=.false.
  add_lines: do i=1,nlines

    input_line=adjustl(abaqus_input(i)) 

    if (input_line(1:27) == '** fds2fem-amplitude-marker') then
      omit_lines=.true.
      cycle add_lines
    end if

    if (.not. omit_lines) then
      write(iochannel(1),'(a)',iostat=ios) trim(abaqus_input(i))
    end if

    if (input_line(1:1) == '*' .and. &
        input_line(2:2) /= '*') then
        
        keyword=get_keyword(input_line)
        if (lowercase(keyword) == 'end assembly') then
              
          write(iochannel(1),'(a)') '** fds2fem-amplitude-marker (do not remove)'
          write(iochannel(1),'(2(a))') '*Include, input=', trim(amplitude_file) 

        end if

    end if
    
    if (omit_lines) omit_lines=.false.

  end do add_lines

  close(unit=iochannel(1))

  deallocate(abaqus_input,stat=ios); call error_allocate(ios) 
 
end subroutine dump_abaqus_amplitudes

subroutine dump_abaqus_boundaries()
!------------------------------------
! Add *Boundary lines to ABAQUS input
!------------------------------------
  use abaqus_reader
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,nlines
  logical :: omit_lines,replace
  character(len=chr80) :: bc_file
  character(len=input_line_length) :: keyword,input_line
  character(len=input_line_length), dimension(:), allocatable :: abaqus_input

  bc_file=trim(fem_input_file) // '.bc'
  open(unit=iochannel(1),file=trim(bc_file),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(bc_file)

  do j=1,nnodes_fem
    write(iochannel(1),'(3(a))') '** Name: Temp-BC-', trim(int2str(j)), ' Type: Temperature'
    write(iochannel(1),'(2(a))') '*Boundary, amplitude=AMP-', trim(int2str(j)) 
    write(iochannel(1),'(2(a))') trim(fem_node_name(j)), ', 11, 11, 1.0'
  end do

  close(unit=iochannel(1))

  write(*,'(t3,5(a))') trim(int2str(nnodes_fem)), &
    ' nodal boundary conditions written in file ', &
    trim(quote(fem_input_file))

  open(unit=iochannel(1),file=trim(fem_input_file),status='old',iostat=ios)
  if (ios /= 0) call error_open_file(fem_input_file)

  !------------------------------------------
  ! Count the number of lines in ABAQUS input
  !------------------------------------------
  nlines=0; replace=.false.
  count_lines: do
    read(iochannel(1),'(a)',iostat=ios) input_line
    select case(ios)
    case (:-1)
      exit count_lines
    case(1:)
      call error_read_file(fem_input_file)
    case default
      if (input_line(1:26) == '** fds2fem-boundary-marker') then
        replace=.true.
      end if
      nlines=nlines+1
    end select
  end do count_lines

  rewind(iochannel(1))
  allocate(abaqus_input(nlines),stat=ios); call error_allocate(ios)

  !---------------------
  ! Read in ABAQUS input
  !---------------------
  read_lines: do i=1,nlines
    read(iochannel(1),'(a)',iostat=ios) abaqus_input(i)
    select case(ios)
    case (:-1)
      exit read_lines
    case(1:)
      call error_read_file(fem_input_file)
    end select
  end do read_lines

  rewind(iochannel(1))

  if (.not. replace) then
    write(*,'(t3,2(a))') 'Adding include entry for boundary conditions in file ', &
      trim(quote(fem_input_file))
  else
    write(*,'(t3,2(a))') 'Replacing old include entry for boundary conditions in file ', &
      trim(quote(fem_input_file))
  end if

  !-----------------------------
  ! Add boundary condition lines
  !-----------------------------
  omit_lines=.false.
  add_lines: do i=1,nlines

    input_line=adjustl(abaqus_input(i)) 

    if (input_line(1:26) == '** fds2fem-boundary-marker') then
      omit_lines=.true.
      cycle add_lines
    end if

    if (input_line(1:1) == '*' .and. &
        input_line(2:2) /= '*') then
        
      keyword=get_keyword(input_line)
      if (lowercase(keyword) == 'end step') then

        write(iochannel(1),'(a)') '** fds2fem-boundary-marker (do not remove)'
        write(iochannel(1),'(2(a))') '*Include, input=', trim(bc_file)

      end if
    end if

    if (.not. omit_lines) then
      write(iochannel(1),'(a)',iostat=ios) trim(abaqus_input(i))
    end if
       
    if (omit_lines) omit_lines=.false.

  end do add_lines

  close(unit=iochannel(1))

  deallocate(abaqus_input,stat=ios); call error_allocate(ios)

end subroutine dump_abaqus_boundaries

subroutine dump_abaqus_cradiate()
!------------------------------------
! Add *Boundary lines to ABAQUS input
!------------------------------------
  use abaqus_reader
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,nlines
  logical :: omit_lines,replace
  character(len=chr80) :: cradiate_file
  character(len=input_line_length) :: keyword,input_line
  character(len=input_line_length), dimension(:), allocatable :: abaqus_input

  cradiate_file=trim(fem_input_file) // '.cradiate'
  open(unit=iochannel(1),file=trim(cradiate_file),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(cradiate_file)

  do j=1,nnodes_fem
    write(iochannel(1),'(3(a))') '** Name: Cradiate-BC-', trim(int2str(j)), ' Type: Radiative heat flux'
    write(iochannel(1),'(2(a))') '*Cradiate, amplitude=AMP-', trim(int2str(j)) 
    write(iochannel(1),'(2(a),es15.7e3,a,es15.7e3)') trim(fem_node_name(j)), ', ', abaqus_node_area(j), ', 1.0, ', abaqus_node_emissivity(j)
  end do

  close(unit=iochannel(1))

  write(*,'(t3,5(a))') trim(int2str(nnodes_fem)), &
    ' nodal boundary conditions (Cradiate) written in file ', &
    trim(quote(fem_input_file))

  open(unit=iochannel(1),file=trim(fem_input_file),status='old',iostat=ios)
  if (ios /= 0) call error_open_file(fem_input_file)

  !------------------------------------------
  ! Count the number of lines in ABAQUS input
  !------------------------------------------
  nlines=0; replace=.false.
  count_lines: do
    read(iochannel(1),'(a)',iostat=ios) input_line
    select case(ios)
    case (:-1)
      exit count_lines
    case(1:)
      call error_read_file(fem_input_file)
    case default
      if (input_line(1:26) == '** fds2fem-cradiate-marker') then
        replace=.true.
      end if
      nlines=nlines+1
    end select
  end do count_lines

  rewind(iochannel(1))
  allocate(abaqus_input(nlines),stat=ios); call error_allocate(ios)

  !---------------------
  ! Read in ABAQUS input
  !---------------------
  read_lines: do i=1,nlines
    read(iochannel(1),'(a)',iostat=ios) abaqus_input(i)
    select case(ios)
    case (:-1)
      exit read_lines
    case(1:)
      call error_read_file(fem_input_file)
    end select
  end do read_lines

  rewind(iochannel(1))

  if (.not. replace) then
    write(*,'(t3,2(a))') 'Adding include entry for boundary conditions (Cradiate) in file ', &
      trim(quote(fem_input_file))
  else
    write(*,'(t3,2(a))') 'Replacing old include entry for boundary conditions (Cradiate) in file ', &
      trim(quote(fem_input_file))
  end if

  !-----------------------------
  ! Add boundary condition lines
  !-----------------------------
  omit_lines=.false.
  add_lines: do i=1,nlines

    input_line=adjustl(abaqus_input(i)) 

    if (input_line(1:26) == '** fds2fem-cradiate-marker') then
      omit_lines=.true.
      cycle add_lines
    end if

    if (input_line(1:1) == '*' .and. &
        input_line(2:2) /= '*') then
        
      keyword=get_keyword(input_line)
      if (lowercase(keyword) == 'end step') then

        write(iochannel(1),'(a)') '** fds2fem-cradiate-marker (do not remove)'
        write(iochannel(1),'(2(a))') '*Include, input=', trim(cradiate_file)

      end if
    end if

    if (.not. omit_lines) then
      write(iochannel(1),'(a)',iostat=ios) trim(abaqus_input(i))
    end if
       
    if (omit_lines) omit_lines=.false.

  end do add_lines

  close(unit=iochannel(1))

  deallocate(abaqus_input,stat=ios); call error_allocate(ios)

end subroutine dump_abaqus_cradiate

subroutine dump_abaqus_cfilm()
!------------------------------------
! Add *Boundary lines to ABAQUS input
!------------------------------------
  use abaqus_reader
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,nlines
  logical :: omit_lines,replace
  character(len=chr80) :: cfilm_file
  character(len=input_line_length) :: keyword,input_line
  character(len=input_line_length), dimension(:), allocatable :: abaqus_input

  cfilm_file=trim(fem_input_file) // '.cfilm'
  open(unit=iochannel(1),file=trim(cfilm_file),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(cfilm_file)

  if (read_hcoeff) then
    ! Heat transfer coefficient provided from FDS
    do j=1,nnodes_fem
      write(iochannel(1),'(3(a))') '** Name: Cfilm-BC-', trim(int2str(j)), ' Type: Convective heat flux'
      write(iochannel(1),'(4(a))') '*Cfilm, amplitude=AMP-', trim(int2str(j)), ', film amplitude=AMP-H-', trim(int2str(j))
      write(iochannel(1),'(2(a),es15.7e3,a)') trim(fem_node_name(j)), ', ', abaqus_node_area(j), ', 1.0, 1.0'
    end do

  else
    ! Constant heat transfer coefficient
    do j=1,nnodes_fem
      write(iochannel(1),'(3(a))') '** Name: Cfilm-BC-', trim(int2str(j)), ' Type: Convective heat flux'
      write(iochannel(1),'(2(a))') '*Cfilm, amplitude=AMP-', trim(int2str(j)) 
      write(iochannel(1),'(2(a),es15.7e3,a,es15.7e3)') trim(fem_node_name(j)), ', ', abaqus_node_area(j), ', 1.0, ', abaqus_node_hcoeff(j)
    end do
  end if

  close(unit=iochannel(1))

  write(*,'(t3,5(a))') trim(int2str(nnodes_fem)), &
    ' nodal boundary conditions (Cfilm) written in file ', &
    trim(quote(fem_input_file))

  open(unit=iochannel(1),file=trim(fem_input_file),status='old',iostat=ios)
  if (ios /= 0) call error_open_file(fem_input_file)

  !------------------------------------------
  ! Count the number of lines in ABAQUS input
  !------------------------------------------
  nlines=0; replace=.false.
  count_lines: do
    read(iochannel(1),'(a)',iostat=ios) input_line
    select case(ios)
    case (:-1)
      exit count_lines
    case(1:)
      call error_read_file(fem_input_file)
    case default
      if (input_line(1:23) == '** fds2fem-cfilm-marker') then
        replace=.true.
      end if
      nlines=nlines+1
    end select
  end do count_lines

  rewind(iochannel(1))
  allocate(abaqus_input(nlines),stat=ios); call error_allocate(ios)

  !---------------------
  ! Read in ABAQUS input
  !---------------------
  read_lines: do i=1,nlines
    read(iochannel(1),'(a)',iostat=ios) abaqus_input(i)
    select case(ios)
    case (:-1)
      exit read_lines
    case(1:)
      call error_read_file(fem_input_file)
    end select
  end do read_lines

  rewind(iochannel(1))

  if (.not. replace) then
    write(*,'(t3,2(a))') 'Adding include entry for boundary conditions (Cfilm) in file ', &
      trim(quote(fem_input_file))
  else
    write(*,'(t3,2(a))') 'Replacing old include entry for boundary conditions (Cfilm) in file ', &
      trim(quote(fem_input_file))
  end if

  !-----------------------------
  ! Add boundary condition lines
  !-----------------------------
  omit_lines=.false.
  add_lines: do i=1,nlines

    input_line=adjustl(abaqus_input(i)) 

    if (input_line(1:23) == '** fds2fem-cfilm-marker') then
      omit_lines=.true.
      cycle add_lines
    end if

    if (input_line(1:1) == '*' .and. &
        input_line(2:2) /= '*') then
        
      keyword=get_keyword(input_line)
      if (lowercase(keyword) == 'end step') then

        write(iochannel(1),'(a)') '** fds2fem-cfilm-marker (do not remove)'
        write(iochannel(1),'(2(a))') '*Include, input=', trim(cfilm_file)

      end if
    end if

    if (.not. omit_lines) then
      write(iochannel(1),'(a)',iostat=ios) trim(abaqus_input(i))
    end if
       
    if (omit_lines) omit_lines=.false.

  end do add_lines

  close(unit=iochannel(1))

  deallocate(abaqus_input,stat=ios); call error_allocate(ios)

end subroutine dump_abaqus_cfilm

end module abaqus_output
