!---------------------------------------------
! Functions and subroutines for FEM data dumps
!---------------------------------------------
module fem_dump

contains

subroutine fem_dump_module()
!---------
! Main sub
!---------
  use fds_head_arrays
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  logical :: something_done

  write(*,'(a)') ''
  write(*,'(a)') 'FEM dump module'

  if (trim(dump_fem_nodes) == 'off' .and. &
      trim(dump_fem_data)  == 'off' .and. &
      trim(dump_fem_model) == 'off') then
    write(*,'(t3,a)') 'No dump requested'
    return
  end if

  if (trim(fem_mode) == 'abaqus') then
    if (trim(fds_chid) == trim(basename(fem_input_file))) then
      write(*,'(a)') 'WARNING: congruent FDS CHID and ABAQUS input file name'
    end if
  end if

  !---------------------
  ! Node coordinate dump
  !---------------------

  something_done=.false.
  if (fem_xyz_available) then
    select case (trim(dump_fem_nodes))
    case ('xyz')
      ! XYZ node dump
      call dump_fem_nodes_xyz()
    case ('vtk')
      ! VTK node dump
      call dump_fem_nodes_vtk()
    end select
    something_done=.true.
  end if

  !----------
  ! Data dump
  !----------

  if (fem_data_available) then
    select case (trim(dump_fem_data))
    case ('sdf')
      ! SDF data dump
      call dump_fem_data_sdf()
    case ('vtk')
      ! VTK data dump
      call dump_fem_data_vtk(1,ubound(fem_time,1))
      if (read_hcoeff) then
        call dump_fem_hcoeff_vtk(1,ubound(fem_time,1))
      end if
    end select
    something_done=.true.
  end if

  !-----------
  ! Model dump
  !-----------

  if (fem_model_available) then
    select case (trim(dump_fem_model))
    case ('vtk')
      ! VTK model dump
      call dump_fem_model_vtk()
    end select
    something_done=.true.
  end if

  if (.not. something_done) then
    write(*,'(t3,a)') 'Nothing to be done'
  end if

  !----------
  ! Otherwise
  !----------

  if (.not. something_done) then
    write(*,'(t3,a)') 'Nothing to be done'
  end if

end subroutine fem_dump_module

!***************
! Auxiliary subs
!***************

subroutine dump_fem_nodes_xyz()
!----------------------------------------
! Dump FEM node coordinates in XYZ format
!----------------------------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,nnodes
  character(len=chr80) :: filename

  nnodes=ubound(fem_xyz,1)

  filename=trim(basename(fem_input_file)) // '_nodes.xyz'
  open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  do i=1,nnodes
    write(iochannel(1),'(3(es15.7e3,1x))') (fem_xyz(i,j),j=1,3)
  end do 

  close(unit=iochannel(1))

  ! Status report
  write(*,'(t3,5(a))') trim(int2str(nnodes)), &
    ' node coordinates written in file ', &
    trim(quote(filename))

end subroutine dump_fem_nodes_xyz

subroutine dump_fem_data_sdf()
!----------------------------
! Dump FEM data in SDF format
!----------------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,ntimes,nnodes
  character(len=chr80) :: filename

  filename=trim(basename(fem_input_file)) // '_nodes.sdf'
  open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  do i=1,ntimes_fem
    write(iochannel(1),'(es15.7e3)',advance='no') fem_time(i)
    do j=1,nnodes_fem
      write(iochannel(1),'(1x,es15.7e3)',advance='no') fem_data(i,j)
    end do
    write(iochannel(1),'(a)') ''
  end do 

  close(unit=iochannel(1))

  ! Status report
  write(*,'(t3,5(a))') trim(int2str(ntimes_fem)), &
    ' time steps of data written in file ', &
    trim(quote(filename))

end subroutine dump_fem_data_sdf

subroutine dump_fem_nodes_vtk(filename)
!----------------------------------------
! Dump FEM node coordinates in VTK format
!----------------------------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,nnodes,ncell,ncell_size
  character(len=chr80), optional :: filename
  character(len=chr80) :: output_file

  nnodes=ubound(fem_xyz,1)
  ncell=nnodes; ncell_size=2*ncell

  if (present(filename)) then
    output_file=trim(filename)
  else
    output_file=trim(basename(fem_input_file)) // '_nodes.vtk'
  end if

  open(unit=iochannel(1),file=trim(output_file),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(output_file)

  write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
  write(iochannel(1),'(a)') 'fds2fem: ABAQUS node coordinate dump'
  write(iochannel(1),'(a)') 'ASCII'
  write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

  write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(ncell)), ' double'
  do i=1,nnodes
    write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (fem_xyz(i,j),j=1,3) 
  end do

  write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', &
    trim(int2str(ncell_size))
  do i=1,nnodes
    write(iochannel(1),'(3(a))') trim(int2str(1)), ' ', trim(int2str(i-1))
  end do

  write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
  do i=1,nnodes
    write(iochannel(1),'(a)') '1'
  end do

  close(unit=iochannel(1))

  write(*,'(t3,5(a))') trim(int2str(nnodes)), &
    ' node coordinates written in file ', &
    trim(quote(output_file))

end subroutine dump_fem_nodes_vtk

subroutine dump_fem_data_vtk(itime_begin,itime_end)
!----------------------------
! Dump FEM data in VTK format
!----------------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,nnodes,ncell,ncell_size,itime
  integer :: itime_begin,itime_end,ntimes
  character(len=chr80) :: filename,ctime,f1,f2

  nnodes=ubound(fem_xyz,1)
  ncell=nnodes; ncell_size=2*ncell

  ! Exception handling
  if (itime_begin < 1 .or. itime_begin > ubound(fem_time,1)) then
    return
  end if

  if (itime_end < 1 .or. itime_end > ubound(fem_time,1)) then
    return
  end if

  if (itime_begin > itime_end) then
    return
  end if

  ! Write loop
  time_loop: do itime=itime_begin,itime_end

    write(ctime,'(i4.4)') itime
    filename=trim(basename(fem_input_file)) // '_data_' // trim(ctime) // '.vtk'
    open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
    if (ios /= 0) call error_open_file(filename)

    write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
    write(iochannel(1),'(a)') 'fds2fem: ABAQUS data dump'
    write(iochannel(1),'(a)') 'ASCII'
    write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

    write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(ncell)), ' double'
    do i=1,nnodes
      write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (fem_xyz(i,j),j=1,3) 
    end do

    write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', trim(int2str(ncell_size))
    do i=1,nnodes
      write(iochannel(1),'(3(a))') trim(int2str(1)), ' ', trim(int2str(i-1))
    end do

    write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
    do i=1,nnodes
      write(iochannel(1),'(a)') '1'
    end do

    write(iochannel(1),'(2(a))') 'CELL_DATA ', trim(int2str(ncell))
    write(iochannel(1),'(a)') 'SCALARS anonymous double 1'
    write(iochannel(1),'(a)') 'LOOKUP_TABLE default'

    do i=1,nnodes
      write(iochannel(1),'(es15.7e3)') fem_data(itime,i)
    end do

    close(unit=iochannel(1))

  end do time_loop

  ! Status report
  ntimes=(itime_end-itime_begin)+1
  f1=trim(basename(fem_input_file)) // '_data_' // trim(int2str(itime_begin)) // '.vtk'
  f2=trim(basename(fem_input_file)) // '_data_' // trim(int2str(itime_end)) // '.vtk'

  write(*,'(t3,5(a))') trim(int2str(ntimes)), ' time steps of data written in files ', &
    trim(quote(f1)), '...', trim(quote(f2)) 

end subroutine dump_fem_data_vtk

subroutine dump_fem_hcoeff_vtk(itime_begin,itime_end)
!----------------------------
! Dump FEM data in VTK format
!----------------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,nnodes,ncell,ncell_size,itime
  integer :: itime_begin,itime_end,ntimes
  character(len=chr80) :: filename,ctime,f1,f2

  nnodes=ubound(fem_xyz,1)
  ncell=nnodes; ncell_size=2*ncell

  ! Exception handling
  if (itime_begin < 1 .or. itime_begin > ubound(fem_time,1)) then
    return
  end if

  if (itime_end < 1 .or. itime_end > ubound(fem_time,1)) then
    return
  end if

  if (itime_begin > itime_end) then
    return
  end if

  ! Write loop
  time_loop: do itime=itime_begin,itime_end

    write(ctime,'(i4.4)') itime
    filename=trim(basename(fem_input_file)) // '_hcoeff_' // trim(ctime) // '.vtk'
    open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
    if (ios /= 0) call error_open_file(filename)

    write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
    write(iochannel(1),'(a)') 'fds2fem: ABAQUS hcoeff dump'
    write(iochannel(1),'(a)') 'ASCII'
    write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

    write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(ncell)), ' double'
    do i=1,nnodes
      write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (fem_xyz(i,j),j=1,3) 
    end do

    write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', trim(int2str(ncell_size))
    do i=1,nnodes
      write(iochannel(1),'(3(a))') trim(int2str(1)), ' ', trim(int2str(i-1))
    end do

    write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
    do i=1,nnodes
      write(iochannel(1),'(a)') '1'
    end do

    write(iochannel(1),'(2(a))') 'CELL_DATA ', trim(int2str(ncell))
    write(iochannel(1),'(a)') 'SCALARS anonymous double 1'
    write(iochannel(1),'(a)') 'LOOKUP_TABLE default'

    do i=1,nnodes
      write(iochannel(1),'(es15.7e3)') fem_hcoeff(itime,i)
    end do

    close(unit=iochannel(1))

  end do time_loop

  ! Status report
  ntimes=(itime_end-itime_begin)+1
  f1=trim(basename(fem_input_file)) // '_hcoeff_' // trim(int2str(itime_begin)) // '.vtk'
  f2=trim(basename(fem_input_file)) // '_hcoeff_' // trim(int2str(itime_end)) // '.vtk'

  write(*,'(t3,5(a))') trim(int2str(ntimes)), ' time steps of hcoeff-data written in files ', &
    trim(quote(f1)), '...', trim(quote(f2)) 

end subroutine dump_fem_hcoeff_vtk

subroutine dump_fem_model_vtk(filename)
!-----------------------------
! Dump FEM model in VTK format
!-----------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,nnodes,ncell,ncell_size
  character(len=chr80), optional :: filename
  character(len=chr80) :: output_file

  nnodes=ubound(model_node_xyz,1)
  ncell=ubound(model_element_int_node,1)
  
  ncell_size=ncell
  do i=1,ncell
   select case (model_element_type(i))
    case (1)
      ncell_size=ncell_size+8
    case (2)
      ncell_size=ncell_size+4
    case (3)
      ncell_size=ncell_size+3
    case (4)
      ncell_size=ncell_size+4
    end select
  end do

  if (present(filename)) then
    output_file=trim(filename)
  else
    output_file=trim(basename(fem_input_file)) // '_model.vtk'
  end if

  open(unit=iochannel(1),file=trim(output_file),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(output_file)

  write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
  write(iochannel(1),'(a)') 'fds2fem: FEM model dump'
  write(iochannel(1),'(a)') 'ASCII'
  write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

  write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(nnodes)), ' double'
  do i=1,nnodes
    write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (model_node_xyz(i,j),j=1,3) 
  end do

  write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', &
    trim(int2str(ncell_size))

  do i=1,ncell
    select case (model_element_type(i))
    case (1)
      write(iochannel(1),'(10(a))') trim(int2str(8)), ' ', (trim(int2str(model_element_int_node(i,j)-1)) // ' ',j=1,8)
    case (2)
      write(iochannel(1),'(10(a))') trim(int2str(4)), ' ', (trim(int2str(model_element_int_node(i,j)-1)) // ' ',j=1,4)
    case (3)
      write(iochannel(1),'(10(a))') trim(int2str(3)), ' ', (trim(int2str(model_element_int_node(i,j)-1)) // ' ',j=1,3)
    case (4)
      write(iochannel(1),'(10(a))') trim(int2str(4)), ' ', (trim(int2str(model_element_int_node(i,j)-1)) // ' ',j=1,4)
    end select
  end do

  write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
  do i=1,ncell
    select case (model_element_type(i))
    case (1)
      write(iochannel(1),'(a)') '12'
    case (2)
      write(iochannel(1),'(a)') '10'
    case (3)
      write(iochannel(1),'(a)') '5'
    case (4)
      write(iochannel(1),'(a)') '9'
    end select
  end do

  close(unit=iochannel(1))

  write(*,'(t3,5(a))') trim(int2str(nnodes)), &
    ' node coordinates written in file ', &
    trim(quote(output_file))

end subroutine dump_fem_model_vtk

end module fem_dump
