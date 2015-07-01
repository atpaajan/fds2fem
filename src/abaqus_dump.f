!------------------------------------------------
! Functions and subroutines for ABAQUS data dumps
!------------------------------------------------
module abaqus_dump

contains

subroutine abaqus_dump_module()
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
  write(*,'(a)') 'ABAQUS dump module'

  !call dump_abaqus_model_faces_vtk()
  !call dump_abaqus_nodal_areas_vtk()

  if (trim(dump_abaqus_nodes) == 'off' .and. &
      trim(dump_abaqus_data)  == 'off' .and. &
      trim(dump_abaqus_model) == 'off') then
    write(*,'(t3,a)') 'No dump requested'
    return
  end if

  if (trim(fds_chid) == trim(basename(fem_input_file))) then
    write(*,'(a)') 'WARNING: congruent FDS CHID and ABAQUS input file name'
  end if

  !---------------------
  ! Node coordinate dump
  !---------------------

  something_done=.false.
  if (abaqus_xyz_available) then
    select case (trim(dump_abaqus_nodes))
    case ('xyz')
      ! XYZ node dump
      call dump_abaqus_nodes_xyz()
    case ('vtk')
      ! VTK node dump
      call dump_abaqus_nodes_vtk()
    end select
    something_done=.true.
  end if

  !----------
  ! Data dump
  !----------

  if (abaqus_data_available) then
    select case (trim(dump_abaqus_data))
    case ('sdf')
      ! SDF data dump
      call dump_abaqus_data_sdf()
    case ('vtk')
      ! VTK data dump
      call dump_abaqus_data_vtk(1,ubound(abaqus_time,1))
      call dump_abaqus_data_faces_vtk(1,ubound(abaqus_time,1))
    end select
    something_done=.true.
  end if

  !-----------
  ! Model dump
  !-----------

  if (abaqus_model_available) then
    select case (trim(dump_abaqus_model))
    case ('vtk')
      ! VTK model dump
      call dump_abaqus_model_vtk()
    end select
    something_done=.true.
  end if

  if (.not. something_done) then
    write(*,'(t3,a)') 'Nothing to be done'
  end if

end subroutine abaqus_dump_module

!***************
! Auxiliary subs
!***************

subroutine dump_abaqus_nodes_xyz()
!-------------------------------------------
! Dump ABAQUS node coordinates in XYZ format
!-------------------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,nnodes
  character(len=chr80) :: filename

  nnodes=ubound(abaqus_xyz,1)

  filename=trim(basename(fem_input_file)) // '_nodes.xyz'
  open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  do i=1,nnodes
    write(iochannel(1),'(3(es15.7e3,1x))') (abaqus_xyz(i,j),j=1,3)
  end do 

  close(unit=iochannel(1))

  ! Status report
  write(*,'(t3,5(a))') trim(int2str(nnodes)), &
    ' node coordinates written in file ', &
    trim(quote(filename))

end subroutine dump_abaqus_nodes_xyz

subroutine dump_abaqus_data_sdf()
!-------------------------------
! Dump ABAQUS data in SDF format
!-------------------------------
  use abaqus_arrays
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

  do i=1,ntimes_abaqus
    write(iochannel(1),'(es15.7e3)',advance='no') abaqus_time(i)
    do j=1,nnodes_abaqus
      write(iochannel(1),'(1x,es15.7e3)',advance='no') abaqus_data(i,j)
    end do
    write(iochannel(1),'(a)') ''
  end do 

  close(unit=iochannel(1))

  ! Status report
  write(*,'(t3,5(a))') trim(int2str(ntimes_abaqus)), &
    ' time steps of data written in file ', &
    trim(quote(filename))

end subroutine dump_abaqus_data_sdf

subroutine dump_abaqus_nodes_vtk(filename)
!-------------------------------------------
! Dump ABAQUS node coordinates in VTK format
!-------------------------------------------
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

  nnodes=ubound(abaqus_xyz,1)
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
    write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (abaqus_xyz(i,j),j=1,3) 
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

end subroutine dump_abaqus_nodes_vtk

subroutine dump_abaqus_data_vtk(itime_begin,itime_end)
!-------------------------------
! Dump ABAQUS data in VTK format
!-------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,nnodes,ncell,ncell_size,itime
  integer :: itime_begin,itime_end,ntimes
  character(len=chr80) :: filename,ctime,f1,f2

  nnodes=ubound(abaqus_xyz,1)
  ncell=nnodes; ncell_size=2*ncell

  ! Exception handling
  if (itime_begin < 1 .or. itime_begin > ubound(abaqus_time,1)) then
    return
  end if

  if (itime_end < 1 .or. itime_end > ubound(abaqus_time,1)) then
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
      write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (abaqus_xyz(i,j),j=1,3) 
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
      write(iochannel(1),'(es15.7e3)') abaqus_data(itime,i)
    end do

    close(unit=iochannel(1))

  end do time_loop

  ! Status report
  ntimes=(itime_end-itime_begin)+1
  f1=trim(basename(fem_input_file)) // '_data_' // trim(int2str(itime_begin)) // '.vtk'
  f2=trim(basename(fem_input_file)) // '_data_' // trim(int2str(itime_end)) // '.vtk'

  write(*,'(t3,5(a))') trim(int2str(ntimes)), ' time steps of data written in files ', &
    trim(quote(f1)), '...', trim(quote(f2)) 

end subroutine dump_abaqus_data_vtk

subroutine dump_abaqus_data_faces_vtk(itime_begin,itime_end)
!-------------------------------
! Dump ABAQUS data in VTK format
!-------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,nnodes,ncell,ncell_size,itime
  integer :: itime_begin,itime_end,ntimes
  character(len=chr80) :: filename,ctime,f1,f2
  real(kind=rk) :: x

  nnodes=ubound(abaqus_xyz,1)
  ncell=ubound(abaqus_face_node,1)
  
  ncell_size=ncell
  do i=1,ncell
    ncell_size=ncell_size+abaqus_face_type(i)
  end do

  ! Exception handling
  if (itime_begin < 1 .or. itime_begin > ubound(abaqus_time,1)) then
    return
  end if

  if (itime_end < 1 .or. itime_end > ubound(abaqus_time,1)) then
    return
  end if

  if (itime_begin > itime_end) then
    return
  end if

  ! Write loop
  time_loop: do itime=itime_begin,itime_end

    write(ctime,'(i4.4)') itime
    filename=trim(basename(fem_input_file)) // '_data_face_' // trim(ctime) // '.vtk'
    open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
    if (ios /= 0) call error_open_file(filename)

    write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
    write(iochannel(1),'(a)') 'fds2fem: ABAQUS data dump'
    write(iochannel(1),'(a)') 'ASCII'
    write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

  write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(nnodes)), ' double'
  do i=1,nnodes
    write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (abaqus_xyz(i,j),j=1,3) 
  end do

  write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', &
    trim(int2str(ncell_size))

  do i=1,ncell
    select case (abaqus_face_type(i))
    case (3)
      write(iochannel(1),'(10(a))') trim(int2str(3)), ' ', (trim(int2str(abaqus_face_node(i,j)-1)) // ' ',j=1,3)
    case (4)
      write(iochannel(1),'(10(a))') trim(int2str(4)), ' ', (trim(int2str(abaqus_face_node(i,j)-1)) // ' ',j=1,4)
    end select
  end do

  write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
  do i=1,ncell
    select case (abaqus_face_type(i))
    case (3)
      write(iochannel(1),'(a)') '5'
    case (4)
      write(iochannel(1),'(a)') '9'
    end select
  end do

  write(iochannel(1),'(2(a))') 'CELL_DATA ', trim(int2str(ncell))
  write(iochannel(1),'(a)') 'SCALARS anonymous double 1'
  write(iochannel(1),'(a)') 'LOOKUP_TABLE default'

  do i=1,ncell
    select case (abaqus_face_type(i))
    case (3)
      x=0.0
      do j=1,3
        x=x+(1.0/3.0)*abaqus_data(itime,abaqus_face_node(i,j))
      end do
      write(iochannel(1),*) x
    case (4)
      x=0.0
      do j=1,4
        x=x+(1.0/4.0)*abaqus_data(itime,abaqus_face_node(i,j))
      end do
      write(iochannel(1),*) x
    end select
  end do

  close(unit=iochannel(1))

  end do time_loop

  ! Status report
  ntimes=(itime_end-itime_begin)+1
  f1=trim(basename(fem_input_file)) // '_data_face_' // trim(int2str(itime_begin)) // '.vtk'
  f2=trim(basename(fem_input_file)) // '_data_face_' // trim(int2str(itime_end)) // '.vtk'

  write(*,'(t3,5(a))') trim(int2str(ntimes)), ' time steps of data written in files ', &
    trim(quote(f1)), '...', trim(quote(f2)) 

end subroutine dump_abaqus_data_faces_vtk

subroutine dump_abaqus_model_vtk(filename)
!-------------------------------------------
! Dump ABAQUS node coordinates in VTK format
!-------------------------------------------
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
  write(iochannel(1),'(a)') 'fds2fem: ABAQUS model dump'
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

end subroutine dump_abaqus_model_vtk

subroutine dump_abaqus_model_faces_vtk(filename)
!-------------------------------------------
! Dump ABAQUS node coordinates in VTK format
!-------------------------------------------
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

  !nnodes=count(model_node_face_mask)
  nnodes=ubound(model_node_xyz,1)
  ncell=ubound(model_face_node,1)
  
  ncell_size=ncell
  do i=1,ncell
    ncell_size=ncell_size+model_face_type(i)
  end do

  if (present(filename)) then
    output_file=trim(filename)
  else
    output_file=trim(basename(fem_input_file)) // '_model_face.vtk'
  end if

  open(unit=iochannel(1),file=trim(output_file),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(output_file)

  write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
  write(iochannel(1),'(a)') 'fds2fem: ABAQUS model dump'
  write(iochannel(1),'(a)') 'ASCII'
  write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

  write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(nnodes)), ' double'
  do i=1,nnodes
    write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (model_node_xyz(i,j),j=1,3) 
  end do

  write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', &
    trim(int2str(ncell_size))

  do i=1,ncell
    select case (model_face_type(i))
    case (3)
      write(iochannel(1),'(10(a))') trim(int2str(3)), ' ', (trim(int2str(model_face_node(i,j)-1)) // ' ',j=1,3)
    case (4)
      write(iochannel(1),'(10(a))') trim(int2str(4)), ' ', (trim(int2str(model_face_node(i,j)-1)) // ' ',j=1,4)
    end select
  end do

  write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
  do i=1,ncell
    select case (model_face_type(i))
    case (3)
      write(iochannel(1),'(a)') '5'
    case (4)
      write(iochannel(1),'(a)') '9'
    end select
  end do

  write(iochannel(1),'(2(a))') 'CELL_DATA ', trim(int2str(ncell))
  write(iochannel(1),'(a)') 'SCALARS anonymous double 1'
  write(iochannel(1),'(a)') 'LOOKUP_TABLE default'

  do i=1,ncell
    write(iochannel(1),*) model_face_area(i) 
  end do

  close(unit=iochannel(1))

  write(*,'(t3,5(a))') trim(int2str(ncell)), &
    ' element faces written in file ', &
    trim(quote(output_file))

end subroutine dump_abaqus_model_faces_vtk

subroutine dump_abaqus_model_nodes_vtk(filename)
!-------------------------------------------
! Dump ABAQUS node coordinates in VTK format
!-------------------------------------------
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
  ncell=nnodes; ncell_size=2*ncell

  if (present(filename)) then
    output_file=trim(filename)
  else
    output_file=trim(basename(fem_input_file)) // '_nodes.vtk'
  end if

  open(unit=iochannel(1),file=trim(output_file),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(output_file)

  write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
  write(iochannel(1),'(a)') 'fds2fem: ABAQUS model dump'
  write(iochannel(1),'(a)') 'ASCII'
  write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

  write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(nnodes)), ' double'
  do i=1,nnodes
    write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (model_node_xyz(i,j),j=1,3) 
  end do

  write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', &
    trim(int2str(ncell_size))

  do i=1,ncell
    write(iochannel(1),'(3(a))') trim(int2str(1)), ' ', trim(int2str(i-1))
  end do

  write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
  do i=1,ncell
      write(iochannel(1),'(a)') '1'
  end do

  write(iochannel(1),'(2(a))') 'CELL_DATA ', trim(int2str(ncell))
  write(iochannel(1),'(a)') 'SCALARS anonymous double 1'
  write(iochannel(1),'(a)') 'LOOKUP_TABLE default'

  do i=1,ncell
    write(iochannel(1),*) model_node_elements(i) 
  end do

  close(unit=iochannel(1))

  write(*,'(t3,5(a))') trim(int2str(nnodes)), &
    ' node coordinates written in file ', &
    trim(quote(output_file))

end subroutine dump_abaqus_model_nodes_vtk

subroutine dump_abaqus_nodal_areas_vtk(filename)
!-------------------------------------------
! Dump ABAQUS node coordinates in VTK format
!-------------------------------------------
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

  nnodes=ubound(abaqus_xyz,1)
  ncell=nnodes; ncell_size=2*ncell

  output_file=trim(basename(fem_input_file)) // '_nodal_areas.vtk'

  open(unit=iochannel(1),file=trim(output_file),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(output_file)

  write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
  write(iochannel(1),'(a)') 'fds2fem: ABAQUS node coordinate dump'
  write(iochannel(1),'(a)') 'ASCII'
  write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

  write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(ncell)), ' double'
  do i=1,nnodes
    write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (abaqus_xyz(i,j),j=1,3) 
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

  write(iochannel(1),'(2(a))') 'CELL_DATA ', trim(int2str(ncell))
  write(iochannel(1),'(a)') 'SCALARS anonymous double 1'
  write(iochannel(1),'(a)') 'LOOKUP_TABLE default'

  do i=1,nnodes
    write(iochannel(1),'(es15.7e3)') abaqus_node_area(i)
  end do

  close(unit=iochannel(1))

  write(*,'(t3,5(a))') trim(int2str(nnodes)), &
    ' node coordinates written in file ', &
    trim(quote(output_file))

end subroutine dump_abaqus_nodal_areas_vtk

subroutine dump_abaqus_node_elements_vtk(filename)
!-------------------------------------------
! Dump ABAQUS node coordinates in VTK format
!-------------------------------------------
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

  nnodes=ubound(abaqus_xyz,1)
  ncell=nnodes; ncell_size=2*ncell

  output_file=trim(basename(fem_input_file)) // '_node_elements.vtk'

  open(unit=iochannel(1),file=trim(output_file),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(output_file)

  write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
  write(iochannel(1),'(a)') 'fds2fem: ABAQUS node coordinate dump'
  write(iochannel(1),'(a)') 'ASCII'
  write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

  write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(ncell)), ' double'
  do i=1,nnodes
    write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (abaqus_xyz(i,j),j=1,3) 
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

  write(iochannel(1),'(2(a))') 'CELL_DATA ', trim(int2str(ncell))
  write(iochannel(1),'(a)') 'SCALARS anonymous double 1'
  write(iochannel(1),'(a)') 'LOOKUP_TABLE default'

  do i=1,nnodes
    write(iochannel(1),'(i3)') abaqus_node_elements(i)
  end do

  close(unit=iochannel(1))

  write(*,'(t3,5(a))') trim(int2str(nnodes)), &
    ' node coordinates written in file ', &
    trim(quote(output_file))

end subroutine dump_abaqus_node_elements_vtk

end module abaqus_dump
