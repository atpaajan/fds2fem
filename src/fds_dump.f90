!---------------------------------------------
! Functions and subroutines for FDS data dumps
!---------------------------------------------
module fds_dump

contains

subroutine fds_dump_module()
!---------
! Main sub
!---------
  use global_constants
  use global_variables
  use mapping_arrays
  implicit none

  logical :: something_done

  write(*,'(a)') ''
  write(*,'(a)') 'FDS dump module'

  if (trim(dump_fds_nodes) == 'off' .and. &
      trim(dump_fds_data)  == 'off' .and. &
      trim(dump_fds_model) == 'off') then
    write(*,'(t3,a)') 'No dump requested'
    return
  end if

  !---------------------
  ! Node coordinate dump
  !---------------------

  something_done=.false.
  if (fds_xyz_available) then
    select case (trim(dump_fds_nodes))
    case ('xyz')
      ! XYZ node dump
      call dump_fds_nodes_xyz() 
    case ('vtk')
      ! VTK node dump
      call dump_fds_nodes_vtk() 
    end select
    something_done=.true.
  end if

  !----------
  ! Data dump
  !----------

  if (fds_data_available) then
    select case (trim(dump_fds_data))
    case ('sdf')
      ! SDF data dump
      call dump_fds_data_sdf()
    case ('vtk')
      ! VTK data dump
      call dump_fds_data_vtk(1,ubound(fds_time,1))
      !call dump_fds_data_faces_vtk(1,ubound(fds_time,1))
      call dump_fds_element_data_vtk(1,ubound(fds_time,1),.true.)
    end select
    something_done=.true.
  end if

  !-----------
  ! Model dump
  !-----------
  
  if (fds_model_available) then
    select case (trim(dump_fds_model))
    case ('vtk')
      ! VTK model dump
      call construct_fds_model()
      call dump_fds_model_vtk()
    end select
    something_done=.true.
  end if

  if (.not. something_done) then
    write(*,'(t3,a)') 'Nothing to be done'
  end if

end subroutine fds_dump_module

!***************
! Auxiliary subs
!***************

subroutine dump_fds_devc_nodes_xyz()
!------------------------------------------
! Dump FDS device coordinates in XYZ format
!------------------------------------------
  use error_messages
  use fds_devc_arrays
  use fds_head_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: j,ios,idevc,ndevc
  character(len=chr80) :: filename

  filename=trim(fds_chid) // "_devc.xyz"
  open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  ndevc=ubound(fds_devc_name,1)

  do idevc=1,ndevc
    write(iochannel(1),'(3(es15.7e3,1x),i2,1x,3(a))') (fds_devc_xyz(idevc,j),j=1,3), &
      fds_devc_ior(idevc), char(34), trim(fds_devc_name(idevc)), char(34)
  end do

  close(unit=iochannel(1))

  ! Status report
  write(*,'(t3,5(a))') trim(int2str(ndevc)), &
    ' FDS device coordinates written in file ', &
    trim(quote(filename))

end subroutine dump_fds_devc_nodes_xyz

subroutine dump_fds_devc_data_sdf
!-----------------------------------
! Dump FDS device data in SDF format
!-----------------------------------
  use error_messages
  use fds_devc_arrays
  use fds_head_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,j,ios,ndevc,ntimes
  character(len=chr80) :: filename

  ntimes=ubound(fds_devc_time,1)
  ndevc=ubound(fds_devc_name,1)

  filename=trim(fds_chid) // "_devc.sdf"
  open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  do i=1,ntimes
    write(iochannel(1),'(es15.7e3)',advance='no') fds_devc_time(i)
    do j=1,ndevc
      write(iochannel(1),'(1x,es15.7e3)',advance='no') fds_devc_data(i,j)
    end do
    write(iochannel(1),'(a)') ''
  end do

  close(unit=iochannel(1))

  ! Status report
  write(*,'(t3,5(a))') trim(int2str(ntimes)), &
    ' time steps of FDS device data written in file ', &
    trim(quote(filename))

end subroutine dump_fds_devc_data_sdf

subroutine dump_fds_devc_nodes_vtk()
!------------------------------------------
! Dump FDS device coordinates in VTK format
!------------------------------------------
  use error_messages
  use fds_devc_arrays
  use fds_head_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,j,ios,ndevc,ncell,ncell_size
  character(len=chr80) :: filename

  ndevc=ubound(fds_devc_name,1)
  ncell=ndevc; ncell_size=2*ncell

  filename=trim(fds_chid) // '_devc.vtk'
  open(unit=iochannel(1),file=filename,status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
  write(iochannel(1),'(a)') 'fds2fem: FDS device coordinate dump'
  write(iochannel(1),'(a)') 'ASCII'
  write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

  write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(ncell)), ' double'
  do i=1,ndevc
    write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (fds_devc_xyz(i,j),j=1,3) 
  end do

  write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', trim(int2str(ncell_size))
  do i=1,ndevc
    write(iochannel(1),'(3(a))') trim(int2str(1)), ' ', trim(int2str(i-1))
  end do

  write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
  do i=1,ndevc
    write(iochannel(1),'(a)') '1'
  end do

  close(unit=iochannel(1))

  write(*,'(t3,a,a,a,a,a)') trim(int2str(ndevc)), &
    ' FDS device coordinates written in file ', &
    trim(quote(filename))

end subroutine dump_fds_devc_nodes_vtk

subroutine dump_fds_devc_data_vtk(itime_begin,itime_end)
!-----------------------------------
! Dump FDS device data in VTK format
!-----------------------------------
  use error_messages
  use fds_devc_arrays
  use fds_head_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,j,ios,ndevc,ncell,ncell_size,itime
  integer :: itime_begin,itime_end,ntimes
  character(len=chr80) :: filename,ctime,f1,f2

  ndevc=ubound(fds_devc_name,1)
  ncell=ndevc; ncell_size=2*ncell

  ! Exception handling
  if (itime_begin < 1 .or. itime_begin > ubound(fds_devc_time,1)) then
    return
  end if

  if (itime_end < 1 .or. itime_end > ubound(fds_devc_time,1)) then
    return
  end if

  if (itime_begin > itime_end) then
    return
  end if

  ! Write loop
  time_loop: do itime=itime_begin,itime_end

    write(ctime,'(i4.4)') itime
    filename=trim(fds_chid) // '_devc_' // trim(ctime) // '.vtk'
    open(unit=iochannel(1),file=filename,status='replace',iostat=ios)
    if (ios /= 0) call error_open_file(filename)

    write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
    write(iochannel(1),'(a)') 'fds2fem: FDS device data dump'
    write(iochannel(1),'(a)') 'ASCII'
    write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

    write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(ncell)), ' double'
    do i=1,ndevc
      write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (fds_devc_xyz(i,j),j=1,3) 
    end do

    write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', trim(int2str(ncell_size))
    do i=1,ndevc
      write(iochannel(1),'(3(a))') trim(int2str(1)), ' ', trim(int2str(i-1))
    end do

    write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
    do i=1,ndevc
      write(iochannel(1),'(a)') '1'
    end do

    write(iochannel(1),'(2(a))') 'CELL_DATA ', trim(int2str(ncell))
    write(iochannel(1),'(a)') 'SCALARS anonymous double 1'
    write(iochannel(1),'(a)') 'LOOKUP_TABLE default'

    do i=1,ndevc
      write(iochannel(1),'(es15.7e3)') fds_devc_data(itime,i)
    end do

    close(unit=iochannel(1))

  end do time_loop

  ! Status report
  ntimes=(itime_end-itime_begin)+1
  f1=trim(fds_chid) // '_devc_' // trim(int2str(itime_begin)) // '.vtk'
  f2=trim(fds_chid) // '_devc_' // trim(int2str(itime_end)) // '.vtk'

  write(*,'(t3,5(a))') trim(int2str(ntimes)), ' time steps of FDS device data written in files ', &
    trim(quote(f1)), '...', trim(quote(f2)) 

end subroutine dump_fds_devc_data_vtk

subroutine dump_fds_bndf_nodes_xyz()
!-------------------------------------------------
! Dump FDS boundary node coordinates in XYZ format
!-------------------------------------------------
  use error_messages
  use fds_bndf_arrays
  use fds_head_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,j,ios,nnode
  character(len=chr80) :: filename,fmt_1

  filename=trim(fds_chid) // "_bndf.xyz"
  open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  nnode=ubound(fds_bndf_xyz,1)

  fmt_1='(2(es15.7e3,1x),es15.7e3,i3)'
  do i=1,nnode
    write(iochannel(1),fmt_1) (fds_bndf_xyz(i,j),j=1,3), fds_bndf_ior(i)
  end do

  close(unit=iochannel(1))

  ! Status report
  write(*,'(t3,5(a))') trim(int2str(nnode)), &
    ' FDS boundary node coordinates written in file ', &
    trim(quote(filename))

end subroutine dump_fds_bndf_nodes_xyz

subroutine dump_fds_bndf_data_sdf
!-------------------------------------
! Dump FDS boundary data in SDF format
!-------------------------------------
  use error_messages
  use fds_bndf_arrays
  use fds_head_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,j,ios,nnodes,ntimes
  character(len=chr80) :: filename

  ntimes=ubound(fds_bndf_time,1)
  nnodes=ubound(fds_bndf_data,2)

  filename=trim(fds_chid) // "_bndf.sdf"
  open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  do i=1,ntimes
    write(iochannel(1),'(es15.7e3)',advance='no') fds_bndf_time(i)
    do j=1,nnodes
      write(iochannel(1),'(1x,es15.7e3)',advance='no') fds_bndf_data(i,j)
    end do
    write(iochannel(1),'(a)') ''
  end do

  close(unit=iochannel(1))

  ! Status report
  write(*,'(t3,5(a))') trim(int2str(ntimes)), &
    ' time steps of FDS boundary data written in file ', &
    trim(quote(filename))

end subroutine dump_fds_bndf_data_sdf

subroutine dump_fds_bndf_nodes_vtk()
!-------------------------------------------------
! Dump FDS boundary node coordinates in VTK format
!-------------------------------------------------
  use error_messages
  use fds_bndf_arrays
  use fds_head_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,j,ios,nnode,ncell,ncell_size
  character(len=chr80) :: filename

  nnode=ubound(fds_bndf_xyz,1)
  ncell=nnode; ncell_size=2*ncell

  filename=trim(fds_chid) // '_bndf.vtk'
  open(unit=iochannel(1),file=filename,status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
  write(iochannel(1),'(a)') 'fds2fem: FDS boundary node coordinate dump'
  write(iochannel(1),'(a)') 'ASCII'
  write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

  write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(nnode)), ' double'
  do i=1,nnode
    write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (fds_bndf_xyz(i,j),j=1,3) 
  end do

  write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', trim(int2str(ncell_size))
  do i=1,nnode
    write(iochannel(1),'(3(a))') trim(int2str(1)), ' ', trim(int2str(i-1))
  end do

  write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
  do i=1,nnode
    write(iochannel(1),'(a)') '1'
  end do

  close(unit=iochannel(1))

  write(*,'(t3,5(a))') trim(int2str(nnode)), &
    ' FDS boundary node coordinates written in file ', &
    trim(quote(filename))

end subroutine dump_fds_bndf_nodes_vtk

subroutine dump_fds_bndf_data_vtk(itime_begin,itime_end)
!-------------------------------------
! Dump FDS boundary data in VTK format
!-------------------------------------
  use error_messages
  use fds_bndf_arrays
  use fds_head_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,j,ios,nnode,ncell,ncell_size,itime
  integer :: itime_begin,itime_end,ntimes
  character(len=chr80) :: filename,ctime,f1,f2

  nnode=ubound(fds_bndf_xyz,1)
  ncell=nnode; ncell_size=2*ncell

  ! Exception handling
  if (itime_begin < 1 .or. itime_begin > ubound(fds_bndf_time,1)) then
    return
  end if

  if (itime_end < 1 .or. itime_end > ubound(fds_bndf_time,1)) then
    return
  end if

  if (itime_begin > itime_end) then
    return
  end if

  ! Write loop
  time_loop: do itime=itime_begin,itime_end

    write(ctime,'(i4.4)') itime
    filename=trim(fds_chid) // '_bndf_' // trim(ctime) // '.vtk'
    open(unit=iochannel(1),file=filename,status='replace',iostat=ios)
    if (ios /= 0) call error_open_file(filename)
    
    write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
    write(iochannel(1),'(a)') 'fds2fem: FDS boundary data dump'
    write(iochannel(1),'(a)') 'ASCII'
    write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

    write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(nnode)), ' double'
    do i=1,nnode
      write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (fds_bndf_xyz(i,j),j=1,3) 
    end do

    write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', trim(int2str(ncell_size))
    do i=1,nnode
      write(iochannel(1),'(3(a))') trim(int2str(1)), ' ', trim(int2str(i-1))
    end do

    write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
    do i=1,nnode
      write(iochannel(1),'(a)') '1'
    end do

    write(iochannel(1),'(2(a))') 'CELL_DATA ', trim(int2str(ncell))
    write(iochannel(1),'(a)') 'SCALARS anonymous double 1'
    write(iochannel(1),'(a)') 'LOOKUP_TABLE default'

    do j=1,nnode
      write(iochannel(1),'(es15.7e3)') fds_bndf_data(itime,j)
    end do

    close(unit=iochannel(1))

  end do time_loop

  ! Status report
  ntimes=(itime_end-itime_begin)+1
  f1=trim(fds_chid) // '_bndf_' // trim(int2str(itime_begin)) // '.vtk'
  f2=trim(fds_chid) // '_bndf_' // trim(int2str(itime_end)) // '.vtk'

  write(*,'(t3,5(a))') trim(int2str(ntimes)), ' time steps of FDS boundary data written in files ', &
    trim(quote(f1)), '...', trim(quote(f2)) 

end subroutine dump_fds_bndf_data_vtk

subroutine dump_fds_nodes_xyz()
!----------------------------------------
! Dump FDS node coordinates in XYZ format
!----------------------------------------
  use error_messages
  use fds_head_arrays
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: ios,j,inode,nnodes
  character(len=chr80) :: filename
  
  nnodes=ubound(fds_xyz,1)

  filename=trim(fds_chid) // "_nodes.xyz"
  open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  do inode=1,nnodes
    write(iochannel(1),'(3(es15.7e3,1x),i2)') (fds_xyz(inode,j),j=1,3), fds_ior(inode)
  end do

  close(unit=iochannel(1))

  ! Status report
  write(*,'(t3,5(a))') trim(int2str(nnodes)), &
    ' node coordinates written in file ', &
    trim(quote(filename))

end subroutine dump_fds_nodes_xyz

subroutine dump_fds_data_sdf
!----------------------------
! Dump FDS data in SDF format
!----------------------------
  use error_messages
  use fds_head_arrays
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,nnodes,ntimes
  character(len=chr80) :: filename

  ntimes=ubound(fds_time,1)
  nnodes=ubound(fds_data,2)

  filename=trim(fds_chid) // "_data.sdf"
  open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  do i=1,ntimes
    write(iochannel(1),'(es15.7e3)',advance='no') fds_time(i)
    do j=1,nnodes
      write(iochannel(1),'(1x,es15.7e3)',advance='no') fds_data(i,j)
    end do
    write(iochannel(1),'(a)') ''
  end do

  close(unit=iochannel(1))

  ! Status report
  write(*,'(t3,5(a))') trim(int2str(ntimes)), &
    ' time steps written in file ', &
    trim(quote(filename))

end subroutine dump_fds_data_sdf

subroutine dump_fds_nodes_vtk
!----------------------------------------
! Dump FDS node coordinates in VTK format
!----------------------------------------
  use error_messages
  use fds_head_arrays
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,nnodes,ncell,ncell_size
  character(len=chr80) :: filename

  nnodes=ubound(fds_xyz,1)
  ncell=nnodes; ncell_size=2*ncell

  filename=trim(fds_chid) // '_nodes.vtk'
  open(unit=iochannel(1),file=filename,status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
  write(iochannel(1),'(a)') 'fds2fem: FDS node coordinate dump'
  write(iochannel(1),'(a)') 'ASCII'
  write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

  write(iochannel(1),'(a,a,a)') 'POINTS ', trim(int2str(ncell)), ' double'
  do i=1,nnodes
    write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (fds_xyz(i,j),j=1,3) 
  end do

  write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', trim(int2str(ncell_size))
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
    trim(quote(filename))

end subroutine dump_fds_nodes_vtk

subroutine dump_fds_data_vtk(itime_begin,itime_end)
!----------------------------
! Dump FDS data in VTK format
!----------------------------
  use error_messages
  use fds_head_arrays
  use global_constants
  use global_variables
  use string_handling
  use mapping_arrays
  implicit none

  integer :: i,j,ios,nnode,ncell,ncell_size,itime
  integer :: itime_begin,itime_end,ntimes
  character(len=chr80) :: filename,ctime,f1,f2

  nnode=ubound(fds_xyz,1)
  ncell=nnode; ncell_size=2*ncell

  ! Exception handling
  if (itime_begin < 1 .or. itime_begin > ubound(fds_time,1)) then
    return
  end if

  if (itime_end < 1 .or. itime_end > ubound(fds_time,1)) then
    return
  end if

  if (itime_begin > itime_end) then
    return
  end if

  ! Write loop
  time_loop: do itime=itime_begin,itime_end

    write(ctime,'(i4.4)') itime
    filename=trim(fds_chid) // '_data_' // trim(ctime) // '.vtk'
    open(unit=iochannel(1),file=filename,status='replace',iostat=ios)
    if (ios /= 0) call error_open_file(filename)

    write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
    write(iochannel(1),'(a)') 'fds2fem: FDS data dump'
    write(iochannel(1),'(a)') 'ASCII'
    write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

    write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(nnode)), ' double'
    do i=1,nnode
      write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (fds_xyz(i,j),j=1,3) 
    end do

    write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', trim(int2str(ncell_size))
    do i=1,nnode
      write(iochannel(1),'(3(a))') trim(int2str(1)), ' ', trim(int2str(i-1))
    end do

    write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
    do i=1,nnode
      write(iochannel(1),'(a)') '1'
    end do

    write(iochannel(1),'(2(a))') 'CELL_DATA ', trim(int2str(ncell))
    write(iochannel(1),'(a)') 'SCALARS anonymous double 1'
    write(iochannel(1),'(a)') 'LOOKUP_TABLE default'

    do j=1,nnode
      write(iochannel(1),'(es15.7e3)') fds_data(itime,j)
    end do

    close(unit=iochannel(1))

  end do time_loop

  ! Status report
  ntimes=(itime_end-itime_begin)+1
  f1=trim(fds_chid) // '_data_' // trim(int2str(itime_begin)) // '.vtk'
  f2=trim(fds_chid) // '_data_' // trim(int2str(itime_end)) // '.vtk'

  write(*,'(t3,5(a))') trim(int2str(ntimes)), ' time steps written in files ', &
    trim(quote(f1)), '...', trim(quote(f2)) 

end subroutine dump_fds_data_vtk

subroutine dump_fds_element_data_vtk(itime_begin,itime_end,average)
!----------------------------------------
! Dump FDS node coordinates in VTK format
!----------------------------------------
  use error_messages
  use fds_head_arrays
  use fds_bndf_arrays
  use global_constants
  use global_variables
  use string_handling
  use mapping_arrays
  implicit none

  integer :: i,j,ios,nnodes,ncell,ncell_size
  integer :: itime,itime_begin,itime_end,ntimes
  character(len=chr80) :: filename,ctime,f1,f2
  logical :: average
  real(kind=rk) :: dsum

  nnodes=ubound(fds_bndf_xyz,1)
  ncell=ubound(fds_bndf_el,1)
  ncell_size=5*ncell

  ! Exception handling
  if (itime_begin < 1 .or. itime_begin > ubound(fds_time,1)) then
    return
  end if

  if (itime_end < 1 .or. itime_end > ubound(fds_time,1)) then
    return
  end if

  if (itime_begin > itime_end) then
    return
  end if

  if (average) then

    write(f1,'(i4.4)') itime_begin; write(f2,'(i4.4)') itime_end
    filename=trim(fds_chid) // '_bndf_data_' // trim(f1) // "_to_" // trim(f2) // '_average.vtk'
    open(unit=iochannel(1),file=filename,status='replace',iostat=ios)
    if (ios /= 0) call error_open_file(filename)

    write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
    write(iochannel(1),'(a)') 'fds2fem: FDS BNDF dump'
    write(iochannel(1),'(a)') 'ASCII'
    write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

    write(iochannel(1),'(a,a,a)') 'POINTS ', trim(int2str(nnodes)), ' double'
    do i=1,nnodes
      write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (fds_bndf_xyz(i,j),j=1,3)
    end do

    write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', trim(int2str(ncell_size))
    do i=1,ncell
      write(iochannel(1),'(6(a))') trim(int2str(4)), ' ', (trim(int2str(fds_bndf_el(i,j)-1)) // ' ',j=1,4)
    end do

    write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
    do i=1,ncell
      write(iochannel(1),'(a)') '9'
    end do

    write(iochannel(1),'(2(a))') 'CELL_DATA ', trim(int2str(ncell))
    write(iochannel(1),'(a)') 'SCALARS anonymous double 1'
    write(iochannel(1),'(a)') 'LOOKUP_TABLE default'

    ntimes=itime_end-itime_begin+1
    do j=1,ncell
      dsum=0.0
      do itime=itime_begin,itime_end
        dsum=dsum+fds_data(itime,fds_bndf_el(j,1))
        dsum=dsum+fds_data(itime,fds_bndf_el(j,2))
        dsum=dsum+fds_data(itime,fds_bndf_el(j,3))
        dsum=dsum+fds_data(itime,fds_bndf_el(j,4))
      end do
      dsum=0.25*dsum/ntimes
      write(iochannel(1),'(es15.7e3)') dsum
    end do

    close(unit=iochannel(1))

  else

  ! Write loop
  time_loop: do itime=itime_begin,itime_end

    write(ctime,'(i4.4)') itime
    filename=trim(fds_chid) // '_bndf_data_' // trim(ctime) // '.vtk'
    open(unit=iochannel(1),file=filename,status='replace',iostat=ios)
    if (ios /= 0) call error_open_file(filename)

    write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
    write(iochannel(1),'(a)') 'fds2fem: FDS BNDF dump'
    write(iochannel(1),'(a)') 'ASCII'
    write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

    write(iochannel(1),'(a,a,a)') 'POINTS ', trim(int2str(nnodes)), ' double'
    do i=1,nnodes
      write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (fds_bndf_xyz(i,j),j=1,3)
    end do

    write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', trim(int2str(ncell_size))
    do i=1,ncell
      write(iochannel(1),'(6(a))') trim(int2str(4)), ' ', (trim(int2str(fds_bndf_el(i,j)-1)) // ' ',j=1,4)
    end do

    write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
    do i=1,ncell
      write(iochannel(1),'(a)') '9'
    end do

    write(iochannel(1),'(2(a))') 'CELL_DATA ', trim(int2str(ncell))
    write(iochannel(1),'(a)') 'SCALARS anonymous double 1'
    write(iochannel(1),'(a)') 'LOOKUP_TABLE default'

    do j=1,ncell
      dsum=0.25*fds_data(itime,fds_bndf_el(j,1))
      dsum=dsum+0.25*fds_data(itime,fds_bndf_el(j,2))
      dsum=dsum+0.25*fds_data(itime,fds_bndf_el(j,3))
      dsum=dsum+0.25*fds_data(itime,fds_bndf_el(j,4))
      write(iochannel(1),'(es15.7e3)') dsum
    end do

    close(unit=iochannel(1))

  end do time_loop

  end if

  ! Status report
  ntimes=(itime_end-itime_begin)+1
  f1=trim(fds_chid) // '_bndf_data_' // trim(int2str(itime_begin)) // '.vtk'
  f2=trim(fds_chid) // '_bndf_data_' // trim(int2str(itime_end)) // '.vtk'

  write(*,'(t3,5(a))') trim(int2str(ntimes)), ' time steps written in files ', &
    trim(quote(f1)), '...', trim(quote(f2)) 

end subroutine dump_fds_element_data_vtk

subroutine dump_fds_model_vtk()
!----------------------------------------
! Dump FDS node coordinates in VTK format
!----------------------------------------
  use error_messages
  use fds_head_arrays
  use fds_obst_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,j,ios,nnodes,ncell,ncell_size
  character(len=chr80) :: filename

  nnodes=ubound(fds_obst_nd,1)
  ncell=ubound(fds_obst_el,1)
  ncell_size=9*ncell

  filename=trim(fds_chid) // '_model.vtk'
  open(unit=iochannel(1),file=filename,status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
  write(iochannel(1),'(a)') 'fds2fem: FDS model dump'
  write(iochannel(1),'(a)') 'ASCII'
  write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

  write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(nnodes)), ' double'
  do i=1,nnodes
    write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (fds_obst_nd(i,j),j=1,3)
  end do

  write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', trim(int2str(ncell_size))
  do i=1,ncell
    write(iochannel(1),'(10(a))') trim(int2str(8)), ' ', (trim(int2str(fds_obst_el(i,j)-1)) // ' ',j=1,8)
  end do

  write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
  do i=1,ncell
    write(iochannel(1),'(a)') '12'
  end do

  close(unit=iochannel(1))

  write(*,'(t3,5(a))') trim(int2str(ncell)), &
    ' obstructions written in file ', &
    trim(quote(filename))

end subroutine dump_fds_model_vtk

subroutine construct_fds_model()
!-----------------------------------------------------------
! Construct an element-based representation of the FDS model
!-----------------------------------------------------------
  use error_messages
  use fds_obst_arrays
  use global_constants
  use global_variables
  implicit none

  integer :: i,ios,n_obst,n_node

  n_obst=ubound(fds_obst_xb,1)
  n_node=8*n_obst

  allocate(fds_obst_nd(n_node,3),stat=ios); call error_allocate(ios)
  allocate(fds_obst_el(n_obst,8),stat=ios); call error_allocate(ios)
  
  do i=1,n_obst
    fds_obst_nd((i-1)*8+1,1)=fds_obst_xb(i,1)
    fds_obst_nd((i-1)*8+1,2)=fds_obst_xb(i,3)
    fds_obst_nd((i-1)*8+1,3)=fds_obst_xb(i,5)
    
    fds_obst_nd((i-1)*8+2,1)=fds_obst_xb(i,2)
    fds_obst_nd((i-1)*8+2,2)=fds_obst_xb(i,3)
    fds_obst_nd((i-1)*8+2,3)=fds_obst_xb(i,5)
    
    fds_obst_nd((i-1)*8+3,1)=fds_obst_xb(i,2)
    fds_obst_nd((i-1)*8+3,2)=fds_obst_xb(i,4)
    fds_obst_nd((i-1)*8+3,3)=fds_obst_xb(i,5)
    
    fds_obst_nd((i-1)*8+4,1)=fds_obst_xb(i,1)
    fds_obst_nd((i-1)*8+4,2)=fds_obst_xb(i,4)
    fds_obst_nd((i-1)*8+4,3)=fds_obst_xb(i,5)
    
    fds_obst_nd((i-1)*8+5,1)=fds_obst_xb(i,1)
    fds_obst_nd((i-1)*8+5,2)=fds_obst_xb(i,3)
    fds_obst_nd((i-1)*8+5,3)=fds_obst_xb(i,6)
    
    fds_obst_nd((i-1)*8+6,1)=fds_obst_xb(i,2)
    fds_obst_nd((i-1)*8+6,2)=fds_obst_xb(i,3)
    fds_obst_nd((i-1)*8+6,3)=fds_obst_xb(i,6)
    
    fds_obst_nd((i-1)*8+7,1)=fds_obst_xb(i,2)
    fds_obst_nd((i-1)*8+7,2)=fds_obst_xb(i,4)
    fds_obst_nd((i-1)*8+7,3)=fds_obst_xb(i,6)
    
    fds_obst_nd((i-1)*8+8,1)=fds_obst_xb(i,1)
    fds_obst_nd((i-1)*8+8,2)=fds_obst_xb(i,4)
    fds_obst_nd((i-1)*8+8,3)=fds_obst_xb(i,6)

    fds_obst_el(i,1)=(i-1)*8+1
    fds_obst_el(i,2)=(i-1)*8+2
    fds_obst_el(i,3)=(i-1)*8+3
    fds_obst_el(i,4)=(i-1)*8+4
    fds_obst_el(i,5)=(i-1)*8+5
    fds_obst_el(i,6)=(i-1)*8+6
    fds_obst_el(i,7)=(i-1)*8+7
    fds_obst_el(i,8)=(i-1)*8+8
  end do

end subroutine construct_fds_model

end module fds_dump
