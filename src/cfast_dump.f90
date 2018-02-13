!-----------------------------------------------
! Functions and subroutines for CFAST data dumps
!-----------------------------------------------
module cfast_dump

  ! Dumps just the targets that are defined in the nset-file

contains

subroutine cfast_dump_module()
!---------
! Main sub
!---------
  use global_constants
  use global_variables
  use mapping_arrays
  implicit none

  logical :: something_done

  write(*,'(a)') ''
  write(*,'(a)') 'CFAST dump module'

  if (trim(dump_cfast_nodes) == 'off' .and. &
      trim(dump_cfast_data)  == 'off') then
    write(*,'(t3,a)') 'No dump requested'
    return
  end if

  !---------------------
  ! Node coordinate dump
  !---------------------

  something_done=.false.
  if (cfast_xyz_available) then
    select case (trim(dump_cfast_nodes))
    case ('xyz')
      ! XYZ node dump
      call dump_cfast_devc_nodes_xyz() 
    case ('vtk')
      ! VTK node dump
      call dump_cfast_devc_nodes_vtk() 
    end select
    something_done=.true.
  end if

  !----------
  ! Data dump
  !----------

  ! Arrays, where just the targets in the nset file:
  ! nnodes_fds                       Number of targets that are actually used in the nset file
  ! fds_id(nnodes_fds)               Target name
  ! fds_time(ntimes_fds)             time axis
  ! fds_idevc(nnodes_fds)            Target running number (in input/_w.csv file)
  ! fds_data(ntimes_fds,nnodes_fds)  Temperature/flux data
  ! fds_ior(nnodes_fds)              +-1,+-2,+-3 (approx. own interpolation)
  
  if (cfast_data_available) then
    select case (trim(dump_cfast_data))
    case ('sdf')
      ! SDF data dump
      call dump_cfast_devc_data_sdf()
    case ('vtk')
      ! VTK data dump
      call dump_cfast_devc_data_vtk(1,ubound(fds_time,1))
    end select
    something_done=.true.
  end if

  if (.not. something_done) then
    write(*,'(t3,a)') 'Nothing to be done'
  end if

end subroutine cfast_dump_module

!***************
! Auxiliary subs
!***************

subroutine dump_cfast_devc_nodes_xyz()
!--------------------------------------------
! Dump CFAST target coordinates in XYZ format
!--------------------------------------------
  use error_messages
  use fds_devc_arrays
  use fds_head_arrays
  use global_constants
  use global_variables
  use string_handling
  Use mapping_arrays
  implicit none

  integer :: j,ios,idevc,ntrg
  character(len=chr80) :: filename

  filename=trim(fds_input_file) // "_trg.xyz"
  open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  ntrg=ubound(fds_idevc,1)

  do idevc=1,ntrg
    write(iochannel(1),'(3(es15.7e3,1x),i2,1x,3(a))') (fds_xyz(idevc,j),j=1,3), &
      fds_ior(idevc), char(34), trim(fds_id(idevc)), char(34)
  end do

  close(unit=iochannel(1))

  ! Status report
  write(*,'(t3,5(a))') trim(int2str(ntrg)), &
    ' CFAST device coordinates written in file ', &
    trim(quote(filename))

end subroutine dump_cfast_devc_nodes_xyz

subroutine dump_cfast_devc_data_sdf
!-------------------------------------
! Dump CFAST target data in SDF format
!-------------------------------------
  use error_messages
  use fds_devc_arrays
  use fds_head_arrays
  use global_constants
  use global_variables
  use string_handling
  use mapping_arrays
  implicit none

  integer :: i,j,ios,ntrg,ntimes
  character(len=chr80) :: filename

  ntimes=ubound(fds_time,1)
  ntrg=ubound(fds_idevc,1)

  filename=trim(fds_input_file) // "_trg.sdf"
  open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  do i=1,ntimes
    write(iochannel(1),'(es15.7e3)',advance='no') fds_time(i)
    do j=1,ntrg
      write(iochannel(1),'(1x,es15.7e3)',advance='no') fds_data(i,j)
    end do
    write(iochannel(1),'(a)') ''
  end do

  close(unit=iochannel(1))

  ! Status report
  write(*,'(t3,5(a))') trim(int2str(ntimes)), &
    ' time steps of CFAST device data written in file ', &
    trim(quote(filename))

end subroutine dump_cfast_devc_data_sdf

subroutine dump_cfast_devc_nodes_vtk()
!--------------------------------------------
! Dump CFAST target coordinates in VTK format
!--------------------------------------------
  use error_messages
  use fds_devc_arrays
  use fds_head_arrays
  use global_constants
  use global_variables
  use string_handling
  use mapping_arrays
  implicit none

  integer :: i,j,ios,ntrg,ncell,ncell_size
  character(len=chr80) :: filename

  ntrg=ubound(fds_idevc,1)
  ncell=ntrg; ncell_size=2*ncell

  filename=trim(fds_input_file) // '_trg.vtk'
  open(unit=iochannel(1),file=filename,status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
  write(iochannel(1),'(a)') 'fds2fem: CFAST target coordinate dump'
  write(iochannel(1),'(a)') 'ASCII'
  write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

  write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(ncell)), ' double'
  do i=1,ntrg
     !write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (fds_devc_xyz(i,j),j=1,3) 
    write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (fds_xyz(i,j),j=1,3) 
  end do

  write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', trim(int2str(ncell_size))
  do i=1,ntrg
    write(iochannel(1),'(3(a))') trim(int2str(1)), ' ', trim(int2str(i-1))
  end do

  write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
  do i=1,ntrg
    write(iochannel(1),'(a)') '1'
  end do

  close(unit=iochannel(1))

  write(*,'(t3,a,a,a,a,a)') trim(int2str(ntrg)), &
    ' CFAST device coordinates written in file ', &
    trim(quote(filename))

end subroutine dump_cfast_devc_nodes_vtk

subroutine dump_cfast_devc_data_vtk(itime_begin,itime_end)
!-------------------------------------
! Dump CFAST target data in VTK format
!-------------------------------------
  use error_messages
  use fds_devc_arrays
  use fds_head_arrays
  use global_constants
  use global_variables
  use string_handling
  use mapping_arrays
  implicit none

  integer :: i,j,ios,ntrg,ncell,ncell_size,itime
  integer :: itime_begin,itime_end,ntimes
  character(len=chr80) :: filename,ctime,f1,f2

  ntrg=ubound(fds_idevc,1)
  ncell=ntrg; ncell_size=2*ncell

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
    filename=trim(fds_input_file) // '_trg_' // trim(ctime) // '.vtk'
    open(unit=iochannel(1),file=filename,status='replace',iostat=ios)
    if (ios /= 0) call error_open_file(filename)

    write(iochannel(1),'(a)') '# vtk DataFile Version 3.0' 
    write(iochannel(1),'(a)') 'fds2fem: CFAST target data dump'
    write(iochannel(1),'(a)') 'ASCII'
    write(iochannel(1),'(a)') 'DATASET UNSTRUCTURED_GRID'

    write(iochannel(1),'(3(a))') 'POINTS ', trim(int2str(ncell)), ' double'
    do i=1,ntrg
      write(iochannel(1),'(2(es15.7e3,1x),es15.7e3)') (fds_xyz(i,j),j=1,3) 
    end do

    write(iochannel(1),'(4(a))') 'CELLS ', trim(int2str(ncell)), ' ', trim(int2str(ncell_size))
    do i=1,ntrg
      write(iochannel(1),'(3(a))') trim(int2str(1)), ' ', trim(int2str(i-1))
    end do

    write(iochannel(1),'(2(a))') 'CELL_TYPES ', trim(int2str(ncell))
    do i=1,ntrg
      write(iochannel(1),'(a)') '1'
    end do

    write(iochannel(1),'(2(a))') 'CELL_DATA ', trim(int2str(ncell))
    write(iochannel(1),'(a)') 'SCALARS anonymous double 1'
    write(iochannel(1),'(a)') 'LOOKUP_TABLE default'

    do i=1,ntrg
      write(iochannel(1),'(es15.7e3)') fds_data(itime,i)
    end do

    close(unit=iochannel(1))

  end do time_loop

  ! Status report
  ntimes=(itime_end-itime_begin)+1
  f1=trim(fds_input_file) // '_trg_' // trim(int2str(itime_begin)) // '.vtk'
  f2=trim(fds_input_file) // '_trg_' // trim(int2str(itime_end)) // '.vtk'

  write(*,'(t3,5(a))') trim(int2str(ntimes)), ' time steps of CFAST target data written in files ', &
    trim(quote(f1)), '...', trim(quote(f2)) 

end subroutine dump_cfast_devc_data_vtk

end module cfast_dump
