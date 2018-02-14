!----------------------------------------
! Miscellaneous functions and subroutines
!----------------------------------------
module miscellaneous

contains

subroutine allocate_channels
!----------------------
! Allocate I/O channels
!----------------------
  use global_constants
  implicit none

  integer :: i,ichannel,nchannel
  logical :: reserved,completed

  nchannel=ubound(iochannel,1)

  ! ubound(iochannel,1) >= 10 expected
  if (nchannel < 10) then
    write(*,'(a)') 'ERROR: in subroutine allocate_channels(): ubound(iochannel,1) < 10'
    stop
  end if

  ichannel=1; completed=.false.
  iochannel_loop: do i=10,127
    inquire(unit=i,opened=reserved)
    if (.not. reserved) then
      if (ichannel == nchannel+1) then
        scratch_channel=i
        completed=.true.
        exit iochannel_loop
      end if
      iochannel(ichannel)=i
      ichannel=ichannel+1
    end if
  end do iochannel_loop
 
  if (.not. completed) then
    write(*,'(a)') 'ERROR: in subroutine allocate_channels(): in allocating I/O channels'
    stop
  end if

end subroutine allocate_channels

subroutine deallocate_excess_fds_arrays()
!-----------------
! Free some memory
!-----------------
  use error_messages
  use fds_bndf_arrays
  use fds_devc_arrays
  use fds_mesh_arrays
  use fds_obst_arrays
  use fds_prop_arrays
  use global_variables
  implicit none

  integer :: ios

  !--------------------
  ! BNDF-related arrays
  !--------------------

  if (allocated(fds_bndf_ior)) then
    deallocate(fds_bndf_ior,stat=ios); call error_allocate(ios)
  end if

  if (allocated(fds_bndf_nnod)) then
    deallocate(fds_bndf_nnod,stat=ios); call error_allocate(ios)
  end if

  if (allocated(fds_bndf_npat)) then
    deallocate(fds_bndf_npat,stat=ios); call error_allocate(ios)
  end if
  
  if (allocated(fds_bndf_ntim)) then
    deallocate(fds_bndf_ntim,stat=ios); call error_allocate(ios)
  end if
  
  if (allocated(fds_bndf_pat)) then
    deallocate(fds_bndf_pat,stat=ios); call error_allocate(ios)
  end if

  if (allocated(fds_bndf_nb)) then
    deallocate(fds_bndf_nb,stat=ios); call error_allocate(ios)
  end if

  if (allocated(fds_bndf_nm)) then
    deallocate(fds_bndf_nm,stat=ios); call error_allocate(ios)
  end if
  
  if (allocated(fds_bndf_np)) then
    deallocate(fds_bndf_np,stat=ios); call error_allocate(ios)
  end if

  if (allocated(fds_bndf_qnty)) then
    deallocate(fds_bndf_qnty,stat=ios); call error_allocate(ios)
  end if
  
  if (allocated(fds_bndf_snam)) then
    deallocate(fds_bndf_snam,stat=ios); call error_allocate(ios)
  end if
  
  if (allocated(fds_bndf_unit)) then
    deallocate(fds_bndf_unit,stat=ios); call error_allocate(ios)
  end if

  if (allocated(fds_bndf_file)) then
    deallocate(fds_bndf_file,stat=ios); call error_allocate(ios)
  end if

  if (allocated(fds_bndf_time)) then
    deallocate(fds_bndf_time,stat=ios); call error_allocate(ios)
  end if
  
  if (trim(dump_fds_model) == 'off' .And. Trim(dump_fds_data) == 'off') then
    if (allocated(fds_bndf_xyz)) then
      deallocate(fds_bndf_xyz,stat=ios); call error_allocate(ios)
    end if
  end if
  
  if (allocated(fds_bndf_data)) then
    deallocate(fds_bndf_data,stat=ios); call error_allocate(ios)
  end if

  !--------------------
  ! DEVC-related arrays
  !--------------------

  !if (allocated(fds_devc_ior)) then
  !  deallocate(fds_devc_ior,stat=ios); call error_allocate(ios)
  !end if

  if (allocated(fds_devc_rows)) then
    deallocate(fds_devc_rows,stat=ios); call error_allocate(ios)
  end if

  if (allocated(fds_devc_cols)) then
    deallocate(fds_devc_cols,stat=ios); call error_allocate(ios)
  end if

  if (allocated(fds_devc_file)) then
    deallocate(fds_devc_file,stat=ios); call error_allocate(ios)
  end if

  !if (allocated(fds_devc_name)) then
  !  deallocate(fds_devc_name,stat=ios); call error_allocate(ios)
  !end if

  if (allocated(fds_devc_name_b)) then
    deallocate(fds_devc_name_b,stat=ios); call error_allocate(ios)
  end if

  if (allocated(fds_devc_qnty)) then
    deallocate(fds_devc_qnty,stat=ios); call error_allocate(ios)
  end if

  if (allocated(fds_devc_unit)) then
    deallocate(fds_devc_unit,stat=ios); call error_allocate(ios)
  end if

  !if (allocated(fds_devc_time)) then
  !  deallocate(fds_devc_time,stat=ios); call error_allocate(ios)
  !end if

  if (allocated(fds_devc_xyz)) then
    deallocate(fds_devc_xyz,stat=ios); call error_allocate(ios)
  end if

  if (allocated(fds_devc_xb)) then
    deallocate(fds_devc_xb,stat=ios); call error_allocate(ios)
  end if

  !if (allocated(fds_devc_data)) then
  !  deallocate(fds_devc_data,stat=ios); call error_allocate(ios)
  !end if

  !--------------------
  ! MESH-related arrays
  !--------------------

  if (allocated(fds_mesh_ijk)) then
    deallocate(fds_mesh_ijk,stat=ios); call error_allocate(ios)
  end if
  
  if (allocated(fds_mesh_xb)) then
    deallocate(fds_mesh_xb,stat=ios); call error_allocate(ios)
  end if

  !--------------------
  ! OBST-related arrays
  !--------------------

  if (trim(dump_fds_model) == 'off') then
    if (allocated(fds_obst_xb)) then
      deallocate(fds_obst_xb,stat=ios); call error_allocate(ios)
    end if

    if (allocated(fds_obst_nd)) then
      deallocate(fds_obst_nd,stat=ios); call error_allocate(ios)
    end if

    if (allocated(fds_obst_el)) then
      deallocate(fds_obst_el,stat=ios); call error_allocate(ios)
    end if
  end if

  !--------------------
  ! PROP-related arrays
  !--------------------

  if (allocated(fds_prop_name)) then
    deallocate(fds_prop_name,stat=ios); call error_allocate(ios)
  end if
  
  if (allocated(fds_prop_qnty)) then
    deallocate(fds_prop_qnty,stat=ios); call error_allocate(ios)
  end if

end subroutine deallocate_excess_fds_arrays

subroutine deallocate_excess_abaqus_arrays()
!-----------------
! Free some memory
!-----------------
  use abaqus_arrays
  use error_messages
  use global_variables
  implicit none

  integer :: ios
  logical, save :: first_call = .true.

  if (first_call) then

    !-------------
    ! Part-related
    !-------------

    if (allocated(part_node_xyz)) then
      deallocate(part_node_xyz,stat=ios); call error_allocate(ios)
    end if

    if (allocated(part_node_number)) then
      deallocate(part_node_number,stat=ios); call error_allocate(ios)
    end if
    
    if (allocated(part_node_name)) then
      deallocate(part_node_name,stat=ios); call error_allocate(ios)
    end if

    if (allocated(part_node_part)) then
      deallocate(part_node_part,stat=ios); call error_allocate(ios)
    end if
    
    if (allocated(part_nset_name)) then
      deallocate(part_nset_name,stat=ios); call error_allocate(ios)
    end if

    if (allocated(part_nset_node)) then
      deallocate(part_nset_node,stat=ios); call error_allocate(ios)
    end if
    
    if (allocated(part_nset_part)) then
      deallocate(part_nset_part,stat=ios); call error_allocate(ios)
    end if

    if (allocated(part_element_number)) then
      deallocate(part_element_number,stat=ios); call error_allocate(ios)
    end if

    if (allocated(part_element_node)) then
      deallocate(part_element_node,stat=ios); call error_allocate(ios)
    end if

    if (allocated(part_element_part)) then
      deallocate(part_element_part,stat=ios); call error_allocate(ios)
    end if

    !-----------------
    ! Instance-related
    !-----------------

    if (allocated(instance_node_xyz)) then
      deallocate(instance_node_xyz,stat=ios); call error_allocate(ios)
    end if

    if (allocated(instance_node_number)) then
      deallocate(instance_node_number,stat=ios); call error_allocate(ios)
    end if
    
    if (allocated(instance_node_name)) then
      deallocate(instance_node_name,stat=ios); call error_allocate(ios)
    end if

    if (allocated(instance_node_instance)) then
      deallocate(instance_node_instance,stat=ios); call error_allocate(ios)
    end if
    
    if (allocated(instance_nset_name)) then
      deallocate(instance_nset_name,stat=ios); call error_allocate(ios)
    end if

    if (allocated(instance_nset_node)) then
      deallocate(instance_nset_node,stat=ios); call error_allocate(ios)
    end if
    
    if (allocated(instance_nset_instance)) then
      deallocate(instance_nset_part,stat=ios); call error_allocate(ios)
    end if

    if (allocated(instance_element_number)) then
      deallocate(instance_element_number,stat=ios); call error_allocate(ios)
    end if

    if (allocated(instance_element_node)) then
      deallocate(instance_element_node,stat=ios); call error_allocate(ios)
    end if

    if (allocated(instance_element_part)) then
      deallocate(instance_element_part,stat=ios); call error_allocate(ios)
    end if
    
    if (allocated(instance_element_instance)) then
      deallocate(instance_element_instance,stat=ios); call error_allocate(ios)
    end if

    first_call=.false.

  else

    !--------------
    ! Model-related
    !--------------

    !if (trim(dump_abaqus_model) == 'off') then
    !  if (allocated(model_node_xyz)) then
    !    deallocate(model_node_xyz,stat=ios); call error_allocate(ios)
    !  end if
    !end if

    if (allocated(model_node_number)) then
      deallocate(model_node_number,stat=ios); call error_allocate(ios)
    end if
    
    if (allocated(model_node_name)) then
      deallocate(model_node_name,stat=ios); call error_allocate(ios)
    end if

   if (allocated(model_node_part)) then
      deallocate(model_node_part,stat=ios); call error_allocate(ios)
    end if
    
    if (allocated(model_node_instance)) then
      deallocate(model_node_instance,stat=ios); call error_allocate(ios)
    end if

    if (allocated(model_node_type)) then
      deallocate(model_node_type,stat=ios); call error_allocate(ios)
    end if
    
    if (allocated(model_node_elements)) then
      deallocate(model_node_elements,stat=ios); call error_allocate(ios)
    end if

    if (allocated(model_node_area)) then
      deallocate(model_node_area,stat=ios); call error_allocate(ios)
    end if
    
    if (allocated(model_nset_name)) then
      deallocate(model_nset_name,stat=ios); call error_allocate(ios)
    end if

    if (allocated(model_nset_part)) then
      deallocate(model_nset_part,stat=ios); call error_allocate(ios)
    end if

    if (allocated(model_nset_instance)) then
      deallocate(model_nset_instance,stat=ios); call error_allocate(ios)
    end if

    if (allocated(model_nset_node)) then
      deallocate(model_nset_node,stat=ios); call error_allocate(ios)
    end if
   
    if (allocated(model_element_part)) then
      deallocate(model_element_part,stat=ios); call error_allocate(ios)
    end if

    if (allocated(model_element_instance)) then
      deallocate(model_element_instance,stat=ios); call error_allocate(ios)
    end if

    if (allocated(model_element_number)) then
      deallocate(model_element_number,stat=ios); call error_allocate(ios)
    end if 

    !if (trim(dump_abaqus_model) == 'off') then
    !  if (allocated(model_element_type)) then
    !    deallocate(model_element_type,stat=ios); call error_allocate(ios)
    !  end if 
    !end if

    if (allocated(model_element_node)) then
      deallocate(model_element_node,stat=ios); call error_allocate(ios)
    end if

    !if (trim(dump_abaqus_model) == 'off') then
    !  if (allocated(model_element_int_node)) then
    !    deallocate(model_element_int_node,stat=ios); call error_allocate(ios)
    !  end if 
    !end if
    
  end if

end subroutine deallocate_excess_abaqus_arrays

subroutine deallocate_excess_ansys_arrays()
!-----------------
! Free some memory
!-----------------
  use error_messages
  use global_variables
  implicit none

  ! Under development

end subroutine deallocate_excess_ansys_arrays

function sec2str(t) result(str)
!-------------------------------
! Return elapsed time in seconds
!-------------------------------
  use global_constants
  implicit none

  real(kind=rk) :: t
  character(len=chr80) :: str

  str=''
  if ((t) < 1.0E1) then 
    ! [0.0 s, 10.0 s)
    write(str,'(f6.3)') (t)
  else if ((t) < 1.0E2) then
    ! [10.0 s, 100.0 s)
    write(str,'(f7.3)') (t)
  else if ((t) < 1.0E3) then
    ! [100.0 s, 1000.0 s)
    write(str,'(f8.3)') (t)
  else if ((t) < 1.0E4) then
    ! [1000.0 s, 10000.0 s)
    write(str,'(f9.3)') (t)
  else if ((t) < 1.0E5) then
    ! [10000.0 s, 100000.0 s)
    write(str,'(f10.3)') (t)
  else
    ! [100000.0 s, Inf)
    write(str,'(es15.7e3)') (t)
  end if

end function sec2str

subroutine progress_bar(x,tab,width)
!-----------------------
! Display a progress bar
!-----------------------
  use global_constants
  use string_handling
  implicit none

  integer :: i,tab,width,cursor,percentage
  integer :: narrow,ndash,nspace
  
  character(len=chr80) :: fmt_1

  real(kind=rk) :: x

  cursor=0
  cursor=nint(x*(width-2))
  percentage=nint(100.0*x)

  ndash=cursor-1
  narrow=cursor
  nspace=width-cursor-2

  write(*,'(2(a))') char(27),'[1A'

  fmt_1='(t' // trim(int2str(tab)) // ',a,$)'
  write(*,fmt_1) '* ['

  if (ndash >= 1) then
  fmt_1='(' // trim(int2str(ndash)) // 'a,$' //')'
  write(*,fmt_1) ('-',i=1,ndash)  
  end if

  if (narrow > 0) then
  write(*,'(a,$)') '>'
  end if

  if (nspace >= 1) then
  fmt_1='(' // trim(int2str(nspace)) // 'a,$' //')'
  write(*,fmt_1) (' ',i=1,nspace)
  end if

  if (percentage < 100) then
    write(*,'(3(a),$)') '] ', trim(int2str(percentage)), '%'
    write(*,'(2(a))') char(27),'[1A'
  else
    write(*,'(3(a))') '] ', trim(int2str(percentage)), '%'
  end if
  
end subroutine progress_bar

subroutine progress_line(x)
!------------------------
! Display a progress line
!------------------------
  use global_constants
  use string_handling
  implicit none

  integer :: percentage
  real(kind=rk) :: x

  percentage=nint(100.0*x)
 
  if (percentage == 0) then
    write(*,'(t6,2(a),$)')  trim(int2str(percentage)), '%'
  else if (percentage < 100) then
    write(*,'(3(a),$)') ', ', trim(int2str(percentage)), '%'
  else
    write(*,'(3(a))') ', ', trim(int2str(percentage)), '%'
  end if 

end subroutine progress_line

end module miscellaneous
