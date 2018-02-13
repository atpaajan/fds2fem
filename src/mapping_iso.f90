!-------------------------------------------
! Functions and subroutines for mesh mapping
! Time-temperature curve inputs
!-------------------------------------------
module iso_mapping

contains

subroutine mapping_iso()
!---------
! Main sub
!---------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: ios

  write(*,'(a)') ''
  write(*,'(a)') 'ISO Mesh mapping module'

  !--------------------------------
  ! Sometimes there's nothing to do
  !--------------------------------

  if (trim(mapping_method) == 'off') then
    write(*,'(t3,a)') 'Nothing to be done'
    fem_data_available=.false.
    return
  end if

  if (.not. fem_xyz_available) then
    write(*,'(t3,a)') 'WARNING: FEM node coordinates unavailable'
    fem_data_available=.false.
    return
  end if

  !--------
  ! Mapping
  !--------

  ntimes_fem=ubound(fds_time,1)
  nnodes_fem=ubound(fem_xyz,1)

  allocate(fem_time(ntimes_fem),stat=ios);            call error_allocate(ios)
  allocate(fem_data(ntimes_fem,nnodes_fem),stat=ios); call error_allocate(ios)
  if (read_hcoeff) then
    allocate(fem_hcoeff(ntimes_fem,nnodes_fem),stat=ios); call error_allocate(ios)
  end if

  fem_time=fds_time; fem_data=0.0

  allocate(fds_mask(nnodes_fds),stat=ios); call error_allocate(ios)
  allocate(fem_mask(nnodes_fem),stat=ios); call error_allocate(ios)

  allocate(abaqus_node_emissivity(nnodes_fem),stat=ios); call error_allocate(ios)
  allocate(abaqus_node_hcoeff(nnodes_fem),stat=ios); call error_allocate(ios)
  abaqus_node_emissivity= emissivity
  abaqus_node_hcoeff=     hcoeff

  !------------------------------
  ! Direct NSET-DEVC connectivity
  !------------------------------
  
  call map_iso_to_nset()
  
  fem_data_available=.true.
  
end subroutine mapping_iso

!***************
! Auxiliary subs
!***************

subroutine map_iso_to_nset()
!---------------------------------------
! Time-temp curve(s) to node set mapping
!---------------------------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none
 
  integer :: i, j, k, k_ttcurve, ios, ishift_begin
  character(len=chr80) :: nset_name
  Real(kind=rk) :: dt, x_tmp
  Real(kind=rk), Dimension(:), Allocatable :: fds_data_tmp

  Write(*,'(t3,3(a))') 'Direct NSET - time-temp curve mapping'

  if (.not. nset_connectivity) then
    write(*,'(a)') 'ERROR: no NSET - time-temp curve connectivity given'
    stop
  end if

  !---------------------------
  ! Assign values to node sets
  !---------------------------

  fem_data = 0.0
  
  ! nnodes_fds: the number of lines in the nset file for time-temp input
  Allocate(fds_data_tmp(ntimes_fds),stat=ios); call error_allocate(ios)

  i_loop: Do i = 1, Ubound(connectivity_table,1)
     nset_name = connectivity_table(i,1)
     k_loop: Do k = 1, nnodes_fds
        If (fds_idevc(k) == connectivity_table_num(i,2)) Then
           k_ttcurve = k
           Exit k_loop
        End If
     End Do k_loop

     If (time_shift_mask(i)) Then
        ! Do the time shift operation/interpolation
        ! The time shift must be positive (Add this check to the input reading part)
        dt           = (iso_tend-iso_tbegin)/Max(iso_ntimes-1,1)
        ishift_begin = Max(0,Int(connectivity_phys_cons(i,3)/dt)) + 1
        fds_data_tmp = fds_data(1,k_ttcurve)
        x_tmp        = (Max(0.0,connectivity_phys_cons(i,3)) - fds_time(ishift_begin))/dt
        j_loop_1: Do j = ishift_begin + 1, ntimes_fds
           fds_data_tmp(j) = fds_data(j-ishift_begin,k_ttcurve) + (1.0-x_tmp)* &
                (fds_data(j-ishift_begin+1,k_ttcurve)-fds_data(j-ishift_begin,k_ttcurve))
        End Do j_loop_1
     Else
        fds_data_tmp = fds_data(:,k_ttcurve)
     End If
     j_loop_2: Do j = 1, nnodes_fem
        If (Trim(fem_nset(j)) == Trim(nset_name)) Then
           abaqus_node_emissivity(j) = connectivity_phys_cons(i,1)
           abaqus_node_hcoeff(j)     = connectivity_phys_cons(i,2)
           fem_data(:,j) = fds_data_tmp(:)
        End If
     End Do j_loop_2
        
  End Do i_loop

  Deallocate(fds_data_tmp)
  
end subroutine map_iso_to_nset

end module iso_mapping
