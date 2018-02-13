!---------------------------------------------------------
! Functions and subroutines for calculating FDS statistics
!---------------------------------------------------------
module fds_stats

contains

subroutine fds_stats_module()
!---------
! Main sub
!---------
  use global_constants
  use global_variables
  implicit none

  write(*,'(a)') ''
  write(*,'(a)') 'FDS statistics module'

  if (fds_statistics) then
    if (fds_xyz_available) then
      call analyze_fds_nodes()
    else
      write(*,'(t3,a)') 'Nothing to be done'
    end if
  else
    write(*,'(t3,a)') 'Nothing to be done'
  end if

end subroutine fds_stats_module

!***************
! Auxiliary subs
!***************

subroutine analyze_fds_nodes()
!----------------------
! FDS node set analysis
!----------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,imindist
  character(len=chr80) :: fmt_1
  real(kind=rk) :: dx,dy,dz,r_sq,r_ave,r_ave_sum,mindist,fds_average

  logical, dimension(:), allocatable :: dmask
  real(kind=rk), dimension(:), allocatable :: distance

  !-------------------
  ! CM of bounding box
  !-------------------

  if (.not. fds_cm_calculated) then
    fds_node_range(1)=minval(fds_xyz(1:nnodes_fds,1))
    fds_node_range(2)=maxval(fds_xyz(1:nnodes_fds,1))
    fds_node_range(3)=minval(fds_xyz(1:nnodes_fds,2)) 
    fds_node_range(4)=maxval(fds_xyz(1:nnodes_fds,2)) 
    fds_node_range(5)=minval(fds_xyz(1:nnodes_fds,3)) 
    fds_node_range(6)=maxval(fds_xyz(1:nnodes_fds,3)) 

    fds_dr(1)=fds_node_range(2)-fds_node_range(1)
    fds_dr(2)=fds_node_range(4)-fds_node_range(3)
    fds_dr(3)=fds_node_range(6)-fds_node_range(5)

    cm_fds_bb(1)=0.5*(fds_node_range(2)+fds_node_range(1))
    cm_fds_bb(2)=0.5*(fds_node_range(4)+fds_node_range(3))
    cm_fds_bb(3)=0.5*(fds_node_range(6)+fds_node_range(5))
  end if

  !-----------------------
  ! CM of node coordinates
  !-----------------------
  
  if (.not. fds_cm_calculated) then
    cm_fds=0.0
    do i=1,nnodes_fds
      cm_fds(1)=cm_fds(1)+fds_xyz(i,1)  
      cm_fds(2)=cm_fds(2)+fds_xyz(i,2)  
      cm_fds(3)=cm_fds(3)+fds_xyz(i,3)  
    end do
    cm_fds=cm_fds/nnodes_fds
  end if

  fds_cm_calculated=.true.

  !----------------
  ! Number of nodes
  !----------------
  
  write(*,'(t3,2(a))') 'Number of nodes:                  ', trim(int2str(nnodes_fds))

  !----------------------------------
  ! Average nearest neighbor distance
  !----------------------------------

  allocate(distance(nnodes_fds),stat=ios);       call error_allocate(ios)
  allocate(dmask(nnodes_fds),stat=ios);          call error_allocate(ios)

  r_ave=0.0
  if (nnodes_fds > 1) then

    r_ave_sum=0.0
    outer_loop: do i=1,nnodes_fds
        
      ! Calculate distance to FDS nodes
      inner_loop: do j=1,nnodes_fds
        if (j /= i) then
          dx=fds_xyz(i,1)-fds_xyz(j,1)
          dy=fds_xyz(i,2)-fds_xyz(j,2)
          dz=fds_xyz(i,3)-fds_xyz(j,3)
          r_sq=dx*dx+dy*dy+dz*dz
          distance(j)=r_sq
        else
          distance(j)=huge(r_sq)
        end if
      end do inner_loop

      ! Find nearest neighbors
      dmask=.true.
      do j=1,mp_n
        mindist=minval(distance,1,dmask)
        imindist=minloc(distance,1,dmask) 
        dmask(imindist)=.false.
        r_ave_sum=r_ave_sum+sqrt(mindist)
      end do

    end do outer_loop

    r_ave=r_ave_sum/(nnodes_fds*mp_n)

  end if
    
  write(*,'(t3,a,es15.7e3)') 'Nearest neighbor distance:      ', r_ave
  
  write(*,'(a)') ''

  fmt_1='(t3,a,3(es15.7e3,a,))'
  write(*,fmt_1) 'Center of mass (nodes):          ', cm_fds(1),    ' ', cm_fds(2),    ' ', cm_fds(3)
  write(*,fmt_1) 'Center of mass (bounding box):   ', cm_fds_bb(1), ' ', cm_fds_bb(2), ' ', cm_fds_bb(3)

  write(*,'(a)') ''

  fmt_1='(t3,a,es15.7e3,a,es15.7e3)'
  write(*,fmt_1) 'Bounds in x-direction:           ', fds_node_range(1), ' ', fds_node_range(2)
  write(*,fmt_1) 'Bounds in y-direction:           ', fds_node_range(3), ' ', fds_node_range(4)
  write(*,fmt_1) 'Bounds in z-direction:           ', fds_node_range(5), ' ', fds_node_range(6)

  !--------------------------------
  ! Minimum and maximum data values
  !--------------------------------

  write(*,'(a)') ''
  
  fmt_1='(t3,a,es15.7e3)'
  write(*,fmt_1) 'Minimum data value:              ', minval(fds_data) 
  write(*,fmt_1) 'Maximum data value:              ', maxval(fds_data) 

  fds_average=0.0
  do i=1,ubound(fds_data,1)
    do j=1,ubound(fds_data,2)
      fds_average=fds_average+fds_data(i,j)
    end do
  end do
  fds_average=fds_average/(ubound(fds_data,1)*ubound(fds_data,2))

  write(*,fmt_1) 'Average data value:              ', fds_average

  !------------------
  ! Deallocate arrays
  !------------------

  deallocate(distance,stat=ios); call error_allocate(ios)
  deallocate(dmask,stat=ios);    call error_allocate(ios)

end subroutine analyze_fds_nodes

subroutine calculate_fds_cm()
!--------------------------
! Calculate CM of FDS nodes
!--------------------------
  use global_constants
  use global_variables
  use mapping_arrays
  implicit none

  integer :: i

  !-------------------
  ! CM of bounding box
  !-------------------

  fds_node_range(1)=minval(fds_xyz(1:nnodes_fds,1))
  fds_node_range(2)=maxval(fds_xyz(1:nnodes_fds,1))
  fds_node_range(3)=minval(fds_xyz(1:nnodes_fds,2)) 
  fds_node_range(4)=maxval(fds_xyz(1:nnodes_fds,2)) 
  fds_node_range(5)=minval(fds_xyz(1:nnodes_fds,3)) 
  fds_node_range(6)=maxval(fds_xyz(1:nnodes_fds,3)) 

  fds_dr(1)=fds_node_range(2)-fds_node_range(1)
  fds_dr(2)=fds_node_range(4)-fds_node_range(3)
  fds_dr(3)=fds_node_range(6)-fds_node_range(5)

  cm_fds_bb(1)=0.5*(fds_node_range(2)+fds_node_range(1))
  cm_fds_bb(2)=0.5*(fds_node_range(4)+fds_node_range(3))
  cm_fds_bb(3)=0.5*(fds_node_range(6)+fds_node_range(5))

  !-----------------------
  ! CM of node coordinates
  !-----------------------

  cm_fds=0.0
  do i=1,nnodes_fds
    cm_fds(1)=cm_fds(1)+fds_xyz(i,1)  
    cm_fds(2)=cm_fds(2)+fds_xyz(i,2)  
    cm_fds(3)=cm_fds(3)+fds_xyz(i,3)  
  end do
  cm_fds=cm_fds/nnodes_fds

  fds_cm_calculated=.true.

end subroutine calculate_fds_cm

end module fds_stats
