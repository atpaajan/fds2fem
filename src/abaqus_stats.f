!------------------------------------------------------------
! Functions and subroutines for calculating ABAQUS statistics
!------------------------------------------------------------
module abaqus_stats

contains

subroutine abaqus_stats_module()
!---------
! Main sub
!---------
  use global_constants
  use global_variables
  implicit none

  write(*,'(a)') ''
  write(*,'(a)') 'ABAQUS statistics module'

  if (abaqus_xyz_available .and. abaqus_statistics) then
    call analyze_abaqus_nodes()
  else
    write(*,'(t3,a)') 'Nothing to be done'
  end if

end subroutine abaqus_stats_module

!***************
! Auxiliary subs
!***************

subroutine analyze_abaqus_nodes()
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
  real(kind=rk) :: dx,dy,dz,r_sq,r_ave,r_ave_sum,mindist,abaqus_average

  logical, dimension(:), allocatable :: dmask
  real(kind=rk), dimension(:), allocatable :: distance

  !-------------------
  ! CM of bounding box
  !-------------------

  abaqus_node_range(1)=minval(abaqus_xyz(1:nnodes_abaqus,1))
  abaqus_node_range(2)=maxval(abaqus_xyz(1:nnodes_abaqus,1))
  abaqus_node_range(3)=minval(abaqus_xyz(1:nnodes_abaqus,2)) 
  abaqus_node_range(4)=maxval(abaqus_xyz(1:nnodes_abaqus,2)) 
  abaqus_node_range(5)=minval(abaqus_xyz(1:nnodes_abaqus,3)) 
  abaqus_node_range(6)=maxval(abaqus_xyz(1:nnodes_abaqus,3)) 

  abaqus_dr(1)=abaqus_node_range(2)-abaqus_node_range(1)
  abaqus_dr(2)=abaqus_node_range(4)-abaqus_node_range(3)
  abaqus_dr(3)=abaqus_node_range(6)-abaqus_node_range(5)

  cm_abaqus_bb(1)=0.5*(abaqus_node_range(2)+abaqus_node_range(1))
  cm_abaqus_bb(2)=0.5*(abaqus_node_range(4)+abaqus_node_range(3))
  cm_abaqus_bb(3)=0.5*(abaqus_node_range(6)+abaqus_node_range(5))

  !-----------------------
  ! CM of node coordinates
  !-----------------------
  
  cm_abaqus=0.0
  do i=1,nnodes_abaqus
    cm_abaqus(1)=cm_abaqus(1)+abaqus_xyz(i,1)  
    cm_abaqus(2)=cm_abaqus(2)+abaqus_xyz(i,2)  
    cm_abaqus(3)=cm_abaqus(3)+abaqus_xyz(i,3)  
  end do
  cm_abaqus=cm_abaqus/nnodes_abaqus

  abaqus_cm_calculated=.true.

  !----------------
  ! Number of nodes
  !----------------
  
  write(*,'(t3,2(a))') 'Number of nodes:                  ', trim(int2str(nnodes_abaqus))

  !----------------------------------
  ! Average nearest neighbor distance
  !----------------------------------

  r_ave=0.0
  if (nnodes_abaqus > 1) then

    allocate(distance(nnodes_abaqus),stat=ios);       call error_allocate(ios)
    allocate(dmask(nnodes_abaqus),stat=ios);          call error_allocate(ios)

    r_ave_sum=0.0
    outer_loop: do i=1,nnodes_abaqus
        
      ! Calculate distance to FDS nodes
      inner_loop: do j=1,nnodes_abaqus
        if (j /= i) then
          dx=abaqus_xyz(i,1)-abaqus_xyz(j,1)
          dy=abaqus_xyz(i,2)-abaqus_xyz(j,2)
          dz=abaqus_xyz(i,3)-abaqus_xyz(j,3)
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

    r_ave=r_ave_sum/(nnodes_abaqus*mp_n)

  end if

  write(*,'(t3,a,es15.7e3)') 'Nearest neighbor distance:      ', r_ave
  
  write(*,'(a)') ''

  fmt_1='(t3,a,3(es15.7e3,a,))'
  write(*,fmt_1) 'Center of mass (nodes):          ', cm_abaqus(1),    ' ', cm_abaqus(2),    ' ', cm_abaqus(3)
  write(*,fmt_1) 'Center of mass (bounding box):   ', cm_abaqus_bb(1), ' ', cm_abaqus_bb(2), ' ', cm_abaqus_bb(3)

  write(*,'(a)') ''

  fmt_1='(t3,a,es15.7e3,a,es15.7e3)'
  write(*,fmt_1) 'Bounds in x-direction:           ', abaqus_node_range(1), ' ', abaqus_node_range(2)
  write(*,fmt_1) 'Bounds in y-direction:           ', abaqus_node_range(3), ' ', abaqus_node_range(4)
  write(*,fmt_1) 'Bounds in z-direction:           ', abaqus_node_range(5), ' ', abaqus_node_range(6)
 
  !--------------------------------
  ! Minimum and maximum data values
  !--------------------------------

  write(*,'(a)') ''
  
  fmt_1='(t3,a,es15.7e3)'
  write(*,fmt_1) 'Minimum data value:              ', minval(abaqus_data) 
  write(*,fmt_1) 'Maximum data value:              ', maxval(abaqus_data) 

  abaqus_average=0.0
  do i=1,ubound(abaqus_data,1)
    do j=1,ubound(abaqus_data,2)
      abaqus_average=abaqus_average+abaqus_data(i,j)
    end do
  end do
  abaqus_average=abaqus_average/(ubound(abaqus_data,1)*ubound(abaqus_data,2))

  write(*,fmt_1) 'Average data value:              ', abaqus_average

  !------------------
  ! Deallocate arrays
  !------------------

  deallocate(distance,stat=ios); call error_allocate(ios)
  deallocate(dmask,stat=ios);    call error_allocate(ios)

end subroutine analyze_abaqus_nodes

subroutine calculate_abaqus_cm()
!-----------------------------
! Calculate CM of ABAQUS nodes
!-----------------------------
  use global_constants
  use global_variables
  use mapping_arrays
  implicit none

  integer :: i

  !-------------------
  ! CM of bounding box
  !-------------------

  abaqus_node_range(1)=minval(abaqus_xyz(1:nnodes_abaqus,1))
  abaqus_node_range(2)=maxval(abaqus_xyz(1:nnodes_abaqus,1))
  abaqus_node_range(3)=minval(abaqus_xyz(1:nnodes_abaqus,2)) 
  abaqus_node_range(4)=maxval(abaqus_xyz(1:nnodes_abaqus,2)) 
  abaqus_node_range(5)=minval(abaqus_xyz(1:nnodes_abaqus,3)) 
  abaqus_node_range(6)=maxval(abaqus_xyz(1:nnodes_abaqus,3)) 

  abaqus_dr(1)=abaqus_node_range(2)-abaqus_node_range(1)
  abaqus_dr(2)=abaqus_node_range(4)-abaqus_node_range(3)
  abaqus_dr(3)=abaqus_node_range(6)-abaqus_node_range(5)

  cm_abaqus_bb(1)=0.5*(abaqus_node_range(2)+abaqus_node_range(1))
  cm_abaqus_bb(2)=0.5*(abaqus_node_range(4)+abaqus_node_range(3))
  cm_abaqus_bb(3)=0.5*(abaqus_node_range(6)+abaqus_node_range(5))

  !-----------------------
  ! CM of node coordinates
  !-----------------------

  cm_abaqus=0.0
  do i=1,nnodes_abaqus
    cm_abaqus(1)=cm_abaqus(1)+abaqus_xyz(i,1)  
    cm_abaqus(2)=cm_abaqus(2)+abaqus_xyz(i,2)  
    cm_abaqus(3)=cm_abaqus(3)+abaqus_xyz(i,3)  
  end do
  cm_abaqus=cm_abaqus/nnodes_abaqus

  abaqus_cm_calculated=.true.

end subroutine calculate_abaqus_cm

end module abaqus_stats
