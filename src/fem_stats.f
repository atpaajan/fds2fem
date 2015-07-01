!------------------------------------------------------------
! Functions and subroutines for calculating ABAQUS statistics
!------------------------------------------------------------
module fem_stats

contains

subroutine fem_stats_module()
!---------
! Main sub
!---------
  use global_constants
  use global_variables
  implicit none

  write(*,'(a)') ''
  write(*,'(a)') 'FEM statistics module'

  if (fem_xyz_available .and. fem_statistics) then
    call analyze_fem_nodes()
  else
    write(*,'(t3,a)') 'Nothing to be done'
  end if

end subroutine fem_stats_module

!***************
! Auxiliary subs
!***************

subroutine analyze_fem_nodes()
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
  real(kind=rk) :: dx,dy,dz,r_sq,r_ave,r_ave_sum,mindist,fem_average

  logical, dimension(:), allocatable :: dmask
  real(kind=rk), dimension(:), allocatable :: distance

  nnodes_fem=ubound(fem_xyz,1)

  !-------------------
  ! CM of bounding box
  !-------------------

  fem_node_range(1)=minval(fem_xyz(1:nnodes_fem,1))
  fem_node_range(2)=maxval(fem_xyz(1:nnodes_fem,1))
  fem_node_range(3)=minval(fem_xyz(1:nnodes_fem,2)) 
  fem_node_range(4)=maxval(fem_xyz(1:nnodes_fem,2)) 
  fem_node_range(5)=minval(fem_xyz(1:nnodes_fem,3)) 
  fem_node_range(6)=maxval(fem_xyz(1:nnodes_fem,3)) 

  fem_dr(1)=fem_node_range(2)-fem_node_range(1)
  fem_dr(2)=fem_node_range(4)-fem_node_range(3)
  fem_dr(3)=fem_node_range(6)-fem_node_range(5)

  cm_fem_bb(1)=0.5*(fem_node_range(2)+fem_node_range(1))
  cm_fem_bb(2)=0.5*(fem_node_range(4)+fem_node_range(3))
  cm_fem_bb(3)=0.5*(fem_node_range(6)+fem_node_range(5))

  !-----------------------
  ! CM of node coordinates
  !-----------------------
  
  cm_fem=0.0
  do i=1,nnodes_fem
    cm_fem(1)=cm_fem(1)+fem_xyz(i,1)  
    cm_fem(2)=cm_fem(2)+fem_xyz(i,2)  
    cm_fem(3)=cm_fem(3)+fem_xyz(i,3)  
  end do
  cm_fem=cm_fem/nnodes_fem

  fem_cm_calculated=.true.

  !----------------
  ! Number of nodes
  !----------------
  
  write(*,'(t3,2(a))') 'Number of nodes:                  ', trim(int2str(nnodes_fem))

  !----------------------------------
  ! Average nearest neighbor distance
  !----------------------------------

  r_ave=0.0
  if (nnodes_fem > 1) then

    allocate(distance(nnodes_fem),stat=ios);       call error_allocate(ios)
    allocate(dmask(nnodes_fem),stat=ios);          call error_allocate(ios)

    r_ave_sum=0.0
    outer_loop: do i=1,nnodes_fem
        
      ! Calculate distance to FDS nodes
      inner_loop: do j=1,nnodes_fem
        if (j /= i) then
          dx=fem_xyz(i,1)-fem_xyz(j,1)
          dy=fem_xyz(i,2)-fem_xyz(j,2)
          dz=fem_xyz(i,3)-fem_xyz(j,3)
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

    r_ave=r_ave_sum/(nnodes_fem*mp_n)

  end if

  write(*,'(t3,a,es15.7e3)') 'Nearest neighbor distance:      ', r_ave
  
  write(*,'(a)') ''

  fmt_1='(t3,a,3(es15.7e3,a,))'
  write(*,fmt_1) 'Center of mass (nodes):          ', cm_fem(1),    ' ', cm_fem(2),    ' ', cm_fem(3)
  write(*,fmt_1) 'Center of mass (bounding box):   ', cm_fem_bb(1), ' ', cm_fem_bb(2), ' ', cm_fem_bb(3)

  write(*,'(a)') ''

  fmt_1='(t3,a,es15.7e3,a,es15.7e3)'
  write(*,fmt_1) 'Bounds in x-direction:           ', fem_node_range(1), ' ', fem_node_range(2)
  write(*,fmt_1) 'Bounds in y-direction:           ', fem_node_range(3), ' ', fem_node_range(4)
  write(*,fmt_1) 'Bounds in z-direction:           ', fem_node_range(5), ' ', fem_node_range(6)
 
  !--------------------------------
  ! Minimum and maximum data values
  !--------------------------------

  write(*,'(a)') ''
  
  fmt_1='(t3,a,es15.7e3)'
  write(*,fmt_1) 'Minimum data value:              ', minval(fem_data) 
  write(*,fmt_1) 'Maximum data value:              ', maxval(fem_data) 

  fem_average=0.0
  do i=1,ubound(fem_data,1)
    do j=1,ubound(fem_data,2)
      fem_average=fem_average+fem_data(i,j)
    end do
  end do
  fem_average=fem_average/(ubound(fem_data,1)*ubound(fem_data,2))

  write(*,fmt_1) 'Average data value:              ', fem_average

  !------------------
  ! Deallocate arrays
  !------------------

  deallocate(distance,stat=ios); call error_allocate(ios)
  deallocate(dmask,stat=ios);    call error_allocate(ios)

end subroutine analyze_fem_nodes

subroutine calculate_fem_cm()
!-----------------------------
! Calculate CM of ABAQUS nodes
!-----------------------------
  use global_constants
  use global_variables
  use mapping_arrays
  implicit none

  integer :: i

  nnodes_fem=ubound(fem_xyz,1)

  !-------------------
  ! CM of bounding box
  !-------------------

  fem_node_range(1)=minval(fem_xyz(1:nnodes_fem,1))
  fem_node_range(2)=maxval(fem_xyz(1:nnodes_fem,1))
  fem_node_range(3)=minval(fem_xyz(1:nnodes_fem,2)) 
  fem_node_range(4)=maxval(fem_xyz(1:nnodes_fem,2)) 
  fem_node_range(5)=minval(fem_xyz(1:nnodes_fem,3)) 
  fem_node_range(6)=maxval(fem_xyz(1:nnodes_fem,3)) 

  fem_dr(1)=fem_node_range(2)-fem_node_range(1)
  fem_dr(2)=fem_node_range(4)-fem_node_range(3)
  fem_dr(3)=fem_node_range(6)-fem_node_range(5)

  cm_fem_bb(1)=0.5*(fem_node_range(2)+fem_node_range(1))
  cm_fem_bb(2)=0.5*(fem_node_range(4)+fem_node_range(3))
  cm_fem_bb(3)=0.5*(fem_node_range(6)+fem_node_range(5))

  !-----------------------
  ! CM of node coordinates
  !-----------------------

  cm_fem=0.0
  do i=1,nnodes_fem
    cm_fem(1)=cm_fem(1)+fem_xyz(i,1)  
    cm_fem(2)=cm_fem(2)+fem_xyz(i,2)  
    cm_fem(3)=cm_fem(3)+fem_xyz(i,3)  
  end do
  cm_fem=cm_fem/nnodes_fem

  fem_cm_calculated=.true.

end subroutine calculate_fem_cm

end module fem_stats
