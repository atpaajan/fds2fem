!---------------------------------------------
! Functions and subroutines for model matching
!---------------------------------------------
module matching

contains

subroutine matching_module()
!---------
! Main sub
!---------
  use global_constants
  use global_variables
  implicit none

  write(*,'(a)') ''
  write(*,'(a)') 'Model matching module'

  if (match_translate .or. match_rotate) then
    
    if (manual_translate) then
      call translate_manual()
    else if (automatic_translate) then
      call translate_automatic()
    end if

    if (manual_rotate) then
      if (fem_xyz_available) then
        call rotate_manual()
      else
        write(*,'(a)') 'Nothing to be done'
      end if
    else if (automatic_rotate) then
      if (fds_xyz_available .and. fem_xyz_available) then
        call rotate_automatic()
      else
        write(*,'(a)') 'Nothing to be done'
      end if
    end if

  else
    write(*,'(t3,a)') 'Nothing to be done'
  end if

end subroutine matching_module

!***************
! Auxiliary subs
!***************

subroutine translate_manual()
!-------------------------------------------------
! Move the user-defined reference points to origin
!-------------------------------------------------
  use global_constants
  use global_variables
  use mapping_arrays
  implicit none

  integer :: i
  real(kind=rk) :: dx,dy,dz,dr

  if (.not. manual_translate) then
    ! Nothing to be done
    return
  end if
  
  write(*,'(t3,a)') 'Manual translation'
  
  if (fds_xyz_available) then

    if (maxval(origin_fds) > eps) then
      ! Translate FDS nodes
      translate_fds: do i=1,nnodes_fds
        fds_xyz(i,1)=fds_xyz(i,1)-origin_fds(1) 
        fds_xyz(i,2)=fds_xyz(i,2)-origin_fds(2)
        fds_xyz(i,3)=fds_xyz(i,3)-origin_fds(3)
      end do translate_fds
    else
      write(*,'(t3,a)') '--< WARNING: small displacement in FDS nodes - skipping translation'
    end if

    origin_fds=0.0
  end if

  if (fem_xyz_available) then

    if (maxval(origin_fem) > eps) then
      ! Translate FDS nodes
      translate_fem: do i=1,nnodes_fem
        fem_xyz(i,1)=fem_xyz(i,1)-origin_fem(1) 
        fem_xyz(i,2)=fem_xyz(i,2)-origin_fem(2)
        fem_xyz(i,3)=fem_xyz(i,3)-origin_fem(3)
      end do translate_fem
    else
      write(*,'(t3,a)') '--< WARNING: small displacement in FEM nodes - skipping translation'
    end if

    origin_fem=0.0
  end if

end subroutine translate_manual

subroutine rotate_manual()
!--------------------------------------
! Rotate  nodes around the z-axis
!--------------------------------------
  use global_constants
  use global_variables
  use mapping_arrays
  use mathematics
  implicit none

  integer :: i

  real(kind=rk) :: phi

  real(kind=rk), dimension(3)   :: v
  real(kind=rk), dimension(3,3) :: rm

  if (.not. manual_rotate) then
    ! Nothing to be done
    return
  end if

  if (.not. fem_xyz_available) then
    ! Nothing to be done
    return
  end if

  write(*,'(t3,a)') 'Manual rotation'

  if (abs(e_alpha) < eps .and. &
      abs(e_beta)  < eps .and. &
      abs(e_gamma) < eps) then 
    ! No need to rotate if all of the angles are very small
    write(*,'(t3,a)') '--< WARNING: small Euler angles - skipping rotation'
    return
  end if

  rm=euler(e_alpha,e_beta,e_gamma)
 
  rotate_loop: do i=1,nnodes_fem
    v(1:3)=fem_xyz(i,1:3)
    v=matmul(rm,v)
    fem_xyz(i,1:3)=v(1:3)
  end do rotate_loop

end subroutine rotate_manual

subroutine translate_automatic
!----------------------------------------------------
! Move the center of mass of both node sets to origin
!----------------------------------------------------
  use fem_stats
  use fds_stats
  use global_constants
  use global_variables
  use mapping_arrays
  implicit none

  integer :: i

  if (.not. automatic_translate) then
    ! Nothing to be done
    return
  end if

  if (.not. fds_xyz_available .and. .not. fem_xyz_available) then
    ! Nothing to be done
    return
  end if

  write(*,'(t3,a)') 'Automatic center of mass adjustment'

  !-------------------
  ! CM of bounding box
  !-------------------

  if (fds_xyz_available .and. .not. fds_cm_calculated) then
    call calculate_fds_cm()
  end if

  if (fem_xyz_available .and. .not. fem_cm_calculated) then
    call calculate_fem_cm()
  end if

  ! * Full output
  if (full_output) then
    write(*,'(a)') ''
    write(*,'(a)') '* Center of mass '
    if (fds_cm_calculated) then 
      write(*,'(a,3x,a,3(es15.7e3,a))') '*', 'FDS:    (', &
        cm_fds(1), ', ',cm_fds(2), ', ', cm_fds(3), ')'
      end if
    if (fem_cm_calculated) then
      write(*,'(a,3x,a,3(es15.7e3,a))') '*', 'FEM: (', &
        cm_fem(1), ', ',cm_fem(2), ', ', cm_fem(3), ')'
      write(*,'(a)') ''
    end if
  end if

  !-------------
  ! Safety check
  !-------------

  if (maxval(cm_fem,1) < eps .and. &
      maxval(cm_fds,1)  < eps) then
    ! No need to rotate if either of the models has only one node 
    write(*,'(t3,a)') '--< WARNING: CM of both models at origin - skipping translation'
    return
  end if

  ! Do something

  !-------------------
  ! Move CMs to origin
  !-------------------

  nnodes_fem=ubound(fem_xyz,1)
  if (abs(maxval(cm_fem,1)) > eps) then
    do i=1,nnodes_fem
      fem_xyz(i,1)=fem_xyz(i,1)-cm_fem(1)
      fem_xyz(i,2)=fem_xyz(i,2)-cm_fem(2)
      fem_xyz(i,3)=fem_xyz(i,3)-cm_fem(3)
    end do
    cm_fem=0.0
  end if

  nnodes_fds=ubound(fds_xyz,1)
  if (abs(maxval(cm_fds,1)) > eps) then
    do i=1,nnodes_fds
      fds_xyz(i,1)=fds_xyz(i,1)-cm_fds(1)
      fds_xyz(i,2)=fds_xyz(i,2)-cm_fds(2)
      fds_xyz(i,3)=fds_xyz(i,3)-cm_fds(3)
    end do
    cm_fds=0.0
  end if

end subroutine translate_automatic

subroutine rotate_automatic()
!------------------------------------------------
! Rotate  node set to find an optimum match
!------------------------------------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use mathematics
  use string_handling
  implicit none

  integer :: i,j,ios,nrows,ncols,nphi_full,nsearch
  integer :: mphi_full_imin,mphi_old_imin,mphi_new_imin

  logical :: looks_bad

  real(kind=rk) :: mphi_full_min,mphi_old_min,mphi_new_min,dphi
  real(kind=rk), dimension(3)   :: v
  real(kind=rk), dimension(3,3) :: rm
  real(kind=rk), dimension(:,:), allocatable :: mphi_full,mphi_old,mphi_new
  real(kind=rk), dimension(:,:), allocatable :: fem_xyz_tmp

  if (.not. automatic_rotate) then
    ! Nothing to be done
    return
  end if

  if (.not. fds_xyz_available .and. .not. fem_xyz_available) then
    ! Nothing to be done
    return
  end if

  write(*,'(t3,a)') 'Automatic azimuthal angle adjustment'
  
  if (nnodes_fds < 2) then
    ! No need to rotate if either of the models has only one node 
    write(*,'(t3,a)') '--< WARNING: not enough FDS nodes - skipping rotation'
    return
  end if

  if (nnodes_fem < 2) then
    ! No need to rotate if either of the models has only one node 
    write(*,'(t3,a)') '--< WARNING: not enough FEM nodes - skipping rotation'
    return
  end if

  !-----------
  ! Parameters
  !-----------

  nphi_full = 8
  nsearch   = 5

  !----------------------------------------------
  ! Make a backup copy of  node coordinates
  !----------------------------------------------

  nrows=ubound(fem_xyz,1)
  ncols=ubound(fem_xyz,2)
  
  allocate(fem_xyz_tmp(nrows,ncols),stat=ios); call error_allocate(ios)
  fem_xyz_tmp=fem_xyz

  !-----------------------
  ! Coarse 360 degree scan
  !-----------------------

  if (nphi_full < 5) then
    write(*,'(a)') 'ERROR: something that should not happen'
    stop
  end if

  allocate(mphi_full(nphi_full,2),stat=ios); call error_allocate(ios)

  mphi_full(1,1)=0.0
  mphi_full(1,2)=measure()

  ! * Full output
  if (full_output) then
    write(*,'(a)') ''
    write(*,'(a)') '* Rotation angle vs. measure '
    write(*,'(a,3x,a,f7.3,a,es15.7e3)') '*', 'Phi= ', 180.0*mphi_full(1,1)/(4.0*atan(1.0)), ', &
      f(Phi)= ', mphi_full(1,2)
  end if

  do i=2,nphi_full
    mphi_full(i,1)=(8.0*atan(1.0)/nphi_full)*(i-1)
    rm=rotate_z(mphi_full(i,1))
    do j=1,nnodes_fem
      v(1:3)=fem_xyz_tmp(j,1:3)
      v=matmul(rm,v)
      fem_xyz(j,1:3)=v(1:3)
    end do
    mphi_full(i,2)=measure()
  end do

  mphi_full_min=minval(mphi_full(1:nphi_full,2),1)
  mphi_full_imin=minloc(mphi_full(1:nphi_full,2),1)

  if (mphi_full_imin /= 1) then
    if (abs(mphi_full(1,2)-mphi_full(mphi_full_imin,2)) < eps) then
      mphi_full_min=mphi_full(1,2)
      mphi_full_imin=1
    end if
  end if

  ! * Full output
  if (full_output) then
    write(*,'(a,3x,a,f7.3,a,es15.7e3,a,f7.3)') '*', 'Phi= ', 180.0*mphi_full(mphi_full_imin,1)/(4.0*atan(1.0)), &
      ', f(Phi)= ', mphi_full(mphi_full_imin,2), ', dPhi= ', (360.0/nphi_full)
  end if

  !-------------------------------------
  ! If there are problems, skip rotation
  !-------------------------------------

  looks_bad=.true.
  check_equal_1: do i=1,nphi_full-1
    if (abs(mphi_full(i,2)-mphi_full(i+1,2)) > eps) then
      looks_bad=.false.
      exit check_equal_1
    end if
  end do check_equal_1

  if (looks_bad) then
    write(*,'(t3,a)') '--< WARNING: problems in finding optimum angle - skipping rotation'
    fem_xyz=fem_xyz_tmp
    return
  end if

  !------------------------
  ! A little bit finer scan
  !------------------------

  allocate(mphi_new(5,2),stat=ios); call error_allocate(ios) 
 
  if (mphi_full_imin == 1) then
    mphi_new(1,1:2)=mphi_full(nphi_full-1,1:2) 
    mphi_new(2,1:2)=mphi_full(nphi_full,1:2) 
    mphi_new(3,1:2)=mphi_full(1,1:2)
    mphi_new(4,1:2)=mphi_full(2,1:2)
    mphi_new(5,1:2)=mphi_full(3,1:2)
  else if (mphi_full_imin == 2) then
    mphi_new(1,1:2)=mphi_full(nphi_full,1:2) 
    mphi_new(2,1:2)=mphi_full(1,1:2) 
    mphi_new(3,1:2)=mphi_full(2,1:2)
    mphi_new(4,1:2)=mphi_full(3,1:2)
    mphi_new(5,1:2)=mphi_full(4,1:2)
  else if (mphi_full_imin == nphi_full) then
    mphi_new(1,1:2)=mphi_full(nphi_full-2,1:2) 
    mphi_new(2,1:2)=mphi_full(nphi_full-1,1:2)
    mphi_new(3,1:2)=mphi_full(nphi_full,1:2)
    mphi_new(4,1:2)=mphi_full(1,1:2)
    mphi_new(5,1:2)=mphi_full(2,1:2)
  else if (mphi_full_imin == nphi_full-1) then
    mphi_new(1,1:2)=mphi_full(nphi_full-3,1:2) 
    mphi_new(2,1:2)=mphi_full(nphi_full-2,1:2)
    mphi_new(3,1:2)=mphi_full(nphi_full-1,1:2)
    mphi_new(4,1:2)=mphi_full(nphi_full,1:2)
    mphi_new(5,1:2)=mphi_full(1,1:2)
  else
    mphi_new(1,1:2)=mphi_full(mphi_full_imin-2,1:2) 
    mphi_new(2,1:2)=mphi_full(mphi_full_imin-1,1:2) 
    mphi_new(3,1:2)=mphi_full(mphi_full_imin,1:2)
    mphi_new(4,1:2)=mphi_full(mphi_full_imin+1,1:2)
    mphi_new(5,1:2)=mphi_full(mphi_full_imin+2,1:2)
  end if

  mphi_new_min  = mphi_new(3,2)
  mphi_new_imin = 3

  !---------------
  ! Iterative scan
  !---------------

  allocate(mphi_old(5,2),stat=ios); call error_allocate(ios) 

  dphi=(360.0/nphi_full)
  search_loop: do i=1,nsearch

    ! New becomes old
    mphi_old_min  = mphi_new_min
    mphi_old_imin = mphi_new_imin
    mphi_old      = mphi_new

    if (mphi_old_imin == 1 .or. mphi_old_imin == 5) then
      ! This should not happen
      write(*,'(a)') 'ERROR: in automatic azimuthal angle adjustment algorithm'
      stop
    else
      mphi_new(1,1:2)=mphi_old(mphi_old_imin-1,1:2) 
      mphi_new(3,1:2)=mphi_old(mphi_old_imin,1:2)
      mphi_new(5,1:2)=mphi_old(mphi_old_imin+1,1:2)
    end if

    ! First new measure
    mphi_new(2,1)=0.5*(mphi_new(1,1)+mphi_new(3,1))

    rm=rotate_z(mphi_new(2,1))
    do j=1,nnodes_fem
      v(1:3)=fem_xyz_tmp(j,1:3)
      v=matmul(rm,v)
      fem_xyz(j,1:3)=v(1:3)
    end do

    mphi_new(2,2)=measure()

    ! Second new measure
    mphi_new(4,1)=0.5*(mphi_new(3,1)+mphi_new(5,1))

    rm=rotate_z(mphi_new(4,1))
    do j=1,nnodes_fem
      v(1:3)=fem_xyz_tmp(j,1:3)
      v=matmul(rm,v)
      fem_xyz(j,1:3)=v(1:3)
    end do

    mphi_new(4,2)=measure()

    mphi_new_min=minval(mphi_new(1:5,2),1)
    mphi_new_imin=minloc(mphi_new(1:5,2),1)

    ! * Full output
    if (full_output) then
      dphi=dphi*0.5
      write(*,'(a,3x,a,f7.3,a,es15.7e3,a,f7.3)') '*', 'Phi= ', 180.0*mphi_new(mphi_new_imin,1)/(4.0*atan(1.0)), &
        ', f(Phi)= ', mphi_new(mphi_new_imin,2), ', dPhi= ', dphi
    end if

  end do search_loop

  !------------------------
  ! Apply selected rotation
  !------------------------

  rm=rotate_z(mphi_new(mphi_new_imin,1))
  do j=1,nnodes_fem
    v(1:3)=fem_xyz_tmp(j,1:3)
    v=matmul(rm,v)
    fem_xyz(j,1:3)=v(1:3)
  end do

end subroutine rotate_automatic

function measure()
!--------------------------------------------------------
! Measure used in optimizing rotation angle around z-axis
!--------------------------------------------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use mathematics
  implicit none
 
  integer :: i,j
  real(kind=rk) :: dr_sq,dx,dy,dz
  real(kind=rk) :: measure,mindist,mindist_sum

  mindist_sum=0.0
  fem_loop: do i=1,nnodes_fem
    mindist=huge(mindist)
    fds_loop: do j=1,nnodes_fds
      dx=fem_xyz(i,1)-fds_xyz(j,1)
      dy=fem_xyz(i,2)-fds_xyz(j,2)
      dr_sq=(dx*dx)+(dy*dy)
      if (dr_sq < mindist) mindist=dr_sq
    end do fds_loop
    mindist_sum=mindist_sum+mindist
  end do fem_loop

  measure=mindist_sum/nnodes_fem

end function measure

end module matching
