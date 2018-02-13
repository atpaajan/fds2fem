!-------------------------------------------
! Functions and subroutines for mesh mapping
!-------------------------------------------
module mapping

contains

subroutine mapping_module()
!---------
! Main sub
!---------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,k,ios
  integer :: nnodes_fds_sub,nnodes_fem_sub

  write(*,'(a)') ''
  write(*,'(a)') 'Mesh mapping module'

  !--------------------------------
  ! Sometimes there's nothing to do
  !--------------------------------

  if (trim(mapping_method) == 'off') then
    write(*,'(t3,a)') 'Nothing to be done'
    fem_data_available=.false.
    return
  end if

  if (.not. fds_xyz_available .And. .Not. cfast_input) then
    write(*,'(t3,a)') 'WARNING: FDS node coordinates unavailable'
    fem_data_available=.false.
    return
  end if

  if (.not. fds_data_available) then
    write(*,'(t3,a)') 'WARNING: FDS data unavailable'
    fem_data_available=.false.
    return  
  end if

  if (.not. fem_xyz_available) then
    write(*,'(t3,a)') 'WARNING: FEM node coordinates unavailable'
    fem_data_available=.false.
    return
  end if

  If (cfast_input .And. .Not. Trim(mapping_method) == 'devc_to_nset') Then
    Write(*,'(a,a)') 'ERROR: CFAST input and mapping method: ',Trim(mapping_method)
    Stop
  End If

  If (cfast_input .And. .Not. nset_connectivity) Then
    Write(*,'(a)') 'ERROR: CFAST input and no nset_connectivity '
    Stop
  End If

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

  if (.not. nset_connectivity) then
    !----------------------------------------------------------
    ! The simplest case: no NSET-DEVC or NSET-BNDF connectivity
    !----------------------------------------------------------

    fds_mask=.true.; fem_mask=.true.

    select case (trim(mapping_method))
    case ('nearest')
      call map_nearest()
    case default
      ! Nothing to be done
    end select

    fem_data_available=.true.

  else
    !------------------------------------
    ! NSET-DEVC or NSET-BNDF connectivity
    !------------------------------------

    if (trim(mapping_method) == 'devc_to_nset') then
      !------------------------------
      ! Direct NSET-DEVC connectivity
      !------------------------------

      call map_devc_to_nset()

      fem_data_available=.true.

    else
      !------------------------------------
      ! General NSET-DEVC/BNDF connectivity
      !------------------------------------
   
      do i=1,ubound(connectivity_table,1)

        !----------------------------------------------------------------------
        ! Create a connectivity table -based logical mask for FDS and FEM nodes
        !----------------------------------------------------------------------

        nnodes_fds_sub=0
        fds_mask=.false.
        fds_mask_loop: do k=1,nnodes_fds
          if (trim(fds_output) == 'devc') then
            !------------
            ! DEVC output
            !------------
            if (sum(connectivity_table_num(i,2:ubound(connectivity_table_num,2))) == 0) then
              fds_mask(k)=.true. 
              nnodes_fds_sub=nnodes_fds_sub+1
              cycle fds_mask_loop
            else
              do j=2,ubound(connectivity_table_num,2)
                if (fds_idevc(k) == connectivity_table_num(i,j)) then
                  fds_mask(k)=.true.
                  nnodes_fds_sub=nnodes_fds_sub+1
                  cycle fds_mask_loop
                end if
              end do
            end if
          else if (trim(fds_output) == 'bndf') then
            !------------
            ! BNDF output
            !------------
            if (sum(connectivity_table_num(i,2:ubound(connectivity_table_num,2))) == 0) then
              fds_mask(k)=.true. 
              nnodes_fds_sub=nnodes_fds_sub+1
              cycle fds_mask_loop
            else
              do j=2,ubound(connectivity_table_num,2)
                if (fds_patch(k) == connectivity_table_num(i,j)) then
                  fds_mask(k)=.true.
                  nnodes_fds_sub=nnodes_fds_sub+1
                  cycle fds_mask_loop
                end if
              end do 
            end if
          end if
        end do fds_mask_loop

        if (nnodes_fds_sub == 0) then
          write(*,'(5(a))') 'ERROR: in NSET connectivity file: no FDS nodes found (file ', &
            trim(quote(nset_input_file)), ', line ', trim(int2str(i)), ')'
            stop
        end if

        nnodes_fem_sub=0
        fem_mask=.false.
        fem_mask_loop: do k=1,nnodes_fem
          if (trim(fem_nset(k)) == trim(connectivity_table(i,1))) then
            fem_mask(k)=.true.
            nnodes_fem_sub=nnodes_fem_sub+1
            cycle fem_mask_loop
          end if
        end do fem_mask_loop
        
        if (nnodes_fem_sub == 0) then
          write(*,'(5(a))') 'ERROR: in NSET connectivity file: no FEM nodes found (file ', &
            trim(quote(nset_input_file)), ', line ', trim(int2str(i)), ')'
            stop
        end if

        select case (trim(mapping_method))
        case ('nearest')
          call map_nearest()
        case default
          ! Nothing to be done
        end select

      end do

      !------------------
      ! Deallocate memory
      !------------------

      deallocate(fds_mask,stat=ios); call error_allocate(ios)
      deallocate(fem_mask,stat=ios); call error_allocate(ios)

      fem_data_available=.true. 

    end if

  end if

end subroutine mapping_module

!***************
! Auxiliary subs
!***************

subroutine map_nearest()
!---------------------------------------------------
! K-nearest mapping with an optional cut-off radius.
! Does also simple radius mapping.
!
! If k-nearest mapping is used then at least one NN
! is used regardless of the actual distance so no 
! orphans are allowed.
!
! If a simple radius mapping is used then an orphan
! treatment should be done, but at this version this
! is not done, the program is just stopped.
!
! Parameters (user given)
!    mp_n: number of nearest neighbors to be used
!    mp_nmx: maximum number of nn points allowed
!    mp_cut: cut-off radius
!    mp_del: search nn points a little bit further
!            to find outliers up to mp_nmx points
!    mp_deg: distance weight power (1/r^deg 
!            distance dependence)
!---------------------------------------------------
  use mapping_arrays
  use error_messages
  use global_constants
  use global_variables
  use string_handling
  implicit none

  logical, save :: first_call=.true.

  integer :: i,j,k,ios,imindist
  integer :: nn_r_cut,nn_dr_cut,count_new,nn_n, nn_nmx
  integer, dimension(:,:), allocatable :: neighbor
  
  logical, dimension(:), allocatable :: dmask
  
  real(kind=rk) :: dx,dy,dz,mindist,w_sum,r_now,weight
  real(kind=rk) :: r_nn,r_sq,r_cut,r_cut_sq,dr_cut_sq
  real(kind=rk), dimension(:), allocatable :: distance
  
  real(kind=rk), dimension(:,:), allocatable :: weight_neighbor

  !---------------
  ! Initialization
  !---------------
  
  nnodes_fds=ubound(fds_xyz,1)
  nnodes_fem=ubound(fem_xyz,1)
  
  if (mp_n < 1) then
    !---------------------
    ! Plain radius mapping
    !---------------------
    nn_n=nnodes_fds
    nn_nmx=nnodes_fds
  else
    !------------------
    ! K-nearest mapping
    !------------------
    nn_n=min(mp_n,nnodes_fds)
    nn_nmx=min(mp_nmx,nnodes_fds)
  end if

  allocate(neighbor(nnodes_fem,nn_nmx),stat=ios);        call error_allocate(ios)
  allocate(weight_neighbor(nnodes_fem,nn_nmx),stat=ios); call error_allocate(ios)
  allocate(distance(nnodes_fds),stat=ios);               call error_allocate(ios)
  allocate(dmask(nnodes_fds),stat=ios);                  call error_allocate(ios)

  !--------------------------
  ! Distance weighted mapping
  !--------------------------

  r_cut=mp_cut
  weight_neighbor=0.0

  if (mp_n < 1) then
    !---------------------
    ! Plain radius mapping
    !---------------------
    r_cut_sq=r_cut**2 
    dr_cut_sq=((1.0+mp_del)*r_cut)**2
    r_nn=r_cut

    if (first_call) then
      write(*,'(t3,3(a))')   'Using plain radius mapping'
      write(*,'(t6,a,f5.3)') 'Radius cut-off:                      ', sqrt(r_cut_sq)
      write(*,'(t6,a,f5.3)') 'Power of distance-weight dependence: ', mp_deg 
    end if

  else
    if (r_cut > epsilon(r_cut)) then
      !-------------------------------
      ! K-nearest mapping with cut-off
      !-------------------------------
      r_cut_sq=r_cut**2 
      dr_cut_sq=((1.0+mp_del)*r_cut)**2
      r_nn=r_cut
    
      if (first_call) then
        write(*,'(t3,3(a))') 'Using k-nearest mapping with a cut-off radius'
        write(*,'(t6,a,f5.3,a,f5.3)') 'Nearest neighbors cut-offs:          ', sqrt(r_cut_sq),', ',sqrt(dr_cut_sq)
      end if
    
    else
      !----------------------------------
      ! K-nearest mapping without cut-off
      !----------------------------------
      r_cut_sq=huge(r_cut_sq)
      dr_cut_sq=huge(dr_cut_sq)
      r_nn=0.0
      
      ! if (.not. mp_del > epsilon(mp_del)) mp_del=0.05
    
      if (first_call) then
        write(*,'(t3,3(a))') 'Using k-nearest mapping without a cut-off radius'
      end if
    
    end if
    
    if (first_call) then
      write(*,'(t6,2(a))')   'Nearest neighbors considered:        ', trim(int2str(nn_n))
      write(*,'(t6,2(a))')   'Nearest neighbors maximum:           ', trim(int2str(nn_nmx))
      write(*,'(t6,a,f5.3)') 'Delta radius factor:                 ', mp_del
      write(*,'(t6,a,f5.3)') 'Power of distance-weight dependence: ', mp_deg 
    end if
  end if
     
  !------------------------
  ! Construct neighbor list
  !------------------------

  ! For each FEM node 
  fem_node_loop: do i=1,nnodes_fem
      
    if (.not. fem_mask(i)) cycle fem_node_loop

    ! Calculate distance (squared) to FDS nodes
    distance=huge(r_sq)
    fds_node_loop: do j=1,nnodes_fds

      if (.not. fds_mask(j)) cycle fds_node_loop

      dx=fem_xyz(i,1)-fds_xyz(j,1)
      dy=fem_xyz(i,2)-fds_xyz(j,2)
      dz=fem_xyz(i,3)-fds_xyz(j,3)
      r_sq=dx*dx+dy*dy+dz*dz
      distance(j)=r_sq

    end do fds_node_loop

    dmask=.true.
    
    nn_r_cut=0
    nn_dr_cut=0

    ! Find nearest neighbors, at most nn_nmx points.
    do j=1,nn_nmx
      ! Count the numbers of points within radius and radius+dr
      mindist=minval(distance,1,dmask)
      imindist=minloc(distance,1,dmask) 
      neighbor(i,j)=imindist
      dmask(imindist)=.false.
      if (distance(imindist) <=  r_cut_sq) nn_r_cut=nn_r_cut+1
      if (distance(imindist) <= dr_cut_sq) nn_dr_cut=nn_dr_cut+1
    end do

    if (mp_n < 1 .and. nn_r_cut == 0) then
      write(*,'(3(a))') 'FEM node ', trim(int2str(i)), &
        ': no FDS nodes inside cut-off radius'
      stop
    end if

    ! Determine radius for neighbor selection
    nn_r_cut=max(1,nn_r_cut)
    nn_dr_cut=max(1,nn_dr_cut)
    if (nn_r_cut > nn_n) then
       ! more than mp_n nn points found. Set the r limit as the distance of the mp_n.
       r_nn=sqrt(distance(neighbor(i,nn_n)))
    else
       ! mp_n or less nn points found. Set the r limit as the last point found 
       ! if no cut-off radius is given. Otherwise use the given cut-off radius.
       r_nn=max(r_nn,sqrt(distance(neighbor(i,nn_r_cut))))
    end if

    w_sum=0.0 ! Sum of weights
    count_new=nn_dr_cut
    if (abs(mp_deg) > epsilon(mp_deg)) then
      ! Non-uniform weight function
       inter_loop: do j=1,max(nn_r_cut,nn_dr_cut)
          r_now=sqrt(distance(neighbor(i,j)))
          if (j > min(nn_n,nn_r_cut) .and. r_now > r_nn*(1.0+mp_del)) then
            count_new=count_new-1
            weight_neighbor(i,j)=0.0
            cycle inter_loop
          end if
          weight=1.0/max(r_now**mp_deg,10.0*tiny(r_now))
          weight_neighbor(i,j)=weight
          w_sum=w_sum+weight 
       end do inter_loop
    else
      ! Uniform weight function, i.e., just a plain arithmetic average.
       ave_loop: do j=1,max(nn_r_cut,nn_dr_cut)
          r_now=sqrt(distance(neighbor(i,j)))
          if (j > min(nn_n,nn_r_cut) .and. r_now > r_nn*(1.0+mp_del)) then
            count_new=count_new-1
            weight_neighbor(i,j)=0.0
            cycle ave_loop
          end if
          weight=1.0
          weight_neighbor(i,j)=weight
          w_sum=w_sum+weight
       end do ave_loop
    end if

    if (w_sum < 10.0*tiny(w_sum)) then
       ! No neighbors within the cut-off radius, use the closest point
       weight_neighbor(i,1)=1.0
    else
       ! Normalize the weight function
       weight_neighbor(i,1:nn_nmx)=weight_neighbor(i,1:nn_nmx)/w_sum
    end if

  end do fem_node_loop

  !-----------------------------------------
  ! Assign neighbor average to each FEM node
  !-----------------------------------------

  if (read_hcoeff) then
    fem_hcoeff_node_loop_2: do i=1,nnodes_fem
      if (.not. fem_mask(i)) cycle fem_hcoeff_node_loop_2
      do j=1,nn_nmx
        do k=1,ntimes_fem
          fem_data(k,i)=fem_data(k,i)+fds_data(k,neighbor(i,j))*weight_neighbor(i,j)
          fem_hcoeff(k,i)=fem_hcoeff(k,i)+fds_hcoeff(k,neighbor(i,j))*weight_neighbor(i,j)
        end do
      end do
    end do fem_hcoeff_node_loop_2

  else
    fem_node_loop_2: do i=1,nnodes_fem
      if (.not. fem_mask(i)) cycle fem_node_loop_2
      do j=1,nn_nmx
        do k=1,ntimes_fem
          fem_data(k,i)=fem_data(k,i)+fds_data(k,neighbor(i,j))*weight_neighbor(i,j)
        end do
      end do
    end do fem_node_loop_2
  end if

  deallocate(neighbor,stat=ios);        call error_allocate(ios)
  deallocate(weight_neighbor,stat=ios); call error_allocate(ios)
  deallocate(distance,stat=ios);        call error_allocate(ios)
  deallocate(dmask,stat=ios);           call error_allocate(ios)

  if (first_call) first_call=.false.

end subroutine map_nearest

subroutine map_devc_to_nset()
!----------------------------------
! FDS:   Device to node set mapping
! CFAST: Target to node set mapping
!----------------------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  Use cfast_arrays
  implicit none
 
  integer :: i,j,k,l,m
  Real(kind=rk) :: cfast_eps
  character(len=chr80) :: nset_name

  If (cfast_input) Then
     write(*,'(t3,3(a))') 'Direct Target-NSET mapping'
  Else
     write(*,'(t3,3(a))') 'Direct DEVC-NSET mapping'
  End If

  If (Trim(fds_output) /= 'devc' .And. .Not. cfast_input) Then
    write(*,'(a)') 'ERROR: DEVC output expected'
    stop
  end if

  if (.not. nset_connectivity) then
     If (cfast_input) Then
        write(*,'(a)') 'ERROR: no NSET-Target connectivity given'
     Else
        write(*,'(a)') 'ERROR: no NSET-DEVC connectivity given'
     End If
     Stop
  end if

  !-----------------------------------------
  ! Assign DEVC/Target averages to node sets
  !-----------------------------------------

  fem_data = 0.0
  i_loop: Do i = 1, Ubound(connectivity_table,1)
     nset_name = Trim(connectivity_table(i,1))

     j_loop: Do j = 1, nnodes_fem
        If (Trim(fem_nset(j)) == Trim(nset_name)) Then
           
           m = 0
           cfast_eps = 0.0
           k_loop: Do k = 1, nnodes_fds

              l_loop: Do l = 2, Ubound(connectivity_table,2)
                 If (fds_idevc(k) == connectivity_table_num(i,l)) Then
                    fem_data(1:ntimes_fem,j) = fem_data(1:ntimes_fem,j) + fds_data(1:ntimes_fds,k)
                    cfast_eps = cfast_eps + cfast_target_epsilon(k)
                    m = m + 1
                 End If
              End Do l_loop
              
           End Do k_loop

           if (m /= 0) then
              fem_data(1:ntimes_fem,j) = fem_data(1:ntimes_fem,j)/m
              If (Trim(fem_mode)=='abaqus' .And. cfast_input) Then
                 abaqus_node_emissivity(j) = cfast_eps/m
              End If
           end if
        
        end if

     end do j_loop

  end do i_loop

end subroutine map_devc_to_nset

end module mapping
