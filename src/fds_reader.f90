!-----------------------------------------------------------------
! Functions and subroutines for reading FDS input and output files
!-----------------------------------------------------------------
module fds_reader

contains

subroutine fds_reader_module()
!---------
! Main sub
!---------
  use fds_bndf_arrays
  use fds_devc_arrays
  use global_constants
  use global_variables
  use mapping_arrays
  use miscellaneous
  use string_handling
  implicit none

  integer :: n_devc_files,n_bndf_files

  fds_xyz_available=.false.
  fds_data_available=.false.
  fds_model_available=.false.

  write(*,'(a)') ''
  write(*,'(a)') 'FDS reader module'

  if (len_trim(fds_input_file) == 0) then 

    !------------------------
    ! No FDS input file given
    !------------------------

    write(*,'(t3,a)') 'Nothing to be done'

  else

    !---------------------
    ! Parse FDS input file
    !---------------------

    call parse_fds_head_namelist()
    call parse_fds_mesh_namelist()
    call parse_fds_trnx_namelist()
    call parse_fds_trny_namelist()
    call parse_fds_trnz_namelist()
    call parse_fds_dump_namelist()
    call parse_fds_obst_namelist()
    call parse_fds_prop_namelist()
    call parse_fds_devc_namelist()
    call parse_fds_bndf_namelist()

    write(*,'(t3,3(a))') 'FDS input file ', trim(quote(fds_input_file)), ' parsed'

    ! Relevant variables and arrays
    !---------------------------------------------------------
    !   fds_chid,         fds_title,         fds_fyi,
    !   fds_mesh_xb(:,:), fds_mesh_ijk(:,:), fds_devc_qnty(:),
    !   fds_devc_name(:), fds_devc_xyz(:,:), fds_devc_ior(:),
    !   fds_bndf_qnty(:), fds_prop_name(:),  fds_prop_qnty(:)
    !---------------------------------------------------------

    !------------------------
    ! Parse BNDF or DEVC file
    !------------------------

    select case (trim(fds_output))
    case ('bndf')
      !------------
      ! BNDF output
      !------------
    
      call check_bndf_exists()
      if (read_hcoeff) then
        call check_hcoeff_bndf_exists()
      end if
      call locate_bndf_files()
      call survey_bndf_files()
      call import_bndf_geom()
      call import_bndf_data()
      call filter_bndf_data()

      ! Relevant arrays
      !---------------------
      !   fds_bndf_ior(:),
      !   fds_bndf_time(:), 
      !   fds_bndf_xyz(:,:),
      !   fds_bndf_nb(:),
      !   fds_bndf_nm(:),
      !   fds_bndf_np(:),
      !   fds_bndf_data(:,:)
      !---------------------

      n_bndf_files=ubound(fds_bndf_file,1)
      if (n_bndf_files == 1) then
        write(*,'(t3,3(a))') 'Boundary file ', &
          trim(quote(fds_bndf_file(1,1))), ' parsed'
      else if (n_bndf_files > 1) then
        write(*,'(t3,5(a))') 'Boundary files ', &
          trim(quote(fds_bndf_file(1,1))), '...', &
          trim(quote(fds_bndf_file(n_bndf_files,1))), ' parsed'
      end if

      call deallocate_excess_fds_arrays()

    case ('devc')
      !------------
      ! DEVC output
      !------------

      call locate_devc_files()
      call import_devc_data()
      call filter_devc_data()

      ! Relevant arrays
      !---------------------
      !   fds_devc_ior(:),
      !   fds_devc_time(:), 
      !   fds_devc_xyz(:,:),
      !   fds_devc_data(:,:)
      !---------------------

      n_devc_files=ubound(fds_devc_file,1)
      if (n_devc_files == 1) then
        write(*,'(t3,3(a))') 'Device file ', &
          trim(quote(fds_devc_file(1))), ' parsed'
      else if (n_devc_files > 1) then
        write(*,'(t3,5(a))') 'Device files ', &
          trim(quote(fds_devc_file(1))), '...', &
          trim(quote(fds_devc_file(n_devc_files))), ' parsed'
      end if
      
      call deallocate_excess_fds_arrays()

    end select
    
    nnodes_fds=ubound(fds_xyz,1)
    if (nnodes_fds >= 1) then
      fds_xyz_available=.true.; fds_data_available=.true.
    else
      fds_xyz_available=.false.; fds_data_available=.false.
    end if

    fds_model_available=.true.

    ! At this point we have
    !----------------------
    !   fds_id(:),
    !   fds_ior(:),
    !   fds_time(:),
    !   fds_idevc(:),
    !   fds_patch(:),
    !   fds_xyz(:,:),
    !   fds_data(:,:)
    !----------------------

  end if

end subroutine fds_reader_module

!***************
! Auxiliary subs
!***************

!******************************************
! Subroutines for parsing an FDS input file
!******************************************

subroutine parse_fds_bndf_namelist()
!---------------------------------------------------
! Parse FDS input file for boundary file information
!
! Fills the array
!   fds_bndf_qnty(:)
!---------------------------------------------------
  use error_messages
  use fds_bndf_arrays
  use fds_head_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none
  
  integer :: i,ios,nbndf

  ! BNDF-namelist group
  logical :: cell_centered,recount_drip
  character(len=30) quantity,part_id,prop_id,spec_id
  character(len=100) :: fyi
  
  namelist /bndf/ cell_centered,fyi,part_id,prop_id,quantity,recount_drip,spec_id

  if (trim(lowercase(fds_output)) /= 'bndf') return

  open(unit=iochannel(1),file=trim(fds_input_file),status='old',iostat=ios) 
  if (ios /= 0) call error_open_file(fds_input_file)

  nbndf=0
  count_bndf: do
    read(iochannel(1),nml=bndf,iostat=ios)
    select case(ios)
    case (:-1)
      exit count_bndf
    case (0)
      nbndf=nbndf+1
    case(1:)
      write(*,'(5(a))') 'ERROR: in reading BNDF-namelist record (file ', &
        trim(quote(fds_input_file)), ', BNDF number ', trim(int2str(nbndf+1)), ')'
      call error_read_file(fds_input_file)
    end select
  end do count_bndf

  if (nbndf > 0) then
    rewind(iochannel(1))
    allocate(fds_bndf_qnty(nbndf),stat=ios); call error_allocate(ios)

    do i=1,nbndf
      quantity=''
      read(iochannel(1),nml=bndf,iostat=ios)
      select case(ios)
      case (:-1)
        exit
      case (0)
        fds_bndf_qnty(i)=trim(quantity) 
      case(1:)
        call error_read_file(fds_input_file)
      end select
    end do
  end if

  close(unit=iochannel(1))

  if (nbndf == 0) then
    write(*,'(2(a))') 'ERROR: no BNDF namelist records found in file ', &
      trim(quote(fds_input_file))
    stop
  end if

end subroutine parse_fds_bndf_namelist

subroutine parse_fds_devc_namelist()
!--------------------------------------------
! Parse FDS input file for device information
!
! Fills the arrays
!   devc_qnty(:),devc_name(:),devc_xyz(:,:),
!   devc_ior(:)
!--------------------------------------------
  use error_messages
  use fds_devc_arrays
  use fds_head_arrays
  use fds_prop_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,j,ios,ndevc
  logical :: xyz_given,xb_given,identical

  ! DEVC-namelist group
  integer :: ior,points,trip_direction,velo_index,virtual_index
  logical :: dry,evacuation,hide_coordinates,initial_state,latch,output,relative,&
             time_averaged
  real(kind=rk) :: bypass_flowrate,conversion_factor,delay,depth,flowrate,rotation,&
                   setpoint,smoothing_factor
  real(kind=rk), dimension(3) :: orientation,xyz
  real(kind=rk), dimension(6) :: xb
  character(len=chr30) :: ctrl_id,devc_id,duct_id,id,matl_id,part_id,prop_id,quantity,&
                       spec_id,statistics,surf_id,units,x_id,y_id,z_id
  character(len=chr30), dimension(2) :: node_id
  character(len=chr100) :: fyi
  
  namelist /devc/ bypass_flowrate,conversion_factor,ctrl_id,delay,depth,devc_id,&
                  dry,duct_id,evacuation,flowrate,fyi,hide_coordinates,id,&
                  initial_state,ior,latch,matl_id,node_id,orientation,output,&
                  part_id,points,prop_id,quantity,relative,rotation,setpoint,&
                  smoothing_factor,spec_id,statistics,surf_id,time_averaged,&
                  trip_direction,units,velo_index,virtual_index,xb,xyz,x_id,&
                  y_id,z_id

  if (trim(lowercase(fds_output)) /= 'devc') return

  open(unit=iochannel(1),file=trim(fds_input_file),status='old',iostat=ios) 
  if (ios /= 0) call error_open_file(fds_input_file)

  ! Count devices
  ndevc=0;  points=0
  count_devc: do
    id=''
    read(iochannel(1),nml=devc,iostat=ios)
    select case(ios)
    case (:-1)
      exit count_devc
    case (0)
      if (points /= 0) then
        !+-+-+-+-+-+-+-+-+-+-+-
        ! Restriction 13.6.2012
        !+-+-+-+-+-+-+-+-+-+-+-
        write(*,'(5(a))') 'ERROR: device arrays are not supported (file ', &
          trim(quote(fds_input_file)), ', DEVC number ', trim(int2str(ndevc+1)), ')'
        stop
      end if
      ndevc=ndevc+1
    case(1:)
      if (len_trim(id) == 0) then
        write(*,'(5(a))') 'ERROR: in reading DEVC-namelist record (file ', &
          trim(quote(fds_input_file)), ', DEVC number ', trim(int2str(ndevc+1)), ')'
      else
        write(*,'(9(a))') 'ERROR: in reading DEVC-namelist record (file ', &
          trim(quote(fds_input_file)), ', DEVC number ', trim(int2str(ndevc+1)), &
          ', DEVC ID ', char(39), trim(id), char(39), ')'
      end if
      call error_read_file(fds_input_file)
    end select
  end do count_devc

  ! Allocate and fill device arrays
  if (ndevc > 0) then
    rewind(iochannel(1))
    allocate(fds_devc_qnty(ndevc),stat=ios);  call error_allocate(ios)
    allocate(fds_devc_name(ndevc),stat=ios);  call error_allocate(ios)
    allocate(fds_devc_xyz(ndevc,3),stat=ios); call error_allocate(ios)
    allocate(fds_devc_xb(ndevc,6),stat=ios);  call error_allocate(ios)
    allocate(fds_devc_ior(ndevc),stat=ios);   call error_allocate(ios)

    do i=1,ndevc
      id=''; quantity=''; prop_id=''; ior=0; statistics=''
      xyz=0.0; xyz(1)=-1.0e10; xb=0.0; xb(1)=-1.0e10
      read(iochannel(1),nml=devc,iostat=ios)
      select case(ios)
      case (:-1)
        exit
      case (0)
        fds_devc_qnty(i)=trim(quantity) 
        fds_devc_name(i)=trim(id)
        fds_devc_xyz(i,1:3)=xyz
        fds_devc_xb(i,1:6)=xb
        fds_devc_ior(i)=ior 
      
        if (numerical_20(fds_devc_name(i))) then
          write(*,'(9(a))') 'ERROR: numerical DEVC ID (file ', &
            trim(quote(fds_input_file)), ', DEVC number ', trim(int2str(i)), &
            ', DEVC ID ', char(39), trim(id), char(39), ')'
          stop
        end if

        !---------------------------------
        ! Check that coordinates are given
        !---------------------------------
        xyz_given=.true.
        if (xyz(1) < -1.0e9) then
          xyz_given=.false.
        end if

        xb_given=.true.
        if (xb(1) < -1.0e9) then
          xb_given=.false.
        end if

        if (xyz_given .and. xb_given) then
          if (len_trim(fds_devc_name(i)) /= 0) then
            write(*,'(5(a))') 'ERROR: coordinates for device ', &
              trim(int2str(i)), ' (', trim(id), ') given using both XYZ and XB'
          else
            write(*,'(3(a))') 'ERROR: coordinates for device ', &
              trim(int2str(i)), ' given using both XYZ and XB'
          end if  
          stop 
        end if

        if ((.not. xyz_given) .and. (.not. xb_given)) then
          if (len_trim(fds_devc_name(i)) /= 0) then
            write(*,'(5(a))') 'ERROR: no coordinates given for device ', &
              trim(int2str(i)), ' (', trim(fds_devc_name(i)), ')'
          else
            write(*,'(2(a))') 'ERROR: no coordinates given for device ', &
              trim(int2str(i))
          end if  
          stop 
        end if

        ! Assign the center of an XB-definition as XYZ
        if (xb_given) then
          fds_devc_xyz(i,1)=0.5*(fds_devc_xb(i,1)+fds_devc_xb(i,2))
          fds_devc_xyz(i,2)=0.5*(fds_devc_xb(i,3)+fds_devc_xb(i,4))
          fds_devc_xyz(i,3)=0.5*(fds_devc_xb(i,5)+fds_devc_xb(i,6))
        end if

        ! Search for quantity in the PROP-namelist record
        if (len_trim(prop_id) /= 0) then
          if (size(fds_prop_name,1) > 0) then
            do j=1,size(fds_prop_name,1)
              if (trim(fds_prop_name(j)) == trim(prop_id)) then
                if (len_trim(fds_prop_qnty(j)) /= 0) then
                  fds_devc_qnty(i)=fds_prop_qnty(j)
                end if
              end if
            end do
          end if
        end if

        if (len_trim(fds_devc_qnty(i)) == 0) then
          if (len_trim(fds_devc_name(i)) /= 0) then
            write(*,'(5(a))') 'ERROR: no QUANTITY associated with device ', &
              trim(int2str(i)), ' (', trim(id), ')'
          else
            write(*,'(2(a))') 'ERROR: no QUANTITY associated with device ', &
              trim(int2str(i))
          end if  
          stop
        end if

      case(1:)
        call error_read_file(fds_input_file)
      end select
    end do
  end if

  close(unit=iochannel(1))

  ! Error handling
  if (ndevc == 0) then
    write(*,'(2(a))') 'ERROR: no DEVC namelist records found in file ', &
      trim(quote(fds_input_file))
    stop
  end if

  ! Check for identical device names
  identical=.false.
  do i=1,ndevc-1
    do j=i+1,ndevc
      if (trim(fds_devc_name(i)) == trim(fds_devc_name(j))) then
        identical=.true.
        if (full_output) then
          write(*,'(t3,5(a))') '* DEVCs ', trim(int2str(i)), &
            & ' and ', trim(int2str(j)), ' have the same ID'
        end if
      end if
    end do
  end do

  if (identical) then
    write(*,'(t3,a)') 'WARNING: identical DEVC IDs found. This might cause problems.'
  end if

end subroutine parse_fds_devc_namelist

subroutine parse_fds_dump_namelist()
!----------------------------------------------
! Parse FDS input file for dump configuration
!
! Assigns values to the variables
!   fds_column_dump_limit,fds_devc_column_limit
!----------------------------------------------
  use error_messages
  use fds_dump_arrays
  use fds_head_arrays
  use global_constants  
  use global_variables
  use string_handling
  implicit none

  integer :: ios,ndump

  integer :: ctrl_column_limit,devc_column_limit,sig_figs,sig_figs_exp,maximum_particles,&
             nframes,plot3d_velo_index
  logical :: debug,flush_file_buffers,smoke3d,timing,mass_file,velocity_error_file,write_xyz,status_files,&
             cutcell_data_file,column_dump_limit
  character(len=chr30) :: smoke3d_quantity,smoke3d_spec_id
  character(len=chr256) :: render_file
  character(len=chr30), dimension(3) :: plot3d_part_id,plot3d_quantity,plot3d_spec_id
  real(kind=eb) :: dt_bnde,dt_bndf,dt_ctrl,dt_devc,dt_devc_line,dt_flush,dt_geom,dt_hrr,dt_isof,dt_mass,&
                   dt_part,dt_pl3d,dt_prof,dt_restart,dt_sl3d,dt_slcf,dt_veg,uvw_timer
  
  namelist /dump/ column_dump_limit,ctrl_column_limit,cutcell_data_file,&
                  debug,devc_column_limit,dt_bnde,dt_bndf,dt_ctrl,dt_devc,dt_devc_line,dt_flush,&
                  dt_geom,dt_hrr,dt_isof,dt_mass,dt_part,dt_pl3d,dt_prof,dt_restart,dt_sl3d,dt_slcf,&
                  dt_veg,flush_file_buffers,mass_file,maximum_particles,nframes,plot3d_part_id,&
                  plot3d_quantity,plot3d_spec_id,plot3d_velo_index,render_file,sig_figs,sig_figs_exp,&
                  smoke3d,smoke3d_quantity,smoke3d_spec_id,status_files,timing,uvw_timer,&
                  velocity_error_file,write_xyz

  column_dump_limit=.true.; devc_column_limit=254
  open(unit=iochannel(1),file=trim(fds_input_file),status='old',iostat=ios) 
  if (ios /= 0) call error_open_file(fds_input_file)

  ndump=0
  dump_loop: do
    read(iochannel(1),nml=dump,iostat=ios)
    select case(ios)
    case (:-1)
      exit dump_loop
    case (0)
      fds_column_dump_limit=column_dump_limit
      fds_devc_column_limit=devc_column_limit
      ndump=ndump+1
    case(1:)
      write(*,'(3(a))') 'ERROR: in reading DUMP-namelist record (file ', &
        trim(quote(fds_input_file)), ')' 
      call error_read_file(fds_input_file)
    end select
  end do dump_loop

  close(unit=iochannel(1))

end subroutine parse_fds_dump_namelist

subroutine parse_fds_head_namelist()
!---------------------------------------------
! Parse FDS input file for CHID, TITLE and FYI
!
! Assigns values to the variables
!   fds_chid,fds_title,fds_fyi
!---------------------------------------------
  use error_messages
  use fds_head_arrays
  use global_constants  
  use global_variables
  use string_handling
  implicit none

  integer :: ios,nhead

  ! HEAD-namelist group
  character(len=chr40) :: chid,title,fyi
  
  namelist /head/ chid,title,fyi

  chid=''; title=''; fyi=''
  open(unit=iochannel(1),file=trim(fds_input_file),status='old',iostat=ios) 
  if (ios /= 0) call error_open_file(fds_input_file)

  ! Read in HEAD-namelist records
  nhead=0
  head_loop: do
    read(iochannel(1),nml=head,iostat=ios)
    select case(ios)
    case (:-1)
      exit head_loop
    case (0)
      fds_chid=trim(chid)
      fds_title=trim(title)
      fds_fyi=trim(fyi)
      nhead=nhead+1
    case(1:)
      write(*,'(3(a))') 'ERROR: in reading HEAD-namelist record (file ', &
        trim(quote(fds_input_file)), ')'
      call error_read_file(fds_input_file)
    end select
  end do head_loop

  close(unit=iochannel(1))

  if (nhead == 0) then
    write(*,'(2(a))') 'ERROR: no HEAD-namelist record found in file ', &
      trim(quote(fds_input_file))
    stop
  end if

  if (nhead > 1) then
    write(*,'(2(a))') 'WARNING: multiple HEAD-namelist records found in file ', &
      trim(quote(fds_input_file))
  end if

  if (len_trim(fds_chid) == 0) then
    write(*,'(2(a))') 'ERROR: no CHID string found in the HEAD-namelist record in file ', &
      trim(quote(fds_input_file))
    stop
  end if


end subroutine parse_fds_head_namelist

subroutine parse_fds_mesh_namelist
!------------------------------------------
! Parse FDS input file for mesh information
!
! Fills the arrays
!   mesh_xb(:,:),mesh_ijk(:,:)
!------------------------------------------
  use error_messages
  use fds_head_arrays
  use fds_mesh_arrays
  use global_constants  
  use global_variables
  use string_handling
  implicit none

  integer :: i,ios,nmesh
 
  ! MESH-namelist group
  integer :: level,mpi_process
  integer, dimension(3) :: ijk,rgb
  logical :: cylindrical,evacuation,evac_humans,synchronize
  character(len=chr25) :: color
  character(len=chr30) :: id,mult_id
  character(len=chr100) :: fyi
  real(kind=rk) :: evac_z_offset
  real(kind=rk), dimension(6) :: xb
  
  namelist /mesh/ color,cylindrical,evacuation,evac_humans,&
                  evac_z_offset,fyi,id,ijk,level,mpi_process,&
                  mult_id,rgb,synchronize,xb

  ijk=0; xb=0.0
  open(unit=iochannel(1),file=trim(fds_input_file),status='old',iostat=ios) 
  if (ios /= 0) call error_open_file(fds_input_file)

  ! How many meshes?
  nmesh=0
  mesh_loop: do
    id=''
    read(iochannel(1),nml=mesh,iostat=ios)
    select case(ios)
    case (:-1)
      exit mesh_loop
    case (0)
      nmesh=nmesh+1
    case(1:)
      if (len_trim(id) == 0) then
        write(*,'(5(a))') 'ERROR: in reading MESH-namelist record (file ', &
          trim(quote(fds_input_file)), ', MESH number ', trim(int2str(nmesh+1)), ')'
      else
        write(*,'(9(a))') 'ERROR: in reading MESH-namelist record (file ', &
          trim(quote(fds_input_file)), ', MESH number ', trim(int2str(nmesh+1)), &
          ', MESH ID ', char(39), trim(id), char(39), ')'
      end if
      call error_read_file(fds_input_file)
    end select
  end do mesh_loop

  !-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  ! Restriction 3.9.2012
  !  Use of multiple meshes is not supported.
  !  This is related to NSET-DEVC/BNDF connectivity.
  !-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

  if (nmesh > 1) then
    write(*,'(a)') 'ERROR: multiple-mesh simulations are not supported'
    !stop
  end if

  rewind(iochannel(1))

  allocate(fds_mesh_ijk(nmesh,3),stat=ios); call error_allocate(ios)
  allocate(fds_mesh_xb(nmesh,6),stat=ios);  call error_allocate(ios)
  fds_mesh_ijk=0; fds_mesh_xb=0

  ! Read in MESH-namelist records
  do i=1,nmesh
    fds_mesh_ijk(i,1:3)=0; fds_mesh_xb(i,1:6)=0.0
    read(iochannel(1),nml=mesh,iostat=ios)
    select case(ios)
    case (:-1)
      exit
    case (0)
      fds_mesh_ijk(i,1:3)=ijk(1:3)
      fds_mesh_xb(i,1:6)=xb(1:6)
    case(1:)
      call error_read_file(fds_input_file)
    end select
  end do

  close(unit=iochannel(1))

  if (nmesh == 0) then
    write(*,'(2(a))') 'ERROR: no MESH-namelist record found in file ', &
      trim(quote(fds_input_file))
    stop
  end if

end subroutine parse_fds_mesh_namelist

subroutine parse_fds_obst_namelist()
!-------------------------------------------------
! Parse FDS input file for obstruction definitions
!
! Fills the arrays
!   mesh_xb(:,:),mesh_ijk(:,:)
!-------------------------------------------------
  use error_messages
  use fds_obst_arrays
  use fds_head_arrays
  use global_constants  
  use global_variables
  use string_handling
  implicit none

  integer :: i,ios,nobst
  
  ! MESH-namelist group
  integer, dimension(3) :: rgb
            
  logical :: sawtooth,thicken,permit_hole,allow_vent,evacuation, removable,bndf_obst,outline,noterrain
  logical, dimension(-3:3) :: bndf_face
             
  character(len=chr25) :: color
  character(len=chr30) :: id,devc_id,prop_id,surf_id,ctrl_id,mult_id
  character(len=chr60) :: mesh_id 
  character(len=chr100) :: fyi

  character(len=chr30), dimension(3) :: surf_ids
  character(len=chr30), dimension(6) :: surf_id6

  real(kind=eb) :: transparency,bulk_density
  real(kind=eb), dimension(3) :: texture_origin
  real(kind=eb), dimension(6) :: xb
  
  namelist /obst/ allow_vent,bndf_face,bndf_obst,bulk_density,color,ctrl_id,devc_id,&
                  evacuation,fyi,id,mesh_id,mult_id,noterrain,outline,permit_hole,prop_id,&
                  removable,rgb,sawtooth,surf_id,surf_id6,surf_ids,texture_origin,thicken,&
                  transparency,xb


  open(unit=iochannel(1),file=trim(fds_input_file),status='old',iostat=ios) 
  if (ios /= 0) call error_open_file(fds_input_file)

  ! How many obstructions?
  nobst=0
  obst_loop: do
    id=''
    read(iochannel(1),nml=obst,iostat=ios)
    select case(ios)
    case (:-1)
      exit obst_loop
    case (0)
      nobst=nobst+1
    case(1:)
      if (len_trim(id) == 0) then
        write(*,'(5(a))') 'ERROR: in reading OBST-namelist record (file ', &
          trim(quote(fds_input_file)), ', OBST number ', trim(int2str(nobst+1)), ')'
      else
        write(*,'(9(a))') 'ERROR: in reading OBST-namelist record (file ', &
          trim(quote(fds_input_file)), ', OBST number ', trim(int2str(nobst+1)), &
          ', OBST ID ', char(39), trim(id), char(39), ')'
      end if
      call error_read_file(fds_input_file)
    end select
  end do obst_loop

  rewind(iochannel(1))

  allocate(fds_obst_xb(nobst,6),stat=ios);  call error_allocate(ios)

  ! Read in OBST-namelist records
  do i=1,nobst
    fds_obst_xb(i,1:6)=0.0
    read(iochannel(1),nml=obst,iostat=ios)
    select case(ios)
    case (:-1)
      exit
    case (0)
      fds_obst_xb(i,1:6)=xb(1:6)
    case(1:)
      call error_read_file(fds_input_file)
    end select
  end do

  close(unit=iochannel(1))

  !if (nobst == 0) then
  !  write(*,'(2(a))') 'ERROR: no OBST-namelist record found in file ', &
  !    trim(quote(fds_input_file))
  !  stop
  !end if

end subroutine parse_fds_obst_namelist

subroutine parse_fds_prop_namelist()
!----------------------------------------
! Parse PROP-namelist record for QUANTITY
!
! Fills the arrays
!   prop_name(:),prop_qnty(:)
!----------------------------------------
  use error_messages
  use fds_head_arrays
  use fds_prop_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none
  
  integer :: i,ios,nprop

  ! PROP-namelist group
  integer :: pdpa_m,pdpa_n,particles_per_second,velocity_component,pdpa_histogram_nbins         
  logical :: pdpa_integrate,pdpa_normalize,pdpa_histogram
  character(len=chr30) :: quantity,part_id,flow_ramp,spray_pattern_table,spec_id,id, &
                          pressure_ramp,spray_pattern_shape
  character(len=chr30), dimension(20) :: smokeview_id,smokeview_parameters
  real(eb) :: activation_obscuration,activation_temperature,alpha_c,alpha_e,beta_c,beta_e, &
              bead_diameter,bead_emissivity,bead_specific_heat,bead_density,bead_h_fixed, &
              conduit_diameter,conduit_thickness,cable_mass_per_length,cable_diameter, &
              cable_jacket_thickness,cable_failure_temperature,c_factor,characteristic_velocity, &
              orifice_diameter,particle_velocity,flow_rate,flow_tau,gauge_temperature, &
              initial_temperature,k_factor,length,offset,operating_pressure,rti,pdpa_start, &
              pdpa_end,pdpa_radius,spray_pattern_mu,spray_pattern_beta,p0
  real(kind=eb), dimension(2) :: pdpa_histogram_limits
  real(kind=eb), dimension(3) :: px
  real(kind=eb), dimension(2,2) :: spray_angle
  real(kind=eb), dimension(3,3) :: pxx

  namelist /prop/ activation_obscuration,activation_temperature,alpha_c,alpha_e, &
                  bead_density,bead_diameter,bead_emissivity,bead_h_fixed,bead_specific_heat, &
                  beta_c,beta_e,cable_diameter,cable_failure_temperature,cable_jacket_thickness, &
                  cable_mass_per_length,characteristic_velocity,conduit_diameter,conduit_thickness, &
                  c_factor,particles_per_second,particle_velocity,flow_ramp,flow_rate,flow_tau, &
                  gauge_temperature,id,initial_temperature,k_factor,length,offset,operating_pressure, &
                  orifice_diameter,p0,part_id,pdpa_end,pdpa_histogram,pdpa_histogram_limits,pdpa_histogram_nbins, &
                  pdpa_integrate,pdpa_m,pdpa_n,pdpa_normalize,pdpa_radius,pdpa_start,pressure_ramp,px,pxx,quantity, &
                  rti,smokeview_id,smokeview_parameters,spec_id,spray_angle,spray_pattern_beta,spray_pattern_mu, &
                  spray_pattern_shape,spray_pattern_table,velocity_component

  open(unit=iochannel(1),file=trim(fds_input_file),status='old',iostat=ios) 
  if (ios /= 0) call error_open_file(fds_input_file)

  nprop=0
  count_prop: do
    id=''
    read(iochannel(1),nml=prop,iostat=ios)
    select case(ios)
    case (:-1)
      exit count_prop
    case (0)
      nprop=nprop+1
    case(1:)
      ! What if the error doesn't come from a PROP-namelist?
      if (len_trim(id) == 0) then
        write(*,'(5(a))') 'ERROR: in reading PROP-namelist record (file ', &
          trim(quote(fds_input_file)), ', PROP number ', trim(int2str(nprop+1)), ')'
      else
        write(*,'(9(a))') 'ERROR: in reading PROP-namelist record (file ', &
          trim(quote(fds_input_file)), ', PROP number ', trim(int2str(nprop+1)), &
          ', PROP ID ', char(39), trim(id), char(39), ')'
      end if
      call error_read_file(fds_input_file)
    end select
  end do count_prop

  ! Allocate and fill property arrays
  if (nprop > 0) then
    rewind(iochannel(1))
    allocate(fds_prop_name(nprop),stat=ios); call error_allocate(ios)
    allocate(fds_prop_qnty(nprop),stat=ios); call error_allocate(ios)

    do i=1,nprop
      id=''; quantity=''
      read(iochannel(1),nml=prop,iostat=ios)
      select case(ios)
      case (:-1)
        exit
      case (0)
        fds_prop_name(i)=trim(id)
        fds_prop_qnty(i)=trim(quantity)
      case(1:)
        call error_read_file(fds_input_file)
      end select
    end do
  end if

  close(unit=iochannel(1))

end subroutine parse_fds_prop_namelist

subroutine parse_fds_trnx_namelist()
!--------------------------------------------------------
! Parse FDS input file for mesh stretching in x-direction
!--------------------------------------------------------
  use error_messages
  use global_constants  
  use global_variables
  use string_handling
  implicit none

  integer :: ios,ntrnx

  ! TRNX-namelist group
  integer :: ideriv,mesh_number
  character(len=chr100) :: fyi
  real(kind=eb) :: pc,cc
  
  namelist /trnx/ cc,fyi,ideriv,mesh_number,pc

  ideriv=0; mesh_number=0; pc=0.0; cc=0.0; fyi=''
  open(unit=iochannel(1),file=trim(fds_input_file),status='old',iostat=ios) 
  if (ios /= 0) call error_open_file(fds_input_file)

  ! Read in TRNX-namelist records
  ntrnx=0
  trnx_loop: do
    read(iochannel(1),nml=trnx,iostat=ios)
    select case(ios)
    case (:-1)
      exit trnx_loop
    case (0)
      ntrnx=ntrnx+1
    case(1:)
      write(*,'(3(a))') 'ERROR: in reading TRNX-namelist record (file ', &
        trim(quote(fds_input_file)), ')'
      call error_read_file(fds_input_file)
    end select
  end do trnx_loop

  close(unit=iochannel(1))

  if (ntrnx /= 0) then
    write(*,'(3(a))') 'ERROR: use of TRNX-namelist record is not supported (file ', &
        trim(quote(fds_input_file)), ')'
    stop
  end if

end subroutine parse_fds_trnx_namelist

subroutine parse_fds_trny_namelist()
!--------------------------------------------------------
! Parse FDS input file for mesh stretching in y-direction
!--------------------------------------------------------
  use error_messages
  use global_constants  
  use global_variables
  use string_handling
  implicit none

  integer :: ios,ntrny

  ! TRNY-namelist group
  integer :: ideriv,mesh_number
  character(len=chr100) :: fyi
  real(kind=eb) :: pc,cc
  
  namelist /trny/ cc,fyi,ideriv,mesh_number,pc

  ideriv=0; mesh_number=0; pc=0.0; cc=0.0; fyi=''
  open(unit=iochannel(1),file=trim(fds_input_file),status='old',iostat=ios) 
  if (ios /= 0) call error_open_file(fds_input_file)

  ! Read in TRNY-namelist records
  ntrny=0
  trny_loop: do
    read(iochannel(1),nml=trny,iostat=ios)
    select case(ios)
    case (:-1)
      exit trny_loop
    case (0)
      ntrny=ntrny+1
    case(1:)
      write(*,'(3(a))') 'ERROR: in reading TRNY-namelist record (file ', &
        trim(quote(fds_input_file)), ')'
      call error_read_file(fds_input_file)
    end select
  end do trny_loop

  close(unit=iochannel(1))

  if (ntrny /= 0) then
    write(*,'(3(a))') 'ERROR: use of TRNY-namelist record is not supported (file ', &
        trim(quote(fds_input_file)), ')'
    stop
  end if

end subroutine parse_fds_trny_namelist

subroutine parse_fds_trnz_namelist()
!--------------------------------------------------------
! Parse FDS input file for mesh stretching in z-direction
!--------------------------------------------------------
  use error_messages
  use global_constants  
  use global_variables
  use string_handling
  implicit none

  integer :: ios,ntrnz

  ! TRNZ-namelist group
  integer :: ideriv,mesh_number
  character(len=chr100) :: fyi
  real(kind=eb) :: pc,cc
  
  namelist /trnz/ cc,fyi,ideriv,mesh_number,pc

  ideriv=0; mesh_number=0; pc=0.0; cc=0.0; fyi=''
  open(unit=iochannel(1),file=trim(fds_input_file),status='old',iostat=ios) 
  if (ios /= 0) call error_open_file(fds_input_file)

  ! Read in TRNZ-namelist records
  ntrnz=0
  trnz_loop: do
    read(iochannel(1),nml=trnz,iostat=ios)
    select case(ios)
    case (:-1)
      exit trnz_loop
    case (0)
      ntrnz=ntrnz+1
    case(1:)
      write(*,'(3(a))') 'ERROR: in reading TRNZ-namelist record (file ', &
        trim(quote(fds_input_file)), ')'
      call error_read_file(fds_input_file)
    end select
  end do trnz_loop

  close(unit=iochannel(1))

  if (ntrnz /= 0) then
    write(*,'(3(a))') 'ERROR: use of TRNZ-namelist record is not supported (file ', &
        trim(quote(fds_input_file)), ')'
    stop
  end if

end subroutine parse_fds_trnz_namelist

!***********************************
! Subroutines for parsing DEVC files
!***********************************

subroutine locate_devc_files()
!----------------------
! Locate FDS DEVC-files
!----------------------
  use error_messages
  use fds_head_arrays
  use fds_devc_arrays
  use fds_dump_arrays
  use global_constants
  use global_variables
  implicit none

  integer :: i,n,ios,ndevc,n_devc_files
  logical :: multiple_files
  character(len=chr80) :: filename

  ! Expectations
  multiple_files=.true.
  if (.not. fds_column_dump_limit) then
    multiple_files=.false.
    n_devc_files=1
  else
    ndevc=ubound(fds_devc_qnty,1)
    if (ndevc <= fds_devc_column_limit) then
      multiple_files=.false.
      n_devc_files=1
    end if
    if (mod(ndevc,fds_devc_column_limit) /= 0) then
      n_devc_files=int(ndevc/fds_devc_column_limit)+1
    else
      n_devc_files=int(ndevc/fds_devc_column_limit)
    end if
  end if

  ! What is found?
  n=0
  if (.not. multiple_files) then
    ! A single file expected
    filename=find_single_devc_file(fds_chid,0)
    if (len_trim(filename) /= 0) then
      n=1; allocate(fds_devc_file(n),stat=ios)
      call error_allocate(ios)
      fds_devc_file(1)=trim(filename)
    end if 
  else
    ! Multiple files expected
    n=1
    loop: do
      filename=find_single_devc_file(fds_chid,n)
      if (len_trim(filename) == 0) exit loop
      n=n+1 
    end do loop
    n=n-1
    if (n /= 0) then
      allocate(fds_devc_file(n),stat=ios)
      call error_allocate(ios)
      do i=1,n
        fds_devc_file(i)=find_single_devc_file(fds_chid,i)
      end do
    end if
  end if

  ! Exception handling
  if (.not. multiple_files) then
    if (n == 0) then
      write(*,'(a)') 'ERROR: FDS DEVC-file not found'
      stop
    end if
  else
    if (n == 0) then
      write(*,'(a)') 'ERROR: FDS DEVC-files not found'
      stop
    else if (n < n_devc_files) then
      write(*,*) n,n_devc_files
      write(*,'(a)') 'ERROR: FDS DEVC-file(s) missing'
      stop
    end if
  end if

end subroutine locate_devc_files

function find_single_devc_file(chid,n) result(filename)
!----------------------------
! Find a single FDS DEVC-file
!----------------------------
  use global_constants
  use string_handling
  implicit none

  integer :: n
  logical :: file_exists
  character(len=chr40) :: chid
  character(len=chr80) :: filename

  if (n <= 0) then
    filename=trim(chid) // '_devc.csv'
  else if (n >= 1) then
    filename=trim(chid) // '_' // trim(int2str(n)) // '_devc.csv'
  end if

  inquire(file=filename,exist=file_exists)
  if (.not. file_exists) filename=''

end function find_single_devc_file

function number_csv_rows(filename) result(rows)
!------------------------------------------
! Find out the number of rows in a CSV-file
!------------------------------------------
  use error_messages
  use global_constants
  implicit none

  integer :: ios,rows
  character(len=chr80) :: filename
  character(len=input_line_length) :: input_line

  open(unit=iochannel(1),file=trim(filename),status='old',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  ! How many rows? 
  rows=0
  rows_loop: do  
    read(iochannel(1),*,iostat=ios) input_line
    select case(ios)
    case (:-1)
      exit rows_loop
    case (0)
      rows=rows+1
    case (1:)
      call error_read_file(filename)
    end select
  end do rows_loop

  close(unit=iochannel(1))

end function number_csv_rows

function number_csv_cols(filename) result(cols)
!---------------------------------------------
! Find out the number of columns in a CSV-file
!
! 26.1.2012 Seems to work correctly
!
! WARNING: If there is only one data column in
! the CSV-file and no commas, the number of
! columns is interpreted as zero.
!
! Also, the number of columns is counted based
! on the first line of the CSV-file.
!---------------------------------------------
  use error_messages
  use global_constants
  implicit none

  integer :: i,ios,cols
  character :: ctmp
  character(len=chr80) :: filename
  character(len=input_line_length) :: input_line

  open(unit=iochannel(1),file=trim(filename),status='old',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  cols=0
  read(iochannel(1),'(a)') input_line
  If (len_trim(input_line) > input_line_length - 1) Then
     Write(*,fmt='(a,a)')'ERROR: Too long line(s) in file ',Trim(filename)
     Stop
  End If
  do i=1,len_trim(input_line)
    read(input_line(i:i),'(a)',iostat=ios) ctmp
    if (ios /= 0) call error_read_file(filename)
    if (ctmp == ',') then
      cols=cols+1
    end if
  end do
  cols=cols+1

  close(unit=iochannel(1))

end function number_csv_cols

subroutine import_devc_data()
!----------------------------------
! Import data from DEVC file(s)
!
! 26.1.2012 Seems to work correctly
!
! WARNING: Device IDs are currently
! not read. They are taken from the
! FDS input file instead.
!
! Should be tested thoroughly.
!----------------------------------
  use error_messages
  use fds_devc_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,j,k,ios,column,ndevc,ndevc_files,ibegin,iend

  character :: ctmp
  character(len=input_line_length) :: input_line

  character(len=chr20), dimension(:), allocatable :: devc_name_tmp
  character(len=chr20), dimension(:), allocatable :: devc_unit_tmp

  real(kind=rk), dimension(:,:), allocatable :: devc_data_tmp 

  ! Count rows and columns in DEVC file(s)
  ndevc_files=ubound(fds_devc_file,1)
  allocate(fds_devc_rows(ndevc_files),stat=ios); call error_allocate(ios)
  allocate(fds_devc_cols(ndevc_files),stat=ios); call error_allocate(ios)
  
  do i=1,ndevc_files
    fds_devc_rows(i)=number_csv_rows(fds_devc_file(i))-2
    fds_devc_cols(i)=number_csv_cols(fds_devc_file(i))
  end do
   
  ! ---> Exception handling <---
  if (ndevc_files > 1) then
    do i=1,ndevc_files-1
      if (fds_devc_rows(i) /= fds_devc_rows(i+1)) then
        write(*,'(a)') 'ERROR: inconsistent number of rows in DEVC files'
        stop
      end if
    end do
  end if
  ! ----------------------------

  ! Allocate memory for time series data
  ndevc=sum(fds_devc_cols(1:ndevc_files))-ndevc_files
  allocate(fds_devc_unit(ndevc),stat=ios); call error_allocate(ios)
  allocate(fds_devc_name_b(ndevc),stat=ios); call error_allocate(ios)
  allocate(fds_devc_time(fds_devc_rows(1)),stat=ios); call error_allocate(ios)
  allocate(fds_devc_data(fds_devc_rows(1),ndevc),stat=ios); call error_allocate(ios)

  ! Read in data from one or more DEVC files
  do i=1,ndevc_files
    open(unit=iochannel(1),file=trim(fds_devc_file(i)),status='old',iostat=ios)
    if (ios /= 0) call error_open_file(fds_devc_file(i))
    
    ! Read units
    allocate(devc_unit_tmp(fds_devc_cols(i)),stat=ios); call error_allocate(ios)
    column=1; devc_unit_tmp=''; j=1
    read(iochannel(1),'(a)',iostat=ios) input_line
    if (ios /= 0) call error_read_file(fds_devc_file(i))
    units_loop: do k=1,len_trim(input_line)
      read(input_line(k:k),'(a)',iostat=ios) ctmp
      if (ctmp == ',') then
        column=column+1; j=1
        cycle units_loop
      end if
      if (j <= chr20) then
        devc_unit_tmp(column)(j:j)=input_line(k:k)
      end if
      j=j+1  
    end do units_loop

    ! ---> Exception handling <---
    if (trim(devc_unit_tmp(1)) /= 's') then
      write(*,'(2(a))') 'WARNING: non-standard time units in DEVC file ', &
        trim(quote(fds_devc_file(i)))
    end if
    ! ---------------------------- 

    ! Read IDs
    allocate(devc_name_tmp(fds_devc_cols(i)),stat=ios); call error_allocate(ios)
    column=1; devc_name_tmp=''; j=1
    read(iochannel(1),'(a)',iostat=ios) input_line
    if (ios /= 0) call error_read_file(fds_devc_file(i))
    ids_loop: do k=1,len_trim(input_line)
      read(input_line(k:k),'(a)',iostat=ios) ctmp
      if (ctmp == ',') then
        column=column+1; j=1
        cycle ids_loop
      end if
      ! Omit quotation marks
      if (ctmp /= char(34) .and. ctmp /= char(39)) then
        if (j <= chr20) then
          devc_name_tmp(column)(j:j)=input_line(k:k)
        end if
        j=j+1  
      end if
    end do ids_loop

    ! Time series
    allocate(devc_data_tmp(fds_devc_rows(i),fds_devc_cols(i)),stat=ios); call error_allocate(ios)
    do j=1,fds_devc_rows(i)
      read(iochannel(1),*) (devc_data_tmp(j,k),k=1,fds_devc_cols(i))
    end do

    if (i == 1) then
      fds_devc_time=devc_data_tmp(1:fds_devc_rows(1),1)
      fds_devc_unit(1:fds_devc_cols(1)-1)=devc_unit_tmp(2:fds_devc_cols(1))
      fds_devc_name_b(1:fds_devc_cols(1)-1)=devc_name_tmp(2:fds_devc_cols(1))
      fds_devc_data(1:fds_devc_rows(1),1:fds_devc_cols(1)-1)=devc_data_tmp(1:fds_devc_rows(1),2:fds_devc_cols(1))
    else
      ! Omit the time vector
      ibegin=2-i+sum(fds_devc_cols(1:i-1)); iend=ibegin+fds_devc_cols(i)-2
      fds_devc_unit(ibegin:iend)=devc_unit_tmp(2:fds_devc_cols(i))
      fds_devc_name_b(ibegin:iend)=devc_name_tmp(2:fds_devc_cols(i))
      fds_devc_data(1:fds_devc_rows(i),ibegin:iend)=devc_data_tmp(1:fds_devc_rows(i),2:fds_devc_cols(i))
    end if

    deallocate(devc_unit_tmp,stat=ios); call error_allocate(ios)
    deallocate(devc_name_tmp,stat=ios); call error_allocate(ios)
    deallocate(devc_data_tmp,stat=ios); call error_allocate(ios)

    close(unit=iochannel(1))
  end do

end subroutine import_devc_data

subroutine filter_devc_data
!---------------------------
! Filter unwanted DEVC data
!---------------------------
  use error_messages
  use fds_devc_arrays
  use fds_head_arrays
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,idevc,ndevc,inode,nrows,ncols
  integer :: ndevc_qnty,ndevc_user,ndevc_both

  integer :: devc_number
  logical :: devc_found

  logical, dimension(:), allocatable :: devc_mask_qnty
  logical, dimension(:), allocatable :: devc_mask_user
  logical, dimension(:), allocatable :: devc_mask_both

  character(len=chr80) :: fds_transfer_quantity

  !---------------------------
  ! Assign a transfer quantity
  !---------------------------

  if (trim(transfer_quantity) == 'wall_temperature') then
    fds_transfer_quantity='wall temperature'
  else if (trim(transfer_quantity) == 'net_heat_flux') then
    fds_transfer_quantity='net heat flux'
  else if (trim(transfer_quantity) == 'adiabatic_surface_temperature') then
    fds_transfer_quantity='adiabatic surface temperature'
  end if

  !-----------------------------------------------
  ! Convert connectivity table into numerical form
  !-----------------------------------------------

  if (nset_connectivity) then

    nrows=ubound(connectivity_table,1); ncols=ubound(connectivity_table,2)
    allocate(connectivity_table_num(nrows,ncols),stat=ios); call error_allocate(ios)

    connectivity_table_num=0 
    row_loop_0: do i=1,ubound(connectivity_table,1)
      col_loop_0: do j=2,ubound(connectivity_table,2)
        if (numerical(connectivity_table(i,j))) then
          devc_number=str2int(connectivity_table(i,j))
          connectivity_table_num(i,j)=devc_number

          ! Exception handling
          if (devc_number > ubound(fds_devc_name,1) .or. devc_number < 1) then
            write(*,'(3(a))') 'ERROR: in NSET-DEVC connectivity table: DEVC ', &
              trim(int2str(devc_number)), ' does not exist'
            stop
          end if
        
        else
          devc_found=.false.
          devc_loop_0: do idevc=1,ubound(fds_devc_xyz,1)
            if (trim(connectivity_table(i,j)) == trim(fds_devc_name(idevc))) then
              connectivity_table_num(i,j)=idevc
              devc_found=.true.
            end if
          end do devc_loop_0

          ! Exception handling
          if (.not. devc_found) then
            write(*,'(3(a))') 'ERROR: in NSET-DEVC connectivity table: DEVC ', &
              trim(quote(connectivity_table(i,j))), ' does not exist'
            stop
          end if

        end if

      end do col_loop_0
    end do row_loop_0

  end if

  !------------------------------------
  ! Generate a list of selected devices
  !------------------------------------

  ndevc=ubound(fds_devc_name,1)
  allocate(devc_mask_qnty(ndevc),stat=ios); call error_allocate(ios)
  allocate(devc_mask_user(ndevc),stat=ios); call error_allocate(ios)
  allocate(devc_mask_both(ndevc),stat=ios); call error_allocate(ios)

  ! Based on measured quantity
  ndevc_qnty=0; devc_mask_qnty=.false.
  devc_loop_1: do idevc=1,ndevc
    if (trim(lowercase_30(fds_devc_qnty(idevc))) == trim(fds_transfer_quantity)) then
      devc_mask_qnty(idevc)=.true.
      ndevc_qnty=ndevc_qnty+1
    end if
  end do devc_loop_1
  
  ! If none are selected
  if (ndevc_qnty == 0) then
    write(*,'(2(a))') 'ERROR: no FDS devices measuring ', trim(quote(uppercase(fds_transfer_quantity)))
    stop
  end if

  if (nset_connectivity) then

    ! Based on user selection
    ndevc_user=0; devc_mask_user=.false.
    devc_loop_2: do idevc=1,ndevc

      row_loop_2: do i=1,ubound(connectivity_table,1)
        col_loop_2: do j=2,ubound(connectivity_table,2)
          if (numerical(connectivity_table(i,j))) then
            if (str2int(connectivity_table(i,j)) == idevc) then
              devc_mask_user(idevc)=.true.
              ndevc_user=ndevc_user+1
            end if
          else
            if (trim(connectivity_table(i,j)) == trim(fds_devc_name(idevc))) then
              devc_mask_user(idevc)=.true.
              ndevc_user=ndevc_user+1
            end if
          end if
        end do col_loop_2
      end do row_loop_2

    end do devc_loop_2

    ! Based on both
    ndevc_both=0; devc_mask_both=.false.
    devc_loop_3: do idevc=1,ndevc
      if (devc_mask_user(idevc)) then
        if (devc_mask_qnty(idevc)) then
          devc_mask_both(idevc)=.true.
          ndevc_both=ndevc_both+1
        else
          write(*,'(5(a))') 'ERROR: Device ', trim(int2str(idevc)), &
            ' (', trim(fds_devc_name(idevc)), ') measures a wrong quantity'
          stop
        end if
      end if 
    end do devc_loop_3

    ! If none are selected (this error message might be unreachable)
    if (ndevc_both == 0) then
      write(*,'(2(a))') 'ERROR: no selected FDS devices measuring ', trim(quote(uppercase(fds_transfer_quantity)))
      stop
    end if

  else

    ! Based on measured quantity
    ndevc_both=0; devc_mask_both=.false.
    devc_loop_4: do idevc=1,ndevc
      if (devc_mask_qnty(idevc)) then 
        devc_mask_both(idevc)=.true.
        ndevc_both=ndevc_both+1
      end if
    end do devc_loop_4

  end if
 
  !-------------------------
  ! Generate FDS data arrays
  !-------------------------

  nnodes_fds=ndevc_both
  ntimes_fds=ubound(fds_devc_time,1)

  allocate(fds_id(nnodes_fds),stat=ios);              call error_allocate(ios)
  allocate(fds_ior(nnodes_fds),stat=ios);             call error_allocate(ios)
  allocate(fds_time(ntimes_fds),stat=ios);            call error_allocate(ios)
  allocate(fds_idevc(nnodes_fds),stat=ios);           call error_allocate(ios)
  allocate(fds_xyz(nnodes_fds,3),stat=ios);           call error_allocate(ios)
  allocate(fds_data(ntimes_fds,nnodes_fds),stat=ios); call error_allocate(ios)

  !-----------------------------------------------------
  ! Create arrays containing device coordinates and data
  !-----------------------------------------------------

  inode=1
  do idevc=1,ndevc
    if (devc_mask_both(idevc)) then
      fds_ior(inode)=fds_devc_ior(idevc)
      fds_id(inode)=fds_devc_name(idevc)
      fds_idevc(inode)=idevc
      fds_xyz(inode,1:3)=fds_devc_xyz(idevc,1:3)
      fds_data(1:ntimes_fds,inode)=fds_devc_data(1:ntimes_fds,idevc)
      inode=inode+1
    end if  
  end do

  fds_time=fds_devc_time

end subroutine filter_devc_data

function devc_within_domain(idevc) result(answer)
!------------------------------------------------------------
! Check that a device is located within the simulation domain
!------------------------------------------------------------
  use error_messages
  use fds_head_arrays
  use fds_mesh_arrays
  use fds_devc_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: idevc,imesh,nmesh
  logical :: answer

  answer=.false.
  nmesh=ubound(fds_mesh_ijk,1)
  mesh_loop: do imesh=1,nmesh
    if (fds_devc_xyz(idevc,1) >= fds_mesh_xb(imesh,1) .and. &
        fds_devc_xyz(idevc,1) <= fds_mesh_xb(imesh,2)) then
      if (fds_devc_xyz(idevc,2) >= fds_mesh_xb(imesh,3) .and. &
          fds_devc_xyz(idevc,2) <= fds_mesh_xb(imesh,4)) then
        if (fds_devc_xyz(idevc,3) >= fds_mesh_xb(imesh,5) .and. &
            fds_devc_xyz(idevc,3) <= fds_mesh_xb(imesh,6)) then
          answer=.true.
        end if
      end if
    end if
    if (answer .eqv. .true.) exit mesh_loop
  end do mesh_loop

  if (.not. answer) then
    if (len_trim(fds_devc_name(idevc)) /= 0) then
      write(*,'(6(a))') 'WARNING: device ', trim(int2str(idevc)), &
        ' (', trim(fds_devc_name(idevc)), ') is not within any mesh'
    else
      write(*,'(3(a))') 'WARNING: device ', trim(int2str(idevc)), &
        ' is not within any mesh'
    end if
  end if

end function devc_within_domain

!***********************************
! Subroutines for parsing BNDF files
!***********************************

subroutine locate_bndf_files()
!----------------------------------------
! Locate FDS BNDF-files
!----------------------------------------
! 26.1.2012 Seems to work correctly
!
! WARNING: Does the file name creator
! work correctly if there are many files?
!----------------------------------------
  use error_messages
  use fds_head_arrays
  use fds_bndf_arrays
  use fds_mesh_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: ios,nbndf,nmesh
  integer :: ibndf,imesh
  character(len=chr80) :: filename
  logical :: file_exists
 
  nbndf=ubound(fds_bndf_qnty,1); nmesh=ubound(fds_mesh_ijk,1)

  allocate(fds_bndf_file(nbndf,nmesh),stat=ios); call error_allocate(ios)

  do ibndf=1,nbndf
    do imesh=1,nmesh
        
      write(filename,'(a,a,g0,a,g0,a)') trim(fds_chid), &
        '_', imesh, '_', ibndf, '.bf'
      
      inquire(file=filename,exist=file_exists)
      if (.not. file_exists) then
        write(*,'(3(a))') 'ERROR: BNDF file ', &
          trim(quote(filename)), ' not found'
        stop
      end if
    
      fds_bndf_file(ibndf,imesh)=trim(filename)
    end do
  end do

end subroutine locate_bndf_files

function bndf_nodes(ibndf,imesh) result(nnode)
!--------------------------------------------
! Find out the number of nodes in a BNDF-file 
!--------------------------------------------
  use error_messages
  use fds_head_arrays
  use fds_bndf_arrays
  use global_constants
  implicit none

  integer :: i,i1,i2,j1,j2,k1,k2,ios
  integer :: ibndf,imesh,npatch,nnode
  character(len=chr30) :: quantity,short_name,units

  open(unit=iochannel(1),file=trim(fds_bndf_file(ibndf,imesh)), &
    status='old',iostat=ios,form='unformatted')
  if (ios /= 0) call error_open_file(fds_bndf_file(ibndf,imesh))

  read(iochannel(1),iostat=ios) quantity
  if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
  read(iochannel(1),iostat=ios) short_name
  if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
  read(iochannel(1),iostat=ios) units
  if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
  read(iochannel(1),iostat=ios) npatch
  if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))

  nnode=0
  do i=1,npatch
    read(iochannel(1),iostat=ios) i1,i2,j1,j2,k1,k2
    if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
    nnode=nnode+(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
  end do

  close(unit=iochannel(1))

end function bndf_nodes

function bndf_patches(ibndf,imesh) result(npatch)
!--------------------------------------------
! Find out the number of nodes in a BNDF-file 
!--------------------------------------------
  use error_messages
  use fds_head_arrays
  use fds_bndf_arrays
  use global_constants
  implicit none

  integer :: ios
  integer :: ibndf,imesh,npatch
  character(len=chr30) :: quantity,short_name,units

  open(unit=iochannel(1),file=trim(fds_bndf_file(ibndf,imesh)), &
    status='old',iostat=ios,form='unformatted')
  if (ios /= 0) call error_open_file(fds_bndf_file(ibndf,imesh))

  read(iochannel(1),iostat=ios) quantity
  if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
  read(iochannel(1),iostat=ios) short_name
  if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
  read(iochannel(1),iostat=ios) units
  if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
  read(iochannel(1),iostat=ios) npatch
  if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))

  close(unit=iochannel(1))

end function bndf_patches

function bndf_elements(ibndf,imesh) result(nelement)
!--------------------------------------------
! Find out the number of nodes in a BNDF-file 
!--------------------------------------------
  use error_messages
  use fds_head_arrays
  use fds_bndf_arrays
  use global_constants
  implicit none

  integer :: i,i1,i2,j1,j2,k1,k2,ios
  integer :: ibndf,imesh,npatch,nelement
  character(len=chr30) :: quantity,short_name,units

  open(unit=iochannel(1),file=trim(fds_bndf_file(ibndf,imesh)), &
    status='old',iostat=ios,form='unformatted')
  if (ios /= 0) call error_open_file(fds_bndf_file(ibndf,imesh))

  read(iochannel(1),iostat=ios) quantity
  if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
  read(iochannel(1),iostat=ios) short_name
  if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
  read(iochannel(1),iostat=ios) units
  if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
  read(iochannel(1),iostat=ios) npatch
  if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))

  nelement=0
  do i=1,npatch
    read(iochannel(1),iostat=ios) i1,i2,j1,j2,k1,k2
    if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
    if (i1 == i2) then
      nelement=nelement+(j2-j1)*(k2-k1)
    else if (j1 == j2) then
      nelement=nelement+(i2-i1)*(k2-k1)
    else if (k1 == k2) then
      nelement=nelement+(i2-i1)*(j2-j1)
    end if
  end do

  close(unit=iochannel(1))

end function bndf_elements

subroutine survey_bndf_files()
!---------------------------------------------------------
! Find out the number of nodes etc. in a set of BNDF-files 
!---------------------------------------------------------
  use error_messages
  use fds_bndf_arrays
  use fds_mesh_arrays
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,i1,i2,j1,j2,k1,k2,ios,nior,nb,nm
  integer :: ibndf,jbndf,imesh,nbndf,nmesh,inode,nsum
  character(len=chr30) :: quantity
  real(kind=fb) :: rtmp

  integer, dimension(:), allocatable :: pat

  nbndf=ubound(fds_bndf_qnty,1); nmesh=ubound(fds_mesh_ijk,1)
  ibndf=selected_bndf()

  allocate(fds_bndf_snam(nmesh),stat=ios); call error_allocate(ios)
  allocate(fds_bndf_unit(nmesh),stat=ios); call error_allocate(ios)
  allocate(fds_bndf_npat(nmesh),stat=ios); call error_allocate(ios)
  allocate(fds_bndf_nnod(nmesh),stat=ios); call error_allocate(ios)
  allocate(fds_bndf_ntim(nmesh),stat=ios); call error_allocate(ios)

  !-------------------------------------------
  ! Boundary file for chosen transfer quantity
  !-------------------------------------------

  mesh_loop: do imesh=1,nmesh

    open(unit=iochannel(1),file=trim(fds_bndf_file(ibndf,imesh)), &
    status='old',iostat=ios,form='unformatted')
    if (ios /= 0) call error_open_file(fds_bndf_file(ibndf,imesh))

    read(iochannel(1),iostat=ios) quantity
    if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
    read(iochannel(1),iostat=ios) fds_bndf_snam(imesh)
    if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
    read(iochannel(1),iostat=ios) fds_bndf_unit(imesh)
    if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
    read(iochannel(1),iostat=ios) fds_bndf_npat(imesh)
    if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))

    deallocate(pat,stat=ios)
    allocate(pat(fds_bndf_npat(imesh)),stat=ios); call error_allocate(ios)

    fds_bndf_nnod(imesh)=0
    patch_loop: do i=1,fds_bndf_npat(imesh)
      read(iochannel(1),iostat=ios) i1,i2,j1,j2,k1,k2,nior,nb,nm
      if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
      fds_bndf_nnod(imesh)=fds_bndf_nnod(imesh)+(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
      pat(i)=(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
    end do patch_loop
    
    nsum=0
    do i=1,fds_bndf_npat(imesh)
      nsum=nsum+pat(i)
    end do

    fds_bndf_ntim(imesh)=0
    time_loop: do
      read(iochannel(1),iostat=ios) rtmp
      select case(ios)
        case (:-1)
          exit time_loop
        case (0)
          fds_bndf_ntim(imesh)=fds_bndf_ntim(imesh)+1
        case (1:)
          call error_read_file(fds_bndf_file(ibndf,imesh))
      end select
      
      node_loop: do i=1,fds_bndf_npat(imesh)
        read(iochannel(1),iostat=ios) (rtmp,inode=1,pat(i))
        if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
      end do node_loop

    end do time_loop

    close(unit=iochannel(1))

  end do mesh_loop

  allocate(fds_bndf_pat(maxval(fds_bndf_npat)))

  !--------------------------------------------
  ! Boundary file for heat transfer coefficient
  !--------------------------------------------

  if (read_hcoeff) then
  
    jbndf=selected_hcoeff_bndf()

    allocate(fds_hcoeff_bndf_snam(nmesh),stat=ios); call error_allocate(ios)
    allocate(fds_hcoeff_bndf_unit(nmesh),stat=ios); call error_allocate(ios)
    allocate(fds_hcoeff_bndf_npat(nmesh),stat=ios); call error_allocate(ios)
    allocate(fds_hcoeff_bndf_nnod(nmesh),stat=ios); call error_allocate(ios)
    allocate(fds_hcoeff_bndf_ntim(nmesh),stat=ios); call error_allocate(ios)

    hcoeff_mesh_loop: do imesh=1,nmesh

      open(unit=iochannel(1),file=trim(fds_bndf_file(jbndf,imesh)), &
      status='old',iostat=ios,form='unformatted')
      if (ios /= 0) call error_open_file(fds_bndf_file(jbndf,imesh))

      read(iochannel(1),iostat=ios) quantity
      if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))
      read(iochannel(1),iostat=ios) fds_hcoeff_bndf_snam(imesh)
      if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))
      read(iochannel(1),iostat=ios) fds_hcoeff_bndf_unit(imesh)
      if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))
      read(iochannel(1),iostat=ios) fds_hcoeff_bndf_npat(imesh)
      if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))

      deallocate(pat,stat=ios)
      allocate(pat(fds_hcoeff_bndf_npat(imesh)),stat=ios); call error_allocate(ios)

      fds_hcoeff_bndf_nnod(imesh)=0
      hcoeff_patch_loop: do i=1,fds_hcoeff_bndf_npat(imesh)
        read(iochannel(1),iostat=ios) i1,i2,j1,j2,k1,k2,nior,nb,nm
        if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))
        fds_hcoeff_bndf_nnod(imesh)=fds_hcoeff_bndf_nnod(imesh)+(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
        pat(i)=(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
      end do hcoeff_patch_loop
      
      nsum=0
      do i=1,fds_hcoeff_bndf_npat(imesh)
        nsum=nsum+pat(i)
      end do

      fds_hcoeff_bndf_ntim(imesh)=0
      hcoeff_time_loop: do
        read(iochannel(1),iostat=ios) rtmp
        select case(ios)
          case (:-1)
            exit hcoeff_time_loop
          case (0)
            fds_hcoeff_bndf_ntim(imesh)=fds_hcoeff_bndf_ntim(imesh)+1
          case (1:)
            call error_read_file(fds_bndf_file(jbndf,imesh))
        end select
        
        hcoeff_node_loop: do i=1,fds_hcoeff_bndf_npat(imesh)
          read(iochannel(1),iostat=ios) (rtmp,inode=1,pat(i))
          if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))
        end do hcoeff_node_loop

      end do hcoeff_time_loop

      close(unit=iochannel(1))

    end do hcoeff_mesh_loop

    allocate(fds_hcoeff_bndf_pat(maxval(fds_bndf_npat)))

  end if

end subroutine survey_bndf_files

subroutine import_bndf_geom()
!------------------------
! Import node coordinates
!------------------------
  use error_messages
  use fds_head_arrays
  use fds_bndf_arrays
  use fds_mesh_arrays
  use global_constants
  use global_variables
  implicit none

  integer :: i,j,k,l,m,i1,i2,j1,j2,k1,k2,ios,nior
  integer :: ibndf,jbndf,imesh,ipatch
  integer :: nmesh,npatch,nnode,nelement,nb,nm

  character(len=chr30) :: quantity,short_name,units
 
  real(kind=rk), dimension(:), allocatable :: x0,x1,dx
  real(kind=rk), dimension(:), allocatable :: y0,y1,dy
  real(kind=rk), dimension(:), allocatable :: z0,z1,dz

  integer, dimension(:,:,:), allocatable :: fds_bndf_xyz_tmp
  integer, dimension(:,:,:), allocatable :: fds_hcoeff_bndf_xyz_tmp

  nmesh=ubound(fds_bndf_file,2)

  !-------------------------
  ! Chosen transfer quantity
  !-------------------------

  ! Total number of nodes, patches and elements
  ibndf=selected_bndf(); nnode=0; nelement=0; npatch=0
  do imesh=1,nmesh
    nnode    = nnode    + bndf_nodes(ibndf,imesh)
    npatch   = npatch   + bndf_patches(ibndf,imesh)
    nelement = nelement + bndf_elements(ibndf,imesh)
  end do

  allocate(fds_bndf_xyz(nnode,3),stat=ios); call error_allocate(ios)
  allocate(fds_bndf_ior(nnode),stat=ios);   call error_allocate(ios)

  allocate(fds_bndf_nb(nnode),stat=ios); call error_allocate(ios)
  allocate(fds_bndf_nm(nnode),stat=ios); call error_allocate(ios)
  allocate(fds_bndf_np(nnode),stat=ios); call error_allocate(ios)

  allocate(x0(nmesh),x1(nmesh),dx(nmesh),stat=ios); call error_allocate(ios)
  allocate(y0(nmesh),y1(nmesh),dy(nmesh),stat=ios); call error_allocate(ios)
  allocate(z0(nmesh),z1(nmesh),dz(nmesh),stat=ios); call error_allocate(ios)

  allocate(fds_bndf_el(nelement,4),stat=ios); call error_allocate(ios)

  ! Discretization of the simulation domain
  do imesh=1,nmesh
    x0(imesh)=fds_mesh_xb(imesh,1); x1(imesh)=fds_mesh_xb(imesh,2)
    y0(imesh)=fds_mesh_xb(imesh,3); y1(imesh)=fds_mesh_xb(imesh,4)
    z0(imesh)=fds_mesh_xb(imesh,5); z1(imesh)=fds_mesh_xb(imesh,6)
    dx(imesh)=(x1(imesh)-x0(imesh))/(fds_mesh_ijk(imesh,1))
    dy(imesh)=(y1(imesh)-y0(imesh))/(fds_mesh_ijk(imesh,2))
    dz(imesh)=(z1(imesh)-z0(imesh))/(fds_mesh_ijk(imesh,3))
  end do

  l=1; m=1
  mesh_loop: do imesh=1,nmesh

    open(unit=iochannel(1),file=trim(fds_bndf_file(ibndf,imesh)), &
      status='old',iostat=ios,form='unformatted')
    if (ios /= 0) call error_open_file(fds_bndf_file(ibndf,imesh))

    read(iochannel(1),iostat=ios) quantity
    if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
    read(iochannel(1),iostat=ios) short_name
    if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
    read(iochannel(1),iostat=ios) units
    if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
    read(iochannel(1),iostat=ios) npatch
    if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))

    patch_loop: do ipatch=1,npatch
      read(iochannel(1),iostat=ios) i1,i2,j1,j2,k1,k2,nior,nb,nm
      if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))

      allocate(fds_bndf_xyz_tmp(i1:i2,j1:j2,k1:k2),stat=ios); call error_allocate(ios)
      fds_bndf_xyz_tmp=0

      do k=k1,k2
        do j=j1,j2
          do i=i1,i2
            ! Node coordinates
            fds_bndf_xyz(l,1)=x0(imesh)+i*dx(imesh)
            fds_bndf_xyz(l,2)=y0(imesh)+j*dy(imesh)
            fds_bndf_xyz(l,3)=z0(imesh)+k*dz(imesh)

            fds_bndf_xyz_tmp(i,j,k)=l

            ! Orientation, boundary number, mesh number and patch number
            fds_bndf_ior(l)=nior
            fds_bndf_nb(l)=nb
            fds_bndf_nm(l)=nm
            fds_bndf_np(l)=ipatch
            
            l=l+1
          end do
        end do
      end do

      ! Elements
      if (i1 == i2) then

        do i=i1,i2
          do j=j1,j2-1
            do k=k1,k2-1
              fds_bndf_el(m,1)=fds_bndf_xyz_tmp(i,j,k)
              fds_bndf_el(m,2)=fds_bndf_xyz_tmp(i,j+1,k)
              fds_bndf_el(m,3)=fds_bndf_xyz_tmp(i,j+1,k+1)
              fds_bndf_el(m,4)=fds_bndf_xyz_tmp(i,j,k+1)
              m=m+1
            end do
          end do
        end do

      else if (j1 == j2) then

        do i=i1,i2-1
          do j=j1,j2
            do k=k1,k2-1
              fds_bndf_el(m,1)=fds_bndf_xyz_tmp(i,j,k)
              fds_bndf_el(m,2)=fds_bndf_xyz_tmp(i+1,j,k)
              fds_bndf_el(m,3)=fds_bndf_xyz_tmp(i+1,j,k+1)
              fds_bndf_el(m,4)=fds_bndf_xyz_tmp(i,j,k+1)
              m=m+1
            end do
          end do
        end do

      else if (k1 == k2) then

        do i=i1,i2-1
          do j=j1,j2-1
            do k=k1,k2
              fds_bndf_el(m,1)=fds_bndf_xyz_tmp(i,j,k)
              fds_bndf_el(m,2)=fds_bndf_xyz_tmp(i+1,j,k)
              fds_bndf_el(m,3)=fds_bndf_xyz_tmp(i+1,j+1,k)
              fds_bndf_el(m,4)=fds_bndf_xyz_tmp(i,j+1,k)
              m=m+1
            end do
          end do
        end do

      end if

      deallocate(fds_bndf_xyz_tmp,stat=ios); call error_allocate(ios)

    end do patch_loop

    close(unit=iochannel(1))
    
  end do mesh_loop

  deallocate(x0,x1,dx,stat=ios); call error_allocate(ios)
  deallocate(y0,y1,dy,stat=ios); call error_allocate(ios)
  deallocate(z0,z1,dz,stat=ios); call error_allocate(ios)

  !--------------------------
  ! Heat transfer coefficient
  !--------------------------

  if (read_hcoeff) then

    ! Total number of nodes, patches and elements
    jbndf=selected_hcoeff_bndf(); nnode=0; nelement=0
    do imesh=1,nmesh
      nnode=nnode+bndf_nodes(jbndf,imesh)
      npatch=npatch+bndf_patches(jbndf,imesh)
      nelement=nelement+bndf_elements(jbndf,imesh)
    end do

    allocate(fds_hcoeff_bndf_xyz(nnode,3),stat=ios); call error_allocate(ios)
    allocate(fds_hcoeff_bndf_ior(nnode),stat=ios);   call error_allocate(ios)

    allocate(fds_hcoeff_bndf_nb(nnode),stat=ios); call error_allocate(ios)
    allocate(fds_hcoeff_bndf_nm(nnode),stat=ios); call error_allocate(ios)
    allocate(fds_hcoeff_bndf_np(nnode),stat=ios); call error_allocate(ios)

    allocate(x0(nmesh),x1(nmesh),dx(nmesh),stat=ios); call error_allocate(ios)
    allocate(y0(nmesh),y1(nmesh),dy(nmesh),stat=ios); call error_allocate(ios)
    allocate(z0(nmesh),z1(nmesh),dz(nmesh),stat=ios); call error_allocate(ios)

    allocate(fds_hcoeff_bndf_el(nelement,4),stat=ios); call error_allocate(ios)

    ! Discretization of the simulation domain
    do imesh=1,nmesh
      x0(imesh)=fds_mesh_xb(imesh,1); x1(imesh)=fds_mesh_xb(imesh,2)
      y0(imesh)=fds_mesh_xb(imesh,3); y1(imesh)=fds_mesh_xb(imesh,4)
      z0(imesh)=fds_mesh_xb(imesh,5); z1(imesh)=fds_mesh_xb(imesh,6)
      dx(imesh)=(x1(imesh)-x0(imesh))/(fds_mesh_ijk(imesh,1))
      dy(imesh)=(y1(imesh)-y0(imesh))/(fds_mesh_ijk(imesh,2))
      dz(imesh)=(z1(imesh)-z0(imesh))/(fds_mesh_ijk(imesh,3))
    end do

    l=1; m=1
    hcoeff_mesh_loop: do imesh=1,nmesh

      open(unit=iochannel(1),file=trim(fds_bndf_file(jbndf,imesh)), &
        status='old',iostat=ios,form='unformatted')
      if (ios /= 0) call error_open_file(fds_bndf_file(jbndf,imesh))

      read(iochannel(1),iostat=ios) quantity
      if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))
      read(iochannel(1),iostat=ios) short_name
      if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))
      read(iochannel(1),iostat=ios) units
      if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))
      read(iochannel(1),iostat=ios) npatch
      if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))

      hcoeff_patch_loop: do ipatch=1,npatch
        read(iochannel(1),iostat=ios) i1,i2,j1,j2,k1,k2,nior,nb,nm
        if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))

        allocate(fds_hcoeff_bndf_xyz_tmp(i1:i2,j1:j2,k1:k2),stat=ios); call error_allocate(ios)

        do k=k1,k2
          do j=j1,j2
            do i=i1,i2
              ! Node coordinates
              fds_hcoeff_bndf_xyz(l,1)=x0(imesh)+i*dx(imesh)
              fds_hcoeff_bndf_xyz(l,2)=y0(imesh)+j*dy(imesh)
              fds_hcoeff_bndf_xyz(l,3)=z0(imesh)+k*dz(imesh)

              fds_hcoeff_bndf_xyz_tmp(i,j,k)=l

              ! Orientation, boundary number, mesh number and patch number
              fds_hcoeff_bndf_ior(l)=nior
              fds_hcoeff_bndf_nb(l)=nb
              fds_hcoeff_bndf_nm(l)=nm
              fds_hcoeff_bndf_np(l)=ipatch
              
              l=l+1
            end do
          end do
        end do

        ! Elements
        if (i1 == i2) then

          do i=i1,i2
            do j=j1,j2-1
              do k=k1,k2-1
                fds_hcoeff_bndf_el(m,1)=fds_hcoeff_bndf_xyz_tmp(i,j,k)
                fds_hcoeff_bndf_el(m,2)=fds_hcoeff_bndf_xyz_tmp(i,j+1,k)
                fds_hcoeff_bndf_el(m,3)=fds_hcoeff_bndf_xyz_tmp(i,j+1,k+1)
                fds_hcoeff_bndf_el(m,4)=fds_hcoeff_bndf_xyz_tmp(i,j,k+1)
                m=m+1
              end do
            end do
          end do

        else if (j1 == j2) then

          do i=i1,i2-1
            do j=j1,j2
              do k=k1,k2-1
                fds_hcoeff_bndf_el(m,1)=fds_hcoeff_bndf_xyz_tmp(i,j,k)
                fds_hcoeff_bndf_el(m,2)=fds_hcoeff_bndf_xyz_tmp(i+1,j,k)
                fds_hcoeff_bndf_el(m,3)=fds_hcoeff_bndf_xyz_tmp(i+1,j,k+1)
                fds_hcoeff_bndf_el(m,4)=fds_hcoeff_bndf_xyz_tmp(i,j,k+1)
                m=m+1
              end do
            end do
          end do

        else if (k1 == k2) then

          do i=i1,i2-1
            do j=j1,j2-1
              do k=k1,k2
                fds_hcoeff_bndf_el(m,1)=fds_hcoeff_bndf_xyz_tmp(i,j,k)
                fds_hcoeff_bndf_el(m,2)=fds_hcoeff_bndf_xyz_tmp(i+1,j,k)
                fds_hcoeff_bndf_el(m,3)=fds_hcoeff_bndf_xyz_tmp(i+1,j+1,k)
                fds_hcoeff_bndf_el(m,4)=fds_hcoeff_bndf_xyz_tmp(i,j+1,k)
                m=m+1
              end do
            end do
          end do

        end if

        deallocate(fds_hcoeff_bndf_xyz_tmp,stat=ios); call error_allocate(ios)

      end do hcoeff_patch_loop

      close(unit=iochannel(1))
      
    end do hcoeff_mesh_loop

    deallocate(x0,x1,dx,stat=ios); call error_allocate(ios)
    deallocate(y0,y1,dy,stat=ios); call error_allocate(ios)
    deallocate(z0,z1,dz,stat=ios); call error_allocate(ios)

  end if

end subroutine import_bndf_geom

subroutine import_bndf_data()
!--------------------------
! Import boundary file data
!--------------------------
  use error_messages
  use fds_head_arrays
  use fds_bndf_arrays
  use fds_mesh_arrays
  use global_constants
  use global_variables
  implicit none

  integer :: i,i1,i2,j1,j2,k1,k2,ios,nior
  integer :: ibndf,jbndf,imesh,ipatch,itime,c,c1,c2
  integer :: nmesh,npatch,nb,nm

  character(len=chr30) :: quantity,short_name,units

  !-------------------------
  ! Chosen transfer quantity
  !-------------------------

  nmesh=ubound(fds_mesh_ijk,1); ibndf=selected_bndf()
  allocate(fds_bndf_time(fds_bndf_ntim(1)))
  allocate(fds_bndf_data(fds_bndf_ntim(1),sum(fds_bndf_nnod,1)))

  mesh_loop: do imesh=1,nmesh

    open(unit=iochannel(1),file=trim(fds_bndf_file(ibndf,imesh)), &
      status='old',iostat=ios,form='unformatted')
    if (ios /= 0) call error_open_file(fds_bndf_file(ibndf,imesh))

    read(iochannel(1),iostat=ios) quantity
    if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
    read(iochannel(1),iostat=ios) short_name
    if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
    read(iochannel(1),iostat=ios) units
    if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
    read(iochannel(1),iostat=ios) npatch
    if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))

    fds_bndf_pat=0
    patch_loop: do ipatch=1,npatch
      read(iochannel(1),iostat=ios) i1,i2,j1,j2,k1,k2,nior,nb,nm
      if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
      fds_bndf_pat(ipatch)=(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
    end do patch_loop

    time_loop: do itime=1,fds_bndf_ntim(imesh)
      read(iochannel(1),iostat=ios) fds_bndf_time(itime)
      if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))

      node_loop: do i=1,npatch
        c1=1+sum(fds_bndf_pat(1:i-1))+sum(fds_bndf_nnod(1:imesh-1))
        c2=sum(fds_bndf_pat(1:i))+sum(fds_bndf_nnod(1:imesh-1))
        read(iochannel(1),iostat=ios) (fds_bndf_data(itime,c),c=c1,c2)
        if (ios /= 0) call error_read_file(fds_bndf_file(ibndf,imesh))
      end do node_loop

    end do time_loop

    close(unit=iochannel(1))
  
  end do mesh_loop

  !--------------------------
  ! Heat transfer coefficient
  !--------------------------

  if (read_hcoeff) then

    nmesh=ubound(fds_mesh_ijk,1); jbndf=selected_hcoeff_bndf()
    allocate(fds_hcoeff_bndf_time(fds_hcoeff_bndf_ntim(1)))
    allocate(fds_hcoeff_bndf_data(fds_hcoeff_bndf_ntim(1),sum(fds_hcoeff_bndf_nnod,1)))

    hcoeff_mesh_loop: do imesh=1,nmesh

      open(unit=iochannel(1),file=trim(fds_bndf_file(jbndf,imesh)), &
        status='old',iostat=ios,form='unformatted')
      if (ios /= 0) call error_open_file(fds_bndf_file(jbndf,imesh))

      read(iochannel(1),iostat=ios) quantity
      if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))
      read(iochannel(1),iostat=ios) short_name
      if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))
      read(iochannel(1),iostat=ios) units
      if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))
      read(iochannel(1),iostat=ios) npatch
      if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))

      fds_hcoeff_bndf_pat=0
      hcoeff_patch_loop: do ipatch=1,npatch
        read(iochannel(1),iostat=ios) i1,i2,j1,j2,k1,k2,nior,nb,nm
        if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))
        fds_hcoeff_bndf_pat(ipatch)=(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
      end do hcoeff_patch_loop

      hcoeff_time_loop: do itime=1,fds_hcoeff_bndf_ntim(imesh)
        read(iochannel(1),iostat=ios) fds_hcoeff_bndf_time(itime)
        if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))

        hcoeff_node_loop: do i=1,npatch
          c1=1+sum(fds_hcoeff_bndf_pat(1:i-1))+sum(fds_hcoeff_bndf_nnod(1:imesh-1))
          c2=sum(fds_hcoeff_bndf_pat(1:i))+sum(fds_hcoeff_bndf_nnod(1:imesh-1))
          read(iochannel(1),iostat=ios) (fds_hcoeff_bndf_data(itime,c),c=c1,c2)
          if (ios /= 0) call error_read_file(fds_bndf_file(jbndf,imesh))
        end do hcoeff_node_loop

      end do hcoeff_time_loop

      close(unit=iochannel(1))
  
    end do hcoeff_mesh_loop

  end if

end subroutine import_bndf_data

subroutine filter_bndf_data
!--------------------------
! Filter unwanted BNDF data
!--------------------------
  use error_messages
  use fds_bndf_arrays
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,i_node,i_node_fds
  integer :: nrows,ncols
  logical :: compatible

  ! Relevant arrays
  !--------------------------
  ! fds_bndf_ior(:),
  ! fds_bndf_time(:), 
  ! fds_bndf_xyz(:,:),
  ! fds_bndf_nb(:),
  ! fds_bndf_nm(:),
  ! fds_bndf_np(:),
  ! fds_bndf_data(:,:)
  !
  ! fds_hcoeff_bndf_ior(:),
  ! fds_hcoeff_bndf_time(:), 
  ! fds_hcoeff_bndf_xyz(:,:),
  ! fds_hcoeff_bndf_nb(:),
  ! fds_hcoeff_bndf_nm(:),
  ! fds_hcoeff_bndf_np(:),
  ! fds_hcoeff_bndf_data(:,:)
  !--------------------------

  !--------------------------------------------------------------- 
  ! Check for full compatibility between the boundary files of the
  ! chosen transfer quantity and heat transfer coefficient
  !---------------------------------------------------------------

  if (read_hcoeff) then
    compatible=.true.
    if (size(fds_bndf_ior) /= size(fds_hcoeff_bndf_ior)) then
      compatible=.false.
    end if
    if (size(fds_bndf_time) /= size(fds_hcoeff_bndf_time)) then
      compatible=.false.
    end if
    if (size(fds_bndf_xyz) /= size(fds_hcoeff_bndf_xyz)) then
      compatible=.false.
    end if
    if (size(fds_bndf_nb) /= size(fds_hcoeff_bndf_nb)) then
      compatible=.false.
    end if
    if (size(fds_bndf_nm) /= size(fds_hcoeff_bndf_nm)) then
      compatible=.false.
    end if
    if (size(fds_bndf_np) /= size(fds_hcoeff_bndf_np)) then
      compatible=.false.
    end if
    if (size(fds_bndf_data) /= size(fds_hcoeff_bndf_data)) then
      compatible=.false.
    end if

    if (.not. compatible) then
      write(*,'(3(a))') 'Error: boundary files for ', trim(transfer_quantity), ' and heat_transfer_coefficient are not compatible'
      stop
    end if

  end if

  !-----------------------------------------------
  ! Convert connectivity table into numerical form
  !-----------------------------------------------

  if (nset_connectivity) then

    nrows=ubound(connectivity_table,1); ncols=ubound(connectivity_table,2)
    allocate(connectivity_table_num(nrows,ncols),stat=ios); call error_allocate(ios)
       
    connectivity_table_num=0 
    row_loop_0: do i=1,ubound(connectivity_table,1)
      col_loop_0: do j=2,ubound(connectivity_table,2)

        if (len_trim(connectivity_table(i,j)) /= 0) then
          if (numerical(connectivity_table(i,j))) then
            connectivity_table_num(i,j)=str2int(connectivity_table(i,j))
            ! Error handling
            if (connectivity_table_num(i,j) <= 0 .or. &
                connectivity_table_num(i,j) > maxval(fds_bndf_np,1)) then
              write(*,'(5(a))') 'ERROR: in NSET-BNDF connectivity table: patch number out of range (file ', &
                trim(quote(nset_input_file)), ', line ', trim(int2str(i)), ')'
              stop
            end if
          else
            write(*,'(a)') 'ERROR: non-numerical entry in NSET-BNDF connectivity table'
            stop
          end if
        end if

      end do col_loop_0
    end do row_loop_0

  end if

  !---------------------
  ! Count selected nodes
  !---------------------

  if (nset_connectivity .and. .not. read_all_fds_data) then

    !-------------------------
    ! Chosen transfer quantity
    !-------------------------

    nnodes_fds=0
    node_loop_1: do i_node=1,ubound(fds_bndf_xyz,1)

      row_loop_1: do i=1,ubound(connectivity_table,1)
        col_loop_1: do j=2,ubound(connectivity_table,2)
          if (fds_bndf_np(i_node) == connectivity_table_num(i,j) .and. &
            fds_bndf_nm(i_node) == 1) then
            !------------------------------------------------------------
            ! N.b. For some reason, only nodes from mesh one are accepted
            !------------------------------------------------------------
              nnodes_fds=nnodes_fds+1
            cycle node_loop_1
          end if
        end do col_loop_1
      end do row_loop_1

    end do node_loop_1

    ntimes_fds=ubound(fds_bndf_time,1)

    if (nnodes_fds <= 0) then
      write(*,'(a)') 'ERROR: no FDS nodes selected'
      stop
    end if

    allocate(fds_ior(nnodes_fds),stat=ios);             call error_allocate(ios)
    allocate(fds_time(ntimes_fds),stat=ios);            call error_allocate(ios)
    allocate(fds_xyz(nnodes_fds,3),stat=ios);           call error_allocate(ios)
    allocate(fds_data(ntimes_fds,nnodes_fds),stat=ios); call error_allocate(ios)
    allocate(fds_patch(nnodes_fds),stat=ios);           call error_allocate(ios)

    !----------------------------------------------------
    ! Create arrays containing nodal coordinates and data
    !----------------------------------------------------
   
    if (read_hcoeff) then
      allocate(fds_hcoeff(ntimes_fds,nnodes_fds),stat=ios); call error_allocate(ios)
    end if

    i_node_fds=1
    node_loop_2: do i_node=1,ubound(fds_bndf_xyz,1)

      row_loop_2: do i=1,ubound(connectivity_table,1)
        col_loop_2: do j=2,ubound(connectivity_table,2)
          if (fds_bndf_np(i_node) == connectivity_table_num(i,j) .and. &
            fds_bndf_nm(i_node) == 1) then

            ! Node coordinates, orientation and patch 
            fds_xyz(i_node_fds,1:3) = fds_bndf_xyz(i_node,1:3)
            fds_ior(i_node_fds)     = fds_bndf_ior(i_node)
            fds_patch(i_node_fds)   = fds_bndf_np(i_node)

            ! Time series data
            fds_data(1:ntimes_fds,i_node_fds)=fds_bndf_data(1:ntimes_fds,i_node)
            if (read_hcoeff) then
              fds_hcoeff(1:ntimes_fds,i_node_fds)=fds_hcoeff_bndf_data(1:ntimes_fds,i_node)
            end if

            i_node_fds=i_node_fds+1

            cycle node_loop_2
          end if
        end do col_loop_2
      end do row_loop_2
 
    end do node_loop_2

  else

    nnodes_fds=ubound(fds_bndf_xyz,1)
    ntimes_fds=ubound(fds_bndf_time,1)

    allocate(fds_ior(nnodes_fds),stat=ios);               call error_allocate(ios)
    allocate(fds_time(ntimes_fds),stat=ios);              call error_allocate(ios)
    allocate(fds_xyz(nnodes_fds,3),stat=ios);             call error_allocate(ios)
    allocate(fds_data(ntimes_fds,nnodes_fds),stat=ios);   call error_allocate(ios)
    allocate(fds_patch(nnodes_fds),stat=ios);             call error_allocate(ios)

    fds_ior   = fds_bndf_ior
    fds_xyz   = fds_bndf_xyz
    fds_data  = fds_bndf_data
    fds_patch = fds_bndf_np

    if (read_hcoeff) then
      allocate(fds_hcoeff(ntimes_fds,nnodes_fds),stat=ios); call error_allocate(ios)
      fds_hcoeff = fds_hcoeff_bndf_data
    end if

  end if

  fds_time = fds_bndf_time

end subroutine filter_bndf_data

function selected_bndf() result(answer)
!---------------------------------------------------
! Return the index of the selected boundary quantity
!---------------------------------------------------
  use error_messages
  use fds_head_arrays
  use fds_bndf_arrays
  use global_constants
  use global_variables
  implicit none

  integer :: ibndf,nbndf,answer
  character(len=chr80) :: transfer_qnty

  if (trim(transfer_quantity) == 'wall_temperature') then
    transfer_qnty='WALL TEMPERATURE'
  else if (trim(transfer_quantity) == 'net_heat_flux') then
    transfer_qnty='NET HEAT FLUX'
  else if (trim(transfer_quantity) == 'adiabatic_surface_temperature') then
    transfer_qnty='ADIABATIC SURFACE TEMPERATURE'
  end if

  answer=0
  nbndf=ubound(fds_bndf_qnty,1)
  do ibndf=1,nbndf
    if (trim(fds_bndf_qnty(ibndf)) == trim(transfer_qnty)) then
      answer=ibndf
    end if
  end do

end function selected_bndf

function selected_hcoeff_bndf() result(answer)
!--------------------------------------------------------------------
! Return the index of the boundary quantity HEAT TRANSFER COEFFICIENT
!--------------------------------------------------------------------
  use error_messages
  use fds_head_arrays
  use fds_bndf_arrays
  use global_constants
  use global_variables
  implicit none

  integer :: ibndf,nbndf,answer
  character(len=chr80) :: transfer_qnty

  transfer_qnty='HEAT TRANSFER COEFFICIENT'

  answer=0
  nbndf=ubound(fds_bndf_qnty,1)
  do ibndf=1,nbndf
    if (trim(fds_bndf_qnty(ibndf)) == trim(transfer_qnty)) then
      answer=ibndf
    end if
  end do

end function selected_hcoeff_bndf

subroutine check_bndf_exists()
!-----------------------------------------------------
! Check that the requested BNDF-namelist record exists
!-----------------------------------------------------
  use global_variables
  use string_handling
  implicit none

  character(len=chr80) :: transfer_qnty

  if (trim(transfer_quantity) == 'wall_temperature') then
    transfer_qnty='WALL TEMPERATURE'
  else if (trim(transfer_quantity) == 'net_heat_flux') then
    transfer_qnty='NET HEAT FLUX'
  else if (trim(transfer_quantity) == 'adiabatic_surface_temperature') then
    transfer_qnty='ADIABATIC SURFACE TEMPERATURE'
  end if

  if (selected_bndf() == 0) then
    write(*,'(4(a))') 'ERROR: BNDF-namelist record for ', &
      trim(quote(transfer_qnty)), ' was not found in file ', trim(quote(fds_input_file))
    stop
  end if

end subroutine check_bndf_exists

subroutine check_hcoeff_bndf_exists()
!-------------------------------------------------------------------------
! Check that the BNDF-namelist record for HEAT TRANSFER COEFFICIENT exists
!-------------------------------------------------------------------------
  use global_variables
  use string_handling
  implicit none

  character(len=chr80) :: transfer_qnty

  transfer_qnty='HEAT TRANSFER COEFFICIENT'

  if (selected_bndf() == 0) then
    write(*,'(4(a))') 'ERROR: BNDF-namelist record for ', &
      trim(quote(transfer_qnty)), ' was not found in file ', trim(quote(fds_input_file))
    stop
  end if

end subroutine check_hcoeff_bndf_exists

end module fds_reader
