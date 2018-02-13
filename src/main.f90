!--------------------------------------------------------------------------
! Program:      fds2fem
! Description:  FDS-ABAQUS/ANSYS coupling tool for fire-structural analysis
!               v2: CFAST and time-temperature curves input added
! Version:      2.0
! Date:         29.1.2018
! Authorq:      Timo Korhonen/Antti Paajanen
! VTT Technical Research Centre of Finland Ltd
!
!--------------------------------------------------------------------------
!
! Fixes, improvements and questions
!
!   - Consistency check for data found in FDS input and output files
!     --* Device units, quantities and names
!     --* Are device coordinates within the simulation domain?
!   - Support for nonuniform meshes
!   - Do things work for overlapping meshes?
!   - How do we handle holes?
!   - How do we handle objects that are created or removed during the simulation?
!   - How to handle hidden BNDF faces?
!   - Does the code handle correctly the use of spaces in e.g. node set names?
!   - Do ABAQUS include-keywords cause any trouble? Probably yes...
!   - What if there are two devices with the same ID?
!   - What if a node belongs to more than one node set in ABAQUS?
!
!--------------------------------------------------------------------------
! Time-temperature curves input
!   - iso_834, ec1-1-1_hc, ec1-1-2_ex, astm_e119 implemented
!   - Only devc_to_nset mapping (devc => time-temp "name")
!   - Always uses nset-connectivity file, e.g.:
!     nset iso_834 epsilon hcoeff time_shift
!                  epsilon hcoeff (time_shift=0s)
!                  time_shift (use default eps&hc)
!   - Hcoeff/epsilon read from config-file (default) or nset-file
!   - Transfer quantities: T_surf, T_ast (same output, because "furnace")
!   - No vtk/xyz etc output for time-temp curve inputs
!
! CFAST zone model input as target outputs: CHID_w.csv file (CHID.in also needed)
!   - Only devc_to_nset mapping (devc => target ID/number)
!   - Always uses nset-connectivity file:
!     Give either target index or target name
!     Many targets can be given on one row => use average (and average emissivity)
!   - Transfer quantities: T_surf, T_ast, q_net (just Ansys, not tested)
!   - Emissivity input read from CFAST input file (CHID.in)
!   - Hcoeff read from config-file
!   - vtk/xyz/etc output possible, tested for vtk
!
!   ToDo: Ansys + Cfast + Tast? How epsilon & hc are transfered?
!         Ansys output: T_ast only with hcoeff read (Cfast: fill the hcoeff time series
!                        with the same constant hcoeff? (Epsilon from CFAST input file)
!         COMPA input, use nset-file (to start with)
!         TARGET input, use xyz ? mikä huone, jos seinän molemmin puolin inputtia
!
!--------------------------------------------------------------------------

program main
  use global_constants
  use global_variables
  use miscellaneous

  use cfg_reader

  use iso_reader
  use cfast_dump
  
  use fds_reader
  use fds_stats
  use fds_dump
  
  use fem_stats
  use fem_dump

  use ansys_reader
  !use ansys_dump
  use ansys_output
  
  use abaqus_reader
  !use abaqus_dump
  use abaqus_output
  
  use matching
  use mapping
  use iso_mapping
  
  implicit none

  real(kind=rk) :: t1,t2

  call cpu_time(t1)

  call allocate_channels()
 
  !=============================
  ! MODULE: Configuration reader
  !=============================

  call cfg_reader_module()

  !===================
  ! MODULE: FDS reader
  !===================

  if (iso_curve .Or. cfast_input) then
     call iso_reader_module()
  else
     call fds_reader_module()
  end if

  ! Relevant new arrays
  !   fds_id(:),
  !   fds_ior(:),
  !   fds_time(:),
  !   fds_idevc(:),
  !   fds_patch(:),
  !   fds_xyz(:,:),
  !   fds_data(:,:)

  !-----------
  ! FEM reader
  !-----------

  select case (trim(fem_mode))
  case ('abaqus')

    !======================
    ! MODULE: ABAQUS reader
    !======================

    call abaqus_reader_module()

    ! Relevant new arrays
    !   abaqus_nset(:)
    !   abaqus_node_name(:),
    !   abaqus_xyz(:,:)

  case ('ansys')

    !=====================
    ! MODULE: ANSYS reader
    !=====================

    call ansys_reader_module()

    ! Relevant new arrays
    !   ansys_node_number(:),
    !   ansys_xyz(:,:)

  case default
    ! Do nothing

  end select

  ! At this point, we need
  !   fem_xyz(:,:)
  !   fem_nset(:)
  !   fem_node(:)

  !=======================
  ! MODULE: Model matching
  !=======================

  call matching_module()

  !=====================
  ! MODULE: Mesh mapping
  !=====================

  If (.not.iso_curve) Call mapping_module() ! FDS and CFAST
  If (     iso_curve) Call mapping_iso()    ! Time-temperature curve

  ! Relevant new arrays
  !   fem_time(:),
  !   fem_data(:,:)

  !=======================
  ! MODULE: FDS statistics
  !=======================

  If (.Not.iso_curve .And. .Not.cfast_input) Call fds_stats_module()

  !==========================
  ! MODULE: ABAQUS statistics
  !==========================

  call fem_stats_module()

  !=================
  ! MODULE: FDS dump
  !=================

  If (.Not.iso_curve .And. .Not.cfast_input) Call fds_dump_module()

  !===================
  ! MODULE: CFAST dump
  !===================

  If (cfast_input) Call cfast_dump_module()

  !================
  ! FEM dump module
  !================

  call fem_dump_module()

  select case (trim(fem_mode))
  case ('abaqus')

    !======================
    ! MODULE: ABAQUS output
    !======================

    call abaqus_output_module()

  case ('ansys')

    !=====================
    ! MODULE: ANSYS output
    !=====================

    call ansys_output_module()

  case default
    ! Do nothing

  end select

  !-------
  ! Lastly
  !-------
  
  call cpu_time(t2)

  write(*,'(a)') ''
  write(*,'(3(a))') 'Ready in ', trim(sec2str(t2-t1)), ' s'

end program main
