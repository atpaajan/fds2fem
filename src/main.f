!--------------------------------------------------------------------------
! Program:      fds2fem
! Description:  FDS-ABAQUS/ANSYS coupling tool for fire-structural analysis
! Version:      1.0
! Date:         06.11.2012
! Author:       VTT Technical Research Centre of Finland
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

program main
  use global_constants
  use global_variables
  use miscellaneous

  use cfg_reader

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

  call fds_reader_module()

  ! Relevant new arrays
  !   fds_id(:),
  !   fds_ior(:),
  !   fds_time(:),
  !   fds_idevc(:),
  !   fds_patch(:),
  !   fds_xyz(:,:),
  !   fds_data(:,:)

  !=======================
  ! MODULE: FDS statistics
  !=======================

 ! call fds_stats_module()

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

  call mapping_module()

  ! Relevant new arrays
  !   fem_time(:),
  !   fem_data(:,:)

  !=======================
  ! MODULE: FDS statistics
  !=======================

  call fds_stats_module()

  !==========================
  ! MODULE: ABAQUS statistics
  !==========================

  call fem_stats_module()

  !=================
  ! MODULE: FDS dump
  !=================

  call fds_dump_module()

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
