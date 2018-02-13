!----------------------------------------------------
! Functions and subroutines for program configuration
!----------------------------------------------------
module cfg_reader

contains

subroutine cfg_reader_module()
!---------
! Main sub
!---------
  implicit none

  !-----------------------------
  ! Parse command-line arguments
  !-----------------------------

  call parse_command_arguments()

  !-------------------
  ! Output header line
  !-------------------

  call display_header()

  !-------------------------------------------
  ! Parse configuration and connectivity files
  !-------------------------------------------

  write(*,'(a)') ''
  write(*,'(a)') 'Configuration reader module'

  call parse_configuration_file()
  call parse_connectivity_file()

end subroutine cfg_reader_module

!***************
! Auxiliary subs
!***************

subroutine display_header()
!-------------------
! Output header line
!-------------------
  use miscellaneous
  implicit none

  character(len=8)  :: date
  character(len=10) :: time
  character(len=80) :: both

  call date_and_time(date,time)

  both=date(7:8) // '.' // date(5:6) // '.' // &
       date(1:4) // ' ' // time(1:2) // ':' // time(3:4)

  write(*,'(11(a))') 'FDS2FEM version 2.0'
  write(*,'(2(a))')  'Run: ', trim(both)

end subroutine display_header

subroutine display_usage(program_name)
!----------------------
! Display help on usage 
!----------------------
  use global_constants
  implicit none

  character(len=chr80) :: program_name

  write(*,'(a,a,a)') 'Usage: ', trim(program_name), ' [option] input_file'
  write(*,'(a)')     'Options:'
  write(*,'(a)')     '  -h, --help         Display this help and exit'
  write(*,'(a)')     '  -v, --version      Print version information'
  write(*,'(a)')     '  -f, --full-output  For debugging purposes'
  write(*,'(a)')     ''
  write(*,'(a)')     'Coupling tool for FDS/CFAST/Time-temp curve - ABAQUS/ANSYS simulations'
  write(*,'(a)')     ''

end subroutine display_usage

subroutine display_version()
!----------------------------
! Display version information
!----------------------------
  use global_constants
  implicit none

  write(*,'(a)') 'FDS2FEM version 2.0'
  write(*,'(a)') 'VTT Technical Research Centre of Finland Ltd, 2018'
  
end subroutine display_version

subroutine parse_command_arguments()
!-----------------------------
! Parse command-line arguments
!-----------------------------
  use global_constants
  use global_variables
  implicit none

  integer :: iargc,nargc
  character(len=chr80) :: argu,program_name

  config_file = ''
  full_output = .false.

  call get_command_argument(0,program_name)
  nargc=command_argument_count()
  if (nargc == 0) then
    call display_usage(program_name)
    stop
  end if
 
  iargc=1
  do while (iargc <= nargc)
    call get_command_argument(iargc,argu)
    argu=trim(adjustl(argu))

    if (argu(1:2) == '--') then

      select case(argu(3:len_trim(argu)))
        case ('help')
          call display_usage(program_name)
          stop
        case ('version')
          call display_version()
          stop        
        case ('full-output')
          full_output=.true.
        case default
          write(*,'(2(a))') 'Unknown option: ', trim(argu)
          stop
      end select

    else if (argu(1:1) == '-') then

        select case(argu(2:len_trim(argu)))
        case ('h')
          call display_usage(program_name)
          stop
        case ('v')
          call display_version()
          stop
        case ('f')
          full_output=.true.
        case default
          write(*,'(2(a))') 'Unknown option: ', trim(argu)
          stop
      end select

    else
      if (len_trim(config_file) == 0) then
        config_file=trim(argu)
      end if
    end if 

    iargc=iargc+1
  end do

  if (len_trim(config_file) == 0) then
    write(*,'(a)') 'ERROR: no configuration file given'
    stop
  end if

end subroutine parse_command_arguments

subroutine parse_configuration_file
!--------------------------
! Configuration file parser
!--------------------------
  use error_messages
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,ios,line_number
  character(len=chr80) :: identifier
  character(len=chr80), dimension(3) :: paramstr
  character(len=input_line_length) :: input_line

  !-------------------
  ! Set default values
  !-------------------

  fds_input_file          = ''
  fem_input_file          = ''
  nset_input_file         = ''

  fds_output              = ''
  transfer_quantity       = ''

  fem_mode                = 'off'
  
  dump_fds_nodes          = 'off'
  dump_fds_data           = 'off'
  dump_fds_model          = 'off'

  dump_cfast_nodes        = 'off'
  dump_cfast_data         = 'off'

  dump_fem_nodes          = 'off'
  dump_fem_data           = 'off'
  dump_fem_model          = 'off'

  mapping_method          = 'nearest'

  mp_n                    = 4
  mp_nmx                  = 8
  mp_cut                  = 0.0
  mp_del                  = 0.1
  mp_deg                  = 2.0

  match_translate         = .false.
  match_rotate            = .false.
  manual_translate        = .false.
  manual_rotate           = .false.
  automatic_translate     = .false.
  automatic_rotate        = .false.

  origin_fds(1:3)         = (/ 0.0, 0.0, 0.0 /)
  origin_fem(1:3)         = (/ 0.0, 0.0, 0.0 /)

  e_alpha                 = 0.0
  e_beta                  = 0.0
  e_gamma                 = 0.0

  fds_statistics          = .false.
  fem_statistics          = .false.

  fds_data_available      = .false.
  fds_xyz_available       = .false.
  cfast_data_available    = .false.
  cfast_xyz_available     = .false.
  fem_data_available      = .false.
  fem_xyz_available       = .false.

  fds_cm_calculated       = .false.
  fds_range_calculated    = .false.
  fem_cm_calculated       = .false.
  fem_range_calculated    = .false.

  read_hcoeff             = .false.
  ansys_ast               = .false.

  hcoeff                  = 25.0
  emissivity              = 0.90

  cfast_input             = .false.
  iso_curve               = .false.
  iso_ntimes              = 361
  iso_tbegin              = 0.0
  iso_tend                = 3600.0
  
  !----------------------------
  ! User given parameter values
  !----------------------------

  open(unit=iochannel(1),file=trim(config_file),status='old',iostat=ios)
  if (ios /= 0) call error_open_file(config_file)

  line_number=0
  parse_loop: do
    read(iochannel(1),'(a)',iostat=ios) input_line
    select case(ios)
    case (:-1)
      exit parse_loop
    case (1:)
      call error_read_file(config_file)
    end select

    line_number=line_number+1
    input_line=adjustl(input_line)

    ! Skip empty lines and comment lines
    if (input_line(1:1) == '#' .or. len_trim(input_line) == 0) cycle parse_loop

    identifier=''; paramstr=''
    read(input_line,*,iostat=ios) identifier,(paramstr(i),i=1,3)

    select case(identifier)

    case ('iso_curve') 
      !------------------------------------------
      ! temperature-time curve (needs nset_input)
      !------------------------------------------
      if (len_trim(paramstr(1)) /= 0) then
        if (trim(lowercase(paramstr(1))) == 'on') then
          iso_curve=.true.
        else if (trim(lowercase(paramstr(1))) == 'off') then
          iso_curve=.false.
        else
          write(*,'(5(a))') 'ERROR: ISO curve has to be either ON or OFF (file ', &
            trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('cfast_input')
      !------------
      ! cfast_input
      !------------
      ! The body of the CFAST input file name
      if (len_trim(paramstr(1)) /= 0) then
        fds_input_file=trim(paramstr(1))
      else
        write(*,'(5(a))') 'ERROR: missing CFAST input file (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if
      cfast_input = .true.

    case ('fds_input')
      !----------
      ! fds_input
      !----------
      if (len_trim(paramstr(1)) /= 0) then
        fds_input_file=trim(paramstr(1))
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('fem_input')
      !----------
      ! fem_input
      !----------
      if (len_trim(paramstr(1)) /= 0) then
        fem_input_file=trim(paramstr(1))
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if
   
    case ('fem_mode')
      !---------
      ! fem_mode
      !---------
      if (len_trim(paramstr(1)) /= 0) then
        fem_mode=trim(paramstr(1))
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      if (trim(fem_mode) /= 'abaqus' .and. &
          trim(fem_mode) /= 'ansys' .and. &
          trim(fem_mode) /= 'off') then
        write(*,'(7(a))') 'ERROR: unknown FEM mode ', trim(quote(fem_mode)), &
            ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
      end if

    case ('nset_input')
      !-----------
      ! nset_input
      !-----------
      if (len_trim(paramstr(1)) /= 0) then
        nset_input_file=trim(paramstr(1))
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('quantity')
      !---------
      ! quantity
      !---------
      if (len_trim(paramstr(1)) /= 0) then
        transfer_quantity=trim(lowercase(paramstr(1)))
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      if (trim(transfer_quantity) /= 'wall_temperature' .and. &
          trim(transfer_quantity) /= 'net_heat_flux' .and. &
          trim(transfer_quantity) /= 'adiabatic_surface_temperature') then
          write(*,'(7(a))') 'ERROR: unknown transfer quantity ', trim(quote(transfer_quantity)), &
            ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('fds_output')
      !-----------
      ! fds_output
      !-----------
      if (len_trim(paramstr(1)) /= 0) then
        fds_output=trim(lowercase(paramstr(1)))
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      if (trim(fds_output) /= 'devc' .and. trim(fds_output) /= 'bndf') then
        write(*,'(7(a))') 'ERROR: unknown FDS output quantity ', trim(quote(fds_output)), &
          ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('dump_fds_nodes')
      !---------------
      ! dump_fds_nodes
      !---------------
      if (len_trim(paramstr(1)) /= 0) then
        dump_fds_nodes=trim(lowercase(paramstr(1)))
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      if (trim(dump_fds_nodes) /= 'off'  .and. &
          trim(dump_fds_nodes) /= 'xyz'  .and. &
          trim(dump_fds_nodes) /= 'vtk') then
        write(*,'(7(a))') 'ERROR: unknown dump format ', trim(quote(dump_fds_nodes)), &
            ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('dump_fds_data')
      !--------------
      ! dump_fds_data
      !--------------
      if (len_trim(paramstr(1)) /= 0) then
        dump_fds_data=trim(lowercase(paramstr(1)))
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      if (trim(dump_fds_data) /= 'off'  .and. &
          trim(dump_fds_data) /= 'sdf'  .and. &
          trim(dump_fds_data) /= 'vtk') then
        write(*,'(7(a))') 'ERROR: unknown dump format ', trim(quote(dump_fds_data)), &
            ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('dump_fds_model')
      !---------------
      ! dump_fds_model
      !---------------
      if (len_trim(paramstr(1)) /= 0) then
        dump_fds_model=trim(lowercase(paramstr(1)))
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      if (trim(dump_fds_model) /= 'off'  .and. &
          trim(dump_fds_model) /= 'vtk') then
        write(*,'(7(a))') 'ERROR: unknown dump format ', trim(quote(dump_fds_data)), &
            ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('dump_cfast_nodes')
      !-----------------
      ! dump_cfast_nodes
      !-----------------
      if (len_trim(paramstr(1)) /= 0) then
        dump_cfast_nodes=trim(lowercase(paramstr(1)))
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      if (trim(dump_cfast_nodes) /= 'off'  .and. &
          trim(dump_cfast_nodes) /= 'xyz'  .and. &
          trim(dump_cfast_nodes) /= 'vtk') then
        write(*,'(7(a))') 'ERROR: unknown dump format ', trim(quote(dump_fds_nodes)), &
            ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('dump_cfast_data')
      !----------------
      ! dump_cfast_data
      !----------------
      if (len_trim(paramstr(1)) /= 0) then
        dump_cfast_data=trim(lowercase(paramstr(1)))
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      if (trim(dump_cfast_data) /= 'off'  .and. &
          trim(dump_cfast_data) /= 'sdf'  .and. &
          trim(dump_cfast_data) /= 'vtk') then
        write(*,'(7(a))') 'ERROR: unknown dump format ', trim(quote(dump_fds_data)), &
            ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('dump_fem_nodes')
      !----------------
      ! dump_fem_nodes
      !----------------
      if (len_trim(paramstr(1)) /= 0) then
        dump_fem_nodes=trim(lowercase(paramstr(1)))
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      if (trim(dump_fem_nodes) /= 'off'  .and. &
          trim(dump_fem_nodes) /= 'xyz'  .and. &
          trim(dump_fem_nodes) /= 'vtk') then
        write(*,'(7(a))') 'ERROR: unknown dump format ', trim(quote(dump_fem_nodes)), &
            ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('dump_fem_data')
      !---------------
      ! dump_fem_data
      !---------------
      if (len_trim(paramstr(1)) /= 0) then
        dump_fem_data=trim(lowercase(paramstr(1)))
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      if (trim(dump_fem_data) /= 'off'  .and. &
          trim(dump_fem_data) /= 'sdf'  .and. &
          trim(dump_fem_data) /= 'vtk') then
        write(*,'(7(a))') 'ERROR: unknown dump format ', trim(quote(dump_fem_data)), &
            ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('dump_fem_model')
      !----------------
      ! dump_fem_model
      !----------------
      if (len_trim(paramstr(1)) /= 0) then
        dump_fem_model=trim(lowercase(paramstr(1)))
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      if (trim(dump_fem_model) /= 'off'  .and. &
          trim(dump_fem_model) /= 'vtk') then
        write(*,'(7(a))') 'ERROR: unknown dump format ', trim(quote(dump_fem_data)), &
            ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('mapping')
      !--------
      ! mapping 
      !--------
      if (len_trim(paramstr(1)) /= 0) then
        mapping_method=trim(lowercase(paramstr(1)))
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      if (trim(mapping_method) /= 'nearest'       .and. &
          trim(mapping_method) /= 'devc_to_nset'  .and. &
          trim(mapping_method) /= 'off') then
        write(*,'(7(a))') 'ERROR: unknown mapping method ', trim(quote(mapping_method)), &
            ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('mp_deg')
      !-------
      ! mp_deg
      !-------
      if (len_trim(paramstr(1)) /= 0) then
        read(paramstr(1),*,iostat=ios) mp_deg
        if (ios /= 0) then
          write(*,'(2(a))') 'ERROR: in reading value for parameter ', trim(quote(identifier))
          stop
        end if
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      if (mp_deg < 0.0) then
        mp_deg = 2.0
        write(*,'(6(a))') 'WARNING: Negative values for mp_deg are not allowed  ', &
          ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
        write(*,'(a,f6.3)') 'WARNING: Default power is used: mp_deg =', mp_deg
      end if
      if (mp_deg > 10.0) then
        mp_deg = 10.0
        write(*,'(6(a))') 'WARNING: Too large power for the distance dependence weights  ', &
          ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
        write(*,'(a,f6.3)') 'WARNING: Power to be used: mp_deg =', mp_deg
      end if

    case ('mp_n')
      !-----
      ! mp_n
      !-----
      if (len_trim(paramstr(1)) /= 0) then
        read(paramstr(1),*,iostat=ios) mp_n
        if (ios /= 0) then
          write(*,'(2(a))') 'ERROR: in reading value for parameter ', trim(quote(identifier))
          stop
        end if
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('mp_nmx')
      !-------
      ! mp_nmx
      !-------
      if (len_trim(paramstr(1)) /= 0) then
        read(paramstr(1),*,iostat=ios) mp_nmx
        if (ios /= 0) then
          write(*,'(2(a))') 'ERROR: in reading value for parameter ', trim(quote(identifier))
          stop
        end if
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('mp_cut')
      !-------
      ! mp_cut
      !-------
      if (len_trim(paramstr(1)) /= 0) then
        read(paramstr(1),*,iostat=ios) mp_cut
        if (ios /= 0) then
          write(*,'(2(a))') 'ERROR: in reading value for parameter ', trim(quote(identifier))
          stop
        end if
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      if (mp_cut < 0.0) then
        mp_cut = 0.0
        write(*,'(6(a))') 'WARNING: Negative cut-off radius is not allowed ', &
          ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
        write(*,'(a,f6.3)') 'WARNING: No cut-off radius is used (mp_cut=0.0)'
      end if

    case ('mp_del')
      !-------
      ! mp_del
      !-------
      if (len_trim(paramstr(1)) /= 0) then
        read(paramstr(1),*,iostat=ios) mp_del
        if (ios /= 0) then
          write(*,'(2(a))') 'ERROR: in reading value for parameter ', trim(quote(identifier))
          stop
        end if
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      if (mp_del < 0.0) then
        mp_del= 0.1
        write(*,'(6(a))') 'WARNING: Cut-off radius delta-parameter < 0.0 ', &
          ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
        write(*,'(a,f6.3)') 'WARNING: Default cut-off radius delta-parameter is used: mp_del =', mp_del
      end if

    case ('origin_fds')
      !-----------
      ! origin_fds
      !-----------
      do i=1,3
        if (len_trim(paramstr(i)) /= 0) then
          read(paramstr(i),*,iostat=ios) origin_fds(i) 
          if (ios /= 0) then
            write(*,'(2(a))') 'ERROR: in reading value for parameter ', trim(quote(identifier))
            stop
          end if
        else
          write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
            ', line ', trim(int2str(line_number)), ')'
          stop
        end if
      end do

    case ('origin_fem')
      !------------
      ! origin_fem
      !------------
      do i=1,3
        if (len_trim(paramstr(i)) /= 0) then
          read(paramstr(i),*,iostat=ios) origin_fem(i) 
          if (ios /= 0) then
            write(*,'(2(a))') 'ERROR: in reading value for parameter ', trim(quote(identifier))
            stop
          end if
        else
          write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
            ', line ', trim(int2str(line_number)), ')'
          stop
        end if
      end do

    case ('euler_alpha')
      !------------
      ! euler_alpha
      !------------
      if (len_trim(paramstr(1)) /= 0) then
        read(paramstr(1),*,iostat=ios) e_alpha
        if (ios /= 0) then
          write(*,'(2(a))') 'ERROR: in reading value for parameter ', trim(quote(identifier))
          stop
        end if

        ! Degrees to radians
        e_alpha=(e_alpha/180.0)*(4.0*atan(1.0))

      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('euler_beta')
      !------------
      ! euler_beta
      !------------
      if (len_trim(paramstr(1)) /= 0) then
        read(paramstr(1),*,iostat=ios) e_beta
        if (ios /= 0) then
          write(*,'(2(a))') 'ERROR: in reading value for parameter ', trim(quote(identifier))
          stop
        end if

        ! Degrees to radians
        e_beta=(e_beta/180.0)*(4.0*atan(1.0))

      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('euler_gamma')
      !------------
      ! euler_beta
      !------------
      if (len_trim(paramstr(1)) /= 0) then
        read(paramstr(1),*,iostat=ios) e_gamma
        if (ios /= 0) then
          write(*,'(2(a))') 'ERROR: in reading value for parameter ', trim(quote(identifier))
          stop
        end if

        ! Degrees to radians
        e_gamma=(e_gamma/180.0)*(4.0*atan(1.0))

      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('match_translate')
      !----------------
      ! match_translate
      !----------------
      if (len_trim(paramstr(1)) /= 0) then
        if (trim(lowercase(paramstr(1))) == 'manual') then
          match_translate     = .true.
          manual_translate    = .true.
          automatic_translate = .false.
        else if (trim(lowercase(paramstr(1))) == 'automatic') then
          match_translate     = .true.
          manual_translate    = .false.
          automatic_translate = .true.
        else if (trim(lowercase(paramstr(1))) == 'off') then
          match_translate     = .false.
          manual_translate    = .false.
          automatic_translate = .false.
        else
          write(*,'(5(a))') 'ERROR: unknown model matching option (file ', &
            trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('match_rotate')
      !-------------
      ! match_rotate
      !-------------
      if (len_trim(paramstr(1)) /= 0) then
        if (trim(lowercase(paramstr(1))) == 'manual') then
          match_rotate     = .true.
          manual_rotate    = .true.
          automatic_rotate = .false.
        else if (trim(lowercase(paramstr(1))) == 'automatic') then
          match_rotate     = .true.
          manual_rotate    = .false.
          automatic_rotate = .true.
        else if (trim(lowercase(paramstr(1))) == 'off') then
          match_rotate     = .false.
          manual_rotate    = .false.
          automatic_rotate = .false.
        else
          write(*,'(5(a))') 'ERROR: unknown model matching option (file ', &
            trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('fds_statistics')
      !---------------
      ! fds_statistics
      !---------------
      if (len_trim(paramstr(1)) /= 0) then
        if (trim(lowercase(paramstr(1))) == 'on') then
          fds_statistics=.true.
        else if (trim(lowercase(paramstr(1))) == 'off') then
          fds_statistics=.false.
        else
          write(*,'(5(a))') 'ERROR: FDS node set analysis has to be either ON or OFF (file ', &
            trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('fem_statistics')
      !----------------
      ! fem_statistics
      !----------------
      if (len_trim(paramstr(1)) /= 0) then
        if (trim(lowercase(paramstr(1))) == 'on') then
          fem_statistics=.true.
        else if (trim(lowercase(paramstr(1))) == 'off') then
          fem_statistics=.false.
        else
          write(*,'(5(a))') 'ERROR: ABAQUS node set analysis has to be either ON or OFF (file ', &
            trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('e_coeff')
      !--------
      ! e_coeff
      !--------
      if (len_trim(paramstr(1)) /= 0) then
        read(paramstr(1),*,iostat=ios) emissivity
        if (ios /= 0) then
          write(*,'(2(a))') 'ERROR: in reading value for parameter ', trim(quote(identifier))
          stop
        end if

        if (emissivity < 0.0 .or. emissivity > 1.0) then
          write(*,'(2(a))') 'ERROR: unphysical emissivity (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        end if 

      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('h_coeff')
      !--------
      ! h_coeff
      !--------
      if (len_trim(paramstr(1)) /= 0) then
        read(paramstr(1),*,iostat=ios) hcoeff
        if (ios /= 0) then
          write(*,'(2(a))') 'ERROR: in reading value for parameter ', trim(quote(identifier))
          stop
        end if

        if (hcoeff < 0.0) then
          !write(*,'(2(a))') 'ERROR: unphysical heat transfer coefficient (file ', trim(quote(config_file)), &
          !', line ', trim(int2str(line_number)), ')'
        end if 

      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('iso_ntimes')
       !-----------
       ! iso_ntimes
       !-----------
      if (len_trim(paramstr(1)) /= 0) then
        read(paramstr(1),*,iostat=ios) iso_ntimes
        if (ios /= 0) then
          write(*,'(2(a))') 'ERROR: in reading value for parameter ', trim(quote(identifier))
          stop
        end if
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('iso_tbegin')
       !-----------
       ! iso_tbegin
       !-----------
      if (len_trim(paramstr(1)) /= 0) then
        read(paramstr(1),*,iostat=ios) iso_tbegin
        if (ios /= 0) then
          write(*,'(2(a))') 'ERROR: in reading value for parameter ', trim(quote(identifier))
          stop
        end if
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case ('iso_tend')
       !---------
       ! iso_tend
       !---------
      if (len_trim(paramstr(1)) /= 0) then
        read(paramstr(1),*,iostat=ios) iso_tend
        if (ios /= 0) then
          write(*,'(2(a))') 'ERROR: in reading value for parameter ', trim(quote(identifier))
          stop
        end if
      else
        write(*,'(5(a))') 'ERROR: missing parameter (file ', trim(quote(config_file)), &
          ', line ', trim(int2str(line_number)), ')'
        stop
      end if

    case default
      write(*,'(7(a))') 'ERROR: unknown identifier ', trim(quote(identifier)), &
        ' (file ', trim(quote(config_file)), ', line ', trim(int2str(line_number)), ')'
      stop

    end select

  end do parse_loop

  close(unit=iochannel(1))

  !------------------------
  ! Input sensibility check
  !------------------------

  ! Inside-radius mapping
  if (mp_n < 1) then
    mp_nmx=mp_n; mp_del=0.0
    if (mp_cut < epsilon(mp_cut)) then
      write(*,'(a,f6.3)') 'ERROR: inside-radius mapping requires a cut-off radius greater than zero'
      stop 
    end if
  end if

  ! Hybrid mapping
  if (mp_n > 0) then
    if (mp_nmx < mp_n) then
      mp_nmx=mp_n
      write(*,'(a)') 'WARNING: mp_nmx < mp_n,  mp_nmx = mp_n is used'
    end if
  end if

  if ((trim(transfer_quantity) == 'adiabatic_surface_temperature' .or. &
      trim(transfer_quantity) == 'net_heat_flux') .and. &
      len_trim(nset_input_file) == 0 .and. trim(fem_mode) == 'abaqus') then
    write(*,'(a)') 'ERROR: an NSET connectivity file is required when transferring ADIABATIC SURFACE TEMPERATURE'
    stop
  end if

  if (trim(transfer_quantity) == 'adiabatic_surface_temperature') then
    if (hcoeff < 0.0) then
      ! If negative heat transfer coefficient is given, read from bndf-file
      read_hcoeff=.true. 
    end if
    if (trim(fem_mode) == 'ansys') then
      ansys_ast=.true.
    end if
  end if

  ! Time-temp curve, some input checks
  if (iso_curve) then
     if (fds_input_file .ne. '') then
        If (cfast_input) Then
           write(*,'(a)') 'ERROR: ISO curve and CFAST input are both selected'
           stop
        Else
           write(*,'(a)') 'ERROR: ISO curve and FDS input are both selected'
           stop
        End If
     end if
     if (trim(transfer_quantity) == 'net_heat_flux') then
        write(*,'(a)') 'ERROR: ISO curve and net_heat_flux as transfer quantity'
        stop
     end if
     if (hcoeff < 0.0) then
        ! If negative heat transfer coefficient is given, read from bndf-file, not for ISO
        write(*,'(a)') 'ERROR: ISO curve and hcoeff < 0'
        stop
     end if
     if (trim(transfer_quantity) == 'adiabatic_surface_temperature') then
        Write(*,'(a)') 'ISO curve and adiabatic surface temperature transfer quantity'
        Write(*,'(a,f6.3)') '   emissivity:           ',emissivity
        Write(*,'(a,f6.3)') '   heat transfer coeff.: ',hcoeff
        If (emissivity < 0.0 .Or. emissivity > 1.0 .Or. hcoeff < 0.0) Then
           Write(*,'(a)') 'ERROR: Unphysical parameter values for emissivity and/or h_coeff'
           Stop
        End If
     End If
     If (Trim(mapping_method) /= 'devc_to_nset') Then
        Write(*,'(a)') 'WARNING: ISO curve and mapping is not devc_to_nset'
        Write(*,'(a)') '         Mapping method is set to devc_to_nset'
     End If
     If (match_translate .Or. manual_translate .Or. automatic_translate) Then
        write(*,'(a)') 'WARNING: ISO curve, no translations done for devc_to_nset mapping'
     End If
     If (match_rotate .Or. manual_rotate .Or. automatic_rotate) Then
        write(*,'(a)') 'WARNING: ISO curve, no rotations done for devc_to_nset mapping'
     End If
     If (iso_tend <= iso_tbegin) Then
        Write(*,'(a)') 'ERROR: ISO curve, tend <= tbegin'
        Stop
     End If
     If (iso_ntimes < 2) Then
        Write(*,'(a)') 'ERROR: ISO curve, too few time points'
        Stop
     End If

     cfast_input         = .false.
     fds_output          = 'devc'
     fds_input_file      = ''
     dump_fds_nodes      = 'off'
     dump_fds_data       = 'off'
     dump_fds_model      = 'off'
     fds_statistics      = .false.
     mapping_method      = 'devc_to_nset'
     match_translate     = .false.
     match_rotate        = .false.
     manual_translate    = .false.
     manual_rotate       = .false.
     automatic_translate = .false.
     automatic_rotate    = .false.
  end if

  ! CFAST input, set some varibles (ToDo: input checks)
  If (cfast_input) Then
     iso_curve           = .false.
     fds_output          = 'devc'
     dump_fds_nodes      = 'off'
     dump_fds_data       = 'off'
     dump_fds_model      = 'off'
     fds_statistics      = .false.
     mapping_method      = 'devc_to_nset'
     match_translate     = .false.
     match_rotate        = .false.
     manual_translate    = .false.
     manual_rotate       = .false.
     automatic_translate = .false.
     automatic_rotate    = .false.
  End If
  
  !-----------------
  ! Ready to proceed
  !-----------------

  write(*,'(t3,2(a))') 'Instructions read from file ', trim(quote(config_file))

end subroutine parse_configuration_file

subroutine parse_connectivity_file()
!---------------------------------------------------------
! Parse input file for NSET-DEVC or NSET-BNDF connectivity 
!---------------------------------------------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,nlines,ncols_min,ncols_max
  integer :: i_ncols_min,i_ncols_max

  integer, dimension(:), allocatable :: ncols

  character(len=chr80) :: fmt_1
  character(len=input_line_length) :: input_line

  !--------------------
  ! Skipping conditions
  !--------------------

  if (len_trim(nset_input_file) == 0) then
    ! Nothing to be done
    nset_connectivity=.false.
    read_all_fds_data=.true.
    return
  end if

  if (len_trim(fds_input_file) == 0 .and. &
      len_trim(fem_input_file) == 0) then
    nset_connectivity=.false.
    read_all_fds_data=.true.
    return
  end if

  if (trim(fem_mode) == 'ansys') then
    write(*,'(a)') 'WARNING: NSET connectivity is not supported in FDS-ANSYS coupling' 
    return
  end if

  !------------------------
  ! Parse connectivity file
  !------------------------

  open(unit=iochannel(1),file=trim(nset_input_file),status='old',iostat=ios)
  if (ios /= 0) call error_open_file(nset_input_file)

  !------------
  ! Count lines
  !------------
  nlines=0
  count_lines: do
  read(iochannel(1),'(a)',iostat=ios) input_line
    select case(ios)
    case (:-1)
      exit count_lines
    case (0)
      nlines=nlines+1
    case(1:)
      call error_read_file(nset_input_file,nlines+1)
    end select
  end do count_lines

  if (nlines == 0) then
     if (.not.iso_curve) then
        If (cfast_input) Then
           write(*,'(2(a))') 'ERROR: no input found in NSET-Target connectivity file ', &
                trim(quote(nset_input_file))
           stop
        Else
           if (trim(fds_output) == 'devc') then
              write(*,'(2(a))') 'ERROR: no input found in NSET-DEVC connectivity file ', &
                   trim(quote(nset_input_file))
              stop
           else if (trim(fds_output) == 'bndf') then
              write(*,'(2(a))') 'ERROR: no input found in NSET-BNDF connectivity file ', &
                   trim(quote(nset_input_file))
              stop
           else
              write(*,'(2(a))') 'ERROR: no input found in NSET connectivity file ', &
                   trim(quote(nset_input_file))
              stop
           end if
        End If
     else
        write(*,'(2(a))') 'ERROR: ISO curve, no input found in NSET connectivity file ', &
             trim(quote(nset_input_file))
        stop
     end if
  end if

  rewind(iochannel(1))

  allocate(ncols(nlines),stat=ios); call error_allocate(ios)

  !--------------
  ! Count columns
  !--------------
  count_columns: do i=1,nlines
    read(iochannel(1),'(a)',iostat=ios) input_line
    select case(ios)
    case (:-1)
      exit count_columns
    case (0)
      ncols(i)=count_sdf_columns(input_line)
    case(1:)
      call error_read_file(nset_input_file,i)
    end select
  end do count_columns

  ncols_min=minval(ncols); i_ncols_min=minloc(ncols,1)
  ncols_max=maxval(ncols); i_ncols_max=maxloc(ncols,1)

  read_all_fds_data=.false.
  if (ncols_min == 1) then
    read_all_fds_data=.true.
  end if

  ! ISO nset: first column fem-node, second column "iso_834" etc time-temp curve.
  ! Additional columns (3-5):
  !                   eps hcoeff time_shift
  ! nset-name iso_834 0.9 25.0   0.0
  !                   eps hcoeff
  ! nset-name iso_834 0.9 25.0  
  !                   time_shift
  ! nset-name iso_834 0.0
  
  rewind(iochannel(1))

  allocate(connectivity_table(nlines,ncols_max),stat=ios); call error_allocate(ios)
  connectivity_table=''

  !--------------------
  ! Read connectivities
  !--------------------
  read_loop: do i=1,nlines
    read(iochannel(1),'(a)',iostat=ios) input_line
    if (ios /= 0) call error_read_file(nset_input_file,i)

    if (ncols(i) == 0) then
      write(*,*) 'ERROR: in NSET connectivity file: empty lines are not allowed'
      stop
    end if

    input_line=adjustl(input_line)

    read(input_line,*,iostat=ios) (connectivity_table(i,j),j=1,ncols(i))
    if (ios /= 0) call error_read_file(nset_input_file,i)

  end do read_loop

  close(unit=iochannel(1))

  write(*,'(t3,2(a))') 'Connectivities read from file ', trim(quote(nset_input_file))

  !------------------
  ! Check sensibility
  !------------------
  do i=1,ubound(connectivity_table,1)
    do j=1,ubound(connectivity_table,1)
      if (i /= j) then
        if (trim(connectivity_table(i,1)) == trim(connectivity_table(j,1))) then
          if (trim(fds_output) == 'devc') then
            write(*,'(5(a))') 'ERROR: duplicate NSET entry in NSET-DEVC connectivity file (file ', &
              trim(quote(nset_input_file)), ', line ', trim(int2str(j)), ')'
            stop
          else if (trim(fds_output) == 'bndf') then   
            write(*,'(5(a))') 'ERROR: duplicate NSET entry in NSET-BNDF connectivity file (file ', &
              trim(quote(nset_input_file)), ', line ', trim(int2str(j)), ')'
            stop
          else
            write(*,'(5(a))') 'ERROR: duplicate NSET entry in NSET connectivity file (file ', &
              trim(quote(nset_input_file)), ', line ', trim(int2str(j)), ')'
            stop
          end if
        end if
      end if
    end do
  end do

  !------------
  ! Full output
  !------------

  if (full_output) then
    
    write(*,'(a)') ''
    write(*,'(a)') '* Connectivity table'

    do i=1,nlines
      write(*,'(a,3x)',advance='no') '*'
      if (ncols(i) == 1) then
        fmt_1=format_char_sdf(ncols(i)+1)
        write(*,fmt_1) trim(connectivity_table(i,1)), '(empty)'
      else
        fmt_1=format_char_sdf(ncols(i))
        write(*,fmt_1) (trim(connectivity_table(i,j)),j=1,ncols(i))
      end if
    end do

  end if

  !------------------
  ! Deallocate memory
  !------------------
 
  deallocate(ncols,stat=ios); call error_allocate(ios)

  nset_connectivity=.true.

end subroutine parse_connectivity_file

end module cfg_reader

