!-----------------------------------------------------------
! Functions and subroutines for time-tempearture curve input
! Functions and subroutines for CFAST input
!-----------------------------------------------------------
Module iso_reader

  Logical :: cfast_old_format

Contains

Subroutine iso_reader_module()
!---------
! Main sub
!---------
  Use fds_bndf_arrays
  Use fds_devc_arrays
  Use global_constants
  Use global_variables
  Use mapping_arrays
  Use miscellaneous
  Use string_handling
  Implicit None

  Write(*,'(a)') ''
  Write(*,'(a)') 'ISO/CFAST reader module'

  If (Len_Trim(nset_input_file) == 0) Then 

    !------------------------
    ! No connectivities given
    !------------------------

    Write(*,'(t3,a)') 'Nothing to be done, no nset input file'

  Else
     
     If (cfast_input) Then
        Call read_cfast_input_file()
        Call import_cfast_data()
        Call filter_cfast_data()
        ! Relevant arrays
        !---------------------------------------------------------
        !   fds_devc_qnty(:), fds_devc_name(:), fds_devc_name_b(:)
        !   fds_devc_unit(:), fds_devc_time(:), fds_devc_data(:,:)
        !---------------------------------------------------------
     Else
        Call import_iso_data()
        Call filter_iso_data()
        ! Relevant arrays
        !---------------------
        !   fds_devc_time(:), 
        !   fds_devc_data(:,:)
     !---------------------
     End If
     
     Call deallocate_excess_fds_arrays()
     If (cfast_input) Then
        fds_xyz_available  =.False.; fds_model_available =.False.
        cfast_xyz_available=.True.;  cfast_data_available=.True.
     Else
        fds_xyz_available  =.False.; fds_data_available  =.False.; fds_model_available=.False.
        cfast_xyz_available=.False.; cfast_data_available=.False.
     End If
     
  End If

End Subroutine iso_reader_module

!***************
! Auxiliary subs
!***************

Function number_csv_rows(filename) Result(rows)
!------------------------------------------
! Find out the number of rows in a CSV-file
!------------------------------------------
  Use error_messages
  Use global_constants
  Implicit None

  Integer :: ios,rows
  Character(len=chr80) :: filename
  Character(len=input_line_length) :: input_line

  Open(unit=iochannel(1),file=Trim(filename),status='old',iostat=ios)
  If (ios /= 0) Call error_open_file(filename)

  ! How many rows? 
  rows=0
  rows_loop: Do  
    Read(iochannel(1),*,iostat=ios) input_line
    Select Case(ios)
    Case (:-1)
      Exit rows_loop
    Case (0)
      rows=rows+1
    Case (1:)
      Call error_read_file(filename)
    End Select
  End Do rows_loop

  Close(unit=iochannel(1))

End Function number_csv_rows

Function number_csv_cols(filename) Result(cols)
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
  Use error_messages
  Use global_constants
  Implicit None

  Integer :: i,ios,cols
  Character :: ctmp
  Character(len=chr80) :: filename
  Character(len=input_line_length) :: input_line

  Open(unit=iochannel(1),file=Trim(filename),status='old',iostat=ios)
  If (ios /= 0) Call error_open_file(filename)

  cols=0
  Read(iochannel(1),'(a)') input_line
  If (len_Trim(input_line) > input_line_length - 1) Then
     Write(*,fmt='(a,a)')'ERROR: Too long line(s) in file ',Trim(filename)
     Stop
  End If
  Do i=1,len_Trim(input_line)
    Read(input_line(i:i),'(a)',iostat=ios) ctmp
    If (ios /= 0) Call error_read_file(filename)
    If (ctmp == ',') Then
      cols=cols+1
    End If
  End Do
  cols=cols+1

  Close(unit=iochannel(1))

End Function number_csv_cols

Function csv_line_column(input_line,icol) Result(column)
!----------------------------------------
! Return the column as a character string
!----------------------------------------
  Use error_messages
  Use global_constants
  Implicit None

  Integer :: i, icol, cols, i_start_col, i_end_col
  Character(len=input_line_length) :: input_line
  Character(len=chr20) :: column
  Character :: ctmp

  cols = 1; i_start_col = 1; i_end_col = 0

  Do i = 1, Len_Trim(input_line)
     Read(input_line(i:i),'(a)') ctmp
     If (ctmp == ',') Then
        If (cols == icol) i_end_col = i - 1
        If (cols == icol -1) i_start_col = i + 1
        cols = cols + 1
     End If
  End Do
  If (i_end_col == 0 ) i_end_col = Len_Trim(input_line)
  
  column = input_line(i_start_col:i_start_col+Min(chr20,i_end_col-i_start_col))
  
End Function csv_line_column

Subroutine read_cfast_input_file()
  !--------------------------------------------------
  ! Read the CFAST input file: TARGET, MATL and COMPA
  ! Find: emissivities, xyz of the targets
  !--------------------------------------------------
  Use error_messages
  Use global_constants
  Use global_variables
  Use string_handling
  Use cfast_arrays
  Implicit None

  Integer :: ios, j, k
  Integer :: n_trg, n_mat, n_com
  ! New Cfast input file format is similar to FDS namelists
  !   TARGET => &DEVC
  !   MATL   => &MATL
  !   COMPA  => &COMP
  Integer :: n_trg_newfmt, n_mat_newfmt, n_com_newfmt
  Logical :: file_exists
  Character(len=chr80) :: filename
  Character(len=input_line_length) :: input_line
  Character(len=chr20) :: input_column
  
  filename = Trim(fds_input_file) // '.in'
  cfast_old_format = .False. ! Default is the new input file format
  Inquire(file=filename, exist=file_exists)
  If (file_exists) Then
     Open(unit=iochannel(1),file=Trim(filename), status='old', iostat=ios)
     If (ios /= 0) Call error_open_file(filename)
     n_trg = 0; n_mat = 0; n_com = 0
     n_trg_newfmt = 0; n_mat_newfmt = 0; n_com_newfmt = 0
     row_loop: Do
        Read(iochannel(1),*,iostat=ios) input_line
        Select Case(ios)
        Case (:-1)
           Exit row_loop
        Case (0)
           If (input_line(1:6) == 'TARGET') n_trg = n_trg + 1
           If (input_line(1:4) == 'MATL')   n_mat = n_mat + 1
           If (input_line(1:5) == 'COMPA')  n_com = n_com + 1
           If (input_line(1:6) == '&DEVC ') n_trg_newfmt = n_trg_newfmt + 1
           If (input_line(1:6) == '&MATL ') n_mat_newfmt = n_mat_newfmt + 1
           If (input_line(1:6) == '&COMP ') n_com_newfmt = n_com_newfmt + 1
        Case (1:)
           Call error_read_file(filename)
        End Select
     End Do row_loop
     Rewind(iochannel(1))
  Else
     Write(*,'(a,a,a)') 'ERROR: CFAST input file ',Trim(filename),' is not found'
     Stop
  End If

  If (n_trg + n_mat + n_com > 0) Then
     ! Old Cfast input file format
     cfast_old_format = .True.
     If (n_trg_newfmt + n_mat_newfmt + n_com_newfmt > 0) Then
        Write(*,'(a,a,a)') 'ERROR: CFAST input file ',Trim(filename), &
             ' has both old and new input file format lines.'
        Stop
     End If
     Write(*,'(t3,a)') 'Old Cfast input file format detected'
  Else
     ! New Cfast input ile format
     Write(*,'(t3,a)') 'New Cfast input file format detected'
     n_trg = n_trg_newfmt
     n_mat = n_mat_newfmt
     n_com = n_com_newfmt
  End If
  
  Allocate(cfast_target_xyz(  3,n_trg),stat=ios); Call error_allocate(ios)
  Allocate(cfast_target_name(   n_trg),stat=ios); Call error_allocate(ios)
  Allocate(cfast_target_epsilon(n_trg),stat=ios); Call error_allocate(ios)

  Allocate(cfast_target_matl(    n_trg),stat=ios); Call error_allocate(ios)
  Allocate(cfast_target_comp(    n_trg),stat=ios); Call error_allocate(ios)
  Allocate(cfast_target_ior(   3,n_trg),stat=ios); Call error_allocate(ios)
  Allocate(cfast_target_compname(n_trg),stat=ios); Call error_allocate(ios)

  Allocate(cfast_comp_xyz(3,n_com),stat=ios); Call error_allocate(ios)
  Allocate(cfast_comp_name( n_com),stat=ios); Call error_allocate(ios)

  Allocate(cfast_matl_epsilon(n_mat),stat=ios); Call error_allocate(ios)
  Allocate(cfast_matl_name(   n_mat),stat=ios); Call error_allocate(ios)

  cfast_target_xyz   = 0.0
  cfast_target_ior   = 0.0
  cfast_target_comp  = 0
  cfast_target_name  = ''
  cfast_target_matl  = ''
  cfast_target_compname = ''
  cfast_comp_name    = ''
  
  If (cfast_old_format) Then
     n_trg = 0; n_mat = 0; n_com = 0
     row_loop_2: Do
        Read(iochannel(1),fmt='(a)',iostat=ios) input_line
        Select Case(ios)
        Case (:-1)
           Exit row_loop_2
        Case (0)
           If (input_line(1:6) == 'TARGET') Then
              !   1  2   3   4     5  6  7  8    9     10       11  12   13
              !     cmp  x   y     z  nx ny nz matl_id                    id
              !TARGET,1,2.2,1.88,2.34,0, 0, 1,CONCRETE,EXPLICIT,PDE,0.5,Targ 1
              n_trg = n_trg + 1
              input_column = csv_line_column(input_line,2) ; Read(input_column,*) cfast_target_comp( n_trg)
              input_column = csv_line_column(input_line,3) ; Read(input_column,*) cfast_target_xyz(1,n_trg)
              input_column = csv_line_column(input_line,4) ; Read(input_column,*) cfast_target_xyz(2,n_trg)
              input_column = csv_line_column(input_line,5) ; Read(input_column,*) cfast_target_xyz(3,n_trg)
              input_column = csv_line_column(input_line,6) ; Read(input_column,*) cfast_target_ior(1,n_trg) ! nx
              input_column = csv_line_column(input_line,7) ; Read(input_column,*) cfast_target_ior(2,n_trg) ! ny
              input_column = csv_line_column(input_line,8) ; Read(input_column,*) cfast_target_ior(3,n_trg) ! nz
              input_column = csv_line_column(input_line,9) ; Read(input_column,fmt='(a)') cfast_target_matl(n_trg)
              input_column = csv_line_column(input_line,13); Read(input_column,fmt='(a)') cfast_target_name(n_trg)
           End If
           If (input_line(1:4) == 'MATL') Then
              !  1   2        3 4    5    6     7   
              !     matl_id   k cp   rho  d    eps
              !MATL,STEEL1/8,48,559,7854,0.003,0.9,"Steel, Plain Carbon (1/8 in)"
              n_mat = n_mat + 1
              input_column = csv_line_column(input_line,7) ; Read(input_column,*) cfast_matl_epsilon( n_mat)
              input_column = csv_line_column(input_line,2) ; Read(input_column,fmt='(a)') cfast_matl_name( n_mat)
           End If
           If (input_line(1:5) == 'COMPA') Then
              !  1     2     3   4  5   6 7 8
              !       id     w   d   h  x y z  ceil     floor     wall 
              !COMPA,Comp 1,3.6,2.4,2.4,0,0,0,CONCRETE,CONCRETE,CONCRETE,50,50,50
              n_com = n_com + 1
              input_column = csv_line_column(input_line,6) ; Read(input_column,*) cfast_comp_xyz(1,n_com)
              input_column = csv_line_column(input_line,7) ; Read(input_column,*) cfast_comp_xyz(2,n_com)
              input_column = csv_line_column(input_line,8) ; Read(input_column,*) cfast_comp_xyz(3,n_com)
           End If
        Case (1:)
           Call error_read_file(filename)
        End Select
     End Do row_loop_2
  Else ! new Cfast input file format
     Call parse_cfast_head_namelist()
     Write(*,'(t3,a,i4)') 'Cfast version: ',cfast_version
     Write(*,'(t3,a,a)') 'Cfast title: ',Trim(cfast_title)
     Call parse_cfast_devc_namelist()
!!$     Do j=1,n_trg
!!$        Write(*,'(t3,a,a)') 'trg: ',Trim(cfast_target_name(j))
!!$        Write(*,'(t3,a,a)') 'mat: ',Trim(cfast_target_matl(j))
!!$        Write(*,'(t3,a,a)') 'cmp: ',Trim(cfast_target_compname(j))
!!$        Write(*,'(t3,a,3f12.6)') 'xyz: ',cfast_target_xyz(1,j),cfast_target_xyz(2,j),cfast_target_xyz(3,j)
!!$        Write(*,'(t3,a,3f12.6)') 'ior: ',cfast_target_ior(1,j),cfast_target_ior(2,j),cfast_target_ior(3,j)
!!$     End do
     Call parse_cfast_matl_namelist()
!!$     Do j=1,n_mat
!!$        Write(*,'(t3,a,a,f12.6)') 'mat, eps: ',Trim(cfast_matl_name(j)), cfast_matl_epsilon(j)
!!$     End do
     Call parse_cfast_comp_namelist()
!!$     Do j=1,n_com
!!$        Write(*,'(t3,a,i4,a,3f12.6)') 'com ',n_com,' xyz: ',cfast_comp_xyz(1,j),cfast_comp_xyz(2,j),cfast_comp_xyz(3,j)
!!$     End Do
     Do j = 1, Ubound(cfast_target_compname,1)
        Do k = 1, Ubound(cfast_comp_name,1)
           If (Trim(cfast_target_compname(j)) == Trim(cfast_comp_name(k))) Then
              cfast_target_comp(j) = k
              Write(*,'(t3,a,a,a,a)') 'trg ',Trim(cfast_target_name(j)), ' is in comp ',cfast_comp_name(k)
           End If
        End Do
        If (cfast_target_comp(j) == 0) Then
           Write(*,'(4(a))') 'ERROR: Target ',Trim(cfast_target_name(j)), ': No comp found: ',&
                Trim(cfast_target_compname(j))
           Call error_read_file(filename)
        End If
     End Do
  End If

  ! Fill the target input data arrays: emissivity and xyz=xyz_comp+xyz_targ
  cfast_target_epsilon = emissivity ! Default
  Do j = 1, Ubound(cfast_target_name,1)
     k = cfast_target_comp(j)
     cfast_target_xyz(:,j) = cfast_target_xyz(:,j) + cfast_comp_xyz(:,k)
     Do k = 1, Ubound(cfast_matl_epsilon,1)
        If(Trim(cfast_target_matl(j)) == Trim(cfast_matl_name(k))) &
             cfast_target_epsilon(j) = cfast_matl_epsilon(k)
     End Do
  End Do

  Deallocate(cfast_target_matl)
  Deallocate(cfast_target_comp)
  Deallocate(cfast_comp_xyz)
  Deallocate(cfast_matl_epsilon)
  Deallocate(cfast_matl_name)
  
End Subroutine read_cfast_input_file

Subroutine parse_cfast_head_namelist()
!---------------------------------------------
! Parse CFAST input file for TITLE and VERSION
!
! Assigns values to the variables
!   cfast_chid,cfast_title
!---------------------------------------------
  Use error_messages
  Use cfast_arrays
  Use global_constants  
  Use global_variables, only : fds_input_file
  Use string_handling
  Implicit None

  Integer :: ios,nhead
  Character(len=chr80) :: filename

  ! HEAD-namelist group
  Namelist /head/ version,title
  Integer :: version
  Character(128) :: title

  filename = Trim(fds_input_file) // '.in'
  version=0; title=''
  Open(unit=iochannel(1),file=Trim(filename),status='old',iostat=ios) 
  If (ios /= 0) Call error_open_file(filename)

  ! Read in HEAD-namelist records
  nhead=0
  head_loop: Do
     Read(iochannel(1),nml=head,iostat=ios)
     Select Case(ios)
     Case (:-1)
        Exit head_loop
     Case (0)
        cfast_version=version
        cfast_title=Trim(title)
        nhead=nhead+1
     Case(1:)
        Write(*,'(3(a))') 'ERROR: in reading HEAD-namelist record (file ', &
             Trim(quote(filename)), ')'
        Call error_read_file(filename)
     End Select
  End Do head_loop

  Close(unit=iochannel(1))

  If (nhead == 0) Then
     Write(*,'(2(a))') 'ERROR: no HEAD-namelist record found in file ', &
          Trim(quote(filename))
     Stop
  End If

  If (nhead > 1) Then
     Write(*,'(2(a))') 'WARNING: multiple HEAD-namelist records found in file ', &
          Trim(quote(filename))
  End If

End Subroutine parse_cfast_head_namelist

Subroutine parse_cfast_devc_namelist()
!---------------------------------------------
! Parse CFAST input file for DEVC namelist (targets)
!
! Assigns values to the variables
!   cfast_target_xyz,_ior,_matl,_name
!---------------------------------------------
  Use error_messages
  Use cfast_arrays
  Use global_constants  
  Use global_variables, only : fds_input_file
  Use string_handling
  Implicit None

  Integer :: ios,ntrg
  Character(len=chr80) :: filename

  ! DEVC-namelist group
  real(kind=eb) :: temperature_depth,rti,setpoint,spray_density
  real(kind=eb),dimension(3) :: location,normal
  real(kind=eb),dimension(2) :: setpoints
  character(64) :: comp_id,id,matl_id
  character(64) :: type
  logical :: adiabatic_target
  real(kind=eb), dimension(2) :: convection_coefficients
  namelist /devc/ comp_id, type, id, temperature_depth, location, matl_id, normal, rti, &
       setpoint, spray_density, setpoints, adiabatic_target, convection_coefficients

  filename = Trim(fds_input_file) // '.in'
  Open(unit=iochannel(1),file=Trim(filename),status='old',iostat=ios) 
  If (ios /= 0) Call error_open_file(filename)

  ! Read in DEVC-namelist records
  comp_id=''; location=0.0; normal=0.0; matl_id=''; id=''
  ntrg=0
  head_loop: Do
     Read(iochannel(1),nml=devc,iostat=ios)
     Select Case(ios)
     Case (:-1)
        Exit head_loop
     Case (0)
        ntrg=ntrg+1
        cfast_target_xyz(1,ntrg) = location(1)
        cfast_target_xyz(2,ntrg) = location(2)
        cfast_target_xyz(3,ntrg) = location(3)
        cfast_target_ior(1,ntrg) = normal(1)
        cfast_target_ior(2,ntrg) = normal(2)
        cfast_target_ior(3,ntrg) = normal(3)
        cfast_target_matl(ntrg)  = Trim(matl_id)
        cfast_target_name(ntrg) = Trim(id)
        cfast_target_compname(ntrg) = Trim(comp_id)
     Case(1:)
        Write(*,'(3(a))') 'ERROR: in reading DEVC-namelist record (file ', &
             Trim(quote(filename)), ')'
        Call error_read_file(filename)
     End Select
  End Do head_loop

  Close(unit=iochannel(1))

  If (ntrg == 0) Then
     Write(*,'(2(a))') 'ERROR: no DEVC-namelist record found in file ', &
          Trim(quote(filename))
     Stop
  End If

End Subroutine parse_cfast_devc_namelist

Subroutine parse_cfast_matl_namelist()
!---------------------------------------------
! Parse CFAST input file for MATL namelist
!
! Assigns values to the variables
!   cfast_matl_epsilon,_name
!---------------------------------------------
  Use error_messages
  Use cfast_arrays
  Use global_constants  
  Use global_variables, only : fds_input_file
  Use string_handling
  Implicit None

  Integer :: ios,nmatl
  Character(len=chr80) :: filename

  ! MATL-namelist group
  Character(64) :: id, material
  Real(kind=eb) :: conductivity, density, emissivity, specific_heat, thickness
  Namelist /matl/ conductivity, density, emissivity, id, material, specific_heat, thickness

  filename = Trim(fds_input_file) // '.in'
  Open(unit=iochannel(1),file=Trim(filename),status='old',iostat=ios) 
  If (ios /= 0) Call error_open_file(filename)

  ! Read in MATL-namelist records
  nmatl=0
  matl_loop: Do
     Read(iochannel(1),nml=matl,iostat=ios)
     Select Case(ios)
     Case (:-1)
        Exit matl_loop
     Case (0)
        nmatl=nmatl+1
        cfast_matl_epsilon(nmatl) = emissivity
        cfast_matl_name(   nmatl) = Trim(id)
     Case(1:)
        Write(*,'(3(a))') 'ERROR: in reading MATL-namelist record (file ', &
             Trim(quote(filename)), ')'
        Call error_read_file(filename)
     End Select
  End Do matl_loop

  Close(unit=iochannel(1))

  If (nmatl == 0) Then
     Write(*,'(2(a))') 'ERROR: no MATL-namelist record found in file ', &
          Trim(quote(filename))
     Stop
  End If

End Subroutine parse_cfast_matl_namelist

Subroutine parse_cfast_comp_namelist()
!---------------------------------------------
! Parse CFAST input file for COMP namelist
!
! Assigns values to the variables
!   cfast_comp_xyz,_name
!---------------------------------------------
  Use error_messages
  Use cfast_arrays
  Use global_constants  
  Use global_variables, only : fds_input_file
  Use string_handling
  Implicit None

  Integer :: ios,ncomp
  Character(len=chr80) :: filename

  ! COMP-namelist group
  integer,dimension(3) :: grid
  real(eb) :: depth, height ,width
  real(eb),dimension(3) :: origin
  real(eb), dimension(200) :: cross_sect_areas, cross_sect_heights
  ! integer, parameter :: mxpts = 200 ! maximum number of data points in a input curve/ramp ! Cfast source
  ! real(eb), dimension(mxpts) :: cross_sect_areas, cross_sect_heights ! Cfast source
  logical :: hall, shaft
  character(64) :: id, ceiling_matl_id, floor_matl_id, wall_matl_id
  Namelist /comp/ cross_sect_areas, cross_sect_heights, depth, grid, hall, height, id, &
       ceiling_matl_id, floor_matl_id, wall_matl_id, origin, shaft, width

  filename = Trim(fds_input_file) // '.in'
  Open(unit=iochannel(1),file=Trim(filename),status='old',iostat=ios) 
  If (ios /= 0) Call error_open_file(filename)

  ! Read in COMP-namelist records
  ncomp=0
  comp_loop: Do
     Read(iochannel(1),nml=comp,iostat=ios)
     Select Case(ios)
     Case (:-1)
        Exit comp_loop
     Case (0)
        ncomp=ncomp+1
        cfast_comp_xyz(1,ncomp) = origin(1)
        cfast_comp_xyz(2,ncomp) = origin(2)
        cfast_comp_xyz(3,ncomp) = origin(3)
        cfast_comp_name( ncomp) = Trim(id)
     Case(1:)
        Write(*,'(3(a))') 'ERROR: in reading COMP-namelist record (file ', &
             Trim(quote(filename)), ')'
        Call error_read_file(filename)
     End Select
  End Do comp_loop

  Close(unit=iochannel(1))

  If (ncomp == 0) Then
     Write(*,'(2(a))') 'ERROR: no COMP-namelist record found in file ', &
          Trim(quote(filename))
     Stop
  End If

End Subroutine parse_cfast_comp_namelist

Subroutine import_cfast_data()
!-----------------------------------------------------------------------------
! Import data from the CFAST "_w.csv" file that contains target output
! Calculate T_ast if that is needed
!  
!  Columns: Old format   (cfast_old_format=T), validation output
!    1st column:    Time(s)
!    2nd-4*Ncomp+1: compatrment boundaries temperatures (_w: wall temperatures)
!    then targets, 15 per one target 
!  
!  Columns: New format, normal output
!    1st column:    Time(s)
!    2nd-4*Ncomp+1: compatrment boundaries temperatures (_w: wall temperatures)
!    then targets, 9 per one target 
!TRGGAST_1       TRGSURT_1       TRGCENT_1       TRGFLXI_1       TRGFLXT_1       TRGFEDG_1       TRGDFEDG_1      TRGFEDH_1       TRGDFEDH_1
! Target Surrounding Gas Temperature      Target Surface Temperature      Target Center Temperature       Target Incident Flux
!                                                                Target Net Flux Target Gas FED  Target GasFED Increment Target Heat FED Target Heat FED Increment
!DoorLeft        DoorLeft        DoorLeft        DoorLeft        DoorLeft        DoorLeft        DoorLeft        DoorLeft        DoorLeft
!C               C               C               KW/m^2          KW/m^2                          
!  
!  Rows: Both old and new formats
!    1st row: TargetQuantityLabel_#, #=target running index
!    2nd row: Target quantity long name
!    3rd row: CompName / TargetName (user can give these in CFAST)
!    4th row: units
!    5th row- time series (first column is the time axis)
!  
!  NOTE: Now just 20 character long target labels are read in
!  NOTE: Now just max. 5096 characters long csv lines read in
!
!-----------------------------------------------------------------------------
  Use error_messages
  Use fds_devc_arrays
  Use global_constants
  Use global_variables
  Use mapping_arrays
  Use string_handling
  Use cfast_arrays
  Implicit None

  Integer :: i,j,k,ios,n_targets,ibegin,iend,column,ndevc_files
  Character(len=chr80) :: filename
  Real(kind=rk) a_ast, b_ast, c_ast, eps_rflux_inc, m_ast, t_ast, alp_ast, bet_ast, gam_ast
  
  Character :: ctmp
  Character(len=input_line_length) :: input_line
  Character(len=chr20) :: quantity_label
  Logical :: cfast_valid_output

  
  Character(len=chr20), Dimension(:), Allocatable :: devc_name_tmp
  Character(len=chr20), Dimension(:), Allocatable :: devc_unit_tmp
  Character(len=chr20), Dimension(:), Allocatable :: devc_label_tmp

  Real(kind=rk), Dimension(:,:), Allocatable :: devc_data_tmp 
  Real(kind=rk), Dimension(:,:), Allocatable :: fds_devc_rflux
  Real(kind=rk), Dimension(:,:), Allocatable :: fds_devc_rloss

  filename = find_single_cfast_file(fds_input_file)  ! fds_input_file = Cfast case name
  If (Trim(filename)=='') Then
     Write(*,'(a,a,a)') 'ERROR: CFAST ',Trim(fds_input_file),'_w.csv file not found'
     Stop
  End If

  ndevc_files = 1 ! Just one _w.csv file for CFAST
  Allocate(fds_devc_file(ndevc_files), stat=ios)
  Call error_allocate(ios)
  fds_devc_file(1) = Trim(filename)
  
  Allocate(fds_devc_rows(ndevc_files),stat=ios); Call error_allocate(ios)
  Allocate(fds_devc_cols(ndevc_files),stat=ios); Call error_allocate(ios)

  Do i = 1, ndevc_files
    fds_devc_rows(i) = number_csv_rows(fds_devc_file(i)) - 4 ! 4 header rows excluded
    fds_devc_cols(i) = number_csv_cols(fds_devc_file(i))
  End Do
   
  ! Error handling
  If (Sum(fds_devc_cols(1:ndevc_files))-ndevc_files == 0) Then
     Write(*,'(2(a))') 'ERROR: no records found in file ', Trim(quote(fds_devc_file(1)))
     Stop
  End If

  If (Trim(transfer_quantity) == 'wall_temperature') Then
     quantity_label = 'TRGSURT' ! New format: 2nd value
  Else If (Trim(transfer_quantity) == 'net_heat_flux') Then
     quantity_label = 'TRGFLXT' ! New format: 5th value
  Else If (Trim(transfer_quantity) == 'adiabatic_surface_temperature') Then
     quantity_label = 'TRGGAST'  ! Gas temp at target location, New format: 1st value
  End If

  ! Read in data from one or more CFAST _w.csf files (now just one file)
  csv_file_loop: Do i = 1, ndevc_files
     Open(unit=iochannel(1),file=Trim(fds_devc_file(i)),status='old',iostat=ios)
     If (ios /= 0) Call error_open_file(fds_devc_file(i))

     ! Read labels
     Allocate(devc_label_tmp(fds_devc_cols(i)),stat=ios); Call error_allocate(ios)
     column = 1 ; devc_label_tmp=''; j = 1
     Read(iochannel(1),'(a)',iostat=ios) input_line
     If (ios /= 0) Call error_read_file(fds_devc_file(i))
     If (len_Trim(input_line) > input_line_length - 1) Then
        Write(*,*)'ERROR: Too long header lines in CFAST target csv file'
        Stop
     End If
     label_loop: Do k = 1, len_Trim(input_line)
        Read(input_line(k:k),'(a)',iostat=ios) ctmp
        If (ctmp == ',') Then
           column = column + 1 ; j = 1
           Cycle label_loop
        End If
        ! Omit quotation marks
        If (ctmp /= Char(34) .And. ctmp /= Char(39)) Then
           If (j <= chr20) Then
              devc_label_tmp(column)(j:j) = input_line(k:k)
           End If
           j = j + 1  
        End If
     End Do label_loop

     ! Read quantity line (not used
     Read(iochannel(1),'(a)',iostat=ios) input_line
     If (ios /= 0) Call error_read_file(fds_devc_file(i))

     ! Read Target names
     Allocate(devc_name_tmp(fds_devc_cols(i)),stat=ios); Call error_allocate(ios)
     column = 1 ; devc_name_tmp=''; j = 1
     Read(iochannel(1),'(a)',iostat=ios) input_line
     If (len_Trim(input_line) > input_line_length - 1) Then
        Write(*,*)'ERROR: Too long header lines in CFAST target csv file'
        Stop
     End If
     If (ios /= 0) Call error_read_file(fds_devc_file(i))
     ids_loop: Do k = 1, len_Trim(input_line)
        Read(input_line(k:k),'(a)',iostat=ios) ctmp
        If (ctmp == ',') Then
           column = column + 1 ; j = 1
           Cycle ids_loop
        End If
        If (j <= chr20) Then
           devc_name_tmp(column)(j:j) = input_line(k:k)
        End If
        j = j + 1  
     End Do ids_loop

     ! Read units
     Allocate(devc_unit_tmp(fds_devc_cols(i)),stat=ios); Call error_allocate(ios)
     column=1; devc_unit_tmp=''; j=1
     Read(iochannel(1),'(a)',iostat=ios) input_line
     If (len_Trim(input_line) > input_line_length - 1) Then
        Write(*,*)'ERROR: Too long header lines in CFAST target csv file'
        Stop
     End If
     If (ios /= 0) Call error_read_file(fds_devc_file(i))
     units_loop: Do k = 1, len_Trim(input_line)
        Read(input_line(k:k),'(a)',iostat=ios) ctmp
        If (ctmp == ',') Then
           column = column + 1 ; j = 1
           Cycle units_loop
        End If
        If (j <= chr20) Then
           devc_unit_tmp(column)(j:j) = input_line(k:k)
        End If
        j = j + 1  
     End Do units_loop

     ! ---> Exception handling <---
     If (Trim(devc_unit_tmp(1)) /= 's') Then
        Write(*,'(2(a))') 'WARNING: non-standard time units in CFAST csv file ', &
             Trim(quote(fds_devc_file(i)))
     End If
     
     ! Time series
     Allocate(devc_data_tmp(fds_devc_rows(i),fds_devc_cols(i)),stat=ios); Call error_allocate(ios)
     Do j = 1, fds_devc_rows(i)
        Read(iochannel(1),*) (devc_data_tmp(j,k),k=1,fds_devc_cols(i))
     End Do

     ! ---------------------------- 
     ! Count the number of different CFAST targets, use _w.csv file information
     n_targets = 0
     cfast_valid_output = .FALSE.
     Do j = 2, fds_devc_cols(1)
        If (devc_label_tmp(j)(1:7) == 'TRGGAST') Then
           n_targets = n_targets + 1
        End If
        If (devc_label_tmp(j)(1:7) == 'TRGFLXR') Then
           cfast_valid_output = .TRUE.
        End If
     End Do
     If (cfast_valid_output)  Then
        Write(*,'(t3,a,a)') 'Cfast validation output detected, file: ',Trim(filename)
     Else
        Write(*,'(t3,a,a)') 'Cfast normal output detected, file: ',Trim(filename)
     End If
     
     If (n_targets /= Ubound(cfast_target_name,1)) Then
        Write(*,*)'ERROR: Different number of targets in CFAST input file vs csv file'
        Stop
     End If
     
     Allocate(fds_devc_name_b(n_targets),               stat=ios); Call error_allocate(ios) ! label
     Allocate(fds_devc_qnty(n_targets),                 stat=ios); Call error_allocate(ios) ! quantity
     Allocate(fds_devc_name(n_targets),                 stat=ios); Call error_allocate(ios) ! target name
     Allocate(fds_devc_unit(n_targets),                 stat=ios); Call error_allocate(ios) ! units
     Allocate(fds_devc_time(fds_devc_rows(1)),          stat=ios); Call error_allocate(ios)
     Allocate(fds_devc_data(fds_devc_rows(1),n_targets),stat=ios); Call error_allocate(ios)
     If (Trim(transfer_quantity) == 'adiabatic_surface_temperature') Then
        Allocate(fds_devc_rflux(fds_devc_rows(1),n_targets),stat=ios); Call error_allocate(ios)
        Allocate(fds_devc_rloss(fds_devc_rows(1),n_targets),stat=ios); Call error_allocate(ios)
     End If

     If (i == 1) Then
        fds_devc_time = devc_data_tmp(1:fds_devc_rows(1),1)
        k = 0
        Do j = 2, fds_devc_cols(1)
           If (devc_label_tmp(j)(1:8) == Trim(Trim(quantity_label)//'_')) Then
              k = k + 1  ! Target index
              fds_devc_unit(k)   = devc_unit_tmp(j)
              fds_devc_name_b(k) = devc_label_tmp(j) ! label
              fds_devc_name(k)   = devc_name_tmp(j)  ! target name
              fds_devc_data(1:fds_devc_rows(1),k) = devc_data_tmp(1:fds_devc_rows(1),j)
              If (Trim(transfer_quantity) == 'adiabatic_surface_temperature' .and. cfast_valid_output) Then
                 fds_devc_rflux(1:fds_devc_rows(1),k) = devc_data_tmp(1:fds_devc_rows(1),j+5)
                 fds_devc_rloss(1:fds_devc_rows(1),k) = devc_data_tmp(1:fds_devc_rows(1),j+10)
                 ! Old format: (actually, validation output)
                 ! T_gas                                            rad.flux                                          rad.loss
                 !   1         2        3         4         5          6        7         8         9         10        11
                 !TRGGAST_1,TRGSURT_1,TRGCENT_1,TRGFLXI_1,TRGFLXT_1,TRGFLXR_1,TRGFLXC_1,TRGFLXF_1,TRGFLXS_1,TRGFLXG_1,TRGFLXRE_1
              End If
              If (Trim(transfer_quantity) == 'adiabatic_surface_temperature' .And. .Not.cfast_valid_output) Then
                 fds_devc_rflux(1:fds_devc_rows(1),k) = devc_data_tmp(1:fds_devc_rows(1),j+3)
                 fds_devc_rloss(1:fds_devc_rows(1),k) = 0.0
                 ! New format: (actually, normal output, not validation output)
                 ! T_gas                        inc.rad.flux
                 !   1         2        3         4         5          6        7         8         9 
                 ! TRGGAST_1 TRGSURT_1 TRGCENT_1 TRGFLXI_1 TRGFLXT_1 TRGFEDG_1 TRGDFEDG_1 TRGFEDH_1 TRGDFEDH_1
              End If
              If (Trim(fds_devc_name(k)) /= Trim(cfast_target_name(k))) Then
                 Write(*,*)'ERROR: Different name of targets in CFAST input file vs csv file'
                 Write(*,*)'ERROR: ',Trim(fds_devc_name(k)),' ',Trim(cfast_target_name(k))
                 Stop
              End If

           End If
        End Do
     Else
        Write(*,'(a)')'ERROR: CFAST Target input can use only one input _w.csv file'
        Stop
        ! Omit the time vector
        ibegin = 2 - i + Sum(fds_devc_cols(1:i-1)) ; iend = ibegin + fds_devc_cols(i) - 2
        fds_devc_unit(ibegin:iend)   = devc_unit_tmp(2:fds_devc_cols(i))
        fds_devc_name_b(ibegin:iend) = devc_label_tmp(2:fds_devc_cols(i))
        fds_devc_name(ibegin:iend)   = devc_name_tmp(2:fds_devc_cols(i))
        fds_devc_data(1:fds_devc_rows(i),ibegin:iend) = devc_data_tmp(1:fds_devc_rows(i),2:fds_devc_cols(i))
     End If

     Deallocate(devc_data_tmp,stat=ios); Call error_allocate(ios)
     Deallocate(devc_unit_tmp,stat=ios); Call error_allocate(ios)
     Deallocate(devc_name_tmp,stat=ios); Call error_allocate(ios)
     Deallocate(devc_label_tmp,stat=ios); Call error_allocate(ios)

     Close(unit=iochannel(1))

     ast_output_if: If (Trim(transfer_quantity) == 'adiabatic_surface_temperature') Then
        ast_target_loop: Do j = 1, n_targets
           ! Stefan-Boltzmann https://physics.nist.gov/cuu/Constants/, Source: 2014 CODATA 5.670367 E-8
           fds_devc_rloss(1,j) = -0.001*cfast_target_epsilon(j)*5.670367E-8*((fds_devc_data(1,j)+273.15)**4) ! t=0 correction, kW/m2
           ast_row_loop: Do k = 1, fds_devc_rows(1)
              ! New format: Target Incident flux (=rad.inc.flux), so this is easy
              ! Below old format: net rad flux and rad loss are in the _w.csv file (validation output):
              eps_rflux_inc = 1000.0*(fds_devc_rflux(k,j) - fds_devc_rloss(k,j)) ! CFAST rloss is negative
              a_ast   = cfast_target_epsilon(j)*5.670367E-8
              b_ast   = hcoeff ! W/m2.K, default, no target specific information available
              c_ast   = -eps_rflux_inc - hcoeff*(fds_devc_data(k,j)+273.15)  ! Kelvin
              alp_ast = (Sqrt(3.0)*Sqrt(27.0*a_ast**2*b_ast**4 - 256.0*a_ast**3*c_ast**3) + 9.0*a_ast*b_ast**2)**(1.0/3.0)
              bet_ast = 4.0*((2.0/3.0)**(1.0/3.0))*c_ast
              gam_ast = (18.0)**(1.0/3.0)*a_ast
              m_ast   = Sqrt(bet_ast/alp_ast + alp_ast/gam_ast)
              t_ast   = 0.5*(-m_ast + Sqrt((2.0*b_ast/(a_ast*m_ast)) - m_ast**2)) - 273.15 ! Celsius
              fds_devc_data(k,j) = t_ast  ! Replase gas temp with ast (celcius)
           End Do ast_row_loop
        End Do ast_target_loop
     End If ast_output_if
        
  End Do csv_file_loop

End Subroutine import_cfast_data

Subroutine filter_cfast_data
!-------------------------------
! Filter unwanted CFAST csv data
!-------------------------------
  Use error_messages
  Use fds_devc_arrays
  Use fds_head_arrays
  Use global_constants
  Use global_variables
  Use mapping_arrays
  Use string_handling
  Use cfast_arrays
  Implicit None

  Integer :: i,j,ios,idevc,ndevc,inode,nrows,ncols
  Integer :: ndevc_user, quantity_label_len

  Integer :: devc_number
  Logical :: devc_found

  Logical, Dimension(:), Allocatable :: devc_mask_user

  Character(len=chr80) :: fds_transfer_quantity
  Character(len=chr20) :: quantity_label

  !---------------------------
  ! Assign a transfer quantity
  !---------------------------

  If (Trim(transfer_quantity) == 'wall_temperature') Then
     fds_transfer_quantity='wall temperature'
     quantity_label = 'TRGSURT'
     quantity_label_len = 7
  Else If (Trim(transfer_quantity) == 'net_heat_flux') Then
     fds_transfer_quantity='net heat flux'
     quantity_label = 'TRGFLXT'
     quantity_label_len = 7
  Else If (Trim(transfer_quantity) == 'adiabatic_surface_temperature') Then
     fds_transfer_quantity='adiabatic surface temperature'
     quantity_label = 'TRGSURT'
     quantity_label_len = 7
  End If

  !-----------------------------------------------
  ! Convert connectivity table into numerical form
  !-----------------------------------------------

  If (nset_connectivity) Then

    nrows = Ubound(connectivity_table,1); ncols = Ubound(connectivity_table,2)
    Allocate(connectivity_table_num(nrows,ncols),stat=ios); Call error_allocate(ios)
    Allocate(connectivity_phys_cons(nrows,3),stat=ios); Call error_allocate(ios)
    Allocate(time_shift_mask(nrows),         stat=ios); Call error_allocate(ios)

    time_shift_mask = .False.
    connectivity_table_num = 0 
    connectivity_phys_cons(:,1)= emissivity ! Default
    connectivity_phys_cons(:,2)= hcoeff     ! Default
    connectivity_phys_cons(:,3)= 0.0        ! Default time shift (s), zero for CFAST
    row_loop_0: Do i = 1, Ubound(connectivity_table,1)
       col_loop_0: Do j = 2, Ubound(connectivity_table,2)
          If (connectivity_table(i,j) == '') Cycle col_loop_0
          If (numerical(connectivity_table(i,j))) Then
             devc_number = str2int(connectivity_table(i,j))
             connectivity_table_num(i,j) = devc_number

             ! Exception handling
             If (devc_number > Ubound(fds_devc_name,1) .Or. devc_number < 1) Then
                Write(*,'(3(a))') 'ERROR: in NSET-Target connectivity table: Target ', &
                     Trim(int2str(devc_number)), ' does not exist, target numbers'
                Stop
             End If
        
          Else
             devc_found = .False.
             devc_loop_0: Do idevc = 1, Ubound(fds_devc_name,1)
                If (Trim(connectivity_table(i,j)) == Trim(fds_devc_name(idevc))) Then
                   connectivity_table_num(i,j) = idevc  ! Target index
                   devc_found = .True.
                End If
             End Do devc_loop_0

             ! Exception handling
             If (.Not. devc_found) Then
                Write(*,'(3(a))') 'ERROR: in NSET-Target connectivity table: Target ', &
                     Trim(quote(connectivity_table(i,j))), ' does not exist, target names'
                Stop
             End If

          End If

      End Do col_loop_0
    End Do row_loop_0

 Else
    Write(*,'(a)')'ERROR: CFAST input needs nset connectivity mode'
    Stop
 End If

  !------------------------------------
  ! Generate a list of selected devices
  !------------------------------------

  ndevc = Ubound(fds_devc_name,1) ! Number of targets in the Cfast _w.csv file
  Allocate(devc_mask_user(ndevc),stat=ios); Call error_allocate(ios)

  If (nset_connectivity) Then

     ! Based on user selection: Are the targets found in the nset file?
     ndevc_user = 0; devc_mask_user = .False.
     devc_loop_2: Do idevc = 1, ndevc

        row_loop_2: Do i = 1, Ubound(connectivity_table,1)
           col_loop_2: Do j = 2, Ubound(connectivity_table,2)
              If (connectivity_table(i,j) == '') Cycle col_loop_2
              If (numerical(connectivity_table(i,j))) Then
                 If (str2int(connectivity_table(i,j)) == idevc) Then
                    devc_mask_user(idevc) = .True. ! This target is used
                    ndevc_user = ndevc_user + 1
                    Cycle devc_loop_2
                 End If
              Else
                 If (Trim(connectivity_table(i,j)) == Trim(fds_devc_name(idevc))) Then
                    devc_mask_user(idevc) = .True.
                    ndevc_user = ndevc_user + 1 ! This target is used
                    Cycle devc_loop_2
                 End If
              End If
           End Do col_loop_2
        End Do row_loop_2
        
     End Do devc_loop_2
     
     ! If none are selected (this error message might be unreachable)
     If (ndevc_user == 0) Then
        Write(*,'(2(a))') 'ERROR: no selected CFAST targets found '
        Stop
     End If

  Else
     Write(*,'(a)')'ERROR: CFAST input and no NSET-Target connectivities given'
     Stop
  End If
 
  !-------------------------
  ! Generate FDS data arrays
  !-------------------------

  nnodes_fds = ndevc_user  ! Number of targets that are actually used in the nset file
  ntimes_fds = Ubound(fds_devc_time,1)
  If (nnodes_fds >= 1) Then
     fds_xyz_available=.False.; fds_data_available=.True.; fds_model_available=.False.
  Else
     fds_xyz_available=.False.; fds_data_available=.False.; fds_model_available=.False.
  End If

  Allocate(fds_id(nnodes_fds),stat=ios);              Call error_allocate(ios)
  Allocate(fds_time(ntimes_fds),stat=ios);            Call error_allocate(ios)
  Allocate(fds_idevc(nnodes_fds),stat=ios);           Call error_allocate(ios)
  Allocate(fds_data(ntimes_fds,nnodes_fds),stat=ios); Call error_allocate(ios)
  Allocate(fds_xyz(nnodes_fds,3),stat=ios);           Call error_allocate(ios)
  Allocate(fds_ior(nnodes_fds),stat=ios);             Call error_allocate(ios)
  fds_ior = 0
  
  !------------------------------
  ! Create arrays containing data
  !------------------------------

  inode = 1
  Do idevc = 1, ndevc  ! Loop over all targets in _w.csv file
    If (devc_mask_user(idevc)) Then
      fds_id(inode)                = fds_devc_name(idevc)  ! Target name
      fds_idevc(inode)             = idevc                 ! Running index of the target
      fds_data(1:ntimes_fds,inode) = fds_devc_data(1:ntimes_fds,idevc)
      fds_xyz(inode,:)             = cfast_target_xyz(:,idevc)
      ! Note: Next is not working generally, it assumes "FDS type" normal vector
      !       nx,ny,nz should be in the x,y, or z direction
      If (Abs(cfast_target_ior(1,idevc)) > 0.5) Then
         fds_ior(inode) = Nint(Sign(1.0_rk,cfast_target_ior(1,idevc)))
      Else If (Abs(cfast_target_ior(2,idevc)) > 0.5) Then
         fds_ior(inode) = Nint(Sign(2.0_rk,cfast_target_ior(2,idevc)))
      Else If (Abs(cfast_target_ior(3,idevc)) > 0.5) Then
         fds_ior(inode) = Nint(Sign(3.0_rk,cfast_target_ior(3,idevc)))
      End If
      inode = inode + 1
    End If  
  End Do
  fds_time = fds_devc_time

End Subroutine filter_cfast_data

Function find_single_cfast_file(chid) Result(filename)
!--------------------------------
! Find a single CFAST _w.csv file
!--------------------------------
  Use global_constants
  Use string_handling
  Implicit None

  Logical :: file_exists
  Character(len=chr40) :: chid
  Character(len=chr80) :: filename

  filename = Trim(chid) // '_w.csv'

  Inquire(file=filename,exist=file_exists)
  If (.Not. file_exists) filename=''

End Function find_single_cfast_file


!***********************************************************
! Subroutines for parsing the information from the Nset file
!***********************************************************

Subroutine import_iso_data()
!----------------------------------------------------------------------
! Import data from time-temperature curves
! (Actually: generate the time-temp curves using mathematical formulas)
!  
! 18.12.2017 implemented
!   iso_834:    ISO 834 curve, "standard curve"
!   ec1-1-2_hc: EC1-1-2 cl3.2.3 Hydrocarbon curve
!   ec1-1-2_ex: EC1-1-2 cl3.2.2 External fire curve
!   astm_e119:  ASTM E 119 curve ("US standard curve")
!----------------------------------------------------------------------
  Use error_messages
  Use fds_devc_arrays
  Use global_constants
  Use global_variables
  Use mapping_arrays
  Use string_handling
  Implicit None

  Integer :: i,j,ios,ndevc,ndevc_files
  Integer :: n_iso, n_hc, n_ex, n_e119
  Real(kind=rk) :: t_min

  ! Count the number of different fire curves
  ndevc_files=0
  n_iso  =0
  n_hc   =0
  n_ex   =0
  n_e119 =0
  Do j=1,Ubound(connectivity_table,1)
     Select Case(Trim(Lowercase(connectivity_table(j,2))))
     Case('iso_834')
        n_iso  = 1
     Case('ec1-1-2_hc')
        n_hc   = 1
     Case('ec1-1-2_ex')
        n_ex   = 1
     Case('astm_e119')
        n_e119 = 1
     End Select
  End Do
  ndevc_files = n_iso + n_hc + n_ex + n_e119
  
  Allocate(fds_devc_rows(ndevc_files),stat=ios); Call error_allocate(ios)
  Allocate(fds_devc_cols(ndevc_files),stat=ios); Call error_allocate(ios)

  Do i=1,ndevc_files
    fds_devc_rows(i)=iso_ntimes
    fds_devc_cols(i)=2  ! time, temperature
  End Do
   
  ! Allocate memory for time series data
  ndevc=Sum(fds_devc_cols(1:ndevc_files))-ndevc_files
  Allocate(fds_devc_name_b(ndevc),               stat=ios); Call error_allocate(ios)
  Allocate(fds_devc_time(fds_devc_rows(1)),      stat=ios); Call error_allocate(ios)
  Allocate(fds_devc_data(fds_devc_rows(1),ndevc),stat=ios); Call error_allocate(ios)

  Do j = 1, iso_ntimes
     fds_devc_time(j)=(j-1)*(iso_tend-iso_tbegin)/Max(iso_ntimes-1,1)
  End Do

  ! ISO 834: T = 345 log10(8t+1) + 20 (T in C and time in minutes)
  ! ASTM E 119 curve: 750*(1-exp(-3.79553*sqrt(t)))+170.41*sqrt(t) + 20 [t]=hour
  ! EC1-1-2 cl3.2.3 Hydrocarbon curve 1080*(1-0.325*exp(-0.167*t)-0.675*exp(-2.5*t)) + 20 [t]=min
  ! EC1-1-2 cl3.2.2 External fire curve 660*(1-0.687*exp(-0.32*t)-0.313*exp(-3.8*t)) + 20 [t]=min
  !                                     temp constant after 22 mins at 660 C
     
  i = 0
  If (n_iso > 0) Then
     i = i + 1
     fds_devc_name_b(i)='iso_834'
     Do j = 1, iso_ntimes
        t_min = fds_devc_time(j)/60.0
        fds_devc_data(j,i)=345.0*Log10(8.0*t_min+1.0) + 20.0
     End Do
  End If
  If (n_hc > 0) Then
     i = i + 1
     fds_devc_name_b(i)='ec1-1-2_hc'
     Do j = 1, iso_ntimes
        t_min = fds_devc_time(j)/60.0
        fds_devc_data(j,i)=1080.0*(1.0-0.325*Exp(-0.167*t_min)-0.675*Exp(-2.5*t_min)) + 20.0
     End Do
  End If
  If (n_ex > 0) Then
     i = i + 1
     fds_devc_name_b(i)='ec1-1-2_ex'
     Do j = 1, iso_ntimes
        t_min = fds_devc_time(j)/60.0
        fds_devc_data(j,i)=660.0*(1.0-0.687*Exp(-0.32*t_min)-0.313*Exp(-3.8*t_min)) + 20.0
     End Do
  End If
  If (n_e119 > 0) Then
     i = i + 1
     fds_devc_name_b(i)='astm_e119'
     Do j = 1, iso_ntimes
        t_min = fds_devc_time(j)/3600.0  ! hours
        fds_devc_data(j,i)=750.0*(1.0-Exp(-3.79553*Sqrt(t_min)))+170.41*Sqrt(t_min) + 20.0
     End Do
  End If

End Subroutine import_iso_data

Subroutine filter_iso_data
!----------------------------------------------------------------------
! Filter unwanted ISO data
! (Actually: Define arrays just for the time-temp curves that are used)
!----------------------------------------------------------------------
  Use error_messages
  Use fds_devc_arrays
  Use fds_head_arrays
  Use global_constants
  Use global_variables
  Use mapping_arrays
  Use string_handling
  Implicit None

  Integer :: i,j,ios,idevc,ndevc,nrows,ncols_i
  Integer :: n_iso, n_hc, n_ex, n_e119

  !-----------------------------------------------
  ! Convert connectivity table into numerical form
  !-----------------------------------------------

  If (nset_connectivity) Then

     n_iso  =0
     n_hc   =0
     n_ex   =0
     n_e119 =0
     nrows = Ubound(connectivity_table,1)

     Do j = 1, nrows
        Select Case(Trim(Lowercase(connectivity_table(j,2))))
        Case('iso_834')
           n_iso  = 1
        Case('ec1-1-2_hc')
           n_hc   = 1
        Case('ec1-1-2_ex')
           n_ex   = 1
        Case('astm_e119')
           n_e119 = 1
        End Select
     End Do
     
    Allocate(connectivity_table_num(nrows,2),stat=ios); Call error_allocate(ios)
    Allocate(connectivity_phys_cons(nrows,3),stat=ios); Call error_allocate(ios)
    Allocate(time_shift_mask(nrows),         stat=ios); Call error_allocate(ios)

    time_shift_mask = .False.
    connectivity_table_num= 0 
    connectivity_phys_cons(:,1)= emissivity ! Default
    connectivity_phys_cons(:,2)= hcoeff     ! Default
    connectivity_phys_cons(:,3)= 0.0        ! Default time shift (s)
    row_loop_0: Do i = 1, nrows
       Select Case(Trim(Lowercase(connectivity_table(i,2))))
       Case('iso_834')
          connectivity_table_num(i,2)= 1
       Case('ec1-1-2_hc')
          connectivity_table_num(i,2)= n_iso + 1
       Case('ec1-1-2_ex')
          connectivity_table_num(i,2)= n_iso + n_hc + 1
       Case('astm_e119')
          connectivity_table_num(i,2)= n_iso + n_hc + n_ex + 1
       Case default
          ! Exception handling
          Write(*,'(3(a))') 'ERROR: in NSET connectivity table: time-temp curve ', &
               Trim(quote(connectivity_table(i,2))), ' is not supported'
          Stop
       End Select
       
       ! Count the columns on this nset connectivity line
       ncols_i = 2  ! defaul for time-temp curve input
       Do j = 3, Ubound(connectivity_table,2)
          If (Trim(connectivity_table(i,j)) /= '') ncols_i = ncols_i + 1
       End Do
       
       Select Case(ncols_i)
       Case(5)
          ! Read all 
          connectivity_phys_cons(i,1)= str2float(connectivity_table(i,3))
          connectivity_phys_cons(i,2)= str2float(connectivity_table(i,4))
          connectivity_phys_cons(i,3)= str2float(connectivity_table(i,5))
          time_shift_mask(i) = .True.
       Case(4)
          ! Read epsilon and h_coeff
          connectivity_phys_cons(i,1)= str2float(connectivity_table(i,3))
          connectivity_phys_cons(i,2)= str2float(connectivity_table(i,4))
       Case(3)
          ! Read time_shift
          connectivity_phys_cons(i,3) = str2float(connectivity_table(i,3))
          time_shift_mask(i) = .True.
       Case(2)
          ! Use default epsilon and h_coeff
       Case Default
          ! Exception handling
          Write(*,'(3(a))') 'ERROR: in NSET connectivity table at row ', &
               Trim(int2str(i)), ', wrong number of columns '
          Stop
       End Select
 
       If (connectivity_phys_cons(i,1) < 0.0 .Or. connectivity_phys_cons(i,1) > 1.0) Then
          Write(*,'(a)') 'ERROR: Time-temp curve and unphysical emissivity in nset file'
          Stop
       End If
       If (connectivity_phys_cons(i,2) < 0.0) Then
          Write(*,'(a)') 'ERROR: Time-temp curve and hcoeff < 0 in nset file'
          Stop
       End If
       If (connectivity_phys_cons(i,3) < 0.0) Then
          Write(*,'(a)') 'ERROR: Time-temp curve and time_shift < 0 in nset file'
          Stop
       End If

    End Do row_loop_0

  End If

  !-------------------------
  ! Generate FDS data arrays
  !-------------------------

  ndevc = Ubound(fds_devc_name_b,1)
  nnodes_fds = ndevc
  ntimes_fds = iso_ntimes
  
  Allocate(fds_id(nnodes_fds),             stat=ios); Call error_allocate(ios)
  Allocate(fds_time(ntimes_fds),           stat=ios); Call error_allocate(ios)
  Allocate(fds_idevc(nnodes_fds),          stat=ios); Call error_allocate(ios)
  Allocate(fds_data(ntimes_fds,nnodes_fds),stat=ios); Call error_allocate(ios)

  !------------------------------
  ! Create arrays containing data
  !------------------------------

  Do idevc = 1, ndevc
     row_loop_3: Do i = 1, nrows
        If (Trim(connectivity_table(i,2)) == Trim(fds_devc_name_b(idevc))) Then
           fds_id(idevc)    = fds_devc_name_b(idevc)
           fds_idevc(idevc) = idevc
           fds_data(1:iso_ntimes,idevc) = fds_devc_data(1:iso_ntimes,idevc)
        End If
     End Do row_loop_3
  End Do

  fds_time = fds_devc_time

End Subroutine filter_iso_data

End Module iso_reader
