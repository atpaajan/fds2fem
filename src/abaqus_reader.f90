!---------------------------------------------------------
! Functions and subroutines for reading ABAQUS input files
!---------------------------------------------------------
module abaqus_reader

contains

subroutine abaqus_reader_module()
!---------
! Main sub
!---------
  use abaqus_arrays
  !use abaqus_dump
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use miscellaneous
  use string_handling
  implicit none

  integer :: i,j,ios
  real(kind=rk) :: t1,t2
  character(len=chr80) :: stmp 

  !---------------
  ! Initialization
  !---------------  

  abaqus_xyz_available=.false.
  abaqus_data_available=.false.
  abaqus_model_available=.false.

  write(*,'(a)') ''
  write(*,'(a)') 'ABAQUS reader module'

  if (len_trim(fem_input_file) == 0) then
    !----------------------
    ! No ABAQUS input given
    !----------------------
    
    write(*,'(t3,a)') 'Nothing to be done'

  else
    !---------------------------
    ! Parse ABAQUS input file(s)
    !---------------------------
    
    open(unit=scratch_channel,status='scratch',iostat=ios)
    if (ios /= 0) call error_open_scratch()

    call cat_input_files()
    call parse_keywords()
   
    write(*,'(t3,a)') 'Reading part definitions'

    write(*,'(t6,a,t23)') '* Nodes'
    call cpu_time(t1); call read_part_nodes(); call cpu_time(t2)
    !write(*,'(3(a))') '(', trim(sec2str(t2-t1)), ' s)'

    write(*,'(t6,a,t23)') '* Node sets'
    call cpu_time(t1); call read_part_nsets(); call cpu_time(t2)
    !write(*,'(3(a))') '(', trim(sec2str(t2-t1)), ' s)'

    write(*,'(t6,a,t23)') '* Elements'
    if (trim(transfer_quantity) == 'adiabatic_surface_temperature' .or. &
      trim(transfer_quantity) == 'net heat flux') then
      call cpu_time(t1); call read_part_elements(); call cpu_time(t2)
      !write(*,'(3(a))') '(', trim(sec2str(t2-t1)), ' s)'
    else
      !write(*,'(a)') '(omitted)'
      !write(*,'(3(a))') '(', trim(sec2str(t2-t1)), ' s)'
    end if

    ! Relevant arrays
    !--------------------------
    !   part_node_xyz(:),
    !   part_node_number(:),
    !   part_node_name(:),
    !   part_node_part(:),
    !   part_nset_name(:),
    !   part_nset_node(:),
    !   part_nset_part(:),
    !   part_element_number(:),
    !   part_element_node(:,:),
    !   part_element_part(:)
    !--------------------------
    
    write(*,'(t3,a)') 'Reading part instance definitions'

    write(*,'(t6,a,t23)') '* Nodes'
    call cpu_time(t1); call read_instance_nodes(); call cpu_time(t2)
    !write(*,'(3(a))') '(', trim(sec2str(t2-t1)), ' s)'

    write(*,'(t6,a,t23)') '* Node sets'
    call cpu_time(t1); call read_instance_nsets(); call cpu_time(t2)
    !write(*,'(3(a))') '(', trim(sec2str(t2-t1)), ' s)'

    write(*,'(t6,a,t23)') '* Elements'
    if (trim(transfer_quantity) == 'adiabatic_surface_temperature' .or. &
      trim(transfer_quantity) == 'net heat flux') then
      call cpu_time(t1); call read_instance_elements(); call cpu_time(t2)
      !write(*,'(3(a))') '(', trim(sec2str(t2-t1)), ' s)'
    else
      !write(*,'(a)') '(omitted)'
      !write(*,'(3(a))') '(', trim(sec2str(t2-t1)), ' s)'
    end if
    
    write(*,'(t3,a)') 'Reading assembly definition'
    write(*,'(t6,a,t23)') '* Node sets'
    call cpu_time(t1); call read_assembly_nsets(); call cpu_time(t2)
    !write(*,'(3(a))') '(', trim(sec2str(t2-t1)), ' s)'
    
    ! Relevant arrays
    !-------------------------------
    !   instance_node_xyz(:),
    !   instance_node_number(:),
    !   instance_node_name(:),
    !   instance_node_instance(:),
    !   instance_nset_name(:),
    !   instance_nset_node(:),
    !   instance_nset_instance(:)
    !   instance_element_number(:),
    !   instance_element_node(:,:),
    !   instance_element_part(:)
    !   instance_element_instance(:)
    !-------------------------------

    close(scratch_channel)
    
    write(*,'(t3,a)') 'Constructing assembly'

    write(*,'(t6,a,t23)') '* Nodes'
    call cpu_time(t1); call inherit_instance_nodes(); call cpu_time(t2)
    !write(*,'(3(a))') '(', trim(sec2str(t2-t1)), ' s)'

    write(*,'(t6,a,t23)') '* Node sets'
    call cpu_time(t1); call inherit_instance_nsets(); call cpu_time(t2) 
    !write(*,'(3(a))') '(', trim(sec2str(t2-t1)), ' s)'

    write(*,'(t6,a,t23)') '* Elements'
    if (trim(transfer_quantity) == 'adiabatic_surface_temperature' .or. &
      trim(transfer_quantity) == 'net heat flux') then
      call cpu_time(t1); call inherit_instance_elements(); call cpu_time(t2)
      !write(*,'(3(a))') '(', trim(sec2str(t2-t1)), ' s)'
    else
      !write(*,'(a)') '(omitted)'
      !write(*,'(3(a))') '(', trim(sec2str(t2-t1)), ' s)'
    end if

    write(*,'(t6,a,t23)') '* Translations'
    call cpu_time(t1); call transl_instance_nodes(); call cpu_time(t2)
    !write(*,'(3(a))') '(', trim(sec2str(t2-t1)), ' s)'

    write(*,'(t6,a,t23)') '* Rotations'
    call cpu_time(t1); call rotate_instance_nodes(); call cpu_time(t2)
    !write(*,'(3(a))') '(', trim(sec2str(t2-t1)), ' s)'

    ! Relevant arrays
    !----------------------------
    !   model_node_xyz(:),
    !   model_node_number(:),
    !   model_node_name(:),
    !   model_node_instance(:),
    !   model_nset_name(:),
    !   model_nset_node(:),
    !   model_nset_instance(:)
    !   model_element_number(:),
    !   model_element_type(:),
    !   model_element_node(:,:),
    !   model_element_instance(:)
    !----------------------------

    ! Check if user-selected node sets exist

    if (nset_connectivity) then
      loop_1: do i=1,ubound(connectivity_table,1)
        loop_2: do j=1,ubound(model_nset_name,1)
          stmp=trim(lowercase(model_nset_instance(j))) // &
            '.' // trim(lowercase(model_nset_name(j)))
          if (trim(lowercase(connectivity_table(i,1))) == trim(stmp)) then
            cycle loop_1
          end if
        end do loop_2
        write(*,'(7(a))') 'ERROR: in NSET connectivity table: NSET ', &
          trim(quote(connectivity_table(i,1))), ' not found (file ', &
          trim(quote(nset_input_file)), ', line ', trim(int2str(i)), ')'
        stop
      end do loop_1
    end if

    if (trim(transfer_quantity) == 'adiabatic_surface_temperature' .or. &
      trim(transfer_quantity) == 'net heat flux') then
      write(*,'(t3,a)') 'Constructing element map'
      call cpu_time(t1); call create_element_node_map(); call cpu_time(t2)
      ! write(*,'(t6,3(a))') 'Ready in ', trim(sec2str(t2-t1)), ' s'

      write(*,'(t3,a)') 'Calculating nodal surface areas'
      call cpu_time(t1); call create_node_surface_map(); call cpu_time(t2)
      ! write(*,'(t6,3(a))') 'Ready in ', trim(sec2str(t2-t1)), ' s'

      abaqus_model_available=.true.
      fem_model_available=.true.

    end if

    ! Now we have
    !----------------------------
    ! model_element_int_node(:,:)
    !----------------------------

    call deallocate_excess_abaqus_arrays()

    call filter_abaqus_data() 
    
    !if (trim(transfer_quantity) == 'adiabatic_surface_temperature' .or. &
    !  trim(transfer_quantity) == 'net heat flux') then
    !  call dump_abaqus_nodal_areas_vtk()
    !end if

    ! At this point we have
    !----------------------
    !   abaqus_nset(:),
    !   abaqus_xyz(:,:),
    !   abaqus_node_name(:),
    !   abaqus_node_number(:)
    !----------------------

    ! ### TEMPORARY SOLUTION ###

    allocate(fem_xyz(ubound(abaqus_xyz,1),3))
    allocate(fem_nset(ubound(abaqus_nset,1)))
    allocate(fem_node_name(ubound(abaqus_node_name,1)))
    allocate(fem_node_number(ubound(abaqus_node_number,1)))
    
    fem_xyz=abaqus_xyz
    fem_nset=abaqus_nset
    fem_node_name=abaqus_node_name
    fem_node_number=abaqus_node_number

    fem_xyz_available=.true.

    deallocate(abaqus_xyz,stat=ios)
    deallocate(abaqus_nset,stat=ios)
    deallocate(abaqus_node_name,stat=ios)
    deallocate(abaqus_node_number,stat=ios)

    ! ### TEMPORARY SOLUTION ###

    call deallocate_excess_abaqus_arrays()
    
    abaqus_xyz_available=.true.

    write(*,'(t3,3(a))') 'ABAQUS input file ', trim(quote(fem_input_file)), ' parsed'

  end if

end subroutine abaqus_reader_module

!***************
! Auxiliary subs
!***************

subroutine cat_input_files()
!-----------------------------------------------------
! Expand include keywords etc. in an ABAQUS input file
!-----------------------------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use string_handling
  implicit none 

  integer :: i,ios,ilevel,nlevel
  logical :: skip_line,skip_next_line,recursion
  character(len=chr80) :: output_file,keyword,paramstr,ctmp1,ctmp2,ctmp3,ctmp4
  character(len=input_line_length) :: input_line

  integer, dimension(:), allocatable :: line_number

  character(len=chr80), dimension(:), allocatable :: input_file  

  !---------------
  ! Initialization
  !---------------

  nlevel=ubound(iochannel,1)
  if (nlevel < 2) then
    write(*,'(a)') 'ERROR: in subroutine cat_input_files(): not enough I/O channels'  
    stop
  end if

  allocate(input_file(nlevel),stat=ios);  call error_allocate(ios)
  allocate(line_number(nlevel),stat=ios); call error_allocate(ios)

  input_file(1)=trim(fem_input_file); input_file(2:nlevel)=''

  open(unit=iochannel(1),file=trim(input_file(1)),status='old',iostat=ios)
  if (ios /= 0) call error_open_file(input_file(1))

  !-------------------
  ! Concatenation loop
  !-------------------

  ilevel=1
  line_number=0
  recursion=.false.
  skip_line=.false.
  skip_next_line=.false.

  loop: do

    line_number(ilevel)=line_number(ilevel)+1

    read(iochannel(ilevel),'(a)',iostat=ios) input_line
    select case(ios)
    case (:-1)
      if (ilevel == 1) then
        ! At the end of the primary input file
        exit loop
      else
        ! At the end of one of the additional input files
        close(unit=iochannel(ilevel))
        ilevel=ilevel-1
        cycle loop
      end if
    case(1:)
      call error_read_file(input_file(ilevel),line_number(ilevel))
    end select
  
    if (skip_next_line) then
      skip_line=.true.
      skip_next_line=.false.
    else
      skip_line=.false.
    end if

    if (skip_line) cycle loop

    input_line=adjustl(input_line)
    if (input_line(1:1) == '*' .and. &
        input_line(2:2) /= '*') then

      !----------
      ! Branching
      !----------
      keyword=get_keyword(input_line)
      if (trim(lowercase(keyword)) == 'include') then

        if (ilevel < (nlevel-1)) then

          ilevel=ilevel+1
          input_file(ilevel)=''
          paramstr='input'; input_file(ilevel)=get_parameter_string(input_line,paramstr)

          do i=1,(ilevel-1)
            if (trim(input_file(ilevel)) == trim(input_file(i))) then
              write(*,'(5(a))') 'WARNING: recursive use of include-keyword encountered (file ', &
                trim(quote(input_file(ilevel-1))), ', line ', trim(int2str(line_number(ilevel-1))), ')'
              recursion=.true.
              ilevel=ilevel-1
            end if
          end do

          if (len_trim(input_file(ilevel)) == 0) then
            ! No input file given as a parameter to the include-keyword
            ilevel=ilevel-1 
          else if (.not. recursion) then
            ! Otherwise, if no recursive use of include-keyword is encountered
            open(unit=iochannel(ilevel),file=trim(input_file(ilevel)),status='old',iostat=ios)
            if (ios /= 0) then
              ilevel=ilevel-1
            end if
            skip_line=.true.
          end if

        end if

      else if (trim(lowercase(keyword)) == 'node') then

        if (ilevel < (nlevel-1)) then

          ilevel=ilevel+1
          input_file(ilevel)=''
          paramstr='input';  input_file(ilevel)=get_parameter_string(input_line,paramstr)
          paramstr='nset';   ctmp1             =get_parameter_string(input_line,paramstr)
          paramstr='system'; ctmp2             =get_parameter_string(input_line,paramstr)

          do i=1,(ilevel-1)
            if (trim(input_file(ilevel)) == trim(input_file(i))) then
              write(*,'(5(a))') 'WARNING: recursive use of node-keyword input-parameter encountered (file ', &
                trim(quote(input_file(ilevel-1))), ', line ', trim(int2str(line_number(ilevel-1))), ')'
              recursion=.true.
              ilevel=ilevel-1
            end if
          end do

          if (len_trim(input_file(ilevel)) == 0) then
            ! No input file given as a parameter to the node-keyword
            ilevel=ilevel-1 
          else if (.not. recursion) then
            ! Otherwise, if no recursive use of input/include-keyword is encountered
            open(unit=iochannel(ilevel),file=trim(input_file(ilevel)),status='old',iostat=ios)
            if (ios /= 0) then
              ilevel=ilevel-1
            end if

            ! Omit the input-parameter from the keyword line
            input_line=trim(input_line(1:5))
            if (len_trim(ctmp1) /= 0) then
              input_line=trim(input_line) // ', nset=' // trim(ctmp1)
            end if 

            if (len_trim(ctmp2) /= 0) then
              input_line=trim(input_line) // ', system=' // trim(ctmp2)
            end if 

          end if

        end if

      else if (trim(lowercase(keyword)) == 'element') then

        if (ilevel < (nlevel-1)) then

          ilevel=ilevel+1
          input_file(ilevel)=''
          paramstr='input';  input_file(ilevel)=get_parameter_string(input_line,paramstr)
          paramstr='type';   ctmp1             =get_parameter_string(input_line,paramstr)
          paramstr='elset';  ctmp2             =get_parameter_string(input_line,paramstr)
          paramstr='file';   ctmp3             =get_parameter_string(input_line,paramstr)
          paramstr='offset'; ctmp4             =get_parameter_string(input_line,paramstr)
          ! paramstr='solid element numbering'; element_sen(i_element)=get_parameter_string(input_line,paramstr)

          do i=1,(ilevel-1)
            if (trim(input_file(ilevel)) == trim(input_file(i))) then
              write(*,'(5(a))') 'WARNING: recursive use of element-keyword input-parameter encountered (file ', &
                trim(quote(input_file(ilevel-1))), ', line ', trim(int2str(line_number(ilevel-1))), ')'
              recursion=.true.
              ilevel=ilevel-1
            end if
          end do

          if (len_trim(input_file(ilevel)) == 0) then
            ! No input file given as a parameter to the element-keyword
            ilevel=ilevel-1 
          else if (.not. recursion) then
            ! Otherwise, if no recursive use of input/include-keyword is encountered
            open(unit=iochannel(ilevel),file=trim(input_file(ilevel)),status='old',iostat=ios)
            if (ios /= 0) then
              ilevel=ilevel-1
            end if

            ! Omit the input-parameter from the keyword line
            input_line=trim(input_line(1:5))
            if (len_trim(ctmp1) /= 0) then
              input_line=trim(input_line) // ', type=' // trim(ctmp1)
            end if 

            if (len_trim(ctmp2) /= 0) then
              input_line=trim(input_line) // ', elset=' // trim(ctmp2)
            end if

            if (len_trim(ctmp3) /= 0) then
              input_line=trim(input_line) // ', file=' // trim(ctmp2)
            end if 

            if (len_trim(ctmp4) /= 0) then
              input_line=trim(input_line) // ', offset=' // trim(ctmp4)
            end if  

          end if

        end if

      end if

    else if (input_line(1:2) == '**') then

      if (input_line(1:27) == '** fds2fem-amplitude-marker') then
        skip_line=.true.; skip_next_line=.true.
      else if (input_line(1:26) == '** fds2fem-boundary-marker') then
        skip_line=.true.; skip_next_line=.true.
      else if (input_line(1:26) == '** fds2fem-cradiate-marker') then
        skip_line=.true.; skip_next_line=.true.
      else if (input_line(1:23) == '** fds2fem-cfilm-marker') then
        skip_line=.true.; skip_next_line=.true.
      else if (input_line(1:26) == '** fds2fem-physical-marker') then
        skip_line=.true.; skip_next_line=.true.
      end if

    end if

    if (.not. skip_line) then
      write(scratch_channel,'(a)',iostat=ios) trim(input_line)
      if (ios /= 0) call error_write_file(output_file)
    end if

  end do loop

  close(unit=iochannel(1))

  rewind(scratch_channel)

  if (.not. recursion) then
    ! Do nothing
  else
    write(*,'(a)') 'ERROR: recursive inclusion of external files is not allowed'
    stop
  end if

end subroutine cat_input_files

subroutine transl_instance_nodes()
!-----------------------------------
! Translate nodes of a part instance
!-----------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use mathematics
  use string_handling
  implicit none

  integer :: i,ios,i_instance
  real(kind=rk) :: dx,dy,dz,dr

  logical, dimension(:), allocatable :: translate

  allocate(translate(n_instance),stat=ios); call error_allocate(ios)

  !------------------------------
  ! Select instances to translate
  !------------------------------

  translate=.false.
  do i_instance=1,n_instance
    dx=instance_translate(i_instance,1)
    dy=instance_translate(i_instance,2)
    dz=instance_translate(i_instance,3)
    dr=sqrt(dx*dx+dy*dy+dz*dz)
    if (abs(dr) > epsilon(dr)) then
      translate(i_instance)=.true.
    else
      translate(i_instance)=.false.
    end if
  end do

  !----------------------
  ! Apply the translation
  !----------------------

  do i=1,ubound(model_node_xyz,1)
    instance_loop: do i_instance=1,n_instance
      if (trim(model_node_instance(i)) == trim(instance_name(i_instance))) then
        if (translate(i_instance)) then

          model_node_xyz(i,1)=model_node_xyz(i,1)+instance_translate(i_instance,1)
          model_node_xyz(i,2)=model_node_xyz(i,2)+instance_translate(i_instance,2)
          model_node_xyz(i,3)=model_node_xyz(i,3)+instance_translate(i_instance,3)
          
          exit instance_loop
        end if
      end if
    end do instance_loop
  end do

end subroutine transl_instance_nodes

subroutine rotate_instance_nodes()
!-----------------------------------
! Translate nodes of a part instance
!-----------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use mathematics
  use string_handling
  implicit none

  integer :: i,ios,i_instance
  real(kind=rk) :: dx,dy,dz,dr,theta

  logical, dimension(:), allocatable :: rotate
  real(kind=rk), dimension(4) :: vtmp
  real(kind=rk), dimension(4,4) :: rotmat
  
  allocate(rotate(n_instance),stat=ios); call error_allocate(ios)

  !---------------------------
  ! Select instances to rotate
  !---------------------------

  rotate=.false.
  do i_instance=1,n_instance
    dx=abs(instance_rotate(i_instance,1)-instance_rotate(i_instance,4))
    dy=abs(instance_rotate(i_instance,2)-instance_rotate(i_instance,5))
    dz=abs(instance_rotate(i_instance,3)-instance_rotate(i_instance,6))
    dr=sqrt(dx*dx+dy*dy+dz*dz)
    theta=instance_rotate(i_instance,7)
    if (abs(dr) > epsilon(dr) .and. abs(theta) > epsilon(theta)) then 
      rotate(i_instance)=.true.
    else
      rotate(i_instance)=.false.
    end if
  end do

  !-------------------
  ! Apply the rotation
  !-------------------

  instance_loop: do i_instance=1,n_instance
    theta=(instance_rotate(i_instance,7)/180.0)*(4.0*atan(1.0))
    rotmat=rotate_arbitrary_axis(instance_rotate(i_instance,1:3), &
      instance_rotate(i_instance,4:6),theta)

    do i=1,ubound(model_node_xyz,1)
      if (trim(model_node_instance(i)) == trim(instance_name(i_instance))) then
        if (rotate(i_instance)) then
            vtmp(1:3)=model_node_xyz(i,1:3); vtmp(4)=1.0; vtmp=matmul(rotmat,vtmp)
            model_node_xyz(i,1:3)=vtmp(1:3)
        end if
      end if
    end do
  end do instance_loop

end subroutine rotate_instance_nodes

subroutine parse_keywords()
!---------------------------------------------------
! Parse ABAQUS input file for relevant keyword lines
!---------------------------------------------------
  use global_constants
  use global_variables
  use abaqus_arrays
  use string_handling
  use error_messages
  implicit none

  integer :: ios,line_number,i,j,nlines,from_instance,n1,n2
  integer :: i_nset,i_ngen,i_ncopy,i_node,i_part,i_nfill,i_nmap,i_system,i_element
  integer :: i_assembly,i_include,i_instance,i_end_assembly,i_end_instance,i_end_part

  logical :: in_part,in_instance,in_assembly,in_node,in_nset,in_element

  character(len=1) :: lastchr
  character(len=chr80) :: paramstr
  character(len=input_line_length) :: keyword,input_line

  in_part=.false.; in_instance=.false.; in_assembly=.false.
  in_node=.false.; in_nset=.false.;     in_element=.false.

  !----------------------------------
  ! Count the number of keyword lines
  !----------------------------------

  n_assembly=0; n_end_assembly=0; n_include=0
  n_instance=0; n_end_instance=0; n_ncopy=0
  n_nfill=0;    n_ngen=0;         n_nmap=0
  n_node=0;     n_nset=0;         n_part=0
  n_end_part=0; n_system=0;       n_element=0

  line_number=0
  count_loop: do

    read(scratch_channel,'(a)',iostat=ios) input_line
    select case(ios)
    case (:-1)
      exit count_loop
    case(1:)
      call error_read_file(fem_input_file)
    end select

    line_number=line_number+1

    input_line=adjustl(input_line) 
    if (input_line(1:1) == '*' .and. &
        input_line(2:2) /= '*') then
      !------------------
      ! On a keyword line
      !------------------

      keyword=get_keyword(input_line)
      select case (lowercase(keyword))
      case ('assembly')
        !----------
        ! *Assembly
        !----------
        if (in_assembly) then
          write(*,'(5(a))') 'ERROR: isolated assembly-keyword encountered (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if
        n_assembly=n_assembly+1
        in_assembly=.true.
      case ('end assembly')
        !--------------
        ! *End assembly
        !--------------
        if (.not. in_assembly) then
          write(*,'(5(a))') 'ERROR: isolated end assembly -keyword encountered (file ', &
          trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if
        n_end_assembly=n_end_assembly+1
        in_assembly=.false.
      case ('element')
        !---------
        ! *Element
        !---------
        n_element=n_element+1
      case ('include')
        !---------
        ! *Include
        !---------
        n_include=n_include+1
      case ('instance')
        !----------
        ! *Instance
        !----------
        if (in_instance) then
          write(*,'(5(a))') 'ERROR: isolated instance-keyword encountered (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if
        n_instance=n_instance+1
        in_instance=.true.
      case ('end instance')
        !--------------
        ! *End instance
        !--------------
        if (.not. in_instance) then
          write(*,'(5(a))') 'ERROR: isolated end instance -keyword encountered (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if
        n_end_instance=n_end_instance+1
        in_instance=.false.
      case ('ncopy')
        !-------
        ! *Ncopy
        !-------
        n_ncopy=n_ncopy+1
      case ('nfill')
        !-------
        ! *Nfill
        !-------
        n_nfill=n_nfill+1
      case ('ngen')
        !------
        ! *Ngen
        !------
        n_ngen=n_ngen+1
      case ('nmap')
        !------        
        ! *Nmap
        !------
        n_nmap=n_nmap+1
      case ('node')
        !------
        ! *Node
        !------
        n_node=n_node+1
      case ('nset')
        !------
        ! *Nset
        !------
        n_nset=n_nset+1
      case ('part')
        !------
        ! *Part
        !------
        if (in_part) then
          write(*,'(5(a))') 'ERROR: isolated part-keyword encountered (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if
        n_part=n_part+1
        in_part=.true.
      case ('end part')
        !----------
        ! *End part
        !----------
        if (.not. in_part) then
            write(*,'(5(a))') 'ERROR: isolated end part -keyword encountered (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if
        n_end_part=n_end_part+1
        in_part=.false.
      case ('system')
        n_system=n_system+1
      case default
        ! Do nothing
      end select
    end if

  end do count_loop

  nlines=line_number
  if (nlines == 0) then
    write(*,'(3(a))') 'ERROR: file ', trim(quote(fem_input_file)), ' is empty'
    stop
  end if

  rewind(scratch_channel)

  !-------------------
  ! Exception handling
  !-------------------

  if (n_assembly == 0) then
    write(*,'(3(a))') 'ERROR: no assembly definition found (file ' , &
      trim(quote(fem_input_file)), ')'
    stop
  else if (n_assembly > 1) then
    write(*,'(3(a))') 'ERROR: more than one assembly definition found (file ' , &
      trim(quote(fem_input_file)), ')'
    stop
  end if

  if (n_instance == 0) then
    write(*,'(3(a))') 'ERROR: no part instance definition found (file ' , &
      trim(quote(fem_input_file)), ')'
    stop
  end if

  if (n_part == 0) then
    write(*,'(3(a))') 'ERROR: no part definition found (file ' , &
      trim(quote(fem_input_file)), ')'
    stop
  end if

  if (in_assembly) then
    write(*,'(3(a))') 'ERROR: unterminated assembly definition found (file ' , &
      trim(quote(fem_input_file)), ')'
    stop
  end if

  if (in_instance) then
    write(*,'(3(a))') 'ERROR: unterminated instance definition(s) found (file ' , &
      trim(quote(fem_input_file)), ')'
    stop
  end if

  if (in_part) then
    write(*,'(3(a))') 'ERROR: unterminated part definition(s) found (file ' , &
      trim(quote(fem_input_file)), ')'
    stop
  end if

  !----------------
  ! Allocate memory
  !----------------

  if (n_assembly > 0) then
    allocate(assembly_name(n_assembly),stat=ios); call error_allocate(ios)
    assembly_name = ''
    
    allocate(assembly_lines(n_assembly,2),stat=ios); call error_allocate(ios)
    assembly_lines = 0
  end if

  if (n_element > 0) then
    allocate(element_type(n_element),stat=ios);   call error_allocate(ios)
    allocate(element_elset(n_element),stat=ios);  call error_allocate(ios)
    allocate(element_file(n_element),stat=ios);   call error_allocate(ios)
    allocate(element_input(n_element),stat=ios);  call error_allocate(ios)
    allocate(element_offset(n_element),stat=ios); call error_allocate(ios)
    allocate(element_sen(n_element),stat=ios);    call error_allocate(ios)
    element_type   = ''
    element_elset  = ''
    element_file   = ''
    element_input  = ''
    element_offset = 100000.0 ! Abaqus default
    element_sen    = ''        ! Abaqus default

    allocate(element_part(n_element),stat=ios);     call error_allocate(ios)
    allocate(element_instance(n_element),stat=ios); call error_allocate(ios)
    allocate(element_assembly(n_element),stat=ios); call error_allocate(ios)
    element_part     = ''
    element_instance = ''
    element_assembly = ''

    allocate(element_lines(n_element,3),stat=ios); call error_allocate(ios)
    allocate(element_elements(n_node),stat=ios);   call error_allocate(ios)
    element_lines    = 0
    element_elements = 0
  end if

  if (n_include > 0) then
    allocate(include_input(n_include),stat=ios);    call error_allocate(ios)
    allocate(include_password(n_include),stat=ios); call error_allocate(ios)
    include_input    = ''
    include_password = ''

    allocate(include_line_number(n_include),stat=ios); call error_allocate(ios)
    include_line_number=0
  end if

  if (n_instance > 0) then
    allocate(instance_name(n_instance),stat=ios);     call error_allocate(ios)
    allocate(instance_part(n_instance),stat=ios);     call error_allocate(ios)
    allocate(instance_instance(n_instance),stat=ios); call error_allocate(ios)
    allocate(instance_library(n_instance),stat=ios);  call error_allocate(ios)
    instance_name     = ''
    instance_part     = ''
    instance_instance = ''
    instance_library  = ''

    allocate(instance_translate(n_instance,3),stat=ios); call error_allocate(ios)
    allocate(instance_rotate(n_instance,7),stat=ios);    call error_allocate(ios)
    instance_translate = 0.0
    instance_rotate    = 0.0
    
    allocate(instance_lines(n_instance,2),stat=ios); call error_allocate(ios)
    instance_lines = 0
  end if

  if (n_ncopy > 0) then
    allocate(ncopy_change_number(n_ncopy),stat=ios); call error_allocate(ios)
    allocate(ncopy_old_set(n_ncopy),stat=ios);       call error_allocate(ios)
    allocate(ncopy_pole(n_ncopy),stat=ios);          call error_allocate(ios)
    allocate(ncopy_reflect(n_ncopy),stat=ios);       call error_allocate(ios)
    allocate(ncopy_shift(n_ncopy),stat=ios);         call error_allocate(ios)
    allocate(ncopy_multiple(n_ncopy),stat=ios);      call error_allocate(ios)
    allocate(ncopy_new_set(n_ncopy),stat=ios);       call error_allocate(ios)
    ncopy_change_number = 0
    ncopy_old_set       = ''
    ncopy_pole          = .false.
    ncopy_reflect       = ''
    ncopy_shift         = .false.
    ncopy_multiple      = 0
    ncopy_new_set       = ''
  end if
  
  if (n_nfill > 0) then
    allocate(nfill_bias(n_nfill),stat=ios);     call error_allocate(ios)
    allocate(nfill_nset(n_nfill),stat=ios);     call error_allocate(ios)
    allocate(nfill_singular(n_nfill),stat=ios); call error_allocate(ios)
    allocate(nfill_two_step(n_nfill),stat=ios); call error_allocate(ios)
    nfill_bias     = 0.0
    nfill_nset     = ''
    nfill_singular = 0
    nfill_two_step = .false.
  end if

  if (n_ngen > 0) then
    allocate(ngen_line(n_ngen),stat=ios);   call error_allocate(ios)
    allocate(ngen_nset(n_ngen),stat=ios);   call error_allocate(ios)
    allocate(ngen_system(n_ngen),stat=ios); call error_allocate(ios)
    ngen_line   = ''
    ngen_nset   = ''
    ngen_system = ''
  end if

  if (n_nmap > 0) then
    allocate(nmap_nset(n_nmap),stat=ios);       call error_allocate(ios)
    allocate(nmap_type(n_nmap),stat=ios);       call error_allocate(ios)
    allocate(nmap_definition(n_nmap),stat=ios); call error_allocate(ios)
    nmap_nset       = ''
    nmap_type       = ''
    nmap_definition = ''
  end if

  if (n_node > 0) then
    allocate(node_input(n_node),stat=ios);  call error_allocate(ios)
    allocate(node_nset(n_node),stat=ios);   call error_allocate(ios)
    allocate(node_system(n_node),stat=ios); call error_allocate(ios)
    node_input  = ''
    node_nset   = ''
    node_system = ''

    allocate(node_part(n_node),stat=ios);     call error_allocate(ios)
    allocate(node_instance(n_node),stat=ios); call error_allocate(ios)
    allocate(node_assembly(n_node),stat=ios); call error_allocate(ios)
    node_part     = ''
    node_instance = ''
    node_assembly = ''

    allocate(node_lines(n_node,3),stat=ios); call error_allocate(ios)
    allocate(node_nodes(n_node),stat=ios);   call error_allocate(ios)
    node_lines = 0
    node_nodes = 0
  end if

  if (n_nset > 0) then
    allocate(nset_nset(n_nset),stat=ios);     call error_allocate(ios)
    allocate(nset_elset(n_nset),stat=ios);    call error_allocate(ios)
    allocate(nset_generate(n_nset),stat=ios); call error_allocate(ios)
    allocate(nset_instance(n_nset),stat=ios); call error_allocate(ios)
    allocate(nset_internal(n_nset),stat=ios); call error_allocate(ios)
    allocate(nset_unsorted(n_nset),stat=ios); call error_allocate(ios)
    nset_nset     = ''
    nset_elset    = ''
    nset_instance = ''
    nset_generate = .false.
    nset_internal = .false.
    nset_unsorted = .false.

    allocate(nset_part(n_nset),stat=ios);    call error_allocate(ios)
    allocate(nset_inst(n_nset),stat=ios);    call error_allocate(ios)
    allocate(nset_assm(n_nset),stat=ios);    call error_allocate(ios)
    allocate(nset_lines(n_nset,3),stat=ios); call error_allocate(ios)
    allocate(nset_nodes(n_nset),stat=ios);   call error_allocate(ios)
    nset_part  = ''
    nset_inst  = ''
    nset_assm  = ''
    nset_lines = 0
    nset_nodes = 0
  end if

  if (n_part > 0) then
    allocate(part_name(n_part),stat=ios); call error_allocate(ios)
    part_name  = ''

    allocate(part_lines(n_part,2),stat=ios); call error_allocate(ios)
    part_lines = 0
  end if

  !--------------------------
  ! Read in the keyword lines
  !--------------------------

  i_assembly=0; i_end_assembly=0; i_include=0
  i_instance=0; i_end_instance=0; i_ncopy=0
  i_nfill=0;    i_ngen=0;         i_nmap=0
  i_node=0;     i_nset=0;         i_part=0
  i_end_part=0; i_system=0;       i_element=0

  in_part=.false.; in_instance=.false.; in_assembly=.false.
  in_node=.false.; in_nset=.false.;     in_element=.false.

  line_number=0
  read_loop: do
 
    read(scratch_channel,'(a)',iostat=ios) input_line
    select case(ios)
    case (:-1)
      exit read_loop
    case(1:)
      call error_read_file(fem_input_file)
    end select

    line_number=line_number+1

    input_line=adjustl(input_line) 
    if (input_line(1:2) == '**') then
      !------------------
      ! On a comment line
      !------------------
      if (in_node)    in_node=.false.
      if (in_nset)    in_nset=.false.
      if (in_element) in_element=.false.

    else if (input_line(1:1) == '*' .and. &
             input_line(2:2) /= '*') then
      !------------------
      ! On a keyword line
      !------------------
      if (in_node)    in_node=.false.
      if (in_nset)    in_nset=.false.
      if (in_element) in_element=.false.

      lastchr=input_line(len_trim(input_line):len_trim(input_line))
      if (lastchr == ',') then
        write(*,'(5(a))') 'ERROR: continuation of keyword lines is not supported (file ', &
          trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
        stop
      end if

      keyword=get_keyword(input_line)
      select case (trim(lowercase(keyword)))
      case ('assembly')
        !----------
        ! *Assembly
        !----------
        i_assembly=i_assembly+1
        paramstr='name'; assembly_name(i_assembly)=get_parameter_string(input_line,paramstr)

        if (in_assembly .or. in_instance .or. in_part) then
          write(*,'(5(a))') 'ERROR: in assembly definition: incorrect level (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if

        if (len_trim(assembly_name(i_assembly)) == 0) then
          write(*,'(5(a))') 'ERROR: in assembly definition: name parameter is required (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if

        assembly_lines(i_assembly,1)=line_number
        in_assembly=.true.
      case ('end assembly')
        !--------------
        ! *End Assembly
        !--------------
        assembly_lines(i_assembly,2)=line_number
        in_assembly=.false.
      case ('element')
        !---------
        ! *Element
        !---------
        i_element=i_element+1
        paramstr='type';   element_type(i_element)  =get_parameter_string(input_line,paramstr)
        paramstr='elset';  element_elset(i_element) =get_parameter_string(input_line,paramstr)
        paramstr='file';   element_file(i_element)  =get_parameter_string(input_line,paramstr)
        paramstr='input';  element_input(i_element) =get_parameter_string(input_line,paramstr)
        paramstr='offset'; element_offset(i_element)=get_parameter_real(input_line,paramstr)
        ! paramstr='solid element numbering'; element_sen(i_element)=get_parameter_string(input_line,paramstr)

        if (in_part) then
          element_part(i_element)=part_name(i_part)
        else if (in_instance) then
          element_instance(i_element)=instance_name(i_instance)
          element_part(i_element)=instance_part(i_instance)

          !+-+-+-+-+-+-+-+-+-+-+
          ! Restriction 1.6.2012
          !+-+-+-+-+-+-+-+-+-+-+

          !write(*,'(5(a))') 'ERROR: in element definition: element definition on part instance level &
          !  is currently not supported (file ', trim(quote(fem_input_file)), ', line ', &
          !  trim(int2str(line_number)), ')'
          !stop

        else if (in_assembly) then
          element_assembly(i_element)=assembly_name(i_assembly)
        else
          write(*,'(5(a))') 'ERROR: in element definition: incorrect level (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if

        ! Check for supported element types
        if (trim(lowercase(element_type(i_element))) /= 'dc3d8'   .and. &
            trim(lowercase(element_type(i_element))) /= 'dcc3d8'  .and. &
            trim(lowercase(element_type(i_element))) /= 'dcc3d8d' .and. &
            trim(lowercase(element_type(i_element))) /= 'dc3d4'   .and. &
            trim(lowercase(element_type(i_element))) /= 'ds3'     .and. &
            trim(lowercase(element_type(i_element))) /= 'ds4') then
          write(*,'(5(a))') 'ERROR: in element definition: unsupported element type (file ', &
              trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if

        if (len_trim(element_file(i_element)) /= 0) then
          write(*,'(5(a))') 'ERROR: in element definition: substructures are currently not supported (file ', &
              trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if

        element_lines(i_element,1)=line_number
        element_lines(i_element,2)=0
        element_lines(i_element,3)=0

        in_element=.true.

      case ('include')
        !---------
        ! *Include
        !---------
        i_include=i_include+1
        paramstr='input';    include_input(i_include)   =get_parameter_string(input_line,paramstr)
        paramstr='password'; include_password(i_include)=get_parameter_string(input_line,paramstr)
        include_line_number(i_include)=line_number

        write(*,'(5(a))') 'ERROR: including external files is not supported (file ', &
          trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
        stop

      case ('instance')
        !----------
        ! *Instance
        !----------
        i_instance=i_instance+1
        paramstr='name';     instance_name(i_instance)    =get_parameter_string(input_line,paramstr)
        paramstr='part';     instance_part(i_instance)    =get_parameter_string(input_line,paramstr)
        paramstr='instance'; instance_instance(i_instance)=get_parameter_string(input_line,paramstr)
        paramstr='library';  instance_library(i_instance) =get_parameter_string(input_line,paramstr)

        if (.not. in_assembly) then
          write(*,'(5(a))') 'ERROR: in instance definition: incorrect level (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if

        if (len_trim(instance_instance(i_instance)) /= 0) then
          write(*,'(5(a))') 'ERROR: in instance definition: imported instances are not supported (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if

        if (len_trim(instance_name(i_instance)) == 0) then
          write(*,'(5(a))') 'ERROR: in instance definition: name parameter is required (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if

        if (len_trim(instance_part(i_instance)) == 0) then
          write(*,'(5(a))') 'ERROR: in instance definition: part parameter is required (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if
       
        from_instance=0 
        instance_lines(i_instance,1)=line_number
        in_instance=.true.

      case ('end instance')
        !--------------
        ! *End instance
        !--------------
        from_instance=0 
        instance_lines(i_instance,2)=line_number
        in_instance=.false.

      case ('ncopy')
        !-------
        ! *Ncopy
        !-------
        i_ncopy=i_ncopy+1
        paramstr='change_number'; ncopy_change_number(i_node)=get_parameter_integer(input_line,paramstr)
        paramstr='old_set';       ncopy_old_set(i_node)      =get_parameter_string(input_line,paramstr)
        paramstr='pole';          ncopy_pole(i_node)         =get_parameter_logical(input_line,paramstr)
        paramstr='reflect';       ncopy_reflect(i_node)      =get_parameter_string(input_line,paramstr)
        paramstr='shift';         ncopy_shift(i_node)        =get_parameter_logical(input_line,paramstr)
        paramstr='multiple';      ncopy_multiple(i_node)     =get_parameter_integer(input_line,paramstr)
        paramstr='new_set';       ncopy_new_set(i_node)      =get_parameter_string(input_line,paramstr)

        if (.not. in_part .and. .not. in_instance) then
          write(*,'(5(a))') 'ERROR: in node definition: incorrect level (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if

        write(*,'(5(a))') 'ERROR: in node definition: ncopy-keyword is currently not supported (file ', &
          trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
        stop

      case ('nfill')
        !-------
        ! *Nfill
        !-------
        i_nfill=i_nfill+1  
        paramstr='bias';     nfill_bias(i_node)    =get_parameter_real(input_line,paramstr)
        paramstr='nset';     nfill_nset(i_node)    =get_parameter_string(input_line,paramstr)
        paramstr='singular'; nfill_singular(i_node)=get_parameter_integer(input_line,paramstr)
        paramstr='two_step'; nfill_two_step(i_node)=get_parameter_logical(input_line,paramstr)

        if (.not. in_part .and. .not. in_instance) then
          write(*,'(5(a))') 'ERROR: in node definition: incorrect level (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if

        write(*,'(5(a))') 'ERROR: in node definition: nfill-keyword is currently not supported (file ', &
          trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
        stop

      case ('ngen')
        !------
        ! *Ngen
        !------
        i_ngen=i_ngen+1
        paramstr='line';   ngen_line(i_ngen)  =get_parameter_string(input_line,paramstr)
        paramstr='nset';   ngen_nset(i_ngen)  =get_parameter_string(input_line,paramstr)
        paramstr='system'; ngen_system(i_ngen)=get_parameter_string(input_line,paramstr)

        if (.not. in_part .and. .not. in_instance) then
          write(*,'(5(a))') 'ERROR: in node definition: incorrect level (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if

        if (len_trim(ngen_system(i_ngen)) /= 0) then
          if (trim(lowercase(ngen_system(i_ngen))) /= 'rc' ) then
            write(*,'(5(a))') 'ERROR: in node definition: only Cartesian coordinate system is supported (file ', &
              trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
            stop
          end if
        end if

        write(*,'(5(a))') 'ERROR: in node definition: ngen-keyword is currently not supported (file ', &
          trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
        stop

      case ('nmap')
        !------
        ! *Nmap
        !------
        i_nmap=i_nmap+1  
        paramstr='nset';       nmap_nset(i_node)      =get_parameter_string(input_line,paramstr)
        paramstr='type';       nmap_type(i_node)      =get_parameter_string(input_line,paramstr)
        paramstr='definition'; nmap_definition(i_node)=get_parameter_string(input_line,paramstr)

      case ('node')
        !------
        ! *Node
        !------
        i_node=i_node+1
        paramstr='input';  node_input(i_node) =get_parameter_string(input_line,paramstr)
        paramstr='nset';   node_nset(i_node)  =get_parameter_string(input_line,paramstr)
        paramstr='system'; node_system(i_node)=get_parameter_string(input_line,paramstr)

        if (in_part) then
          node_part(i_node)=part_name(i_part)
        else if (in_instance) then
          node_instance(i_node)=instance_name(i_instance)
          node_part(i_node)=instance_part(i_instance)

          !+-+-+-+-+-+-+-+-+-+-+
          ! Restriction 1.6.2012
          !+-+-+-+-+-+-+-+-+-+-+

          !write(*,'(5(a))') 'ERROR: in node definition: node definition on part instance level &
          !  is currently not supported (file ', trim(quote(fem_input_file)), ', line ', &
          !  trim(int2str(line_number)), ')'
          !stop

        else if (in_assembly) then
          node_assembly(i_node)=assembly_name(i_assembly)
        else
          write(*,'(5(a))') 'ERROR: in node definition: incorrect level (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if

        if (len_trim(node_system(i_node)) /= 0) then
          if (trim(lowercase(node_system(i_node))) /= 'rc') then
            write(*,'(5(a))') 'ERROR: in node definition: only Cartesian coordinate system is supported (file ', &
              trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
            stop
          end if
        end if

        node_lines(i_node,1)=line_number
        node_lines(i_node,2)=0
        node_lines(i_node,3)=0

        in_node=.true.

      case ('nset')
        !------
        ! *Nset
        !------
        i_nset=i_nset+1
        paramstr='nset';     nset_nset(i_nset)    =get_parameter_string(input_line,paramstr)
        paramstr='elset';    nset_elset(i_nset)   =get_parameter_string(input_line,paramstr)
        paramstr='instance'; nset_instance(i_nset)=get_parameter_string(input_line,paramstr)
        paramstr='generate'; nset_generate(i_nset)=get_parameter_logical(input_line,paramstr)
        paramstr='internal'; nset_internal(i_nset)=get_parameter_logical(input_line,paramstr)
        paramstr='unsorted'; nset_unsorted(i_nset)=get_parameter_logical(input_line,paramstr)
       
        nset_lines(i_nset,1)=line_number
        nset_lines(i_nset,2)=0
        nset_lines(i_nset,3)=0

        if (.not. in_part .and. .not. in_instance .and. .not. in_assembly) then
          write(*,'(5(a))') 'ERROR: in node set definition: incorrect level (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if

        if (in_part .and. .not. in_instance) then
          nset_part(i_nset)=part_name(i_part)
        end if

        if (in_instance .and. .not. in_part) then
          nset_inst(i_nset)=instance_name(i_instance)
          nset_part(i_nset)=instance_part(i_instance)
        end if

        if (in_assembly .and. .not. in_instance) then
          nset_inst(i_nset)=nset_instance(i_nset)
          nset_assm(i_nset)=assembly_name(i_assembly)
        end if

        in_nset=.true.

      case ('part')
        !------
        ! *Part
        !------
        i_part=i_part+1
        paramstr='name'; part_name(i_part)=get_parameter_string(input_line,paramstr)

        if (in_assembly .or. in_instance .or. in_part) then
          write(*,'(5(a))') 'ERROR: in part definition: incorrect level (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if

        if (len_trim(part_name(i_part)) == 0) then
          write(*,'(5(a))') 'ERROR: in part definition: name parameter is required (file ', &
            trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
          stop
        end if

        part_lines(i_part,1)=line_number

        in_part=.true.

      case ('end part')
        part_lines(i_part,2)=line_number
        in_part=.false.

      case ('system')
        i_system=i_system+1

      case default
        ! Do nothing

      end select

    else

      !---------------------------------
      ! Not on a comment or keyword line
      !---------------------------------
      if (in_node) then
        !--------------------
        ! On a node data line
        !--------------------
        if (node_lines(i_node,2) == 0) then
          node_lines(i_node,2)=line_number
          node_lines(i_node,3)=line_number
        else
          node_lines(i_node,3)=line_number
        end if
      end if

      if (in_nset) then
        !---------------------
        ! On an nset data line
        !---------------------
        if (nset_lines(i_nset,2) == 0) then
          nset_lines(i_nset,2)=line_number
          nset_lines(i_nset,3)=line_number
        else
          nset_lines(i_nset,3)=line_number
        end if

        ! Is the generate-parameter user?
        if (.not. nset_generate(i_nset)) then
          nset_nodes(i_nset)=nset_nodes(i_nset)+number_integer_cols(input_line)
        else
          ! Parse an nset data line for the number of nodes
          if (number_integer_cols(input_line) /= 3) then
            write(*,'(5(a))') 'ERROR: in nset definition: invalid data (file ', &
              trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
            stop
          end if
          
          read(input_line,*,iostat=ios) n1,n2,i
          if (ios /= 0) call error_read_file(fem_input_file)

          ! Exception handling
          ! Assumption: mod(n2-n1,i) is not evaluated if i < 1
          if (n2 < n1  .or. n1 < 1 .or. n2 < 1 .or. i < 1 .or. mod(n2-n1,i) /= 0) then
            write(*,'(5(a))') 'ERROR: in nset definition: invalid data (file ', &
              trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
            stop
          end if

          nset_nodes(i_nset)=nset_nodes(i_nset)+(n2-n1)/i+1
        end if
      end if

      if (in_element) then
        !------------------------
        ! On an element data line
        !------------------------
        if (element_lines(i_element,2) == 0) then
          element_lines(i_element,2)=line_number
          element_lines(i_element,3)=line_number
        else
          element_lines(i_element,3)=line_number
        end if
      end if

      if (in_instance) then
        !-------------------------
        ! On an instance data line 
        !-------------------------
        
        ! Read instance translations and rotations
        from_instance=from_instance+1
        if (from_instance == 1) then
          if (.not. in_node .and. .not. in_nset) then
            ! Translation
            read(input_line,*,iostat=ios) (instance_translate(i_instance,j),j=1,3)
            if (ios /= 0) then
            write(*,'(5(a))') 'ERROR: in reading instance definition (file ', &
              trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
            end if
          end if
        else if (from_instance == 2) then
          if (.not. in_node .and. .not. in_nset) then
            ! Rotation
            read(input_line,*,iostat=ios) (instance_rotate(i_instance,j),j=1,7)
            if (ios /= 0) then
            write(*,'(5(a))') 'ERROR: in reading instance definition (file ', &
              trim(quote(fem_input_file)), ', line ', trim(int2str(line_number)), ')'
            end if
          end if
        end if
      end if 

    end if 

  end do read_loop

  !close(unit=iochannel(1))
  rewind(scratch_channel)

  !----------------
  ! Post-processing
  !----------------
  do i_node=1,n_node
    node_nodes(i_node)=node_lines(i_node,3)-node_lines(i_node,2)+1
  end do

  do i_element=1,n_element
    element_elements(i_element)=element_lines(i_element,3)-element_lines(i_element,2)+1
  end do

end subroutine parse_keywords

subroutine read_part_nodes
!---------------------------------------
! Read nodes defined in part definitions
!---------------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: j,ios,line_number
  integer :: i_node,i_part_node

  character(len=input_line_length) :: input_line

  !-------------
  ! Preparations
  !-------------

  n_part_nodes=0
  do i_node=1,n_node
    if (len_trim(node_part(i_node)) /= 0 .and. len_trim(node_instance(i_node)) == 0) then
      n_part_nodes=n_part_nodes+node_nodes(i_node)
    end if
  end do

  allocate(part_node_xyz(n_part_nodes,3),stat=ios);  call error_allocate(ios)
  allocate(part_node_number(n_part_nodes),stat=ios); call error_allocate(ios)
  allocate(part_node_name(n_part_nodes),stat=ios);   call error_allocate(ios)
  allocate(part_node_part(n_part_nodes),stat=ios);   call error_allocate(ios)

  !--------
  ! Actions
  !--------

  i_part_node=1
  node_loop: do i_node=1,n_node

    line_number=0
    read_loop: do
    
      read(scratch_channel,'(a)',iostat=ios) input_line
      select case(ios)
      case (:-1)
        exit read_loop
      case(1:)
        call error_read_file(fem_input_file)
      end select

      line_number=line_number+1

      if (line_number <= node_lines(i_node,3) .and. &
          line_number >= node_lines(i_node,2)) then

        if (len_trim(node_part(i_node)) /= 0 .and. len_trim(node_instance(i_node)) == 0) then
          read(input_line,*,iostat=ios) part_node_number(i_part_node), &
            (part_node_xyz(i_part_node,j),j=1,3)
          if (ios /= 0) call error_read_file(fem_input_file)
          
          part_node_part(i_part_node)=trim(node_part(i_node))
          part_node_name(i_part_node)=trim(node_part(i_node)) &
            // '.' // trim(int2str(part_node_number(i_part_node)))

          i_part_node=i_part_node+1
        end if

      end if

    end do read_loop

    rewind(scratch_channel)

  end do node_loop

end subroutine read_part_nodes

subroutine read_part_nsets
!-------------------------------------------
! Read node sets defined in part definitions
!-------------------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,j,n1,n2,ios,line_number
  integer :: i_nset,i_part_nset,ncols

  character(len=input_line_length) :: input_line

  !-------------
  ! Preparations
  !-------------
  n_part_nset_nodes=0
  do i_nset=1,n_nset
    if (len_trim(nset_part(i_nset)) /= 0 .and. len_trim(nset_inst(i_nset)) == 0) then
      n_part_nset_nodes=n_part_nset_nodes+nset_nodes(i_nset)
    end if
  end do

  allocate(part_nset_name(n_part_nset_nodes),stat=ios); call error_allocate(ios)
  allocate(part_nset_node(n_part_nset_nodes),stat=ios); call error_allocate(ios)
  allocate(part_nset_part(n_part_nset_nodes),stat=ios); call error_allocate(ios)

  !--------
  ! Actions
  !--------
  i_part_nset=1
  nset_loop: do i_nset=1,n_nset

    line_number=0
    read_loop: do
    
      read(scratch_channel,'(a)',iostat=ios) input_line
      select case(ios)
      case (:-1)
        exit read_loop
      case(1:)
        call error_read_file(fem_input_file)
      end select

      line_number=line_number+1

      if (line_number <= nset_lines(i_nset,3) .and. &
          line_number >= nset_lines(i_nset,2)) then

        if (len_trim(nset_part(i_nset)) /= 0 .and. len_trim(nset_inst(i_nset)) == 0) then

          if (.not. nset_generate(i_nset)) then
            ncols=number_integer_cols(input_line)
            part_nset_name(i_part_nset:(i_part_nset+ncols-1))=nset_nset(i_nset)
            part_nset_part(i_part_nset:(i_part_nset+ncols-1))=nset_part(i_nset)
         
            read(input_line,*,iostat=ios) (part_nset_node(i),i=i_part_nset,(i_part_nset+ncols-1))
            if (ios /= 0) call error_read_file(fem_input_file)

            i_part_nset=i_part_nset+ncols
          else
            read(input_line,*,iostat=ios) n1,n2,i
            if (ios /= 0) call error_read_file(fem_input_file)

            ncols=(n2-n1)/i+1 

            part_nset_name(i_part_nset:(i_part_nset+ncols-1))=nset_nset(i_nset)
            part_nset_part(i_part_nset:(i_part_nset+ncols-1))=nset_part(i_nset)

            do j=i_part_nset,(i_part_nset+ncols-1)
              part_nset_node(j)=n1+(j-i_part_nset)*i
            end do

            i_part_nset=i_part_nset+ncols
          end if
        end if

      end if

    end do read_loop

    rewind(scratch_channel)

  end do nset_loop

end subroutine read_part_nsets

subroutine read_part_elements()
!------------------------------------------
! Read elements defined in part definitions
!------------------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: j,ios,line_number
  integer :: i_element,i_part_element

  character(len=chr80) :: etype
  character(len=input_line_length) :: input_line

  !-------------
  ! Preparations
  !-------------

  n_part_elements=0
  do i_element=1,n_element
    if (len_trim(element_part(i_element)) /= 0 .and. len_trim(element_instance(i_element)) == 0) then
      n_part_elements=n_part_elements+element_elements(i_element)
    end if
  end do

  allocate(part_element_number(n_part_elements),stat=ios); call error_allocate(ios)
  allocate(part_element_type(n_part_elements),stat=ios);   call error_allocate(ios)
  allocate(part_element_node(n_part_elements,8),stat=ios); call error_allocate(ios)
  allocate(part_element_part(n_part_elements),stat=ios);   call error_allocate(ios)

  !--------
  ! Actions
  !--------

  i_part_element=1
  element_loop: do i_element=1,n_element

    line_number=0
    read_loop: do
    
      read(scratch_channel,'(a)',iostat=ios) input_line
      select case(ios)
      case (:-1)
        exit read_loop
      case(1:)
        call error_read_file(fem_input_file)
      end select

      line_number=line_number+1

      if (line_number <= element_lines(i_element,3) .and. &
          line_number >= element_lines(i_element,2)) then

        if (len_trim(element_part(i_element)) /= 0 .and. len_trim(element_instance(i_element)) == 0) then
          etype=trim(lowercase(element_type(i_element)))
          
          if (trim(etype) == 'dc3d8'  .or. &
              trim(etype) == 'dcc3d8' .or. &
              trim(etype) == 'dcc3d8d') then
            !----------------------
            ! Linear brick elements
            !----------------------
            read(input_line,*,iostat=ios) part_element_number(i_part_element), &
              (part_element_node(i_part_element,j),j=1,8)
            if (ios /= 0) call error_read_file(fem_input_file)
            part_element_type(i_part_element)=1

          else if (trim(etype) == 'dc3d4') then
            !---------------------------
            ! Linear tetrahedral element
            !---------------------------
            read(input_line,*,iostat=ios) part_element_number(i_part_element), &
              (part_element_node(i_part_element,j),j=1,4)
            if (ios /= 0) call error_read_file(fem_input_file)
            part_element_type(i_part_element)=2

          else if (trim(etype) == 'ds3') then
            !--------------------------------
            ! Linear triangular shell element
            !--------------------------------
            read(input_line,*,iostat=ios) part_element_number(i_part_element), &
              (part_element_node(i_part_element,j),j=1,3)
            if (ios /= 0) call error_read_file(fem_input_file)
            part_element_type(i_part_element)=3

          else if (trim(etype) == 'ds4') then
            !-----------------------------------
            ! Linear quadrilateral shell element
            !-----------------------------------
           read(input_line,*,iostat=ios) part_element_number(i_part_element), &
              (part_element_node(i_part_element,j),j=1,4)
            if (ios /= 0) call error_read_file(fem_input_file)
            part_element_type(i_part_element)=4

          end if
            
          part_element_part(i_part_element)=trim(element_part(i_element))
          
          i_part_element=i_part_element+1
        end if

      end if

    end do read_loop

    rewind(scratch_channel)

  end do element_loop

end subroutine read_part_elements

subroutine read_instance_nodes
!------------------------------------------------
! Read nodes defined in part instance definitions
!------------------------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: j,ios,line_number
  integer :: i_node,i_instance_node

  character(len=input_line_length) :: input_line

  !-------------
  ! Preparations
  !-------------

  n_instance_nodes=0
  do i_node=1,n_node
    if (len_trim(node_instance(i_node)) /= 0) then
      n_instance_nodes=n_instance_nodes+node_nodes(i_node)
    end if
  end do

  allocate(instance_node_xyz(n_instance_nodes,3),stat=ios);    call error_allocate(ios)
  allocate(instance_node_number(n_instance_nodes),stat=ios);   call error_allocate(ios)
  allocate(instance_node_name(n_instance_nodes),stat=ios);     call error_allocate(ios)
  allocate(instance_node_instance(n_instance_nodes),stat=ios); call error_allocate(ios)
  allocate(instance_node_part(n_instance_nodes),stat=ios);     call error_allocate(ios)

  !--------
  ! Actions
  !--------

  i_instance_node=1
  node_loop: do i_node=1,n_node

    line_number=0
    read_loop: do
    
      read(scratch_channel,'(a)',iostat=ios) input_line
      select case(ios)
      case (:-1)
        exit read_loop
      case(1:)
        call error_read_file(fem_input_file)
      end select

      line_number=line_number+1

      if (line_number <= node_lines(i_node,3) .and. &
          line_number >= node_lines(i_node,2)) then

        if (len_trim(node_instance(i_node)) /= 0) then
          read(input_line,*,iostat=ios) instance_node_number(i_instance_node), &
            (instance_node_xyz(i_instance_node,j),j=1,3)
          if (ios /= 0) call error_read_file(fem_input_file)
          
          instance_node_instance(i_instance_node)=trim(node_instance(i_node))
          instance_node_part(i_instance_node)=trim(node_part(i_node))
          instance_node_name(i_instance_node)=trim(node_instance(i_node)) &
            // '.' // trim(int2str(instance_node_number(i_instance_node)))

          i_instance_node=i_instance_node+1
        end if

      end if

    end do read_loop

    rewind(scratch_channel)

  end do node_loop

end subroutine read_instance_nodes

subroutine read_instance_nsets
!-----------------------------------------------
! Read node sets defined in instance definitions
!-----------------------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,j,n1,n2,ios,line_number
  integer :: i_nset,i_instance_nset,ncols

  character(len=input_line_length) :: input_line

  !-------------
  ! Preparations
  !-------------
  n_instance_nset_nodes=0
  do i_nset=1,n_nset
    if (len_trim(nset_inst(i_nset)) /= 0) then
      n_instance_nset_nodes=n_instance_nset_nodes+nset_nodes(i_nset)
    end if
  end do

  allocate(instance_nset_name(n_instance_nset_nodes),stat=ios);     call error_allocate(ios)
  allocate(instance_nset_node(n_instance_nset_nodes),stat=ios);     call error_allocate(ios)
  allocate(instance_nset_instance(n_instance_nset_nodes),stat=ios); call error_allocate(ios)
  allocate(instance_nset_part(n_instance_nset_nodes),stat=ios);     call error_allocate(ios)

  !--------
  ! Actions
  !--------
  i_instance_nset=1
  nset_loop: do i_nset=1,n_nset

    open(unit=iochannel(1),file=trim(fem_input_file),status='old',iostat=ios)
    if (ios /= 0) call error_open_file(fem_input_file)

    line_number=0
    read_loop: do
    
      read(iochannel(1),'(a)',iostat=ios) input_line
      select case(ios)
      case (:-1)
        exit read_loop
      case(1:)
        call error_read_file(fem_input_file)
      end select

      line_number=line_number+1

      if (line_number <= nset_lines(i_nset,3) .and. &
          line_number >= nset_lines(i_nset,2)) then

        if (len_trim(nset_inst(i_nset)) /= 0) then
          if (.not. nset_generate(i_nset)) then
            ncols=number_integer_cols(input_line)
            instance_nset_name(i_instance_nset:(i_instance_nset+ncols-1))=nset_nset(i_nset)
            instance_nset_instance(i_instance_nset:(i_instance_nset+ncols-1))=nset_inst(i_nset)
            instance_nset_part(i_instance_nset:(i_instance_nset+ncols-1))=nset_part(i_nset)
            
            read(input_line,*,iostat=ios) (instance_nset_node(i),i=i_instance_nset,(i_instance_nset+ncols-1))
            if (ios /= 0) call error_read_file(fem_input_file)
            
            i_instance_nset=i_instance_nset+ncols
          else
            read(input_line,*,iostat=ios) n1,n2,i
            if (ios /= 0) call error_read_file(fem_input_file)

            ncols=(n2-n1)/i+1

            instance_nset_name(i_instance_nset:(i_instance_nset+ncols-1))=nset_nset(i_nset)
            instance_nset_instance(i_instance_nset:(i_instance_nset+ncols-1))=nset_inst(i_nset)
            instance_nset_part(i_instance_nset:(i_instance_nset+ncols-1))=nset_part(i_nset)

            do j=i_instance_nset,(i_instance_nset+ncols-1)
              instance_nset_node(j)=n1+(j-i_instance_nset)*i
            end do

            i_instance_nset=i_instance_nset+ncols
          end if
        end if

      end if

    end do read_loop

    close(unit=iochannel(1))

  end do nset_loop

end subroutine read_instance_nsets

subroutine read_instance_elements()
!---------------------------------------------------
! Read elements defined in instance instance definitions
!---------------------------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: j,ios,line_number
  integer :: i_element,i_instance_element

  character(len=chr80) :: etype
  character(len=input_line_length) :: input_line

  !-------------
  ! Preparations
  !-------------

  n_instance_elements=0
  do i_element=1,n_element
    if (len_trim(element_instance(i_element)) /= 0 .and. len_trim(element_instance(i_element)) == 0) then
      n_instance_elements=n_instance_elements+element_elements(i_element)
    end if
  end do

  allocate(instance_element_number(n_instance_elements),stat=ios);   call error_allocate(ios)
  allocate(instance_element_type(n_instance_elements),stat=ios);     call error_allocate(ios)
  allocate(instance_element_node(n_instance_elements,8),stat=ios);   call error_allocate(ios)
  allocate(instance_element_part(n_instance_elements),stat=ios);     call error_allocate(ios)
  allocate(instance_element_instance(n_instance_elements),stat=ios); call error_allocate(ios)

  !--------
  ! Actions
  !--------

  i_instance_element=1
  element_loop: do i_element=1,n_element

    line_number=0
    read_loop: do
    
      read(scratch_channel,'(a)',iostat=ios) input_line
      select case(ios)
      case (:-1)
        exit read_loop
      case(1:)
        call error_read_file(fem_input_file)
      end select

      line_number=line_number+1

      if (line_number <= element_lines(i_element,3) .and. &
          line_number >= element_lines(i_element,2)) then

        if (len_trim(element_instance(i_element)) /= 0 .and. len_trim(element_instance(i_element)) == 0) then
          etype=trim(lowercase(element_type(i_element)))
          
          if (trim(etype) == 'dc3d8'  .or. &
              trim(etype) == 'dcc3d8' .or. &
              trim(etype) == 'dcc3d8d') then
            !----------------------
            ! Linear brick elements
            !----------------------
            read(input_line,*,iostat=ios) instance_element_number(i_instance_element), &
              (instance_element_node(i_instance_element,j),j=1,8)
            if (ios /= 0) call error_read_file(fem_input_file)
            instance_element_type(i_instance_element)=1

          else if (trim(etype) == 'dc3d4') then
            !---------------------------
            ! Linear tetrahedral element
            !---------------------------
            read(input_line,*,iostat=ios) instance_element_number(i_instance_element), &
              (instance_element_node(i_instance_element,j),j=1,4)
            if (ios /= 0) call error_read_file(fem_input_file)
            instance_element_type(i_instance_element)=2

          else if (trim(etype) == 'ds3') then
            !--------------------------------
            ! Linear triangular shell element
            !--------------------------------
            read(input_line,*,iostat=ios) instance_element_number(i_instance_element), &
              (instance_element_node(i_instance_element,j),j=1,3)
            if (ios /= 0) call error_read_file(fem_input_file)
            instance_element_type(i_instance_element)=3

          else if (trim(etype) == 'ds4') then
            !-----------------------------------
            ! Linear quadrilateral shell element
            !-----------------------------------
           read(input_line,*,iostat=ios) instance_element_number(i_instance_element), &
              (instance_element_node(i_instance_element,j),j=1,4)
            if (ios /= 0) call error_read_file(fem_input_file)
            instance_element_type(i_instance_element)=4

          end if

          instance_element_part(i_instance_element)=trim(element_part(i_element))
          instance_element_instance(i_instance_element)=trim(element_instance(i_element))
          
          i_instance_element=i_instance_element+1
        end if

      end if

    end do read_loop

    rewind(scratch_channel)

  end do element_loop

end subroutine read_instance_elements

subroutine read_assembly_nsets
!------------------------------------------------------------------
! Read node sets defined outside part and part instance definitions
!------------------------------------------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use string_handling
  implicit none

  integer :: i,j,n1,n2,ios,line_number
  integer :: i_nset,i_assembly_nset,ncols
  integer :: i_instance

  character(len=input_line_length) :: input_line

  !-------------
  ! Preparations
  !-------------
  n_assembly_nset_nodes=0
  do i_nset=1,n_nset
    if (len_trim(nset_part(i_nset)) == 0 .and. &
        len_trim(nset_inst(i_nset)) == 0 .and. &
        len_trim(nset_assm(i_nset)) /= 0) then
      n_assembly_nset_nodes=n_assembly_nset_nodes+nset_nodes(i_nset)
    end if
  end do

  allocate(assembly_nset_name(n_assembly_nset_nodes),stat=ios);     call error_allocate(ios)
  allocate(assembly_nset_node(n_assembly_nset_nodes),stat=ios);     call error_allocate(ios)
  allocate(assembly_nset_instance(n_assembly_nset_nodes),stat=ios); call error_allocate(ios)
  allocate(assembly_nset_part(n_assembly_nset_nodes),stat=ios);     call error_allocate(ios)

  assembly_nset_name     = ''
  assembly_nset_node     = 0
  assembly_nset_instance = ''
  assembly_nset_part     = ''

  !--------
  ! Actions
  !--------
  i_assembly_nset=1
  nset_loop: do i_nset=1,n_nset

    open(unit=iochannel(1),file=trim(fem_input_file),status='old',iostat=ios)
    if (ios /= 0) call error_open_file(fem_input_file)

    line_number=0
    read_loop: do
    
      read(iochannel(1),'(a)',iostat=ios) input_line
      select case(ios)
      case (:-1)
        exit read_loop
      case(1:)
        call error_read_file(fem_input_file)
      end select

      line_number=line_number+1

      if (line_number <= nset_lines(i_nset,3) .and. &
          line_number >= nset_lines(i_nset,2)) then

        if (len_trim(nset_part(i_nset)) == 0 .and. &
            len_trim(nset_inst(i_nset)) == 0 .and. &
            len_trim(nset_assm(i_nset)) /= 0) then
          if (.not. nset_generate(i_nset)) then
            ncols=number_integer_cols(input_line)
            
            assembly_nset_name(i_assembly_nset:(i_assembly_nset+ncols-1))=nset_nset(i_nset)
            assembly_nset_instance(i_assembly_nset:(i_assembly_nset+ncols-1))=nset_inst(i_nset)

            do i_instance=1,n_instance
              if (trim(instance_name(i_instance)) == trim(nset_inst(i_nset))) then
                assembly_nset_part(i_assembly_nset:(i_assembly_nset+ncols-1))=instance_part(i_instance)
              end if
            end do
            
            read(input_line,*,iostat=ios) (assembly_nset_node(i),i=i_assembly_nset,(i_assembly_nset+ncols-1))
            if (ios /= 0) call error_read_file(fem_input_file)
            
            i_assembly_nset=i_assembly_nset+ncols
          else
            read(input_line,*,iostat=ios) n1,n2,i
            if (ios /= 0) call error_read_file(fem_input_file)

            ncols=(n2-n1)/i+1

            assembly_nset_name(i_assembly_nset:(i_assembly_nset+ncols-1))=nset_nset(i_nset)
            assembly_nset_instance(i_assembly_nset:(i_assembly_nset+ncols-1))=nset_inst(i_nset)

            do i_instance=1,n_instance
              if (trim(instance_name(i_instance)) == trim(nset_inst(i_nset))) then
                assembly_nset_part(i_assembly_nset:(i_assembly_nset+ncols-1))=instance_part(i_instance)
              end if
            end do

            do j=i_assembly_nset,(i_assembly_nset+ncols-1)
              assembly_nset_node(j)=n1+(j-i_assembly_nset)*i
            end do

            i_assembly_nset=i_assembly_nset+ncols
          end if
        end if

      end if

    end do read_loop

    close(unit=iochannel(1))

  end do nset_loop

end subroutine read_assembly_nsets

subroutine inherit_instance_nodes()
!----------------------------------------------
! Instances inherit node definitions from parts
!----------------------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use string_handling
  implicit none
  
  integer :: i,ios,i_instance,inode_model,nnodes_model
  character(len=chr80) :: stmp80

  !----------------------------------------------
  ! Count the number of nodes in the ABAQUS model
  !----------------------------------------------

  ! Nodes defined on part instance level
  nnodes_model=ubound(instance_node_part,1)

  ! Nodes inherited from part level
  do i_instance=1,n_instance
    do i=1,ubound(part_node_part,1)
      if (trim(instance_part(i_instance)) == trim(part_node_part(i))) then
        nnodes_model=nnodes_model+1
      end if
    end do
  end do

  !--------------------------------------
  ! Allocate memory and initialize arrays
  !--------------------------------------

  allocate(model_node_name(nnodes_model),stat=ios);     call error_allocate(ios)
  allocate(model_node_part(nnodes_model),stat=ios);     call error_allocate(ios)
  allocate(model_node_instance(nnodes_model),stat=ios); call error_allocate(ios)
  allocate(model_node_number(nnodes_model),stat=ios);   call error_allocate(ios)
  allocate(model_node_xyz(nnodes_model,3),stat=ios);    call error_allocate(ios)

  model_node_name     = ''
  model_node_part     = ''
  model_node_instance = ''
  model_node_number   = 0
  model_node_xyz      = 0.0

  !--------------------------------------------------------------
  ! Copy data from nodes defined on part and part instance levels
  !--------------------------------------------------------------

  inode_model=1
  do i_instance=1,n_instance

    do i=1,ubound(instance_node_instance,1)
      if (trim(instance_node_instance(i)) == trim(instance_name(i_instance))) then

        model_node_name(inode_model)       = instance_node_name(i)
        model_node_part(inode_model)       = instance_node_part(i)
        model_node_instance(inode_model)   = instance_node_instance(i)
        model_node_number(inode_model)     = instance_node_number(i)
        model_node_xyz(inode_model,1:3)    = instance_node_xyz(i,1:3)
        
        inode_model=inode_model+1

      end if
    end do

    do i=1,ubound(part_node_part,1)

      if (trim(part_node_part(i)) == trim(instance_part(i_instance))) then
       
        stmp80=trim(instance_name(i_instance)) // '.' // int2str(part_node_number(i)) 

        model_node_name(inode_model)       = stmp80
        model_node_part(inode_model)       = part_node_part(i)
        model_node_instance(inode_model)   = instance_name(i_instance)
        model_node_number(inode_model)     = part_node_number(i)
        model_node_xyz(inode_model,1:3)    = part_node_xyz(i,1:3)

        inode_model=inode_model+1

      end if

    end do

  end do

  if (nnodes_model /= (inode_model-1)) then
    write(*,'(a)') 'ERROR: in subroutine inherit_instance_nodes()'
    stop
  end if

end subroutine inherit_instance_nodes

subroutine inherit_instance_nsets()
!--------------------------------------------
! Instances can inherit node sets from parts, 
! they can also be defined on assembly level
!--------------------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use string_handling
  implicit none
  
  integer :: i,ios,i_instance,inode_model,nnodes_model

  !---------------------------------------------------------------------
  ! Count the number of nodes belonging to node sets in the ABAQUS model
  !---------------------------------------------------------------------

  ! Node set nodes defined on part instance level
  nnodes_model=ubound(instance_nset_node,1)

  ! Node set nodes inherited from part level
  do i_instance=1,n_instance
    do i=1,ubound(part_nset_part,1)
      if (trim(instance_part(i_instance)) == trim(part_nset_part(i))) then
        nnodes_model=nnodes_model+1
      end if
    end do
  end do

  ! Node set nodes defined on assembly level, but not on part instance level
  nnodes_model=nnodes_model+ubound(assembly_nset_node,1)

  !--------------------------------------
  ! Allocate memory and initialize arrays
  !--------------------------------------

  allocate(model_nset_name(nnodes_model),stat=ios);     call error_allocate(ios)
  allocate(model_nset_part(nnodes_model),stat=ios);     call error_allocate(ios)
  allocate(model_nset_instance(nnodes_model),stat=ios); call error_allocate(ios)
  allocate(model_nset_node(nnodes_model),stat=ios);     call error_allocate(ios)

  model_nset_name     = ''
  model_nset_part     = ''
  model_nset_instance = ''
  model_nset_node     = 0

  !-----------------------------------------------------------------------
  ! Copy data from node set nodes defined on part and part instance levels
  !-----------------------------------------------------------------------

  inode_model=1
  do i_instance=1,n_instance

    do i=1,ubound(instance_nset_instance,1)
      if (trim(instance_nset_instance(i)) == trim(instance_name(i_instance))) then

        model_nset_name(inode_model)     = instance_nset_name(i)
        model_nset_part(inode_model)     = instance_nset_part(i)
        model_nset_instance(inode_model) = instance_nset_instance(i)
        model_nset_node(inode_model)     = instance_nset_node(i)
        
        inode_model=inode_model+1

      end if
    end do

    do i=1,ubound(part_nset_part,1)

      if (trim(part_nset_part(i)) == trim(instance_part(i_instance))) then
       
        model_nset_name(inode_model)     = part_nset_name(i)
        model_nset_part(inode_model)     = part_nset_part(i)
        model_nset_instance(inode_model) = instance_name(i_instance)
        model_nset_node(inode_model)     = part_nset_node(i)

        inode_model=inode_model+1

      end if

    end do

    do i=1,ubound(assembly_nset_instance,1)
      if (trim(assembly_nset_instance(i)) == trim(instance_name(i_instance))) then

        model_nset_name(inode_model)     = assembly_nset_name(i)
        model_nset_part(inode_model)     = assembly_nset_part(i)
        model_nset_instance(inode_model) = assembly_nset_instance(i)
        model_nset_node(inode_model)     = assembly_nset_node(i)
        
        inode_model=inode_model+1

      end if
    end do

  end do

  if (nnodes_model /= (inode_model-1)) then
    write(*,'(a)') 'ERROR: in subroutine inherit_instance_nsets()'
    stop
  end if

end subroutine inherit_instance_nsets

subroutine inherit_instance_elements()
!-------------------------------------------------
! Instances inherit element definitions from parts
!-------------------------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use string_handling
  implicit none
  
  integer :: i,ios,i_instance,ielement_model,nelements_model

  !-------------------------------------------------
  ! Count the number of elements in the ABAQUS model
  !-------------------------------------------------

  ! Elements defined on part instance level
  nelements_model=ubound(instance_element_instance,1)

  ! Elements inherited from part level
  do i_instance=1,n_instance
    do i=1,ubound(part_element_part,1)
      if (trim(instance_part(i_instance)) == trim(part_element_part(i))) then
        nelements_model=nelements_model+1
      end if
    end do
  end do

  !--------------------------------------
  ! Allocate memory and initialize arrays
  !--------------------------------------

  allocate(model_element_part(nelements_model),stat=ios);     call error_allocate(ios)
  allocate(model_element_instance(nelements_model),stat=ios); call error_allocate(ios)
  allocate(model_element_number(nelements_model),stat=ios);   call error_allocate(ios)
  allocate(model_element_type(nelements_model),stat=ios);   call error_allocate(ios)
  allocate(model_element_node(nelements_model,8),stat=ios);   call error_allocate(ios)

  model_element_part     = ''
  model_element_instance = ''
  model_element_number   = 0
  model_element_type     = 0
  model_element_node     = 0

  !-----------------------------------------------------------------
  ! Copy data from elements defined on part and part instance levels
  !-----------------------------------------------------------------

  ielement_model=1
  do i_instance=1,n_instance

    do i=1,ubound(instance_element_instance,1)
      if (trim(instance_element_instance(i)) == trim(instance_name(i_instance))) then

        model_element_part(ielement_model)       = instance_element_part(i)
        model_element_instance(ielement_model)   = instance_element_instance(i)
        model_element_number(ielement_model)     = instance_element_number(i)
        model_element_type(ielement_model)       = instance_element_type(i)
        model_element_node(ielement_model,1:8)   = instance_element_node(i,1:8)
        
        ielement_model=ielement_model+1

      end if
    end do

    do i=1,ubound(part_element_part,1)

      if (trim(part_element_part(i)) == trim(instance_part(i_instance))) then
       
        model_element_part(ielement_model)       = part_element_part(i)
        model_element_instance(ielement_model)   = instance_name(i_instance)
        model_element_number(ielement_model)     = part_element_number(i)
        model_element_type(ielement_model)       = part_element_type(i)
        model_element_node(ielement_model,1:8)   = part_element_node(i,1:8)

        ielement_model=ielement_model+1

      end if

    end do

  end do

  if (nelements_model /= (ielement_model-1)) then
    write(*,'(a)') 'ERROR: in subroutine inherit_instance_elements()'
    stop
  end if

end subroutine inherit_instance_elements

subroutine create_element_node_map()
!---------------------------------
! Create a global element-node map
!---------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use miscellaneous
  use string_handling
  implicit none
  
  integer :: i,j,k,ios,n_model_node,n_model_element
  logical :: found
  integer :: ifrac
  real(kind=rk) :: rfrac

  n_model_node=ubound(model_node_number,1)
  n_model_element=ubound(model_element_node,1)
  
  allocate(model_element_int_node(n_model_element,8),stat=ios); call error_allocate(ios)
  model_element_int_node=0

  ! Relevant arrays
  !----------------------------
  !   model_node_xyz(:),
  !   model_node_number(:),
  !   model_node_name(:),
  !   model_node_instance(:),
  !   model_nset_name(:),  
  !   model_nset_node(:),
  !   model_nset_instance(:)
  !   model_element_number(:),
  !   model_element_node(:,:),
  !   model_element_instance(:)
  !----------------------------

  ! Progress bar
  ifrac=int(n_model_element/10)
  rfrac=0.0

  element_loop: do i=1,n_model_element

    ! Progress bar
    if (i == 1 .or. mod(i,ifrac) == 0) then
      rfrac=real(i,kind=rk)/n_model_element
      call progress_line(rfrac)
    end if  

    node_loop_1: do j=1,8
      found=.false.
      if (model_element_node(i,j) == 0) cycle element_loop
      
      node_loop_2: do k=1,n_model_node 
        if (trim(model_element_instance(i)) == trim(model_node_instance(k))) then
          if (model_element_node(i,j) == model_node_number(k)) then
            model_element_int_node(i,j)=k
            found=.true.
            cycle node_loop_1
          end if
        end if
      end do node_loop_2

    end do node_loop_1
    
  end do element_loop

end subroutine create_element_node_map

subroutine create_node_surface_map()
!------------------------------------------------
! Calculate nodal areas for selected ABAQUS nodes
!------------------------------------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use mathematics
  use miscellaneous
  use string_handling
  implicit none

  integer :: i,j,k,l,m,n,ios,n_model_node,n_nodes_face
  integer :: n_model_nset,nnodes_patch,common_nodes
  integer :: n_model_element,n_faces,node_order,n_faces_max
  integer :: type_a,type_b,ifrac,i_face,n_face

  integer, dimension(8) :: node_a,node_b

  integer, dimension(:), allocatable :: node_sources

  logical :: found

  logical, dimension(6) :: face_a,face_b

  logical, dimension(3) :: fc3
  logical, dimension(4) :: fc4
  logical, dimension(8) :: fm

  logical, dimension(6) :: face
  logical, dimension(6) :: eface
  logical, dimension(6) :: nface

  logical, dimension(4,6,8) :: mask

  real(kind=rk) :: t_area,d_area,one_third,one_fourth,rfrac
  real(kind=rk), dimension(3) :: p1,p2,p3,p4
  real(kind=rk), dimension(3) :: p5,p6,p7,p8

  character(len=chr80) :: stmp

  logical, dimension(:),   allocatable :: model_nset_mask
  logical, dimension(:),   allocatable :: model_node_mask
  logical, dimension(:),   allocatable :: model_element_mask
  logical, dimension(:,:), allocatable :: model_element_node_mask
  logical, dimension(:,:), allocatable :: model_element_face_mask
  logical, dimension(:,:), allocatable :: model_element_face_mask_2
 
  !-------------
  ! Preparations
  !-------------

  one_third=1.0/3.0; one_fourth=1.0/4.0

  if (len_trim(nset_input_file) == 0) then
    !----------------------------------------
    ! Patch-node set correspondence not given
    !----------------------------------------
    return
  end if
  
  !------------------------------------
  ! Patch-node set correspondence given
  !------------------------------------

  ! Relevant arrays

  !   model_node_xyz(:),
  !   model_node_number(:),
  !   model_node_name(:),
  !   model_node_instance(:),
  !   model_nset_name(:),
  !   model_nset_node(:),
  !   model_nset_instance(:)
  !   model_element_number(:),
  !   model_element_type(:),
  !   model_element_node(:,:),
  !   model_element_instance(:)
  !   connectivity_table(:,:)

  !-------------------
  ! Create mask arrays
  !-------------------

  ! Node sets
  n_model_nset=ubound(model_nset_name,1)
  allocate(model_nset_mask(n_model_nset),stat=ios); call error_allocate(ios)
  model_nset_mask=.false.

  nnodes_abaqus=0
  count_loop: do i=1,n_model_nset
    nset_loop_0: do k=1,ubound(connectivity_table,1)

      stmp=trim(lowercase(model_nset_instance(i))) // '.' // trim(lowercase(model_nset_name(i)))
      if (trim(stmp) == trim(lowercase(connectivity_table(k,1)))) then
        nnodes_abaqus=nnodes_abaqus+1
        model_nset_mask(i)=.true.
        exit nset_loop_0
      end if
    
    end do nset_loop_0
  end do count_loop
  
  ! Nodes
  n_model_node=ubound(model_node_name,1)
  allocate(model_node_mask(n_model_node),stat=ios); call error_allocate(ios)
  model_node_mask=.false.

  do i=1,n_model_nset
    if (model_nset_mask(i) .eqv. .true.) then

      node_loop_0: do j=1,n_model_node
        if (trim(model_node_instance(j)) == trim(model_nset_instance(i))) then
          if (model_node_number(j) == model_nset_node(i)) then
            model_node_mask(j) = .true.
          end if
        end if
      end do node_loop_0

    end if
  end do
  
  ! Elements
  n_model_element=ubound(model_element_number,1)
  allocate(model_element_mask(n_model_element),stat=ios); call error_allocate(ios)
  model_element_mask=.false.
  allocate(model_element_node_mask(n_model_element,8),stat=ios); call error_allocate(ios)
  model_element_node_mask=.false.

  element_loop_0: do i=1,n_model_element
    column_loop_0: do j=1,8
      if (model_element_int_node(i,j) /= 0) then
        if (model_node_mask(model_element_int_node(i,j)) .eqv. .true.) then
          model_element_mask(i)=.true.
          model_element_node_mask(i,j)=.true.
        end if
      end if
    end do column_loop_0

  end do element_loop_0

  ! Identify inter-element surfaces
  allocate(model_element_face_mask(n_model_element,6),stat=ios); call error_allocate(ios)
  model_element_face_mask=.true.
  
  allocate(model_element_face_mask_2(n_model_element,6),stat=ios); call error_allocate(ios)
  model_element_face_mask_2=.false.

  mask=.false.

  mask(1,1,1:4)=.true.
  mask(1,2,5:8)=.true.
  mask(1,3,1:2)=.true.
  mask(1,3,5:6)=.true.
  mask(1,4,2:3)=.true.
  mask(1,4,6:7)=.true.
  mask(1,5,3:4)=.true.
  mask(1,5,7:8)=.true.
  mask(1,6,1)  =.true.
  mask(1,6,4:5)=.true.
  mask(1,6,8)  =.true.
  mask(2,1,1:3)=.true.
  mask(2,2,1:2)=.true.
  mask(2,2,4)  =.true.
  mask(2,3,2:4)=.true.
  mask(2,4,1)  =.true.
  mask(2,4,3:4)=.true.
  mask(3,1,1:3)=.true.
  mask(4,1,1:4)=.true.

  ! Progress bar
  ifrac=int((n_model_node+n_model_element)/10)
  rfrac=0.0

  element_loop_0a: do i=1,n_model_element

    ! Progress bar
    if (i == 1 .or. mod(i,ifrac) == 0) then
      rfrac=real(i,kind=rk)/(n_model_element+n_model_node)
      call progress_line(rfrac)
    end if  

    if (.not. model_element_mask(i)) cycle element_loop_0a
  
    type_a=model_element_type(i)
    face_a=.false.; node_a=0
    select case(type_a)
    case (1)
      ! Brick
      face_a(1:6)=.true.
      node_a(1:8)=model_element_int_node(i,1:8)
      n_nodes_face=4
    case (2)
      ! Tetrahedron
      face_a(1:4)=.true.
      node_a(1:4)=model_element_int_node(i,1:4)
      n_nodes_face=3
    case (3)
      ! Triangle
      face_a(1)=.true.
      node_a(1:3)=model_element_int_node(i,1:3)
      n_nodes_face=3
    case (4)
      ! Quadrangle
      face_a(1)=.true.
      node_a(1:4)=model_element_int_node(i,1:4)
      n_nodes_face=4
    end select  

    element_loop_0b: do j=1,n_model_element
      if (i == j) cycle element_loop_0b

      type_b=model_element_type(j)
      face_b=.false.; node_b=0
      select case(type_b)
      case (1)
        ! Brick
        face_b(1:6)=.true.
        node_b(1:8)=model_element_int_node(j,1:8)
      case (2)
        ! Tetrahedron
        face_b(1:4)=.true.
        node_b(1:4)=model_element_int_node(j,1:4)
      case (3)
        ! Triangle
        face_b(1)=.true.
        node_b(1:3)=model_element_int_node(j,1:3)
      case (4)
        ! Quadrangle
        face_b(1)=.true.
        node_b(1:4)=model_element_int_node(j,1:4)
      end select  
      
      ! Find out if these elements share a common face

      if (type_a == type_b) then
        ! Brick vs. brick

        face_loop_a1: do k=1,6
          if (.not. face_a(k)) cycle face_loop_a1
          if (model_element_face_mask(i,k) .eqv. .false.) cycle face_loop_a1

          face_loop_b1: do l=1,6
            if (.not. face_b(l)) cycle face_loop_b1
            if (model_element_face_mask(j,l) .eqv. .false.) cycle face_loop_b1

            common_nodes=0
            outer1: do m=1,8
              if (.not. mask(type_a,k,m)) cycle outer1

              inner1: do n=1,8
                if (.not. mask(type_b,l,n)) cycle inner1

                if (node_a(m) == node_b(n)) then
                  common_nodes=common_nodes+1
                end if

              end do inner1

            end do outer1
    
            if (common_nodes == n_nodes_face) then
              model_element_face_mask(i,k)=.false.
              model_element_face_mask(j,l)=.false.
            end if

          end do face_loop_b1
        end do face_loop_a1

      else if ((type_a == 1 .and. type_b == 2) .or. &
               (type_a == 2 .and. type_b == 1)) then
        ! Brick vs. tetrahedron

        face_loop_a2: do k=1,6
          if (.not. face_a(k)) cycle face_loop_a2
          if (model_element_face_mask(i,k) .eqv. .false.) cycle face_loop_a2

          face_loop_b2: do l=1,6
            if (.not. face_b(l)) cycle face_loop_b2
            if (model_element_face_mask(j,l) .eqv. .false.) cycle face_loop_b2

            common_nodes=0
            outer2: do m=1,8
              if (.not. mask(type_a,k,m)) cycle outer2

              inner2: do n=1,8
                if (.not. mask(type_b,l,n)) cycle inner2

                if (node_a(m) == node_b(n)) then
                  common_nodes=common_nodes+1
                end if

              end do inner2

            end do outer2
    
            if (common_nodes == 3) then
              model_element_face_mask(i,k)=.false.
              model_element_face_mask(j,l)=.false.
            end if

          end do face_loop_b2
        end do face_loop_a2

      end if

    end do element_loop_0b
  end do element_loop_0a

  !----------------------
  ! Calculate nodal areas
  !----------------------

  allocate(node_sources(n_model_node))
  node_sources=0

  allocate(model_node_area(n_model_node),stat=ios); call error_allocate(ios)
  model_node_area=0.0

  ! Loop through all nodes
  t_area=0.0
  node_loop_1: do i=1,n_model_node

    ! Progress bar
    if (mod((i+n_model_element),ifrac) == 0) then
      rfrac=(real(i,kind=rk)+n_model_element)/(n_model_node+n_model_element)
      call progress_line(rfrac)
    end if  

    ! If this node belongs to a selected node set
    if (model_node_mask(i) .eqv. .true.) then

      n_faces_max=0

      ! Loop through all elements
      element_loop_1: do j=1,n_model_element
        
        ! If this element has one or more nodes that
        ! belong to a selected node set
        if (model_element_mask(j) .eqv. .true.) then

          found=.false.
          column_loop_1: do k=1,8
            if (model_element_int_node(j,k) == i) then
              found=.true.; node_order=k
              exit column_loop_1
            end if
          end do column_loop_1

          ! If the current node doesn't belong to the current element
          if (.not. found) cycle element_loop_1

          ! How many nodes from this element belong to a selected node set
          nnodes_patch=count(model_element_node_mask(j,1:8))
          
          p1(1:3)=model_node_xyz(model_element_int_node(j,1),1:3)
          p2(1:3)=model_node_xyz(model_element_int_node(j,2),1:3)
          p3(1:3)=model_node_xyz(model_element_int_node(j,3),1:3)
          p4(1:3)=model_node_xyz(model_element_int_node(j,4),1:3)
          p5(1:3)=model_node_xyz(model_element_int_node(j,5),1:3)
          p6(1:3)=model_node_xyz(model_element_int_node(j,6),1:3)
          p7(1:3)=model_node_xyz(model_element_int_node(j,7),1:3)
          p8(1:3)=model_node_xyz(model_element_int_node(j,8),1:3)
          
          !-----------------------
          ! Selected element faces
          !-----------------------
          eface(1:6)=.false.
          select case(model_element_type(j))
          case (1)
            ! Brick
            if (nnodes_patch >= 4) then
              fm=model_element_node_mask(j,1:8)

              fc4(1)=fm(1); fc4(2)=fm(2); fc4(3)=fm(3); fc4(4)=fm(4)
              if (all(fc4)) eface(1)=.true.
              fc4(1)=fm(5); fc4(2)=fm(6); fc4(3)=fm(7); fc4(4)=fm(8)
              if (all(fc4)) eface(2)=.true.
              fc4(1)=fm(1); fc4(2)=fm(2); fc4(3)=fm(5); fc4(4)=fm(6)
              if (all(fc4)) eface(3)=.true.
              fc4(1)=fm(2); fc4(2)=fm(3); fc4(3)=fm(6); fc4(4)=fm(7)
              if (all(fc4)) eface(4)=.true.
              fc4(1)=fm(3); fc4(2)=fm(4); fc4(3)=fm(7); fc4(4)=fm(8)
              if (all(fc4)) eface(5)=.true.
              fc4(1)=fm(1); fc4(2)=fm(4); fc4(3)=fm(5); fc4(4)=fm(8)
              if (all(fc4)) eface(6)=.true.

            end if

          case (2)
            ! Tetrahedron
            if (nnodes_patch >= 3) then
              fm=model_element_node_mask(j,1:8)

              fc3(1)=fm(1); fc3(2)=fm(2); fc3(3)=fm(3)
              if (all(fc3)) eface(1)=.true.
              fc3(1)=fm(1); fc3(2)=fm(2); fc3(3)=fm(4)
              if (all(fc3)) eface(2)=.true.
              fc3(1)=fm(2); fc3(2)=fm(3); fc3(3)=fm(4)
              if (all(fc3)) eface(3)=.true.
              fc3(1)=fm(1); fc3(2)=fm(3); fc3(3)=fm(4)
              if (all(fc3)) eface(4)=.true.
              
            end if

          case (3)
            ! Triangle
              fm=model_element_node_mask(j,1:8)

              fc3(1)=fm(1); fc3(2)=fm(2); fc3(3)=fm(3)
              if (all(fc3)) eface(1)=.true.
              
          case (4)
            ! Quadrangle
              fm=model_element_node_mask(j,1:8)

              fc4(1)=fm(1); fc4(2)=fm(2); fc4(3)=fm(3); fc4(4)=fm(4)
              if (all(fc4)) eface(1)=.true.

          end select

          !--------------------------------
          ! Faces contributing to this node
          !--------------------------------
          nface(1:6)=.false.
          select case(model_element_type(j))
          case (1)
            ! Brick
            select case (node_order)
            case (1)
              nface(1)=.true.; nface(3)=.true.; nface(6)=.true.
            case (2)
              nface(1)=.true.; nface(3)=.true.; nface(4)=.true.
            case (3)
              nface(1)=.true.; nface(4)=.true.; nface(5)=.true.
            case (4)
              nface(1)=.true.; nface(5)=.true.; nface(6)=.true.
            case (5)
              nface(2)=.true.; nface(3)=.true.; nface(6)=.true.
            case (6)
              nface(2)=.true.; nface(3)=.true.; nface(4)=.true.
            case (7)
              nface(2)=.true.; nface(4)=.true.; nface(5)=.true.
            case (8)
              nface(2)=.true.; nface(5)=.true.; nface(6)=.true.
            end select

          case (2)
            ! Tetrahedron
            select case (node_order)
            case (1)
              nface(1)=.true.; nface(2)=.true.; nface(4)=.true.
            case (2)
              nface(1)=.true.; nface(2)=.true.; nface(3)=.true.
            case (3)
              nface(1)=.true.; nface(3)=.true.; nface(4)=.true.
            case (4)
              nface(2)=.true.; nface(3)=.true.; nface(4)=.true. 
            end select

          case (3)
            ! Triangle
            nface(1)=.true. 
          case (4)
            ! Quadrangle
            nface(1)=.true.

          end select

          !----------------
          ! Combined effect
          !----------------

          face(1)=((eface(1) .and. nface(1)) .and. model_element_face_mask(j,1))
          face(2)=((eface(2) .and. nface(2)) .and. model_element_face_mask(j,2))
          face(3)=((eface(3) .and. nface(3)) .and. model_element_face_mask(j,3))
          face(4)=((eface(4) .and. nface(4)) .and. model_element_face_mask(j,4))
          face(5)=((eface(5) .and. nface(5)) .and. model_element_face_mask(j,5))
          face(6)=((eface(6) .and. nface(6)) .and. model_element_face_mask(j,6))

          model_element_face_mask_2(j,1)=((eface(1) .and. model_element_face_mask(j,1)))
          model_element_face_mask_2(j,2)=((eface(2) .and. model_element_face_mask(j,2)))
          model_element_face_mask_2(j,3)=((eface(3) .and. model_element_face_mask(j,3)))
          model_element_face_mask_2(j,4)=((eface(4) .and. model_element_face_mask(j,4)))
          model_element_face_mask_2(j,5)=((eface(5) .and. model_element_face_mask(j,5)))
          model_element_face_mask_2(j,6)=((eface(6) .and. model_element_face_mask(j,6)))
          
          n_faces=count(face)
          if (n_faces > n_faces_max) n_faces_max=n_faces

          !-----------------------------------------------
          ! Determine element's contribution to nodal area
          !-----------------------------------------------

          d_area=0.0
          select case (nnodes_patch)
          case (0)
            !------------------
            ! Should not happen
            !------------------
            write(*,'(a)') 'WARNING: from subroutine create_node_surface(): type 0'

          case (1)
            !---------
            ! One node
            !---------
            ! No area contribution from this element
          case (2)
            !----------
            ! Two nodes
            !----------

            ! No area contribution from this element

          case (3)
            !------------
            ! Three nodes
            !------------
            select case (model_element_type(j))
            case (1)
              ! Brick (no area contribution)

            case (2)
              ! Tetrahedral
              if (face(1)) d_area=d_area+one_third*area_triangle(p1,p2,p3)
              if (face(2)) d_area=d_area+one_third*area_triangle(p1,p2,p4)
              if (face(3)) d_area=d_area+one_third*area_triangle(p2,p3,p4)
              if (face(4)) d_area=d_area+one_third*area_triangle(p1,p3,p4)

            case (3)
              ! Triangular
              if (face(1)) d_area=d_area+one_third*area_triangle(p1,p2,p3)

            case (4)
              ! Quadrilateral (no area contribution)  

            end select

          case (4)
            !-----------
            ! Four nodes
            !-----------
            select case (model_element_type(j))
            case (1)
              ! Brick              
              if (face(1)) d_area=d_area+one_fourth*area_quadrangle(p1,p2,p3,p4)
              if (face(2)) d_area=d_area+one_fourth*area_quadrangle(p5,p6,p7,p8)
              if (face(3)) d_area=d_area+one_fourth*area_quadrangle(p1,p2,p6,p5)
              if (face(4)) d_area=d_area+one_fourth*area_quadrangle(p2,p3,p7,p6)
              if (face(5)) d_area=d_area+one_fourth*area_quadrangle(p3,p4,p8,p7)
              if (face(6)) d_area=d_area+one_fourth*area_quadrangle(p1,p4,p8,p5)

            case (2)
              ! Tetrahedral

              if (face(1)) d_area=d_area+one_third*area_triangle(p1,p2,p3)
              if (face(2)) d_area=d_area+one_third*area_triangle(p1,p2,p4)
              if (face(3)) d_area=d_area+one_third*area_triangle(p2,p3,p4)
              if (face(4)) d_area=d_area+one_third*area_triangle(p1,p3,p4)

            case (3)
              ! Triangular
              write(*,'(a)') 'WARNING: from subroutine create_node_surface(): type 4.3'

            case (4)
              ! Quadrilateral
               if (face(1)) d_area=d_area+one_fourth*area_quadrangle(p1,p2,p3,p4)

            end select

          case (5)
            !-----------
            ! Five nodes
            !-----------
            select case (model_element_type(j))
            case (1)
              ! Brick
              if (face(1)) d_area=d_area+one_fourth*area_quadrangle(p1,p2,p3,p4)
              if (face(2)) d_area=d_area+one_fourth*area_quadrangle(p5,p6,p7,p8)
              if (face(3)) d_area=d_area+one_fourth*area_quadrangle(p1,p2,p6,p5)
              if (face(4)) d_area=d_area+one_fourth*area_quadrangle(p2,p3,p7,p6)
              if (face(5)) d_area=d_area+one_fourth*area_quadrangle(p3,p4,p8,p7)
              if (face(6)) d_area=d_area+one_fourth*area_quadrangle(p1,p4,p8,p5)

            case default
              write(*,'(a)') 'WARNING: from subroutine create_node_surface(): type 5'
            end select

          case (6)
            !----------
            ! Six nodes
            !----------
            select case (model_element_type(j))
            case (1)
              ! Brick

              if (face(1)) d_area=d_area+one_fourth*area_quadrangle(p1,p2,p3,p4)
              if (face(2)) d_area=d_area+one_fourth*area_quadrangle(p5,p6,p7,p8)
              if (face(3)) d_area=d_area+one_fourth*area_quadrangle(p1,p2,p6,p5)
              if (face(4)) d_area=d_area+one_fourth*area_quadrangle(p2,p3,p7,p6)
              if (face(5)) d_area=d_area+one_fourth*area_quadrangle(p3,p4,p8,p7)
              if (face(6)) d_area=d_area+one_fourth*area_quadrangle(p1,p4,p8,p5)

            case default
              write(*,'(a)') 'WARNING: from subroutine create_node_surface(): type 6'
            end select

          case (7)
            !------------
            ! Seven nodes
            !------------
            ! Three quadrangles

            select case (model_element_type(j))
            case (1)
              ! Brick
              if (face(1)) d_area=d_area+one_fourth*area_quadrangle(p1,p2,p3,p4)
              if (face(2)) d_area=d_area+one_fourth*area_quadrangle(p5,p6,p7,p8)
              if (face(3)) d_area=d_area+one_fourth*area_quadrangle(p1,p2,p6,p5)
              if (face(4)) d_area=d_area+one_fourth*area_quadrangle(p2,p3,p7,p6)
              if (face(5)) d_area=d_area+one_fourth*area_quadrangle(p3,p4,p8,p7)
              if (face(6)) d_area=d_area+one_fourth*area_quadrangle(p1,p4,p8,p5)

            case default
              write(*,'(a)') 'WARNING: from subroutine create_node_surface(): type 7'
            end select

          case (8)
            !------------
            ! Eight nodes
            !------------
            select case (model_element_type(j))
            case (1)
              ! Brick

              ! N.b. This needs to be fixed

              if (face(1)) d_area=d_area+one_fourth*area_quadrangle(p1,p2,p3,p4)
              if (face(2)) d_area=d_area+one_fourth*area_quadrangle(p5,p6,p7,p8)
              if (face(3)) d_area=d_area+one_fourth*area_quadrangle(p1,p2,p6,p5)
              if (face(4)) d_area=d_area+one_fourth*area_quadrangle(p2,p3,p7,p6)
              if (face(5)) d_area=d_area+one_fourth*area_quadrangle(p3,p4,p8,p7)
              if (face(6)) d_area=d_area+one_fourth*area_quadrangle(p1,p4,p8,p5)

            case default
              write(*,'(a)') 'WARNING: from subroutine create_node_surface(): type 8'
            end select

          end select

          model_node_area(i)=model_node_area(i)+d_area
          t_area=t_area+d_area

        end if      

      end do element_loop_1

    end if

  end do node_loop_1

  !-----------------------------
  ! Calculate element face areas
  !-----------------------------

  n_face=0
  do i=1,n_model_element
    n_face=n_face+count(model_element_face_mask_2(i,1:6))
  end do

  if (n_face <= 0) stop

  allocate(model_face_node(n_face,4),stat=ios);          call error_allocate(ios)
  allocate(model_face_type(n_face),stat=ios);            call error_allocate(ios)
  allocate(model_face_area(n_face),stat=ios);            call error_allocate(ios)
  allocate(model_node_face_mask(n_model_node),stat=ios); call error_allocate(ios)

  model_face_node=0
  model_face_type=0
  model_face_area=0.0
  model_node_face_mask=.false.

  i_face=1
  do i=1,n_model_element
    if (.not. any(model_element_face_mask_2(i,1:6))) cycle
    
    p1(1:3)=model_node_xyz(model_element_int_node(i,1),1:3)
    p2(1:3)=model_node_xyz(model_element_int_node(i,2),1:3)
    p3(1:3)=model_node_xyz(model_element_int_node(i,3),1:3)
    p4(1:3)=model_node_xyz(model_element_int_node(i,4),1:3)
    p5(1:3)=model_node_xyz(model_element_int_node(i,5),1:3)
    p6(1:3)=model_node_xyz(model_element_int_node(i,6),1:3)
    p7(1:3)=model_node_xyz(model_element_int_node(i,7),1:3)
    p8(1:3)=model_node_xyz(model_element_int_node(i,8),1:3)

    select case (model_element_type(i))
    case (1)
      ! Brick
      
      if (model_element_face_mask_2(i,1)) then
        model_face_type(i_face)  =4
        model_face_node(i_face,1)=model_element_int_node(i,1)
        model_face_node(i_face,2)=model_element_int_node(i,2)
        model_face_node(i_face,3)=model_element_int_node(i,3)
        model_face_node(i_face,4)=model_element_int_node(i,4)

        model_face_area(i_face)=area_quadrangle(p1,p2,p3,p4)

        model_node_face_mask(model_element_int_node(i,1))=.true.
        model_node_face_mask(model_element_int_node(i,2))=.true.
        model_node_face_mask(model_element_int_node(i,3))=.true.
        model_node_face_mask(model_element_int_node(i,4))=.true.

        i_face=i_face+1
      end if

      if (model_element_face_mask_2(i,2)) then
        model_face_type(i_face)  =4
        model_face_node(i_face,1)=model_element_int_node(i,5)
        model_face_node(i_face,2)=model_element_int_node(i,6)
        model_face_node(i_face,3)=model_element_int_node(i,7)
        model_face_node(i_face,4)=model_element_int_node(i,8)
        
        model_face_area(i_face)=area_quadrangle(p5,p6,p7,p8)
        
        model_node_face_mask(model_element_int_node(i,5))=.true.
        model_node_face_mask(model_element_int_node(i,6))=.true.
        model_node_face_mask(model_element_int_node(i,7))=.true.
        model_node_face_mask(model_element_int_node(i,8))=.true.
        
        i_face=i_face+1
      end if

      if (model_element_face_mask_2(i,3)) then
        model_face_type(i_face)  =4
        model_face_node(i_face,1)=model_element_int_node(i,1)
        model_face_node(i_face,2)=model_element_int_node(i,2)
        model_face_node(i_face,3)=model_element_int_node(i,6)
        model_face_node(i_face,4)=model_element_int_node(i,5)
        
        model_face_area(i_face)=area_quadrangle(p1,p2,p6,p5)
        
        model_node_face_mask(model_element_int_node(i,1))=.true.
        model_node_face_mask(model_element_int_node(i,2))=.true.
        model_node_face_mask(model_element_int_node(i,6))=.true.
        model_node_face_mask(model_element_int_node(i,5))=.true.
        
        i_face=i_face+1
      end if

      if (model_element_face_mask_2(i,4)) then
        model_face_type(i_face)  =4
        model_face_node(i_face,1)=model_element_int_node(i,2)
        model_face_node(i_face,2)=model_element_int_node(i,3)
        model_face_node(i_face,3)=model_element_int_node(i,7)
        model_face_node(i_face,4)=model_element_int_node(i,6)
        
        model_face_area(i_face)=area_quadrangle(p2,p3,p7,p6)
        
        model_node_face_mask(model_element_int_node(i,2))=.true.
        model_node_face_mask(model_element_int_node(i,3))=.true.
        model_node_face_mask(model_element_int_node(i,7))=.true.
        model_node_face_mask(model_element_int_node(i,6))=.true.
        
        i_face=i_face+1
      end if

      if (model_element_face_mask_2(i,5)) then
        model_face_type(i_face)  =4
        model_face_node(i_face,1)=model_element_int_node(i,3)
        model_face_node(i_face,2)=model_element_int_node(i,4)
        model_face_node(i_face,3)=model_element_int_node(i,8)
        model_face_node(i_face,4)=model_element_int_node(i,7)
        
        model_face_area(i_face)=area_quadrangle(p3,p4,p8,p7)
        
        model_node_face_mask(model_element_int_node(i,3))=.true.
        model_node_face_mask(model_element_int_node(i,4))=.true.
        model_node_face_mask(model_element_int_node(i,8))=.true.
        model_node_face_mask(model_element_int_node(i,7))=.true.
        
        i_face=i_face+1
      end if

      if (model_element_face_mask_2(i,6)) then
        model_face_type(i_face)  =4
        model_face_node(i_face,1)=model_element_int_node(i,1)
        model_face_node(i_face,2)=model_element_int_node(i,4)
        model_face_node(i_face,3)=model_element_int_node(i,8)
        model_face_node(i_face,4)=model_element_int_node(i,5)
        
        model_face_area(i_face)=area_quadrangle(p1,p4,p8,p5)
        
        model_node_face_mask(model_element_int_node(i,1))=.true.
        model_node_face_mask(model_element_int_node(i,4))=.true.
        model_node_face_mask(model_element_int_node(i,8))=.true.
        model_node_face_mask(model_element_int_node(i,5))=.true.
        
        i_face=i_face+1
      end if 

    case (2)
      ! Tetrahedral
      if (model_element_face_mask_2(i,1)) then
        model_face_type(i_face)  =3
        model_face_node(i_face,1)=model_element_int_node(i,1)
        model_face_node(i_face,2)=model_element_int_node(i,2)
        model_face_node(i_face,3)=model_element_int_node(i,3)
        
        model_face_area(i_face)=area_triangle(p1,p2,p3)

        model_node_face_mask(model_element_int_node(i,1))=.true.
        model_node_face_mask(model_element_int_node(i,2))=.true.
        model_node_face_mask(model_element_int_node(i,3))=.true.

        i_face=i_face+1
      end if

      if (model_element_face_mask_2(i,2)) then
        model_face_type(i_face)  =3
        model_face_node(i_face,1)=model_element_int_node(i,1)
        model_face_node(i_face,2)=model_element_int_node(i,2)
        model_face_node(i_face,3)=model_element_int_node(i,4)
        
        model_face_area(i_face)=area_triangle(p1,p2,p4)

        model_node_face_mask(model_element_int_node(i,1))=.true.
        model_node_face_mask(model_element_int_node(i,2))=.true.
        model_node_face_mask(model_element_int_node(i,4))=.true.

        i_face=i_face+1
      end if

      if (model_element_face_mask_2(i,3)) then
        model_face_type(i_face)  =3
        model_face_node(i_face,1)=model_element_int_node(i,2)
        model_face_node(i_face,2)=model_element_int_node(i,3)
        model_face_node(i_face,3)=model_element_int_node(i,4)
        
        model_face_area(i_face)=area_triangle(p2,p3,p4)

        model_node_face_mask(model_element_int_node(i,2))=.true.
        model_node_face_mask(model_element_int_node(i,3))=.true.
        model_node_face_mask(model_element_int_node(i,4))=.true.

        i_face=i_face+1
      end if

      if (model_element_face_mask_2(i,4)) then
        model_face_type(i_face)  =3
        model_face_node(i_face,1)=model_element_int_node(i,1)
        model_face_node(i_face,2)=model_element_int_node(i,3)
        model_face_node(i_face,3)=model_element_int_node(i,4)
        
        model_face_area(i_face)=area_triangle(p1,p3,p4)

        model_node_face_mask(model_element_int_node(i,1))=.true.
        model_node_face_mask(model_element_int_node(i,3))=.true.
        model_node_face_mask(model_element_int_node(i,4))=.true.

        i_face=i_face+1
      end if

    case (3)
      ! Triangular
      if (model_element_face_mask_2(i,1)) then
        model_face_type(i_face)  =3
        model_face_node(i_face,1)=model_element_int_node(i,1)
        model_face_node(i_face,2)=model_element_int_node(i,2)
        model_face_node(i_face,3)=model_element_int_node(i,3)
        
        model_face_area(i_face)=area_triangle(p1,p2,p3)

        model_node_face_mask(model_element_int_node(i,1))=.true.
        model_node_face_mask(model_element_int_node(i,2))=.true.
        model_node_face_mask(model_element_int_node(i,3))=.true.

        i_face=i_face+1
      end if

    case (4)
      ! Quadrilateral

      if (model_element_face_mask_2(i,1)) then
        model_face_type(i_face)  =4
        model_face_node(i_face,1)=model_element_int_node(i,1)
        model_face_node(i_face,2)=model_element_int_node(i,2)
        model_face_node(i_face,3)=model_element_int_node(i,3)
        model_face_node(i_face,4)=model_element_int_node(i,4)

        model_face_area(i_face)=area_quadrangle(p1,p2,p3,p4)

        model_node_face_mask(model_element_int_node(i,1))=.true.
        model_node_face_mask(model_element_int_node(i,2))=.true.
        model_node_face_mask(model_element_int_node(i,3))=.true.
        model_node_face_mask(model_element_int_node(i,4))=.true.

        i_face=i_face+1
      end if

    end select
  end do

end subroutine create_node_surface_map

subroutine filter_abaqus_data()
!-------------------
! Filter ABAQUS data
!-------------------
  use abaqus_arrays
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,k,ios,inode_abaqus
  integer :: n_model_nset
  character(len=chr80) :: stmp

  logical, dimension(:), allocatable :: model_nset_mask
  
  integer, dimension(:), allocatable :: node_map_1,node_map_2
  logical, dimension(:), allocatable :: node_map_mask

  if (len_trim(nset_input_file) == 0) then
    !----------------------------------------
    ! Patch-node set correspondence not given
    !----------------------------------------

    nnodes_abaqus=ubound(model_node_xyz,1)
    ntimes_abaqus=ntimes_fds

    allocate(abaqus_xyz(nnodes_abaqus,3),stat=ios);       call error_allocate(ios)
    allocate(abaqus_node_name(nnodes_abaqus),stat=ios);   call error_allocate(ios)
    allocate(abaqus_nset(nnodes_abaqus),stat=ios);        call error_allocate(ios)
    allocate(abaqus_node_number(nnodes_abaqus),stat=ios); call error_allocate(ios)

    abaqus_xyz        = model_node_xyz
    abaqus_node_name  = model_node_name
    abaqus_nset       = ''
    abaqus_node_number = 0

    return
  end if
  
  !------------------------------------
  ! Patch-node set correspondence given
  !------------------------------------

  ! Relevant arrays
  !--------------------------
  !   model_node_xyz(:),
  !   model_node_number(:),
  !   model_node_name(:),
  !   model_node_instance(:),
  !   model_nset_name(:),
  !   model_nset_node(:),
  !   model_nset_instance(:)
  !--------------------------

!!$  allocate(abaqus_face_type(ubound(model_face_type,1)),stat=ios);   call error_allocate(ios)
!!$  allocate(abaqus_face_node(ubound(model_face_node,1),4),stat=ios); call error_allocate(ios)
!!$  abaqus_face_type=model_face_type; abaqus_face_node=0

  ! Relevant arrays
  !--------------------------
  !   connectivity_table(:,:)
  !--------------------------
  
  n_model_nset=ubound(model_nset_name,1)
  if (.not. allocated(model_nset_mask)) then
    allocate(model_nset_mask(n_model_nset),stat=ios); call error_allocate(ios)
  end if
  model_nset_mask=.false.

  nnodes_abaqus=0
  count_loop: do i=1,ubound(model_nset_name,1)
    nset_loop_1: do k=1,ubound(connectivity_table,1)

      stmp=trim(lowercase(model_nset_instance(i))) // '.' // trim(lowercase(model_nset_name(i)))
      if (trim(stmp) == trim(lowercase(connectivity_table(k,1)))) then
        nnodes_abaqus=nnodes_abaqus+1
        model_nset_mask(i)=.true.
        exit nset_loop_1
      end if
    
    end do nset_loop_1
  end do count_loop

  if (nnodes_abaqus <= 0) then
    write(*,'(t3,a)') 'ERROR: no ABAQUS nodes selected'
    stop
  end if

  allocate(abaqus_xyz(nnodes_abaqus,3),stat=ios);         call error_allocate(ios)
  allocate(abaqus_node_name(nnodes_abaqus),stat=ios);     call error_allocate(ios)
  allocate(abaqus_node_number(nnodes_abaqus),stat=ios);   call error_allocate(ios)
  allocate(abaqus_nset(nnodes_abaqus),stat=ios);          call error_allocate(ios)
  allocate(abaqus_node_area(nnodes_abaqus),stat=ios);     call error_allocate(ios)

  allocate(node_map_1(nnodes_abaqus))
  allocate(node_map_2(ubound(model_node_xyz,1)))
  allocate(node_map_mask(nnodes_abaqus))
  node_map_1=0; node_map_2=0; node_map_mask=.true.

  abaqus_xyz       = 0.0
  abaqus_node_name = ''
  abaqus_nset      = ''
  abaqus_node_area = 0.0
  abaqus_node_number = 0

  if (trim(transfer_quantity) == 'adiabatic_surface_temperature' .or. &
    trim(transfer_quantity) == 'net heat flux') then

     inode_abaqus=1
     do i=1,ubound(model_nset_name,1)
        if (model_nset_mask(i) .eqv. .true.) then
           
           node_loop: do j=1,ubound(model_node_name,1)
              if (trim(model_node_instance(j)) == trim(model_nset_instance(i))) then
                 if (model_node_number(j) == model_nset_node(i)) then
                    
                    abaqus_xyz(inode_abaqus,1:3)       = model_node_xyz(j,1:3)
                    abaqus_node_name(inode_abaqus)     = model_node_name(j)
                    abaqus_nset(inode_abaqus)          = trim(model_nset_instance(i)) // '.' // trim(model_nset_name(i))
                    abaqus_node_area(inode_abaqus)     = model_node_area(j)
                    
                    node_map_1(inode_abaqus)=j
                    node_map_2(j)=inode_abaqus
                    
                    inode_abaqus=inode_abaqus+1
                 end if
              end if
           end do node_loop
           
        end if
     end do

  else

     inode_abaqus=1
     do i=1,ubound(model_nset_name,1)
        if (model_nset_mask(i) .eqv. .true.) then
           
           node_loop_2: do j=1,ubound(model_node_name,1)
              if (trim(model_node_instance(j)) == trim(model_nset_instance(i))) then
                 if (model_node_number(j) == model_nset_node(i)) then
                    
                    abaqus_xyz(inode_abaqus,1:3)       = model_node_xyz(j,1:3)
                    abaqus_node_name(inode_abaqus)     = model_node_name(j)
                    abaqus_nset(inode_abaqus)          = trim(model_nset_instance(i)) // '.' // trim(model_nset_name(i))
                    
                    inode_abaqus=inode_abaqus+1
                 end if
              end if
           end do node_loop_2
           
        end if
     end do

  end if

  !--------------------------
  ! Update element face table
  !--------------------------

  if (trim(transfer_quantity) == 'wall_temperature') then
    return
  end if

  allocate(abaqus_face_type(ubound(model_face_type,1)),stat=ios);   call error_allocate(ios)
  allocate(abaqus_face_node(ubound(model_face_node,1),4),stat=ios); call error_allocate(ios)
  abaqus_face_type=model_face_type; abaqus_face_node=0
  do i=1,ubound(abaqus_face_type,1)
    select case (abaqus_face_type(i))
    case (3)
      do j=1,3
        abaqus_face_node(i,j)=node_map_2(model_face_node(i,j))
      end do
    case (4)
      do j=1,4
       abaqus_face_node(i,j)=node_map_2(model_face_node(i,j))
      end do
    end select
  end do

end subroutine filter_abaqus_data

function get_keyword(input_line) result(keyword)
!--------------------------------------
! Parse ABAQUS input line for a keyword
!--------------------------------------
  use global_constants
  implicit none

  integer :: i,j
  character(len=input_line_length) :: keyword,input_line

  keyword=''
  if (input_line(1:1) /= '*') return

  if (input_line(2:2) == '*') then
    keyword='comment'
    return
  end if

  i=1
  do j=2,len_trim(input_line)
    if (input_line(j:j) /= ',') then
      keyword(i:i)=input_line(j:j)     
      i=i+1
    else
      exit
    end if
  end do

end function get_keyword

function get_parameter_string(input_line,paramstr) result(string)
!---------------------------------------------------------------
! Parse input line for a parameter with a character string value
!---------------------------------------------------------------
  use global_constants
  implicit none

  integer :: i,j,iparam
  character(len=chr80) :: string,paramstr
  character(len=input_line_length) :: input_line,stmp
  
  string=''; iparam=0; stmp=''; i=1
  do j=2,len_trim(input_line)
    if (input_line(j:j) /= ',') then
      stmp(i:i)=input_line(j:j)     
      i=i+1
    else
      if (len_trim(get_string(stmp,paramstr)) /= 0) string=get_string(stmp,paramstr)
      iparam=iparam+1
      stmp=''
      i=1
    end if
  end do

  if (len_trim(get_string(stmp,paramstr)) /= 0) string=get_string(stmp,paramstr)

contains

  function get_string(strin,paramstr) result(string)
    use global_constants
    implicit none

    logical :: found
    character(len=chr80) :: stmp,string,paramstr   
    character(len=input_line_length) :: strin

    stmp=trim(adjustl(strin))
    
!Timo:    found=compstr_initial(stmp(1:len_trim(paramstr)),paramstr)
    found=compstr_initial(stmp,paramstr)
    if (found .eqv. .true.) then
      stmp=adjustl(stmp(len_trim(paramstr)+1:len_trim(stmp)))
      if (stmp(1:1) == '=') then
        string=trim(adjustl(stmp(2:len_trim(stmp))))
      else
        string=''
      end if
    else
      string=''
    end if

  end function get_string

end function get_parameter_string

function get_parameter_logical(input_line,paramstr) result(answer)
!------------------------------------------------------
! Parse input line for a parameter with a logical value
!------------------------------------------------------
  use global_constants
  implicit none

  integer :: i,j,iparam
  logical :: answer
  character(len=chr80) :: paramstr
  character(len=input_line_length) :: input_line,stmp
  
  answer=.false.; iparam=0; stmp=''; i=1
  do j=2,len_trim(input_line)
    if (input_line(j:j) /= ',') then
      stmp(i:i)=input_line(j:j)     
      i=i+1
    else
      if (compstr(stmp,paramstr) .eqv. .true.) then
        answer=.true.
        return
      end if
      iparam=iparam+1
      stmp=''
      i=1
    end if
  end do

  if (answer .eqv. .false.) answer=compstr(stmp,paramstr)

end function get_parameter_logical

function get_parameter_integer(input_line,paramstr) result(n)
!---------------------------------------------------------------
! Parse input line for a parameter with a character string value
!---------------------------------------------------------------
  use global_constants
  implicit none

  integer :: i,j,n,iparam
  character(len=chr80) :: string,paramstr
  character(len=input_line_length) :: input_line,stmp
  
  string=''; iparam=0; stmp=''; i=1
  do j=2,len_trim(input_line)
    if (input_line(j:j) /= ',') then
      stmp(i:i)=input_line(j:j)     
      i=i+1
    else
      if (get_integer(stmp,paramstr) /= 0) n=get_integer(stmp,paramstr) 
      iparam=iparam+1
      stmp=''
      i=1
    end if
  end do

  if (get_integer(stmp,paramstr) == 0) n=get_integer(stmp,paramstr) 

contains

  function get_integer(strin,paramstr) result(n)
    use global_constants
    implicit none

    integer :: n,ios
    logical :: found
    character(len=chr80) :: stmp,string,paramstr   
    character(len=input_line_length) :: strin

    stmp=trim(adjustl(strin))
    
!Timo:    found=compstr_initial(stmp(1:len_trim(paramstr)),paramstr)
    found=compstr_initial(stmp,paramstr)
    if (found .eqv. .true.) then
      stmp=adjustl(stmp(len_trim(paramstr)+1:len_trim(stmp)))
      if (stmp(1:1) == '=') then
        read(stmp(2:len_trim(stmp)),*,iostat=ios) n
        if (ios /= 0) then
          n=0
        end if
      else
        n=0
      end if
    else
      n=0
    end if

  end function get_integer

end function get_parameter_integer

function get_parameter_real(input_line,paramstr) result(x)
!---------------------------------------------------------------
! Parse input line for a parameter with a character string value
!---------------------------------------------------------------
  use global_constants
  implicit none

  integer :: i,j,iparam
  real(kind=rk) :: x
  character(len=chr80) :: string,paramstr
  character(len=input_line_length) :: input_line,stmp
  
  string=''; iparam=0; stmp=''; i=1
  do j=2,len_trim(input_line)
    if (input_line(j:j) /= ',') then
      stmp(i:i)=input_line(j:j)     
      i=i+1
    else
      if (get_real(stmp,paramstr) /= 0.0) x=get_real(stmp,paramstr) 
      iparam=iparam+1
      stmp=''
      i=1
    end if
  end do

  if (get_real(stmp,paramstr) == 0.0) x=get_real(stmp,paramstr) 

contains

  function get_real(strin,paramstr) result(x)
    use global_constants
    implicit none

    integer :: ios
    logical :: found
    real(kind=rk) :: x
    character(len=chr80) :: stmp,string,paramstr   
    character(len=input_line_length) :: strin

    stmp=trim(adjustl(strin))
    
!Timo:    found=compstr_initial(stmp(1:len_trim(paramstr)),paramstr)
    found=compstr_initial(stmp,paramstr)
    if (found .eqv. .true.) then
      stmp=adjustl(stmp(len_trim(paramstr)+1:len_trim(stmp)))
      if (stmp(1:1) == '=') then
        read(stmp(2:len_trim(stmp)),*,iostat=ios) x
        if (ios /= 0) then
          x=0.0
        end if
      else
        x=0.0
      end if
    else
      x=0.0
    end if

  end function get_real

end function get_parameter_real

function compstr(str1,str2) result(match)
!-------------------------------------------------
! Compare two character strings (case insensitive)
!-------------------------------------------------
  use global_constants

  integer :: i
  logical :: match
  character(len=chr80) :: str1,str2

  str1=trim(adjustl(str1))
  str2=trim(adjustl(str2))

  if (len_trim(str1) /= len_trim(str2)) then
    match=.false.
    return
  end if

  match=.true.
  do i=1,len_trim(str2)
    if (iachar(str1(i:i)) == iachar(str2(i:i)) .or. &
        iachar(str1(i:i)) == (iachar(str2(i:i))+32) .or. &
        iachar(str1(i:i)) == (iachar(str2(i:i))-32)) then
      ! Do nothing
    else
      match=.false.
    end if
  end do

end function compstr

function compstr_initial(str1,str2) result(match)
!------------------------------------------------------------------
! Compare the beginning of two character strings (case insensitive)
!------------------------------------------------------------------
  use global_constants

  integer :: i
  logical :: match
  character(len=chr80) :: str1,str2

  str1=trim(adjustl(str1))
  str2=trim(adjustl(str2))

  match=.true.
  do i=1,len_trim(str2)
    if (iachar(str1(i:i)) == iachar(str2(i:i)) .or. &
        iachar(str1(i:i)) == (iachar(str2(i:i))+32) .or. &
        iachar(str1(i:i)) == (iachar(str2(i:i))-32)) then
      ! Do nothing
    else
      match=.false.
    end if
  end do

end function compstr_initial

function uppercase(str1) result(str2)
!----------------------------------------------
! Change lowercase letters to uppercase letters
!----------------------------------------------
  use global_constants

  integer :: i
  character(len=chr80) :: str1,str2

  str2=''
  do i=1,len_trim(str1)
    if (iachar(str1(i:i)) >= 97 .and. iachar(str1(i:i)) <= 122) then
      str2(i:i)=char(iachar(str1(i:i))-32)
    else
      str2(i:i)=str1(i:i)
    end if
  end do

end function uppercase

function number_integer_cols(input_line) result(ncols)
!-----------------------------------------------------------------------
! Find out the number of integer columns in a comma-separated input line
!-----------------------------------------------------------------------
  use error_messages
  use global_constants
  implicit none

  integer :: i,j,n,ios,ncols
  character(len=input_line_length) :: input_line,stmp
  
  stmp=''; i=1; ncols=0
  do j=1,len_trim(input_line)
    if (input_line(j:j) /= ',') then
      stmp(i:i)=input_line(j:j)     
      i=i+1
    else
      read(stmp(1:i-1),*,iostat=ios) n
      if (ios /= 0) then
        write(*,'(a)') 'ERROR: integer data expected' 
        stop
      end if
      ncols=ncols+1
      stmp=''
      i=1
    end if
  end do

  j=len_trim(input_line)
  if (input_line(j:j) /= ',') then
    read(stmp(1:i-1),*,iostat=ios) n
    if (ios /= 0) then
      write(*,'(a)') 'ERROR: integer data expected' 
      stop
    end if
    ncols=ncols+1
    stmp=''
    i=1
  end if

end function number_integer_cols

subroutine read_surface_nodes()
!--------------------------
! Read ABAQUS surface nodes
!--------------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  implicit none

  integer :: i,j,ios,nnode
  character(len=input_line_length) :: input_line

  open(unit=iochannel(1),file=trim(fem_input_file),status='old',iostat=ios)
  if (ios /= 0) call error_open_file(fem_input_file)

  !------------
  ! Count nodes
  !------------

  nnode=0
  count_loop: do
    
    read(iochannel(1),'(a)',iostat=ios) input_line
    select case(ios)
    case (:-1)
      exit count_loop
    case (0)
      nnode=nnode+1
    case (1:)
      call error_read_file(fem_input_file,(nnode+1))
    end select

  end do count_loop

  rewind(iochannel(1))

  !--------------
  ! Read in nodes
  !--------------

  allocate(abaqus_xyz(nnode,3),stat=ios);       call error_allocate(ios)
  allocate(abaqus_node_number(nnode),stat=ios); call error_allocate(ios)
  
  abaqus_xyz=0.0; abaqus_node_number=0

  read_loop: do i=1,nnode
    read(iochannel(1),'(a)',iostat=ios) input_line
    if (ios /= 0) call error_read_file(fem_input_file,(i+1)) 

    read(input_line,*,iostat=ios) abaqus_node_number(i), (abaqus_xyz(i,j),j=1,3)
    if (ios /= 0) then 
      read(input_line,*,iostat=ios) abaqus_node_number(i), (abaqus_xyz(i,j),j=1,2)
      if (ios /= 0) call error_read_file(fem_input_file,(i+1)) 
    end if

  end do read_loop

  close(unit=iochannel(1))

end subroutine read_surface_nodes

end module abaqus_reader
