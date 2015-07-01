!--------------------------------------------------------
! Functions and subroutines for reading ANSYS input files
!--------------------------------------------------------
module ansys_reader

contains

subroutine ansys_reader_module()
!---------
! Main sub
!---------
  use ansys_arrays
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

  ansys_xyz_available=.false.
  ansys_data_available=.false.
  ansys_model_available=.false.

  write(*,'(a)') ''
  write(*,'(a)') 'ANSYS reader module'

  if (len_trim(fem_input_file) == 0) then
    !---------------------
    ! No ANSYS input given
    !---------------------

    write(*,'(t3,a)') 'Nothing to be done' 

  else
    !--------------------------
    ! Parse ANSYS input file(s)
    !--------------------------
  
    if (trim(transfer_quantity) == 'adiabatic_surface_temperature') then
      call read_surface_nodes_ast()
    else
      call read_surface_nodes()
    end if
    fem_xyz_available=.true.

    write(*,'(t3,3(a))') 'ANSYS input file ', trim(quote(fem_input_file)), ' parsed'

  end if

end subroutine ansys_reader_module

!***************
! Auxiliary subs
!***************

subroutine read_surface_nodes()
!-------------------------
! Read ANSYS surface nodes
!-------------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,itmp,nnodes,nelements
  logical :: count_nodes,count_elements
  character(len=input_line_length) :: input_line

  open(unit=iochannel(1),file=trim(fem_input_file),status='old',iostat=ios)
  if (ios /= 0) call error_open_file(fem_input_file)

  !-------------------------
  ! Count nodes and elements
  !-------------------------

  count_nodes=.true.
  count_elements=.false.

  nnodes=0; nelements=0
  count_loop: do
    
    read(iochannel(1),'(a)',iostat=ios) input_line
    select case(ios)
    case (:-1)
      exit count_loop
    case (0)
      if (input_line(1:10) == '# Elements') then
        count_nodes=.false.
        count_elements=.true.
        cycle count_loop
      end if
      if (input_line(1:7) == '# Nodes') then
        count_nodes=.true.
        count_elements=.false.
        cycle count_loop
      end if
      if (count_nodes) then
        nnodes=nnodes+1
      else if (count_elements) then
        nelements=nelements+1
      end if
    case (1:)
      call error_read_file(fem_input_file,(nnodes+1))
    end select

  end do count_loop

  rewind(iochannel(1))

  !--------------
  ! Read in nodes
  !--------------

  nnodes_fem=nnodes

  if (nnodes_fem > 0) then

    allocate(fem_xyz(nnodes_fem,3),stat=ios);       call error_allocate(ios)
    allocate(fem_node_number(nnodes_fem),stat=ios); call error_allocate(ios)
    
    ! Omit comment line
    read(iochannel(1),'(a)',iostat=ios) input_line

    fem_xyz=0.0; fem_node_number=0

    node_loop_1: do i=1,nnodes_fem
      read(iochannel(1),'(a)',iostat=ios) input_line
      if (ios /= 0) call error_read_file(fem_input_file,(i+1)) 

      read(input_line,*,iostat=ios) fem_node_number(i), (fem_xyz(i,j),j=1,3)
      if (ios /= 0) then 
        read(input_line,*,iostat=ios) fem_node_number(i), (fem_xyz(i,j),j=1,2)
        if (ios /= 0) then
          read(input_line,*,iostat=ios) fem_node_number(i), (fem_xyz(i,j),j=1,1)
          if (ios /= 0) then
            read(input_line,*,iostat=ios) fem_node_number(i)
            if (ios /= 0) call error_read_file(fem_input_file,(i+1)) 
          end if
        end if
      end if
    
    end do node_loop_1

    ! Generate node order map
    allocate(fem_node_order(maxval(fem_node_number)),stat=ios); call error_allocate(ios)

    node_loop_2: do i=1,nnodes_fem
      fem_node_order(fem_node_number(i))=i
    end do node_loop_2

    if (nnodes_fem == 1) then
      write(*,'(t3,3(a))') 'Read in ', trim(int2str(nnodes_fem)), ' node'
    else
      write(*,'(t3,3(a))') 'Read in ', trim(int2str(nnodes_fem)), ' nodes'
    end if

  end if

  !-----------------
  ! Read in elements
  !-----------------

  nelements_fem=nelements

  if (nelements_fem > 0) then

    allocate(fem_element_node(nelements_fem,4),stat=ios); call error_allocate(ios)
    allocate(fem_element_number(nelements_fem),stat=ios); call error_allocate(ios)

    ! Omit comment line
    read(iochannel(1),'(a)',iostat=ios) input_line

    fem_element_node=0; fem_element_number=0

    element_loop: do i=1,nelements_fem
      read(iochannel(1),'(a)',iostat=ios) input_line
      if (ios /= 0) call error_read_file(fem_input_file,(nnodes_fem+i+2)) 

      read(input_line,*,iostat=ios) (fem_element_node(i,j),j=1,4), &
        (itmp,j=1,9), fem_element_number(i)
      if (ios /= 0) call error_read_file(fem_input_file,(nnodes_fem+i+2)) 
    end do element_loop
 
    if (nelements_fem == 1) then
      write(*,'(t3,3(a))') 'Read in ', trim(int2str(nelements_fem)), ' elements'
    else
      write(*,'(t3,3(a))') 'Read in ', trim(int2str(nelements_fem)), ' elements'
    end if

  end if

  close(unit=iochannel(1))

end subroutine read_surface_nodes

subroutine read_surface_nodes_ast()
!-------------------------------------
! Read ANSYS surface nodes in AST-mode
!-------------------------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios,itmp,nnodes,nelements
  logical :: count_nodes,count_elements
  character(len=input_line_length) :: input_line

  open(unit=iochannel(1),file=trim(fem_input_file),status='old',iostat=ios)
  if (ios /= 0) call error_open_file(fem_input_file)

  !-------------------------
  ! Count nodes and elements
  !-------------------------

  count_nodes=.true.
  count_elements=.false.

  nnodes=0; nelements=0
  count_loop: do
    
    read(iochannel(1),'(a)',iostat=ios) input_line
    select case(ios)
    case (:-1)
      exit count_loop
    case (0)
      if (input_line(1:10) == '# Elements') then
        count_nodes=.false.
        count_elements=.true.
        cycle count_loop
      end if
      if (input_line(1:7) == '# Nodes') then
        count_nodes=.true.
        count_elements=.false.
        cycle count_loop
      end if
      if (count_nodes) then
        nnodes=nnodes+1
      else if (count_elements) then
        nelements=nelements+1
      end if
    case (1:)
      call error_read_file(fem_input_file,(nnodes+1))
    end select

  end do count_loop

  rewind(iochannel(1))

  !--------------
  ! Read in nodes
  !--------------

  nnodes_fem=nnodes

  if (nnodes_fem > 0) then

    allocate(fem_xyz(nnodes_fem,3),stat=ios);       call error_allocate(ios)
    allocate(fem_node_number(nnodes_fem),stat=ios); call error_allocate(ios)
    
    ! Omit comment line
    read(iochannel(1),'(a)',iostat=ios) input_line

    fem_xyz=0.0; fem_node_number=0

    node_loop_1: do i=1,nnodes_fem
      read(iochannel(1),'(a)',iostat=ios) input_line
      if (ios /= 0) call error_read_file(fem_input_file,(i+1)) 

      read(input_line,*,iostat=ios) fem_node_number(i), (fem_xyz(i,j),j=1,3)
      if (ios /= 0) then 
        read(input_line,*,iostat=ios) fem_node_number(i), (fem_xyz(i,j),j=1,2)
        if (ios /= 0) then
          read(input_line,*,iostat=ios) fem_node_number(i), (fem_xyz(i,j),j=1,1)
          if (ios /= 0) then
            read(input_line,*,iostat=ios) fem_node_number(i)
            if (ios /= 0) call error_read_file(fem_input_file,(i+1)) 
          end if
        end if
      end if
    
    end do node_loop_1

    ! Generate node order map
    allocate(fem_node_order(maxval(fem_node_number)),stat=ios); call error_allocate(ios)

    node_loop_2: do i=1,nnodes_fem
      fem_node_order(fem_node_number(i))=i
    end do node_loop_2

    if (nnodes_fem == 1) then
      write(*,'(t3,3(a))') 'Read in ', trim(int2str(nnodes_fem)), ' node'
    else
      write(*,'(t3,3(a))') 'Read in ', trim(int2str(nnodes_fem)), ' nodes'
    end if

  end if

  !-----------------
  ! Read in elements
  !-----------------

  nelements_fem=nelements

  if (nelements_fem > 0) then

    allocate(fem_element_node(nelements_fem,4),stat=ios); call error_allocate(ios)
    allocate(fem_element_ast_node(nelements_fem),stat=ios); call error_allocate(ios)
    allocate(fem_element_number(nelements_fem),stat=ios); call error_allocate(ios)

    ! Omit comment line
    read(iochannel(1),'(a)',iostat=ios) input_line

    fem_element_node=0; fem_element_number=0

    element_loop: do i=1,nelements_fem
      read(iochannel(1),'(a)',iostat=ios) input_line
      if (ios /= 0) call error_read_file(fem_input_file,(nnodes_fem+i+2)) 

      read(input_line,'(14(i6))',iostat=ios) (fem_element_node(i,j),j=1,4), &
        fem_element_ast_node(i), (itmp,j=1,8), fem_element_number(i)
      if (ios /= 0) call error_read_file(fem_input_file,(nnodes_fem+i+2)) 
    end do element_loop
 
    if (nelements_fem == 1) then
      write(*,'(t3,3(a))') 'Read in ', trim(int2str(nelements_fem)), ' elements'
    else
      write(*,'(t3,3(a))') 'Read in ', trim(int2str(nelements_fem)), ' elements'
    end if

  end if

  close(unit=iochannel(1))

end subroutine read_surface_nodes_ast

end module ansys_reader
