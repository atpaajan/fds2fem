!-------------------------------------------
! Functions and subroutines for ANSYS output
!-------------------------------------------
module ansys_output

contains

subroutine ansys_output_module()
!---------
! Main sub
!---------
  use global_constants
  use global_variables
  implicit none

  write(*,'(a)') ''
  write(*,'(a)') 'ANSYS output module'

  if (fem_data_available) then
    select case(trim(transfer_quantity))
    case ('wall_temperature')
      call dump_nodal_temperature_load()
    case ('net_heat_flux')
      call dump_nodal_hflux_load()
    case ('adiabatic_surface_temperature')
      if (read_hcoeff) then
        call dump_nodal_ast_load()
      end if
    case default
      write(*,'(t3,a)') 'Nothing to be done'
    end select
  else
    write(*,'(t3,a)') 'Nothing to be done'
  end if

end subroutine ansys_output_module

!***************
! Auxiliary subs
!***************

subroutine dump_nodal_temperature_load()
!--------------------------------
! Output surface temperature load
!--------------------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,ios
  character(len=chr80) :: filename

  filename=trim(basename(fem_input_file)) // '_load.txt'
  open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  ntimes_ansys=ubound(fem_data,1)

  time_loop: do i=1,ntimes_ansys,10
    write(iochannel(1),'(2(a))') '! Load step ', trim(int2str(i))
    if (i == 1) then
      write(iochannel(1),'(a,es15.7e3)') 'TIME, ', fem_time(i)+1.0E-7
    else
      write(iochannel(1),'(a,es15.7e3)') 'TIME, ', fem_time(i)
    end if

    node_loop: do j=1,nnodes_fem
      write(iochannel(1),'(3(a),es15.7e3,a)') 'D, ', trim(int2str(fem_node_number(j))), &
        ', TEMP, ', fem_data(i,j), ', 0.0'
    end do node_loop
    
    write(iochannel(1),'(a)') 'SOLVE'

  end do time_loop
    
  write(iochannel(1),'(a)') 'FINISH'

  close(unit=iochannel(1))

  write(*,'(t3,5(a))') 'Temperature load history written in file ', &
    trim(quote(filename))

end subroutine dump_nodal_temperature_load

subroutine dump_nodal_hflux_load()
!--------------------------
! Output net heat flux load
!--------------------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,k,ios
  character(len=chr80) :: filename

  filename=trim(basename(fem_input_file)) // '_load.txt'
  open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  ntimes_ansys=ubound(fem_data,1)

  time_loop: do i=1,ntimes_ansys,4
    write(iochannel(1),'(2(a))') '! Load step ', trim(int2str(i))
    if (i == 1) then
      write(iochannel(1),'(a,es15.7e3)') 'TIME, ', fem_time(i)+1.0E-7
    else
      write(iochannel(1),'(a,es15.7e3)') 'TIME, ', fem_time(i)
    end if

    element_loop: do j=1,nelements_fem
      write(iochannel(1),'(3(a),3(es15.7e3,","),es15.7e3)') 'SFE,', trim(int2str(fem_element_number(j))), &
        ',1,HFLUX,,', (1000.0*fem_data(i,fem_node_order(fem_element_node(j,k))),k=1,4)
    end do element_loop
    
    write(iochannel(1),'(a)') 'SOLVE'

  end do time_loop
    
  write(iochannel(1),'(a)') 'FINISH'

  close(unit=iochannel(1))

  write(*,'(t3,5(a))') 'Net heat flux load history written in file ', &
    trim(quote(filename))

end subroutine dump_nodal_hflux_load

subroutine dump_nodal_ast_load()
!----------------
! Output AST load
!----------------
  use error_messages
  use global_constants
  use global_variables
  use mapping_arrays
  use string_handling
  implicit none

  integer :: i,j,k,ios
  character(len=chr80) :: filename

  filename=trim(basename(fem_input_file)) // '_load.txt'
  open(unit=iochannel(1),file=trim(filename),status='replace',iostat=ios)
  if (ios /= 0) call error_open_file(filename)

  ntimes_ansys=ubound(fem_data,1)

  time_loop: do i=1,ntimes_ansys,1
    write(iochannel(1),'(2(a))') '! Load step ', trim(int2str(i))
    if (i == 1) then
      write(iochannel(1),'(a,es15.7e3)') 'TIME, ', fem_time(i)+1.0E-7
    else
      write(iochannel(1),'(a,es15.7e3)') 'TIME, ', fem_time(i)
    end if

    element_loop_1: do j=1,nelements_fem
      write(iochannel(1),'(3(a),es15.7e3,a)') 'D, ', trim(int2str(fem_element_ast_node(j))), &
        ', TEMP, ', fem_data(i,fem_node_order(fem_element_ast_node(j))), ', 0.0'
    end do element_loop_1   

    element_loop: do j=1,nelements_fem
      write(iochannel(1),'(3(a),3(es15.7e3,","),es15.7e3)') 'SFE,', trim(int2str(fem_element_number(j))), &
        ',1,CONV,1,', (fem_hcoeff(i,fem_node_order(fem_element_node(j,k))),k=1,4)
    end do element_loop
    
    write(iochannel(1),'(a)') 'SOLVE'

  end do time_loop
    
  write(iochannel(1),'(a)') 'FINISH'

  close(unit=iochannel(1))

  write(*,'(t3,5(a))') 'AST load history written in file ', &
    trim(quote(filename))

end subroutine dump_nodal_ast_load

end module ansys_output
