!----------------------------------
! Global constants, variables, etc.
!----------------------------------
module global_constants
!-----------------
! Global constants
!-----------------
  implicit none

  ! I/O variables
  integer :: scratch_channel
  integer, dimension(10) :: iochannel

  ! Floating point precision
  integer, parameter :: rk = selected_real_kind(10,40)
  integer, parameter :: fb = selected_real_kind(6)
  integer, parameter :: eb = selected_real_kind(12)

  ! General purpose
  integer, parameter :: input_line_length = 5096

  ! String length
  integer, parameter :: chr20  = 20
  integer, parameter :: chr25  = 25
  integer, parameter :: chr30  = 30
  integer, parameter :: chr40  = 40
  integer, parameter :: chr60  = 60
  integer, parameter :: chr80  = 80
  integer, parameter :: chr82  = 82
  integer, parameter :: chr100 = 100
  integer, parameter :: chr256 = 256

  ! Other
  real(kind=rk), parameter :: eps=1.0e-10

end module global_constants

module global_variables
!-----------------
! Global variables
!-----------------
  use global_constants
  implicit none

  ! Command-line flags
  logical :: full_output
 
  ! Other variables 
  character(len=chr80) :: config_file

  character(len=chr80) :: fds_input_file
  character(len=chr80) :: fds_output

  character(len=chr80) :: fem_input_file
  character(len=chr80) :: fem_mode
  
  character(len=chr80) :: nset_input_file
  character(len=chr80) :: transfer_quantity
  
  character(len=chr80) :: dump_fds_nodes
  character(len=chr80) :: dump_fds_data
  character(len=chr80) :: dump_fds_model

  character(len=chr80) :: dump_cfast_nodes
  character(len=chr80) :: dump_cfast_data

  character(len=chr80) :: dump_fem_nodes
  character(len=chr80) :: dump_fem_data
  character(len=chr80) :: dump_fem_model

  character(len=chr80) :: mapping_method

  logical :: nset_connectivity
  logical :: read_all_fds_data

  logical :: iso_curve, cfast_input

  ! Mapping
  integer       :: mp_n
  integer       :: mp_nmin
  integer       :: mp_nmx
  real(kind=rk) :: mp_cut
  real(kind=rk) :: mp_del
  real(kind=rk) :: mp_deg

  ! Model matching
  logical :: match_translate
  logical :: match_rotate
  logical :: manual_translate
  logical :: manual_rotate
  logical :: automatic_translate
  logical :: automatic_rotate

  real(kind=rk), dimension(3) :: origin_fds
  real(kind=rk), dimension(3) :: origin_fem

  real(kind=rk) :: e_alpha
  real(kind=rk) :: e_beta
  real(kind=rk) :: e_gamma

  ! Other
  logical :: fds_cm_calculated
  logical :: fds_range_calculated
  logical :: fem_cm_calculated
  logical :: fem_range_calculated

  logical :: fds_xyz_available
  logical :: fds_data_available
  logical :: fds_model_available

  logical :: cfast_xyz_available
  logical :: cfast_data_available
  
  logical :: fem_xyz_available
  logical :: fem_data_available
  logical :: fem_model_available

  logical :: ansys_xyz_available
  logical :: ansys_data_available
  logical :: ansys_model_available

  logical :: abaqus_xyz_available
  logical :: abaqus_data_available
  logical :: abaqus_model_available
  
  logical :: fds_statistics
  logical :: fem_statistics

  logical :: read_hcoeff
  logical :: ansys_ast

  real(kind=rk) :: hcoeff
  real(kind=rk) :: emissivity

  integer :: iso_ntimes
  real(kind=rk) :: iso_tbegin, iso_tend

end module global_variables

module fds_head_arrays
!------------------------
! FDS HEAD-related arrays
!------------------------
  use global_constants
  implicit none

  character(len=chr40) :: fds_chid
  character(len=chr40) :: fds_title
  character(len=chr40) :: fds_fyi
  
end module fds_head_arrays

module fds_mesh_arrays
!------------------------
! FDS MESH-related arrays
!------------------------
  use global_constants

  integer, dimension(:,:), allocatable :: fds_mesh_ijk
  real(kind=rk), dimension(:,:), allocatable :: fds_mesh_xb

end module fds_mesh_arrays

module fds_obst_arrays
!------------------------
! FDS OBST-related arrays
!------------------------
  use global_constants
  implicit none
  
  integer, dimension(:,:),       allocatable :: fds_obst_el
  real(kind=fb), dimension(:,:), allocatable :: fds_obst_xb
  real(kind=fb), dimension(:,:), allocatable :: fds_obst_nd

end module fds_obst_arrays

module fds_prop_arrays
!------------------------
! FDS PROP-related arrays
!------------------------
  use global_constants
  implicit none

  character(len=chr30), dimension(:), allocatable :: fds_prop_name
  character(len=chr30), dimension(:), allocatable :: fds_prop_qnty

end module fds_prop_arrays

module fds_dump_arrays
!------------------------
! FDS DUMP-related arrays
!------------------------
  implicit none

  logical :: fds_column_dump_limit
  integer :: fds_devc_column_limit

end module fds_dump_arrays

module fds_devc_arrays
!------------------------
! FDS DEVC-related arrays
!------------------------
  use global_constants
  implicit none

  integer, dimension(:), allocatable :: fds_devc_ior
  integer, dimension(:), allocatable :: fds_devc_rows
  integer, dimension(:), allocatable :: fds_devc_cols

  character(len=chr80), dimension(:), allocatable :: fds_devc_file
  character(len=chr20), dimension(:), allocatable :: fds_devc_name
  character(len=chr20), dimension(:), allocatable :: fds_devc_name_b
  character(len=chr30), dimension(:), allocatable :: fds_devc_qnty
  character(len=chr30), dimension(:), allocatable :: fds_devc_unit

  real(kind=rk), dimension(:),   allocatable :: fds_devc_time
  real(kind=rk), dimension(:,:), allocatable :: fds_devc_xyz
  real(kind=rk), dimension(:,:), allocatable :: fds_devc_xb
  real(kind=rk), dimension(:,:), allocatable :: fds_devc_data
  
end module fds_devc_arrays

module fds_bndf_arrays
!------------------------
! FDS BNDF-related arrays
!-----------------------
  use global_constants
  implicit none
  
  !-------------------------
  ! Chosen transfer quantity
  !-------------------------

  integer, dimension(:), allocatable :: fds_bndf_ior
  integer, dimension(:), allocatable :: fds_bndf_nnod
  integer, dimension(:), allocatable :: fds_bndf_npat
  integer, dimension(:), allocatable :: fds_bndf_ntim
  integer, dimension(:), allocatable :: fds_bndf_pat

  integer, dimension(:),   allocatable :: fds_bndf_nb
  integer, dimension(:),   allocatable :: fds_bndf_nm
  integer, dimension(:),   allocatable :: fds_bndf_np
  integer, dimension(:,:), allocatable :: fds_bndf_el

  character(len=chr30), dimension(:),   allocatable :: fds_bndf_qnty
  character(len=chr30), dimension(:),   allocatable :: fds_bndf_snam
  character(len=chr30), dimension(:),   allocatable :: fds_bndf_unit
  character(len=chr80), dimension(:,:), allocatable :: fds_bndf_file

  real(kind=fb), dimension(:),   allocatable :: fds_bndf_time
  real(kind=fb), dimension(:,:), allocatable :: fds_bndf_xyz
  real(kind=fb), dimension(:,:), allocatable :: fds_bndf_data

  !--------------------------
  ! Heat transfer coefficient
  !--------------------------

  integer, dimension(:), allocatable :: fds_hcoeff_bndf_ior
  integer, dimension(:), allocatable :: fds_hcoeff_bndf_nnod
  integer, dimension(:), allocatable :: fds_hcoeff_bndf_npat
  integer, dimension(:), allocatable :: fds_hcoeff_bndf_ntim
  integer, dimension(:), allocatable :: fds_hcoeff_bndf_pat

  integer, dimension(:),   allocatable :: fds_hcoeff_bndf_nb
  integer, dimension(:),   allocatable :: fds_hcoeff_bndf_nm
  integer, dimension(:),   allocatable :: fds_hcoeff_bndf_np
  integer, dimension(:,:), allocatable :: fds_hcoeff_bndf_el

  character(len=chr30), dimension(:),   allocatable :: fds_hcoeff_bndf_qnty
  character(len=chr30), dimension(:),   allocatable :: fds_hcoeff_bndf_snam
  character(len=chr30), dimension(:),   allocatable :: fds_hcoeff_bndf_unit
  character(len=chr80), dimension(:,:), allocatable :: fds_hcoeff_bndf_file

  real(kind=fb), dimension(:),   allocatable :: fds_hcoeff_bndf_time
  real(kind=fb), dimension(:,:), allocatable :: fds_hcoeff_bndf_xyz
  real(kind=fb), dimension(:,:), allocatable :: fds_hcoeff_bndf_data

end module fds_bndf_arrays

module mapping_arrays
!-------------------------------
! Arrays related to mesh mapping
!-------------------------------
  use global_constants
  implicit none

  integer :: nnodes_fds
  integer :: ntimes_fds
  integer :: nnodes_fem
  integer :: ntimes_fem
  integer :: nnodes_ansys
  integer :: ntimes_ansys
  integer :: nnodes_abaqus
  integer :: ntimes_abaqus
  integer :: nelements_fem

  integer, dimension(:), allocatable :: fds_idevc
  integer, dimension(:), allocatable :: fds_patch
  integer, dimension(:), allocatable :: fds_selected_devc
  integer, dimension(:), allocatable :: fds_selected_patch
  
  integer, dimension(:), allocatable :: fds_submap
  integer, dimension(:), allocatable :: abaqus_submap
  
  integer, dimension(:),   allocatable :: ansys_node_number

  integer, dimension(:),   allocatable :: fem_node_number
  integer, dimension(:),   allocatable :: fem_node_order

  integer, dimension(:),   allocatable :: abaqus_node_number
  integer, dimension(:),   allocatable :: abaqus_node_elements
  integer, dimension(:),   allocatable :: abaqus_face_type
  integer, dimension(:,:), allocatable :: abaqus_face_node

  integer, dimension(:),   allocatable :: ansys_element_number
  integer, dimension(:,:), allocatable :: ansys_element_node

  integer, dimension(:,:), allocatable :: connectivity_table_num
  real(kind=rk), dimension(:,:), allocatable :: connectivity_phys_cons
  logical, dimension(:), allocatable :: time_shift_mask
  
  integer, dimension(:,:), allocatable :: abaqus_element_node

  logical, dimension(:), allocatable :: fds_mask
  logical, dimension(:), allocatable :: fem_mask
  
  character(len=chr80), dimension(:), allocatable :: fem_nset

  character(len=chr80), dimension(:), allocatable :: fds_id
  character(len=chr80), dimension(:), allocatable :: abaqus_nset
  character(len=chr80), dimension(:), allocatable :: abaqus_selected_nset
  character(len=chr80), dimension(:), allocatable :: abaqus_node_name
  
  character(len=chr80), dimension(:), allocatable :: fem_node_name

  character(len=chr80), dimension(:,:), allocatable :: connectivity_table

  real(kind=rk), dimension(6) :: fds_node_range
  real(kind=rk), dimension(6) :: fem_node_range
  real(kind=rk), dimension(6) :: abaqus_node_range
  real(kind=rk), dimension(3) :: node_ratio

  real(kind=rk), dimension(3) :: cm_fds,cm_fds_bb,fds_dr
  real(kind=rk), dimension(3) :: cm_fem,cm_fem_bb,fem_dr
  
  integer,       dimension(:),   allocatable :: fds_ior
  real(kind=rk), dimension(:),   allocatable :: fds_time
  real(kind=rk), dimension(:,:), allocatable :: fds_xyz
  real(kind=rk), dimension(:,:), allocatable :: fds_xyz_sub
  real(kind=rk), dimension(:,:), allocatable :: fds_data
  real(kind=rk), dimension(:,:), allocatable :: fds_data_sub
  real(kind=rk), dimension(:,:), allocatable :: fds_hcoeff
  real(kind=rk), dimension(:,:), allocatable :: fds_hcoeff_sub
  
  real(kind=rk), dimension(:), allocatable :: abaqus_time
  real(kind=rk), dimension(:), allocatable :: abaqus_node_area
  real(kind=rk), dimension(:), allocatable :: abaqus_node_hcoeff
  real(kind=rk), dimension(:), allocatable :: abaqus_node_emissivity

  real(kind=rk), dimension(:,:), allocatable :: ansys_xyz

  real(kind=rk), dimension(:,:), allocatable :: abaqus_xyz
  real(kind=rk), dimension(:,:), allocatable :: abaqus_data
  
  real(kind=rk), dimension(:), allocatable :: fem_time

  real(kind=rk), dimension(:,:), allocatable :: fem_xyz
  real(kind=rk), dimension(:,:), allocatable :: fem_xyz_sub
  real(kind=rk), dimension(:,:), allocatable :: fem_data
  real(kind=rk), dimension(:,:), allocatable :: fem_data_sub
  real(kind=rk), dimension(:,:), allocatable :: fem_hcoeff
  real(kind=rk), dimension(:,:), allocatable :: fem_hcoeff_sub

  integer, dimension(:),   allocatable :: fem_element_number
  integer, dimension(:),   allocatable :: fem_element_ast_node
  integer, dimension(:,:), allocatable :: fem_element_node

end module mapping_arrays

module ansys_arrays
!---------------------
! ANSYS-related arrays
!---------------------
  use global_constants
  implicit none

end module ansys_arrays

module abaqus_arrays
!----------------------
! ABAQUS-related arrays
!----------------------
  use global_constants
  implicit none

  logical, dimension(:),   allocatable :: model_node_face_mask

  integer, dimension(:),   allocatable :: abaqus_inode

  integer, dimension(:),   allocatable :: part_node_number
  integer, dimension(:),   allocatable :: part_nset_node
  integer, dimension(:),   allocatable :: part_element_number
  integer, dimension(:),   allocatable :: part_element_type
  integer, dimension(:,:), allocatable :: part_element_node

  integer, dimension(:),   allocatable :: instance_node_number
  integer, dimension(:),   allocatable :: instance_nset_node
  integer, dimension(:),   allocatable :: instance_element_number
  integer, dimension(:),   allocatable :: instance_element_type
  integer, dimension(:,:), allocatable :: instance_element_node

  integer, dimension(:),   allocatable :: assembly_nset_node

  integer, dimension(:),   allocatable :: model_node_number
  integer, dimension(:),   allocatable :: model_node_type
  integer, dimension(:),   allocatable :: model_node_elements
  integer, dimension(:),   allocatable :: model_nset_node
  integer, dimension(:),   allocatable :: model_element_number
  integer, dimension(:),   allocatable :: model_element_type
  integer, dimension(:,:), allocatable :: model_element_node
  integer, dimension(:,:), allocatable :: model_element_int_node
  integer, dimension(:),   allocatable :: model_face_type
  integer, dimension(:,:), allocatable :: model_face_node
  
  character(len=chr80), dimension(:), allocatable :: part_node_name
  character(len=chr80), dimension(:), allocatable :: part_node_part
  character(len=chr80), dimension(:), allocatable :: part_nset_name
  character(len=chr80), dimension(:), allocatable :: part_nset_part
  character(len=chr80), dimension(:), allocatable :: part_element_part

  character(len=chr80), dimension(:), allocatable :: instance_node_name
  character(len=chr80), dimension(:), allocatable :: instance_node_part
  character(len=chr80), dimension(:), allocatable :: instance_node_instance
  character(len=chr80), dimension(:), allocatable :: instance_nset_name
  character(len=chr80), dimension(:), allocatable :: instance_nset_part
  character(len=chr80), dimension(:), allocatable :: instance_nset_instance
  character(len=chr80), dimension(:), allocatable :: instance_element_part
  character(len=chr80), dimension(:), allocatable :: instance_element_instance
  character(len=chr80), dimension(:), allocatable :: assembly_nset_name
  character(len=chr80), dimension(:), allocatable :: assembly_nset_part
  character(len=chr80), dimension(:), allocatable :: assembly_nset_instance
  
  character(len=chr80), dimension(:), allocatable :: model_node_name
  character(len=chr80), dimension(:), allocatable :: model_node_part
  character(len=chr80), dimension(:), allocatable :: model_node_instance
  
  character(len=chr80), dimension(:), allocatable :: model_nset_name
  character(len=chr80), dimension(:), allocatable :: model_nset_part
  character(len=chr80), dimension(:), allocatable :: model_nset_instance

  character(len=chr80), dimension(:), allocatable :: model_element_part
  character(len=chr80), dimension(:), allocatable :: model_element_instance

  real(kind=rk), dimension(:,:), allocatable :: part_node_xyz
  real(kind=rk), dimension(:,:), allocatable :: instance_node_xyz
  real(kind=rk), dimension(:,:), allocatable :: model_node_xyz
  real(kind=rk), dimension(:),   allocatable :: model_node_area
  real(kind=rk), dimension(:),   allocatable :: model_face_area
  real(kind=rk), dimension(:,:), allocatable :: abaqus_node

  !-------------------------------------
  ! Keyword-related variables and arrays
  !-------------------------------------

  integer :: n_nset,n_ngen,n_ncopy,n_node,n_part,n_end_part,n_nfill,n_nmap,n_system
  integer :: n_assembly,n_end_assembly,n_include,n_instance,n_end_instance,n_element
  integer :: n_part_nodes,n_part_nset_nodes,n_instance_nset_nodes,n_part_elements
  integer :: n_instance_nodes,n_instance_elements,n_assembly_nset_nodes

  ! ASSEMBLY-keyword
  character(len=chr80), dimension(:), allocatable :: assembly_name

  ! NSET-keyword
  logical, dimension(:), allocatable :: nset_generate
  logical, dimension(:), allocatable :: nset_internal
  logical, dimension(:), allocatable :: nset_unsorted

  character(len=chr80), dimension(:), allocatable :: nset_nset
  character(len=chr80), dimension(:), allocatable :: nset_elset
  character(len=chr80), dimension(:), allocatable :: nset_instance

  character(len=chr80), dimension(:), allocatable :: nset_part
  character(len=chr80), dimension(:), allocatable :: nset_inst
  character(len=chr80), dimension(:), allocatable :: nset_assm

  ! NODE-keyword
  character(len=chr80), dimension(:), allocatable :: node_input
  character(len=chr80), dimension(:), allocatable :: node_nset
  character(len=chr80), dimension(:), allocatable :: node_system

  character(len=chr80), dimension(:), allocatable :: node_part
  character(len=chr80), dimension(:), allocatable :: node_instance
  character(len=chr80), dimension(:), allocatable :: node_assembly

  ! PART-keyword
  character(len=chr80), dimension(:), allocatable :: part_name

  ! NGEN-keyword
  character(len=chr80), dimension(:), allocatable :: ngen_line
  character(len=chr80), dimension(:), allocatable :: ngen_nset
  character(len=chr80), dimension(:), allocatable :: ngen_system

  ! NCOPY-keyword
  integer, dimension(:), allocatable :: ncopy_change_number
  integer, dimension(:), allocatable :: ncopy_multiple

  logical, dimension(:), allocatable :: ncopy_pole
  logical, dimension(:), allocatable :: ncopy_shift

  character(len=chr80), dimension(:), allocatable :: ncopy_old_set
  character(len=chr80), dimension(:), allocatable :: ncopy_reflect
  character(len=chr80), dimension(:), allocatable :: ncopy_new_set

  ! NFILL-keyword
  integer, dimension(:), allocatable :: nfill_singular
  logical, dimension(:), allocatable :: nfill_two_step
  character(len=chr80), dimension(:), allocatable :: nfill_nset
  real(kind=rk), dimension(:), allocatable :: nfill_bias

  ! NMAP-keyword
  character(len=chr80), dimension(:), allocatable :: nmap_nset
  character(len=chr80), dimension(:), allocatable :: nmap_type
  character(len=chr80), dimension(:), allocatable :: nmap_definition

  ! INCLUDE-keyword
  integer, dimension(:), allocatable :: include_line_number

  character(len=chr80), dimension(:), allocatable :: include_input
  character(len=chr80), dimension(:), allocatable :: include_password

  ! INSTANCE-keyword
  character(len=chr80), dimension(:), allocatable :: instance_name
  character(len=chr80), dimension(:), allocatable :: instance_part
  character(len=chr80), dimension(:), allocatable :: instance_instance
  character(len=chr80), dimension(:), allocatable :: instance_library

  real(kind=rk), dimension(:,:), allocatable :: instance_translate
  real(kind=rk), dimension(:,:), allocatable :: instance_rotate

  ! ELEMENT-keyword
  character(len=chr80), dimension(:), allocatable :: element_type
  character(len=chr80), dimension(:), allocatable :: element_elset
  character(len=chr80), dimension(:), allocatable :: element_file
  character(len=chr80), dimension(:), allocatable :: element_input
  character(len=chr80), dimension(:), allocatable :: element_sen

  real(kind=rk), dimension(:), allocatable :: element_offset

  character(len=chr80), dimension(:), allocatable :: element_part
  character(len=chr80), dimension(:), allocatable :: element_instance
  character(len=chr80), dimension(:), allocatable :: element_assembly

  !------
  ! Files
  !------

  character(len=chr80), dimension(:), allocatable :: part_files
  character(len=chr80), dimension(:), allocatable :: instance_files
  character(len=chr80), dimension(:), allocatable :: assembly_files

  !------
  ! Other
  !------
  
  integer, dimension(:), allocatable :: node_nodes
  integer, dimension(:), allocatable :: nset_nodes
  integer, dimension(:), allocatable :: element_elements

  integer, dimension(:,:), allocatable :: node_lines
  integer, dimension(:,:), allocatable :: nset_lines
  integer, dimension(:,:), allocatable :: element_lines

  integer, dimension(:,:), allocatable :: part_lines
  integer, dimension(:,:), allocatable :: instance_lines
  integer, dimension(:,:), allocatable :: assembly_lines
  
end module abaqus_arrays

Module cfast_arrays
!---------------------
! CFAST-related arrays
!---------------------
  Use global_constants
  Implicit none

  Real(kind=rk), Dimension(:,:), Allocatable :: cfast_target_xyz, cfast_comp_xyz, cfast_target_ior
  Real(kind=rk), Dimension(:), Allocatable :: cfast_matl_epsilon, cfast_target_epsilon
  Character(len=chr20), Dimension(:), Allocatable :: cfast_target_name, cfast_target_matl, cfast_matl_name
  Integer, Dimension(:), Allocatable :: cfast_target_comp

End Module cfast_arrays
