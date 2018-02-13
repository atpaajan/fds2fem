!--------------------------------------------------------------
! Functions and subroutines for various mathematical operations
!--------------------------------------------------------------
module mathematics

contains

function euler(e_alpha,e_beta,e_gamma) result(rm)
!----------------------------------------------------------------
! Create a rotation matrix for a rotation defined by Euler angles
!----------------------------------------------------------------
  use global_constants
  implicit none

  real(kind=rk) :: e_alpha,e_beta,e_gamma
  real(kind=rk), dimension(3,3) :: rm

  ! Temporary fix...
  e_alpha=-e_alpha; e_beta=-e_beta; e_gamma=-e_gamma

  rm(1,1)= cos(e_gamma)*cos(e_alpha)-cos(e_beta)*sin(e_alpha)*sin(e_gamma)
  rm(1,2)= cos(e_gamma)*sin(e_alpha)+cos(e_beta)*cos(e_alpha)*sin(e_gamma)
  rm(1,3)= sin(e_gamma)*sin(e_beta)
  rm(2,1)=-sin(e_gamma)*cos(e_alpha)-cos(e_beta)*sin(e_alpha)*cos(e_gamma)
  rm(2,2)=-sin(e_gamma)*sin(e_alpha)+cos(e_beta)*cos(e_alpha)*cos(e_gamma)
  rm(2,3)= cos(e_gamma)*sin(e_beta)
  rm(3,1)= sin(e_beta)*sin(e_alpha)
  rm(3,2)=-sin(e_beta)*cos(e_alpha)
  rm(3,3)= cos(e_beta)

end function euler

function rotate_arbitrary_axis(p1,p2,theta) result(ra)
!-------------------------------------------
! Create a rotation matrix for a rotation
! around an arbitrary axis defined by points
! P1 and P2, and angle theta
!-------------------------------------------
  use global_constants
  implicit none

  real(kind=rk) :: theta,rtmp,st,ct
  real(kind=rk), dimension(3) :: p1,p2,r

  real(kind=rk), dimension(4,4) :: t,t_inv
  real(kind=rk), dimension(4,4) :: rx,rx_inv
  real(kind=rk), dimension(4,4) :: ry,ry_inv
  real(kind=rk), dimension(4,4) :: rz
  
  real(kind=rk), dimension(4,4) :: ra

  t(1,1)=1.0; t(1,2)=0.0; t(1,3)=0.0; t(1,4)=-p1(1)
  t(2,1)=0.0; t(2,2)=1.0; t(2,3)=0.0; t(2,4)=-p1(2)
  t(3,1)=0.0; t(3,2)=0.0; t(3,3)=1.0; t(3,4)=-p1(3)
  t(4,1)=0.0; t(4,2)=0.0; t(4,3)=0.0; t(4,4)=1.0

  t_inv(1,1)=1.0; t_inv(1,2)=0.0; t_inv(1,3)=0.0; t_inv(1,4)=p1(1)
  t_inv(2,1)=0.0; t_inv(2,2)=1.0; t_inv(2,3)=0.0; t_inv(2,4)=p1(2)
  t_inv(3,1)=0.0; t_inv(3,2)=0.0; t_inv(3,3)=1.0; t_inv(3,4)=p1(3)
  t_inv(4,1)=0.0; t_inv(4,2)=0.0; t_inv(4,3)=0.0; t_inv(4,4)=1.0
  
  r(1)=p2(1)-p1(1); r(2)=p2(2)-p1(2); r(3)=p2(3)-p1(3) 
  rtmp=sqrt(r(2)*r(2)+r(3)*r(3))
  If (sqrt(r(1)**2+r(2)**2+r(3)**2) > eps) Then
     r = r/(sqrt(r(1)**2+r(2)**2+r(3)**2))
  Else
     r = 0.0
  End If
  !Timo:!$  r=r/(sqrt(r(1)**2+r(2)**2+r(3)**2))
!Timo:!$  rtmp=sqrt(r(2)*r(2)+r(3)*r(3))
 
  if (rtmp > eps) then
    rx(1,1)=1.0; rx(1,2)=0.0;       rx(1,3)=0.0;        rx(1,4)=0.0
    rx(2,1)=0.0; rx(2,2)=r(3)/rtmp; rx(2,3)=-r(2)/rtmp; rx(2,4)=0.0
    rx(3,1)=0.0; rx(3,2)=r(2)/rtmp; rx(3,3)=r(3)/rtmp;  rx(3,4)=0.0
    rx(4,1)=0.0; rx(4,2)=0.0;       rx(4,3)=0.0;        rx(4,4)=1.0
  else
    rx(1,1)=1.0; rx(1,2)=0.0; rx(1,3)=0.0; rx(1,4)=0.0
    rx(2,1)=0.0; rx(2,2)=1.0; rx(2,3)=0.0; rx(2,4)=0.0
    rx(3,1)=0.0; rx(3,2)=0.0; rx(3,3)=1.0; rx(3,4)=0.0
    rx(4,1)=0.0; rx(4,2)=0.0; rx(4,3)=0.0; rx(4,4)=1.0 
  end if

  if (rtmp > eps) then
    rx_inv(1,1)=1.0; rx_inv(1,2)=0.0;        rx_inv(1,3)=0.0;       rx_inv(1,4)=0.0
    rx_inv(2,1)=0.0; rx_inv(2,2)=r(3)/rtmp;  rx_inv(2,3)=r(2)/rtmp; rx_inv(2,4)=0.0
    rx_inv(3,1)=0.0; rx_inv(3,2)=-r(2)/rtmp; rx_inv(3,3)=r(3)/rtmp; rx_inv(3,4)=0.0
    rx_inv(4,1)=0.0; rx_inv(4,2)=0.0;        rx_inv(4,3)=0.0;       rx_inv(4,4)=1.0
  else
    rx_inv(1,1)=1.0; rx_inv(1,2)=0.0; rx_inv(1,3)=0.0; rx_inv(1,4)=0.0
    rx_inv(2,1)=0.0; rx_inv(2,2)=1.0; rx_inv(2,3)=0.0; rx_inv(2,4)=0.0
    rx_inv(3,1)=0.0; rx_inv(3,2)=0.0; rx_inv(3,3)=1.0; rx_inv(3,4)=0.0
    rx_inv(4,1)=0.0; rx_inv(4,2)=0.0; rx_inv(4,3)=0.0; rx_inv(4,4)=1.0 
  end if
  
  ry(1,1)=rtmp; ry(1,2)=0.0; ry(1,3)=-r(1); ry(1,4)=0.0
  ry(2,1)=0.0;  ry(2,2)=1.0; ry(2,3)=0.0;   ry(2,4)=0.0
  ry(3,1)=r(1); ry(3,2)=0.0; ry(3,3)=rtmp;  ry(3,4)=0.0
  ry(4,1)=0.0;  ry(4,2)=0.0; ry(4,3)=0.0;   ry(4,4)=1.0

  ry_inv(1,1)=rtmp;  ry_inv(1,2)=0.0; ry_inv(1,3)=r(1); ry_inv(1,4)=0.0
  ry_inv(2,1)=0.0;   ry_inv(2,2)=1.0; ry_inv(2,3)=0.0;  ry_inv(2,4)=0.0
  ry_inv(3,1)=-r(1); ry_inv(3,2)=0.0; ry_inv(3,3)=rtmp; ry_inv(3,4)=0.0
  ry_inv(4,1)=0.0;   ry_inv(4,2)=0.0; ry_inv(4,3)=0.0;  ry_inv(4,4)=1.0

  ct=cos(theta); st=sin(theta)

  rz(1,1)=ct;  rz(1,2)=-st; rz(1,3)=0.0; rz(1,4)=0.0
  rz(2,1)=st;  rz(2,2)=ct;  rz(2,3)=0.0; rz(2,4)=0.0
  rz(3,1)=0.0; rz(3,2)=0.0; rz(3,3)=1.0; rz(3,4)=0.0
  rz(4,1)=0.0; rz(4,2)=0.0; rz(4,3)=0.0; rz(4,4)=1.0

  ra=matmul(rx,t); ra=matmul(ry,ra); ra=matmul(rz,ra)
  ra=matmul(ry_inv,ra); ra=matmul(rx_inv,ra); ra=matmul(t_inv,ra)

end function rotate_arbitrary_axis

function rotate_x(theta) result(rx)
!-----------------------
! Rotation around x-axis
!-----------------------
  use global_constants
  implicit none

  real(kind=rk) :: theta,st,ct
  real(kind=rk), dimension(3,3) :: rx

  st=sin(theta); ct=cos(theta)

  rx(1,1)=1.0; rx(1,2)=0.0; rx(1,3)=0.0
  rx(2,1)=0.0; rx(2,2)=ct;  rx(2,3)=-st
  rx(3,1)=0.0; rx(3,2)=st;  rx(3,3)=ct

end function rotate_x

function rotate_y(theta) result(ry)
!-----------------------
! Rotation around y-axis
!-----------------------
  use global_constants
  implicit none

  real(kind=rk) :: theta,st,ct
  real(kind=rk), dimension(3,3) :: ry

  st=sin(theta); ct=cos(theta)

  ry(1,1)=ct;  ry(1,2)=0.0; ry(1,3)=st
  ry(2,1)=0.0; ry(2,2)=1.0; ry(2,3)=0.0
  ry(3,1)=-st; ry(3,2)=0.0; ry(3,3)=ct

end function rotate_y

function rotate_z(theta) result(rz)
!-----------------------
! Rotation around z-axis
!-----------------------
  use global_constants
  implicit none

  real(kind=rk) :: theta,st,ct
  real(kind=rk), dimension(3,3) :: rz

  st=sin(theta); ct=cos(theta)

  rz(1,1)=ct;  rz(1,2)=-st; rz(1,3)=0.0
  rz(2,1)=st;  rz(2,2)=ct;  rz(2,3)=0.0
  rz(3,1)=0.0; rz(3,2)=0.0; rz(3,3)=1.0

end function rotate_z

function area_triangle(p1,p2,p3) result(area)
!-----------------------------------------
! Calculate the area of a general triangle
!-----------------------------------------
  use global_constants
  implicit none

  real(kind=rk) :: area
  real(kind=rk), dimension(3) :: v1,v2,p1,p2,p3
  
  v1(1)=p2(1)-p1(1); v1(2)=p2(2)-p1(2); v1(3)=p2(3)-p1(3)
  v2(1)=p3(1)-p1(1); v2(2)=p3(2)-p1(2); v2(3)=p3(3)-p1(3)

  area=0.5*norm(cross(v1,v2))

end function area_triangle

function area_quadrangle(p1,p2,p3,p4) result(area)
!-------------------------------------------
! Calculate the area of a general quadrangle
!-------------------------------------------
  use global_constants
  implicit none

  real(kind=rk) :: area
  real(kind=rk), dimension(3) :: v1,v2,p1,p2,p3,p4

  v1(1)=p3(1)-p1(1); v1(2)=p3(2)-p1(2); v1(3)=p3(3)-p1(3)
  v2(1)=p4(1)-p2(1); v2(2)=p4(2)-p2(2); v2(3)=p4(3)-p2(3)

  area=0.5*norm(cross(v1,v2))

end function area_quadrangle

function det(m) result(ans)
!-----------------------------------------------
! Calculate the determinant of a 3x3 real matrix
!-----------------------------------------------
  use global_constants
  implicit none

  real(kind=rk) :: part_1,part_2,part_3,ans
  real(kind=rk), dimension(3,3) :: m

  part_1=m(1,1)*m(2,2)*m(3,3)-m(1,1)*m(2,3)*m(3,2)
  part_2=m(1,2)*m(2,3)*m(3,1)-m(1,2)*m(2,1)*m(3,3)
  part_3=m(1,3)*m(2,1)*m(3,2)-m(1,3)*m(2,2)*m(3,1)

  ans=part_1+part_2+part_3

end function det

function cross(v1,v2) result(v3)
!---------------------
! Vector cross product
!---------------------
  use global_constants
  implicit none

  real(kind=rk), dimension(3) :: v1,v2,v3

  v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
  v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
  v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

end function cross

function norm(v) result(ans)
!------------
! Vector norm
!------------
  use global_constants
  implicit none

  real(kind=rk) :: ans
  real(kind=rk), dimension(3) :: v

  ans=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))

end function

end module mathematics
