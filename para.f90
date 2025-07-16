  module para
  implicit none

  real*8 N_a,Pi,Lambda,e,k,epsilon
  integer Max_N_ion
  integer lx,ly,lz                        !geometry size
  parameter(Pi=3.1415926535897933D+0)
  parameter(e=1.602E-19)
  parameter(Lambda=0.25)             !D3Q7 1/4
  parameter(k=1.380622E-23 )       !Boltzmann constant
  parameter(epsilon=6.95e-10 )      !C^2/J m  !Permittivity of a fluid
  parameter(N_a=6.02e23  )           !Avogadro number
  parameter(Max_N_ion=20  )       !Maximum number of ion kinds

  character*80 path
  character*80 read_in_structure

  integer time,numbers,tmaxcirc,Ncirc
  integer,allocatable,dimension(:,:,:)::  obst             !structure
  integer x_m,y_m,z_m  !Monitor Point
  integer ios,boundary_check
  integer MRT,PNP_FICK
  real*8 scale,sigma
  real*8 steperror_diff,fsumdiff_all ! Error of Cycles
  real*8,allocatable,dimension(:)::  dsigma
  real*8,allocatable,dimension(:)::  n_equ_e,e_x,e_y,e_z
  integer N_ion !Number of ions
  real*8,allocatable,dimension(:):: valence         !Ion algebraic valence !!################################
  real*8 temp,Ex
  real*8,allocatable,dimension(:):: x_phy,y_phy,z_phy
  real*8,allocatable,dimension(:,:,:):: Porosity,S_Factor_x,S_Factor_y,S_Factor_z

!***********Nernst Planck Equation*****************
  integer F_N 
  real*8,allocatable,dimension(:)::  Di  !Diffusivity
  real*8,allocatable,dimension(:,:,:,:)::  Diff_x,Diff_y,Diff_z  !Diffusiv
  
  real*8,allocatable,dimension(:)::  c_sqrt_diff,deitat_diff
  real*8,allocatable,dimension(:)::  tao_d,omega_diff
  real*8,allocatable,dimension(:,:,:,:)::  tao_x,tao_y,tao_z

  real*8,allocatable,dimension(:,:,:,:):: different_x_c,different_y_c,different_z_c
  real*8,allocatable,dimension(:,:,:):: different_x_p,different_y_p,different_z_p
  
  ! The concentration of innital & inlet & outlet
  real*8,allocatable,dimension(:):: concent_inital,concent_inlet,concent_outlet 
  real*8  concent_in,concent_out

  real*8 t_d_0,t_d_1

  real*8,allocatable,dimension(:,:,:,:,:):: n_equ_diff,node_diff,n_hlp_diff  !####scource_diff
  real*8,allocatable,dimension(:,:,:,:):: Old_diff,cenc_diff

  real*8 m_sum
  
  !***********Other Equation*****************

  real*8 steperror_F,fsumF_all ! Error of Cycles
  real*8,allocatable,dimension(:,:,:,:)::  K_x,K_y,K_z 
  
  real*8,allocatable,dimension(:)::  c_sqrt_F,deitat_F
  real*8,allocatable,dimension(:)::  tao_F,omega_F
  real*8,allocatable,dimension(:,:,:,:)::  tao_F_x,tao_F_y,tao_F_z

  real*8,allocatable,dimension(:):: F_inital,F_inlet,F_outlet 

  real*8,allocatable,dimension(:,:,:,:,:):: n_equ_F,node_F,n_hlp_F  !####scource_diff
  real*8,allocatable,dimension(:,:,:,:):: Old_F,F_M
  real*8,allocatable,dimension(:,:,:):: u_x,u_y,u_z
	end module
