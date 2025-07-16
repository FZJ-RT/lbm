      subroutine dat
      use para
      implicit none

      integer x,y,z,i,n
      integer If_openmp,Threads_number
      real*8  taoc
      real*8 Max_c_diff,Ne_inf,Max_D
      character*80 filename,mk_dir

      allocate(Di(Max_N_ion),concent_inlet(Max_N_ion),&
        concent_outlet(Max_N_ion),valence(Max_N_ion),concent_inital(Max_N_ion))
      allocate(F_inlet(Max_N_ion),F_outlet(Max_N_ion))

      open (1, file='inputtest.h',Form='formatted', iostat=ios)
      if (ios /= 0 ) stop "Error opening file 1"
      read(1,*)
      read(1,*) scale          ! dx= Lattice size
      read(1,*)
      read(1,'(I10)') lx              !Domain size
      read(1,*)
      read(1,'(I10)') ly
      read(1,*)
      read(1,'(I10)') lz
      read(1,*)
      read(1,'(I10)') N_ion      !Number of ions
      read(1,*)
      read(1,'(F10.2)') temp       !Temperature
      read(1,*)
      read(1,*) sigma       !Dp
      read(1,*)  
      read(1,'(I5)') If_openmp !  IF OPENMP
      read(1,'(I5)') Threads_number !Parallel computing cores number
      read(1,'(a80)') Path
		  read(1,'(a80)') read_in_structure
      read(1,*)
      read(1,*) (Di(i),i=1,N_ion)      ! Diffusivity
      read(1,*)
      read(1,*) (concent_inlet(i),i=1,N_ion)    !Concentration of inlet
      read(1,*) (concent_outlet(i),i=1,N_ion) !Concentration of outlet
      read(1,*)
      read(1,*) (concent_inital(i),i=1,N_ion) !Concentration of innital
      read(1,*)
      read(1,*) (valence(i),i=1,N_ion) !Ion valence
      read(1,*)
      read(1,'(I5)') F_N
      read(1,*)
      read(1,*) (F_inlet(i),i=1,F_N)    !inlet boundary
      read(1,*) (F_outlet(i),i=1,F_N)   !outlet boundary
      read(1,*)
      read(1,*)
      read(1,'(I5)') MRT !  IF MRT
      read(1,*)
      read(1,'(I5)') PNP_FICK !  IF PNP
      read(1,*)
      read(1,'(I10)') tmaxcirc !  Total steps for simulation
      read(1,*)
      read(1,'(I5)') Ncirc ! Number of Outputs 
      close(unit=1, iostat=ios)
      if ( ios /= 0 ) stop "Error closing file unit 1"

!**************   Electric neutrality *****************
      !For inlet boundary condition 
      Ne_inf=0
      do i = 1, N_ion
        Ne_inf=Ne_inf+concent_inlet(i)*valence(i)
      end do
      if ( Ne_inf ==0 ) then
        write(*,*) "Electric neutrality at inlet is satisfied"
      else
        write(*,*) "ERROR! Electric non-neutrality at inlet"
        stop
      end if

      !For outlet boundary condition
      Ne_inf=0
      do i = 1, N_ion
        Ne_inf=Ne_inf+concent_outlet(i)*valence(i)
      end do
      if ( Ne_inf ==0 ) then
        write(*,*) "Electric neutrality at outlet is satisfied"
      else
        write(*,*) "ERROR! Electric non-neutrality at outlet"
        stop
      end if

      !For initial condition
      Ne_inf=0
      do i = 1, N_ion
        Ne_inf=Ne_inf+concent_inital(i)*valence(i)
      end do
      if ( Ne_inf ==0 ) then
        write(*,*) "Electric neutrality of inital C is satisfied"
      else
        write(*,*) "ERROR! Non-electricneutrality of inital C"
        stop
      end if
!**************         Ion bulk concentration          *********************
      do i = 1, N_ion
        if ( valence(i)>0 ) then
          Ne_inf=Ne_inf+concent_inital(i)
        end if
      end do
!********************************************************************
      write(*,*) 'Read data from file and shown as follows:'
      write(*,'(" lx=                      ", i3)')lx
      write(*,'(" ly=                      ", i3)')ly
      write(*,'(" lz=                      ", i3)')lz
      write(*,'(" Lattice size=            ", E10.3)')scale
      write(*,'(" Temperature=             ", F10.2)')temp
      write(*,'(" Precip diff D=          ", E10.3)')sigma
      write(*,'(" Ion concentration:       ", E10.3)')Ne_inf
      write(*,'(" MRT or SRT (1 MRT, 0 SRT):", i2)') MRT
      write(*,'(" PNP or FICK (1 PNP):      ", i2)') PNP_FICK
      write(*,'(" The number of Ion kinds:  ", i3)')N_ion
      write(*,'(" The number of field:  ", i3)')F_N
      write(*,*) " Diffusion coefficient:"
      do i=1,N_ion
        write(*,"(1XE10.3,$)") Di(i)
      end do
      write(*,*)
      write(*,*) " Concent_inlet:"
      do i=1,N_ion
        write(*,"(1XE10.3,$)") concent_inlet(i)
      end do
      write(*,*)
      write(*,*) " Concent_outlet:"
      do i=1,N_ion
        write(*,"(1XE10.3,$)") concent_outlet(i)
      end do
      write(*,*)
      write(*,*) " concent_inital:"
      do i=1,N_ion
        write(*,"(1XE10.3,$)") concent_inital(i)
      end do
      write(*,*)
      write(*,*) "Other field boundary:"
      if(F_N .gt. 0)then
        write(*,*) (F_inlet(i),i=1,F_N)    !inlet boundary
        write(*,*) (F_outlet(i),i=1,F_N)   !outlet boundary
      end if
      write(*,*)
      write(*,*) "Ion valence"
      do i=1,N_ion
        write(*,"(1XE10.1,$)") valence(i)
      end do
      write(*,*) '-----------------------------------------------------'
      write(*,*) '  '

! *****Make the new directory********************************************

      mk_dir='mkdir '//trim(path)
      call system(adjustl(trim(mk_dir)))

!  *********************************************************************
!	     WRITE PARAMETERS
!  *********************************************************************

	    filename=trim(path)//'/Write_Parameter.dat'
      open(2, file=filename, iostat=ios)
      if (ios /= 0 ) stop "Error opening file 2"
      write(2,'(" lx=                      ", i3)')lx
      write(2,'(" ly=                      ", i3)')ly
      write(2,'(" lz=                      ", i3)')lz
      write(2,'(" Lattice size=            ", E10.3)')scale
      write(2,'(" Temperature=             ", F10.2)')temp
      write(2,'(" Precip diff D=          ", E10.3)')sigma
      write(2,'(" MRT or SRT (1 MRT, 0 SRT):", i2)') MRT
      write(2,'(" PNP or FICK (1 PNP):      ", i2)') PNP_FICK
      write(2,'(" Ion concentration:       ", E10.3)')Ne_inf
      write(2,*)
      write(2,'(" The number of Ion kinds:  ", i3)')N_ion
      write(2,'(" The number of field:  ", i3)')F_N
      write(2,*) " Diffusion coefficient:"
      do i=1,N_ion
        write(2,"(1XE10.3,$)") Di(i)
      end do
      write(2,*)
      write(2,*) " Concent_inlet:"
      do i=1,N_ion
        write(2,"(1XE10.3,$)") concent_inlet(i)
      end do
      write(2,*)
      write(2,*) " Concent_outlet:"
      do i=1,N_ion
        write(2,"(1XE10.3,$)") concent_outlet(i)
      end do
      write(2,*)
      write(2,*) " Concent_inital:"
      do i=1,N_ion
        write(2,"(1XE10.3,$)") concent_inital(i)
      end do      
      write(2,*)
      write(2,*) "Ion valence"
      do i=1,N_ion
        write(2,"(1XE10.1,$)") valence(i)
      end do
      write(2,*)
      write(2,*) "Other field boundary:"
      if(F_N .gt. 0)then
        write(2,"(1XE10.3,$)") (F_inlet(i),i=1,F_N)    !inlet boundary
        write(2,"(1XE10.3,$)") (F_outlet(i),i=1,F_N)   !outlet boundary
      end if
      write(2,*)
      write(2,*) '  '
      close(2, iostat=ios)
      if (ios/= 0 ) stop "Error closing file unit 2"

! ***** Open Multi-Processing ******************************

        if ( If_openmp .eq. 1 ) then
          call OMP_SET_NUM_THREADS(Threads_number)
        end if

! **** Define the size of the matrix***************************
       allocate(obst(lx,ly,lz),Porosity(lx,ly,lz),S_Factor_x(lx,ly,lz),S_Factor_y(lx,ly,lz),S_Factor_z(lx,ly,lz))
       allocate(e_x(0:6),e_y(0:6),e_z(0:6))
       allocate(dsigma(lx))
       allocate(c_sqrt_diff(N_ion),deitat_diff(N_ion),omega_diff(N_ion))
       allocate(n_equ_diff(0:6,lx,ly,lz,N_ion),node_diff(0:6,lx,ly,lz,N_ion),n_hlp_diff(0:6,lx,ly,lz,N_ion))
       allocate(cenc_diff(lx,ly,lz,N_ion),Old_diff(lx,ly,lz,N_ion))
       allocate(tao_d(N_ion),tao_x(N_ion,lx,ly,lz),tao_y(N_ion,lx,ly,lz),tao_z(N_ion,lx,ly,lz))
       allocate(Diff_x(N_ion,lx,ly,lz),Diff_y(N_ion,lx,ly,lz),Diff_z(N_ion,lx,ly,lz))
       allocate(n_equ_e(0:6))
       allocate(x_phy(lx),y_phy(ly),z_phy(lz))
       allocate(different_x_c(N_ion,lx,ly,lz),different_y_c(N_ion,lx,ly,lz),different_z_c(N_ion,lx,ly,lz))
       allocate(different_x_p(lx,ly,lz),different_y_p(lx,ly,lz),different_z_p(lx,ly,lz))

       allocate(c_sqrt_F(F_N),deitat_F(F_N),omega_F(F_N))
       allocate(n_equ_F(0:6,lx,ly,lz,F_N),node_F(0:6,lx,ly,lz,F_N),n_hlp_F(0:6,lx,ly,lz,F_N))
       allocate(F_M(lx,ly,lz,F_N),Old_F(lx,ly,lz,F_N))
       allocate(tao_F(F_N),tao_F_x(F_N,lx,ly,lz),tao_F_y(F_N,lx,ly,lz),tao_F_z(F_N,lx,ly,lz))
       allocate(K_x(F_N,lx,ly,lz),K_y(F_N,lx,ly,lz),K_z(F_N,lx,ly,lz))
       allocate(u_x(lx,ly,lz),u_y(lx,ly,lz),u_z(lx,ly,lz))

! ******************************************************
!      Cycles paremeters Initialization
! ******************************************************
      steperror_diff = 1.d0
      steperror_F = 1.d0
      obst = 0.d0
      Old_diff = 0d0
      Old_F = 0.d0
      cenc_diff  = 0d0
      F_M = 0.d0
      different_x_c = 0.d0
      different_y_c = 0.d0
      different_z_c = 0.d0
      different_x_p = 0.d0
      different_y_p = 0.d0
      different_z_p = 0.d0
      Diff_x = 0d0
      Diff_y = 0d0
      Diff_z = 0d0
      K_x = 0d0
      K_y = 0d0
      K_z = 0d0      
      Porosity = 1d0
      S_Factor_x = 1d0
      S_Factor_y = 1d0
      S_Factor_z = 1d0
      dsigma=0d0
      u_x=0.d0
      u_y=0.d0
      u_z=0.d0

! ******************************************************
!      The constant of distribution of mass and energy
! ******************************************************
!      Equilibrium distribution
       t_d_0 = 1.d0 / 4.d0
       t_d_1 = 1.d0 / 8.d0

       n_equ_e( 0)= t_d_0

       n_equ_e( 1)= t_d_1
       n_equ_e( 2)= t_d_1
       n_equ_e( 3)= t_d_1
       n_equ_e( 4)= t_d_1
       n_equ_e( 5)= t_d_1
       n_equ_e( 6)= t_d_1

! ****** Find middle of domain  *******************************
      x_m=int(0.5*lx)
      y_m=int(0.5*ly)
      z_m=int(0.5*lz)
      if(x_m == 0) x_m=1
      if(y_m == 0) y_m=1
      if(z_m == 0) z_m=1
      write(*,'(" Monitoring Point:      ", i3,i3,i3)') x_m, y_m, z_m
! ******        velocity    ************************************
!   The direction of the velocity are numbered as follows:
!		For D3Q7
!                z ^ 5    4
!                  |     /
!                  |   /
!                  | /       x
!       2 -------- 0 ------->1
!                / |
!              /   |
!         3  / y   |
!                  6
       e_x=(/0.,1.,	-1., 0.,  0.,	0.,	 0./)

       e_y=(/0.,0.,	 0., 1., -1.,	0.,	 0./)

       e_z=(/0.,0.,	 0., 0.,  0.,	1.,	-1./)

! ****Get structures **********************************
      call GetStruc
      call Transport_parameter
      !if(sigma .NE. 0d0)then
      !  call charge_domain
      !endif
! **************************************************************
!                    Model paremeters
! **************************************************************
!************* Diffusion *****************      
      print*,'Parameters for diffusion:'
      Max_c_diff=0.d0
      do i = 1, N_ion
        tao_d(i)=1.0d0  !relaxation time of diffusion
        c_sqrt_diff(i)=Di(i)/minval(S_Factor_x)/(tao_d(i)-0.5d0)/Lambda/scale
      end do
      Max_c_diff = 0.1d0*maxval(c_sqrt_diff)
      print*,'lattice speed : ',Max_c_diff
      deitat_diff(:)= scale/Max_c_diff
      print*,'Diffusion time step : ',deitat_diff(1)

      do i= 1,N_ion
        c_sqrt_diff(i) = Max_c_diff
        do x = 1, lx
          do y = 1, ly
            do z = 1, lz
              tao_x(i,x,y,z) = Diff_x(i,x,y,z)/Max_c_diff/Lambda/scale+0.5
              tao_y(i,x,y,z) = Diff_y(i,x,y,z)/Max_c_diff/Lambda/scale+0.5
              tao_z(i,x,y,z) = Diff_z(i,x,y,z)/Max_c_diff/Lambda/scale+0.5
            end do
          end do
        end do
        tao_d(i)      = Di(i)/minval(S_Factor_x)/Max_c_diff/Lambda/scale+0.5d0
        deitat_diff(i)= scale/c_sqrt_diff(i)
        omega_diff(i) = 1.d0/tao_d(i)
      end do
      print*, 'Minimum Tau',minval(tao_x)
      print*, 'Maximum Tau',maxval(tao_x)

!************* Flow *****************
      if(F_N .gt. 0)then 
        print*,'Parameters for flow:'
        tao_F(1)=1.0d0  !relaxation time of diffusion
        c_sqrt_F(1)=maxval(K_x(1,:,:,:))/(tao_F(1)-0.5d0)/Lambda/scale
        print*,'lattice speed : ',c_sqrt_F(1)
        deitat_F(:)= scale/c_sqrt_F(1)
        print*,'Flow time step : ',deitat_F(1)
        do x = 1, lx
          do y = 1, ly
            do z = 1, lz
              tao_F_x(1,x,y,z) = K_x(1,x,y,z)/c_sqrt_F(1)/Lambda/scale+0.5
              tao_F_y(1,x,y,z) = K_y(1,x,y,z)/c_sqrt_F(1)/Lambda/scale+0.5
              tao_F_z(1,x,y,z) = K_z(1,x,y,z)/c_sqrt_F(1)/Lambda/scale+0.5
            end do
          end do
        end do
        omega_F(1) = 1.d0/tao_F(1)
        print*, 'Minimum Tau',minval(tao_F_x(1,:,:,:))
        print*, 'Maximum Tau',maxval(tao_F_x(1,:,:,:))
      end if


!****************************************************************
      do x=1,lx
        x_phy(x)=(x-0.5d0)*scale
      enddo

      do z=1,lz
        z_phy(z)=(z-dble(lz)/2-0.5d0)*scale
!       z_phy(z)=(z-0.5d0)*scale
      enddo

      do y=1,ly
        y_phy(y)=(y-dble(ly)/2-0.5d0)*scale
!       y_phy(y)=(y-0.5d0)*scale
      enddo

      return
      end



! **************************************************************

      subroutine GetStruc
      use para
      implicit none
      integer  x,y,z
      character*80 filename

! *********     The domain and boundary           **************

        do y = 1, ly
          do x = 1, lx
            do   z = 1, lz
                obst(x,y,z) =0
            end do
          end do
        end do

! *************  Got structure from file here  (Diabled) *******
      if(1.eq.1) then
99001   FORMAT(1000I2)
        filename=read_in_structure
	      open(9,file=read_in_structure)

        do z=1,lz
          do y=ly,1,-1
           read(9,99001) (obst(x,y,z),x=1,lx)
          enddo
        enddo
        close(9)
      endif

! **************** Generate convex plate structure  (Enabled)  ************
      if(1.eq.0) then
        do y=1,ly
	        do x= x_m-15,x_m+15
	           obst(x,y,1:30)=1
             obst(x,y,(lz-29):lz)=1
	        end do
	      end do
	    end if
! **************** Generate parallel plates (Enabled)  ************
      if(1.eq.0) then
        do y=1,ly
      	  do x= 1,lx
      	    obst(x,y,1)=1
            obst(x,y,lz)=1
      	  end do
      	end do
      end if

      end subroutine GetStruc