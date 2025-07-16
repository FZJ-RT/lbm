      subroutine Transport_parameter
      use para
      implicit none

      call GetPorosity

      call GetDiffusivity

      call GetPermeability

      return
      end

!*******************************************************

    subroutine GetPorosity
      use para
      implicit none      
      integer  x,y,z
      character*80 filename,read_in_porosity
      real*8 a_a,a_n



! *************  Got Porosity from file here (Diabled) *******
      read_in_porosity='Porosity.txt'
      if(1.eq.0) then
99101   FORMAT(1000F10.2)
        filename=read_in_porosity
        open(10,file=read_in_porosity)

        do z=1,lz
          do y=1,ly
           read(10,99101) (obst(x,y,z),x=1,lx)
          enddo
        enddo
        close(10)
      endif
!***********************************************************
      if(1.eq.1) then
        do z=1,lz
          do y=1,ly
            do x=1,lx
              if(obst(x,y,z) == 2)then
                Porosity(x,y,z) = sigma
                obst(x,y,z) = 0
              endif
            enddo
          enddo
        enddo
      end if

!***********************************************************
      if(1.eq.0) then
        Porosity=0.5d0
        Porosity(1:x_m,:,:)=1.0d0
        Porosity((x_m+1):lx,:,:)=0.1d0
      end if

! **************Percolation threshold*****************
! Percolation threshold here is 0.01 porosity
      do y = 1, ly
        do x = 1, lx
          do   z = 1, lz
            if(Porosity(x,y,z) .lt. 0.01d0 )then
              obst(x,y,z) = 1
            end if
          end do
        end do
      end do
    end subroutine GetPorosity

!*****************************************************
!*****************************************************

    subroutine GetDiffusivity
      use para
      implicit none      
      integer  x,y,z 
      real*8 a_a,a_n 

! Linear relationship of structure factor with porosity
      !S_Factor_x=1.d0/Porosity
      !S_Factor_y=1.d0/Porosity
      !S_Factor_z=1.d0/Porosity


! Archie's law : structure factor with porosity
      a_a = 1d0
      a_n = 1.0d0
      do y = 1, ly
        do x = 1, lx
          do   z = 1, lz
            S_Factor_x(x,y,z) = 1.d0/(a_a*Porosity(x,y,z)**a_n)
            S_Factor_y(x,y,z) = 1.d0/(a_a*Porosity(x,y,z)**a_n)
            S_Factor_z(x,y,z) = 1.d0/(a_a*Porosity(x,y,z)**a_n)
          end do
        end do
      end do

      print*, 'Minimum structure factor',minval(S_Factor_x)
      print*, 'Maximum structure factor',maxval(S_Factor_x)


! Diffusivity matrix
      do y = 1, ly
        do x = 1, lx
          do  z = 1, lz
            Diff_x(1:N_ion,x,y,z)=Di(1:N_ion)/S_Factor_x(x,y,z)
            Diff_y(1:N_ion,x,y,z)=Di(1:N_ion)/S_Factor_y(x,y,z)
            Diff_z(1:N_ion,x,y,z)=Di(1:N_ion)/S_Factor_z(x,y,z)
          end do
        end do
      end do
      print*, 'Minimum D',minval(Diff_x)
      print*, 'Maximum D',maxval(Diff_x)

    end subroutine GetDiffusivity

!*****************************************************
!*****************************************************

    subroutine GetPermeability
      use para
      implicit none      
      integer  x,y,z 
      real*8 a_a,a_n,viscosity
      real*8 k0 
      k0=1.0e-8 !
      viscosity=1.0E-6 !Dynamic viscosity m2s-1
      K_x=1.d0
      K_y=1.d0
      K_z=1.d0
! Archie's law : structure factor with porosity
      a_a = 1d0
      a_n = 1.0d0
      do y = 1, ly
        do x = 1, lx
          do   z = 1, lz
            S_Factor_x(x,y,z) = 1.d0/(a_a*Porosity(x,y,z)**a_n)
            S_Factor_y(x,y,z) = 1.d0/(a_a*Porosity(x,y,z)**a_n)
            S_Factor_z(x,y,z) = 1.d0/(a_a*Porosity(x,y,z)**a_n)
          end do
        end do
      end do
      K_x(1,:,:,:)=k0/S_Factor_x/viscosity
      K_y(1,:,:,:)=k0/S_Factor_y/viscosity
      K_z(1,:,:,:)=k0/S_Factor_z/viscosity
      print*, 'Minimum K',minval(K_x)
      print*, 'Maximum K',maxval(K_x)

    end subroutine GetPermeability