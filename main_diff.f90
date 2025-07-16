      subroutine diff_main
      use para
      implicit none
      
      call diff_propagate
      
      call diff_wall_boundary
      
      call diff_concentration_boundary
      
      call diff_relax

      end subroutine diff_main

!************************************************
!                       Propagate step
!************************************************
      subroutine diff_propagate
      use para
      implicit none

      integer  x,y,z,x_e,x_w,y_n,y_s,z_f,z_b,n
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(x,y,z,n)
      !$OMP DO
      do n = 1,N_ion
        do  x = 1, lx
          do y = 1, ly
            do z = 1, lz
              y_n = mod(y,ly) + 1
              x_e = mod(x,lx) + 1
              z_f = mod(z,lz) + 1
              y_s = ly - mod(ly + 1 - y, ly)
              x_w = lx - mod(lx + 1 - x, lx)
              z_b = lz - mod(lz + 1 - z, lz)

              n_hlp_diff(0 ,x  ,y,  z ,n) =  node_diff(0,x,y,z       ,n)

              n_hlp_diff(1 ,x ,y,   z ,n) =  node_diff(1,x_w,y  ,z  ,n)
              n_hlp_diff(2 ,x ,y,   z ,n) =  node_diff(2,x_e,y  ,z   ,n)
              n_hlp_diff(3 ,x  ,y,  z ,n) =  node_diff(3,x  ,y_s,z  ,n )
              n_hlp_diff(4 ,x  ,y,  z ,n) =  node_diff(4,x  ,y_n,z  ,n )
              n_hlp_diff(5 ,x  ,y,  z ,n) =  node_diff(5,x  ,y  ,z_b ,n)
              n_hlp_diff(6 ,x  ,y,  z ,n) =  node_diff(6,x  ,y  ,z_f ,n)

            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL
      end subroutine diff_propagate
!***********************************************************************************************  
!************************************************
!                     Boundary step
!************************************************
      subroutine diff_wall_boundary

      use para
      implicit none

! *********************************************
      integer x,y,z,x_e,x_w,y_n,y_s,z_f,z_b,n
      real*8  n_equ

!*********************************************
!     Bounce back  boundary
!*********************************************
      if (boundary_check .eq. 1 ) then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(x,y,z,n)
      !$OMP DO
        do n=1, N_ion
          do  x = 1, lx
            do  y =1, ly
              do  z =1,lz
                y_n = mod(y,ly) + 1
                x_e = mod(x,lx) + 1
                z_f = mod(z,lz) + 1
                y_s = ly - mod(ly + 1 - y, ly)
                x_w = lx - mod(lx + 1 - x, lx)
                z_b = lz - mod(lz + 1 - z, lz)
                ! *** 1   ******************************
                if( obst(x  ,y ,z  ).eq.0  .and. obst(x_e,y ,z  )>=1 )then
                  n_hlp_diff(2,x,y,z,n) =  n_hlp_diff(1,x_e  ,y,z,n)
                end if
                ! *** 2   ***************************!***
                if( obst(x  ,y ,z  ).eq.0  .and.  obst(x_w,y ,z  )>=1 )then
                  n_hlp_diff(1,x,y,z,n) =  n_hlp_diff(2,x_w ,y  ,z  ,n)
                end if
                ! *** 3   ******************************
                if( obst(x  ,y ,z ).eq.0  .and. obst(x  ,y_n ,z )>=1  )then
                  n_hlp_diff(4,x,y,z,n) =  n_hlp_diff(3,x  ,y_n, z,n)
                end if
                ! *** 4   ******************************
                if( obst(x  ,y  ,z ).eq.0  .and.  obst(x  ,y_s,z )>=1 )then
                  n_hlp_diff(3,x,y,z,n) =  n_hlp_diff(4,x  ,y_s, z,n)
                end if
                ! *** 5   ******************************
                if( obst(x ,y  ,z ).eq.0  .and.  obst(x  ,y  ,z_f )>=1 )then
                  n_hlp_diff(6,x,y,z,n) =  n_hlp_diff(5,x ,y,z_f,n)
                end if
                ! *** 6   ******************************
                if( obst(x ,y  ,z ).eq.0  .and. obst(x  ,y  ,z_b )>=1  )then
                  n_hlp_diff(5,x,y,z,n) =  n_hlp_diff(6,x  ,y  ,z_b,n)
                end if
              end do
            end do
          end do
        end do
      !$OMP END DO    
      !$OMP END PARALLEL
      end if
!************* up / down ***************
! if bounce-back
    if ( 1 .eq. 1 ) then
      do n=1,N_ion
        y=ly
        do x=1,lx
          do z=1,lz
            if ( obst(x,y,z)==0 ) then
              n_hlp_diff(4,x,y,z,n)=node_diff(3,x,y,z,n)
            end if
          end do
        end do

        y=1
        do x=1,lx
          do z=1,lz
            if ( obst(x,y,z)==0 ) then
              n_hlp_diff(3,x,y,z,n)=node_diff(4,x,y,z,n)
              !Kinetic Reactive model
              !n_hlp_diff(3,x,y,z,n)= node_diff(4,x,y,z,n) & 
              !+2.d0/Lambda*t_d_1*Diff_y(n,x,y,z)*deitat_diff(n)/scale* & 
              !(-0.001d0/Diff_y(n,x,y,z)*(cenc_diff(x,y,z,n)-1.5d0*concent_inlet(n)))
            end if
          end do
        end do
      end do
    end if
  end subroutine diff_wall_boundary
!*********************************************************
!                     Inlet & Outlet
!*********************************************************
  subroutine diff_concentration_boundary
      use para
      implicit none
      integer y,z,n
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(y,z,n)
!$OMP DO  
      do n=1,N_ion
! *    concentration inlet
        do  y = 1, ly
          do  z = 1, lz
            if(obst(1,y,z).eq.0)then
              if(concent_inlet(n).lt.0)then
                n_hlp_diff(1,1,y,z,n)=node_diff(1,1,y,z,n)
              else
                n_hlp_diff(1,1,y,z,n)=2*t_d_1*concent_inlet(n)-node_diff(2,1,y,z,n)
              end if
            end if
          end do
        end do
      end do

!$OMP END DO    
!$OMP END PARALLEL
!   n_hlp_diff(1,1,:,:,:)=node_diff(1,1,:,:,:)
!   n_hlp_diff(1,1,45:55,:,:)=2*t_d_1*concent_inlet(1)-node_diff(2,1,45:55,:,:)

! *    concentration outlets
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(y,z,n)
!$OMP DO 
      do n=1,N_ion
        do  y = 1, ly
          do  z = 1, lz
            if(obst(lx,y,z).eq.0)then
              if(concent_outlet(n).lt.0)then
                n_hlp_diff (2,lx,y,z,n) =  node_diff(2,lx,y,z,n)
              else
                n_hlp_diff(2,lx,y,z,n)=2*t_d_1*concent_outlet(n)-node_diff(1,lx,y,z,n)
              end if
            end if
          end do
        end do
      end do
!$OMP END DO    
!$OMP END PARALLEL
    end subroutine diff_concentration_boundary
!***********************************************************************************************  
!**********************************************************
!             Collision step
!**********************************************************

      subroutine diff_relax
      use para
      implicit none
      integer  x,y,z, i,j,n
      real*8   conc,tempa,S,S_Cx,S_Cy,S_Cz
      real*8   b1,b2
      real*8   S_C_x,S_C_y,S_C_z
      real*8   Qmat(0:6)

      Qmat =0.d0
      S=0.d0


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(x,y,z,n,i)
!$OMP DO 
      do n=1,N_ion
        do  x = 1, lx
          do  y = 1, ly
            do  z = 1,lz
              if(obst(x,y,z).eq.0 )then
! ****  ******  computing macroscopic concentration  **************
                cenc_diff(x,y,z,n)=0.d0
                do  i=0,6
                  cenc_diff(x,y,z,n) = cenc_diff(x,y,z,n) + n_hlp_diff(i,x,y,z,n)
                end do
              end if
            end do
          end do
        enddo
      end do
!$OMP END DO    
!$OMP END PARALLEL
!***********************************************************************************************  
!!*************Equilibrium***************
!!*************************************
      if(PNP_FICK.eq.0)then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(x,y,z,n,i)
!$OMP DO 
        do n=1,N_ion
          do  x = 1, lx
            do  y = 1, ly
              do  z = 1,lz
                if(obst(x,y,z).eq.0 )then
                  do i = 0,6
                    n_equ_diff( i,x,y,z,n) = cenc_diff(x,y,z,n) * n_equ_e(i)
                  end do
                end if
              end do
            end do
          enddo
        end do
!$OMP END DO    
!$OMP END PARALLEL
      end if



!!*************Collision***************
!!******SRT mode******
      if(MRT.eq.0)then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(x,y,z,n,i)
!$OMP DO 
        do n=1,N_ion
          do  x = 1, lx
            do  y = 1, ly
              do  z = 1,lz
                if(obst(x,y,z).eq.0 )then
                  do i = 0,6
                    Qmat(i)=n_equ_diff(i,x,y,z,n) - n_hlp_diff(i,x,y,z,n) 
                    node_diff(i,x,y,z,n) = n_hlp_diff(i,x,y,z,n) + 1.d0/tao_x(n,x,y,z)*Qmat(i) + S*n_equ_e(i)*deitat_diff(n) 
                  end do
                end if
              end do
            end do
          end do
        end do
!$OMP END DO    
!$OMP END PARALLEL
      end if


      end subroutine diff_relax
