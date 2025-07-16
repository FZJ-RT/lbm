      subroutine F_main
      use para
      implicit none
      integer i
      real*8 Pe_dx,u_chr,Cr,Pe
      do i=1,10000000
      
        call F_propagate
      
        call F_wall_boundary
      
        call F_boundary
      
        call F_relax

        if (steperror_F < 1e-6) then
            exit
        end if
        if (mod(i,5000)==0 ) then
          call F_error
          write(*,*) 'Step error :', steperror_F
        end if
      end do
      
      write(*,*) 'Step error :', steperror_F
      call F_results
      u_chr=maxval(u_x)
      Pe_dx=u_chr*scale/(maxval(Diff_x) + minval(Diff_x))*2.d0
      Cr=u_chr/c_sqrt_diff(1)
      Pe=Pe_dx*ly

      write(*,*)
      write(*,'(a40,E11.2)') 'Mesh Peclet Number (Must < 10): ',Pe_dx
      write(*,'(a40,E11.2)') 'Courant Number (Must < 1): ',Cr
      write(*,'(a40,E11.2)') 'Peclet Number: ',Pe
      write(*,*)
      end subroutine F_main

!************************************************
!                       Propagate step
!************************************************

      subroutine F_propagate
      use para
      implicit none

      integer  x,y,z,x_e,x_w,y_n,y_s,z_f,z_b,n
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(x,y,z,n)
      !$OMP DO
      do n = 1,F_N
        do  x = 1, lx
          do y = 1, ly
            do z = 1, lz
              y_n = mod(y,ly) + 1
              x_e = mod(x,lx) + 1
              z_f = mod(z,lz) + 1
              y_s = ly - mod(ly + 1 - y, ly)
              x_w = lx - mod(lx + 1 - x, lx)
              z_b = lz - mod(lz + 1 - z, lz)

              n_hlp_F(0 ,x  ,y,  z ,n) =  node_F(0,x,y,z       ,n)

              n_hlp_F(1 ,x ,y,   z ,n) =  node_F(1,x_w,y  ,z  ,n)
              n_hlp_F(2 ,x ,y,   z ,n) =  node_F(2,x_e,y  ,z   ,n)
              n_hlp_F(3 ,x  ,y,  z ,n) =  node_F(3,x  ,y_s,z  ,n )
              n_hlp_F(4 ,x  ,y,  z ,n) =  node_F(4,x  ,y_n,z  ,n )
              n_hlp_F(5 ,x  ,y,  z ,n) =  node_F(5,x  ,y  ,z_b ,n)
              n_hlp_F(6 ,x  ,y,  z ,n) =  node_F(6,x  ,y  ,z_f ,n)

            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL
      end subroutine F_propagate
!***********************************************************************************************  
!************************************************
!                     Boundary step
!************************************************
      
      subroutine F_wall_boundary

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
        do n=1, F_N
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
                  n_hlp_F(2,x,y,z,n) =  n_hlp_F(1,x_e  ,y,z,n)
                end if
                ! *** 2   ***************************!***
                if( obst(x  ,y ,z  ).eq.0  .and.  obst(x_w,y ,z  )>=1 )then
                  n_hlp_F(1,x,y,z,n) =  n_hlp_F(2,x_w ,y  ,z  ,n)
                end if
                ! *** 3   ******************************
                if( obst(x  ,y ,z ).eq.0  .and. obst(x  ,y_n ,z )>=1  )then
                  n_hlp_F(4,x,y,z,n) =  n_hlp_F(3,x  ,y_n, z,n)
                end if
                ! *** 4   ******************************
                if( obst(x  ,y  ,z ).eq.0  .and.  obst(x  ,y_s,z )>=1 )then
                  n_hlp_F(3,x,y,z,n) =  n_hlp_F(4,x  ,y_s, z,n)
                end if
                ! *** 5   ******************************
                if( obst(x ,y  ,z ).eq.0  .and.  obst(x  ,y  ,z_f )>=1 )then
                  n_hlp_F(6,x,y,z,n) =  n_hlp_F(5,x ,y,z_f,n)
                end if
                ! *** 6   ******************************
                if( obst(x ,y  ,z ).eq.0  .and. obst(x  ,y  ,z_b )>=1  )then
                  n_hlp_F(5,x,y,z,n) =  n_hlp_F(6,x  ,y  ,z_b,n)
                end if
              end do
            end do
          end do
        end do
      !$OMP END DO    
      !$OMP END PARALLEL
      end if

  end subroutine F_wall_boundary
!*********************************************************
!                     Inlet & Outlet
!*********************************************************
  
  subroutine F_boundary
      use para
      implicit none
      integer x,y,z,n
      n=1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(y,z)
!$OMP DO  
! *     inlet
        do  y = 1, ly
          do  z = 1, lz
            if(obst(1,y,z).eq.0)then
              if(F_inlet(n).lt.0.d0)then
                n_hlp_F(1,1,y,z,n)= node_F(1,1,y,z,n) +2.d0/Lambda*t_d_1*deitat_F(n)/scale*(-F_inlet(n))
              else
                n_hlp_F(1,1,y,z,n)=2*t_d_1*F_inlet(n)-node_F(2,1,y,z,n)
              end if
            end if
          end do
        end do
!$OMP END DO    
!$OMP END PARALLEL

! *     outlets
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(y,z)
!$OMP DO 
        do  y = 1, ly
          do  z = 1, lz
            if(obst(lx,y,z).eq.0)then
              if(F_outlet(n).lt.0.d0)then
                n_hlp_F(2,lx,y,z,n)= node_F(2,lx,y,z,n) +2.d0/Lambda*t_d_1*deitat_F(n)/scale*(-F_outlet(n))
              else
                n_hlp_F(2,lx,y,z,n)=2*t_d_1*F_outlet(n)-node_F(1,lx,y,z,n)
              end if
            end if
          end do
        end do
!$OMP END DO    
!$OMP END PARALLEL

      if(F_N .gt. 1)then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(y,z,n)
!$OMP DO  
      do n=2,F_N
! *     inlet
        do  y = 1, ly
          do  z = 1, lz
            if(obst(1,y,z).eq.0)then
              if(F_inlet(n).lt.0.d0)then
                n_hlp_F(1,1,y,z,n)=node_F(1,1,y,z,n)
              else
                n_hlp_F(1,1,y,z,n)=2*t_d_1*F_inlet(n)-node_F(2,1,y,z,n)
              end if
            end if
          end do
        end do
      end do

!$OMP END DO    
!$OMP END PARALLEL

! *     outlets
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(y,z,n)
!$OMP DO 
      do n=2,F_N
        do  y = 1, ly
          do  z = 1, lz
            if(obst(lx,y,z).eq.0)then
              if(F_outlet(n).lt.0.d0)then
                n_hlp_F (2,lx,y,z,n) =  node_F(2,lx,y,z,n)
              else
                n_hlp_F(2,lx,y,z,n)=2*t_d_1*F_outlet(n)-node_F(1,lx,y,z,n)
              end if
            end if
          end do
        end do
      end do
!$OMP END DO    
!$OMP END PARALLEL
      end if

!************* up / down ***************
! if bounce-back
    if ( 1 .eq. 0 ) then
      do n=1,F_N
        y=ly
        do x=1,lx
          do z=1,lz
            if ( obst(x,y,z)==0 ) then
              n_hlp_F(4,x,y,z,n)=node_F(3,x,y,z,n)
            end if
          end do
        end do

        y=1
        do x=1,lx
          do z=1,lz
            if ( obst(x,y,z)==0 ) then
              n_hlp_F(3,x,y,z,n)=node_F(4,x,y,z,n)
            end if
          end do
        end do
      end do
    end if

    end subroutine F_boundary

!***********************************************************************************************  
!**********************************************************
!             Collision step
!**********************************************************

      subroutine F_relax
      use para
      implicit none
      integer  x,y,z, i,j,n
      real*8   conc,tempa
      real*8   Qmat(0:6)

      Qmat =0.d0



!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(x,y,z,n,i)
!$OMP DO 
      do n=1,F_N
        do  x = 1, lx
          do  y = 1, ly
            do  z = 1,lz
              if(obst(x,y,z).eq.0 )then
! ****  ******  computing macroscopic concentration  **************
                F_M(x,y,z,n)=0.d0
                do  i=0,6
                  F_M(x,y,z,n) = F_M(x,y,z,n) + n_hlp_F(i,x,y,z,n)
                end do
              end if
            end do
          end do
        enddo
      end do
!$OMP END DO    
!$OMP END PARALLEL
      
      u_x=0.d0
      u_y=0.d0
      u_z=0.d0
      if ( F_N .gt. 0 )then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(x,y,z,i)
!$OMP DO 
        do  x = 1, lx
          do  y = 1, ly
            do  z = 1,lz
              if(obst(x,y,z).eq.0 )then
! ****  ******  computing velocity  **************
                do  i=0,6
                  u_x(x,y,z) = u_x(x,y,z)+(tao_F_x(1,x,y,z)-0.5d0)/tao_F_x(1,x,y,z)*e_x(i)*n_hlp_F(i,x,y,z,1)*c_sqrt_F(1)
                  u_y(x,y,z) = u_y(x,y,z)+(tao_F_y(1,x,y,z)-0.5d0)/tao_F_y(1,x,y,z)*e_y(i)*n_hlp_F(i,x,y,z,1)*c_sqrt_F(1) 
                  u_z(x,y,z) = u_z(x,y,z)+(tao_F_z(1,x,y,z)-0.5d0)/tao_F_z(1,x,y,z)*e_z(i)*n_hlp_F(i,x,y,z,1)*c_sqrt_F(1) 
                end do
              end if
            end do
          end do
        enddo
!$OMP END DO    
!$OMP END PARALLEL
      end if
!***********************************************************************************************  
!!*************Equilibrium***************
!!*************************************

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(x,y,z,n,i)
!$OMP DO 
        do n=1,F_N
          do  x = 1, lx
            do  y = 1, ly
              do  z = 1,lz
                if(obst(x,y,z).eq.0 )then
                  do i = 0,6
                    n_equ_F( i,x,y,z,n) = F_M(x,y,z,n) * n_equ_e(i)
                  end do
                end if
              end do
            end do
          enddo
        end do
!$OMP END DO    
!$OMP END PARALLEL


!!*************Collision***************
!!******SRT mode******
      if(MRT.eq.0)then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(x,y,z,n,i)
!$OMP DO 
        do n=1,F_N
          do  x = 1, lx
            do  y = 1, ly
              do  z = 1,lz
                if(obst(x,y,z).eq.0 )then
                  do i = 0,6
                    Qmat(i)=n_equ_F(i,x,y,z,n) - n_hlp_F(i,x,y,z,n) 
                    node_F(i,x,y,z,n) = n_hlp_F(i,x,y,z,n) + 1.d0/tao_F_x(n,x,y,z)*Qmat(i)
                  end do
                end if
              end do
            end do
          end do
        end do
!$OMP END DO    
!$OMP END PARALLEL
      end if



!******MRT mode******
      if(MRT.eq.1)then

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(x,y,z,n,i)
!$OMP DO 
        do n=1,F_N
          do  x = 1, lx
            do  y = 1, ly
              do  z = 1,lz
                if(obst(x,y,z).eq.0 )then
                  node_F(0,x,y,z,n)= n_hlp_F(0,x,y,z,n) + (n_equ_F(0,x,y,z,n) - n_hlp_F(0,x,y,z,n)) 
                  node_F(1,x,y,z,n)= n_hlp_F(1,x,y,z,n) + ( 1.d0/tao_F_x(n,x,y,z)+1.d0)*0.5d0*(n_equ_F(1,x,y,z,n) &
                    - n_hlp_F(1,x,y,z,n))+(-1.d0/tao_F_x(n,x,y,z)+1.d0)*0.5d0*(n_equ_F(2,x,y,z,n) - n_hlp_F(2,x,y,z,n))
                  node_F(2,x,y,z,n)= n_hlp_F(2,x,y,z,n) + (-1.d0/tao_F_x(n,x,y,z)+1.d0)*0.5d0*(n_equ_F(1,x,y,z,n) &
                    - n_hlp_F(1,x,y,z,n))+( 1.d0/tao_F_x(n,x,y,z)+1.d0)*0.5d0*(n_equ_F(2,x,y,z,n) - n_hlp_F(2,x,y,z,n)) 
                  node_F(3,x,y,z,n)= n_hlp_F(3,x,y,z,n) + ( 1.d0/tao_F_y(n,x,y,z)+1.d0)*0.5d0*(n_equ_F(3,x,y,z,n) &
                    - n_hlp_F(3,x,y,z,n))+(-1.d0/tao_F_y(n,x,y,z)+1.d0)*0.5d0*(n_equ_F(4,x,y,z,n) - n_hlp_F(4,x,y,z,n)) 
                  node_F(4,x,y,z,n)= n_hlp_F(4,x,y,z,n) + (-1.d0/tao_F_y(n,x,y,z)+1.d0)*0.5d0*(n_equ_F(3,x,y,z,n) &
                    - n_hlp_F(3,x,y,z,n))+( 1.d0/tao_F_y(n,x,y,z)+1.d0)*0.5d0*(n_equ_F(4,x,y,z,n) - n_hlp_F(4,x,y,z,n)) 
                  node_F(5,x,y,z,n)= n_hlp_F(5,x,y,z,n) + ( 1.d0/tao_F_z(n,x,y,z)+1.d0)*0.5d0*(n_equ_F(5,x,y,z,n) &
                    - n_hlp_F(5,x,y,z,n))+(-1.d0/tao_F_z(n,x,y,z)+1.d0)*0.5d0*(n_equ_F(6,x,y,z,n) - n_hlp_F(6,x,y,z,n)) 
                  node_F(6,x,y,z,n)= n_hlp_F(6,x,y,z,n) + (-1.d0/tao_F_z(n,x,y,z)+1.d0)*0.5d0*(n_equ_F(5,x,y,z,n) &
                    - n_hlp_F(5,x,y,z,n))+( 1.d0/tao_F_z(n,x,y,z)+1.d0)*0.5d0*(n_equ_F(6,x,y,z,n) - n_hlp_F(6,x,y,z,n)) 
                end if
              end do
            end do
          enddo
        end do
!$OMP END DO    
!$OMP END PARALLEL
      end if

      end subroutine F_relax

!***********************************************
!        Computing Error

      subroutine  F_error
      use para
      implicit none
      integer x,y,z,i
      real*8 F_sum

      F_sum  = 0.d0
      y=y_m

        do x=1,lx
          do z=1,lz
            if ( obst(x,y_m,z).eq.0 .and. &
              F_sum < abs(Old_F(x,y_m,z,1)-F_M(x,y_m,z,1))/abs(F_M(x,y_m,z,1))) then
              F_sum= abs(Old_F(x,y_m,z,1)-F_M(x,y_m,z,1))/abs(F_M(x,y_m,z,1))
            end if
          end do
        enddo
        steperror_F= F_sum
        Old_F=F_M

      end subroutine  F_error
