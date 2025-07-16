      subroutine initial
      use para
      implicit none
      integer  x,y,z,j,n
      real*8 cenc
      do n=1,N_ion
        
        do x = 1, lx
          cenc=concent_inital(n)
          do y = 1, ly
            do z=  1, lz
              if ( obst(x,y,z).eq.0 ) then
                do j=0,6
                  node_diff(j ,x,y,z, n)  =  n_equ_e(j) * cenc
                end do
              end if
            end do
          end do
        end do
      end do

      do n=1,F_N
        do x = 1, lx
          do y = 1, ly
            do z=  1, lz
              if ( obst(x,y,z).eq.0 ) then
                do j=0,6
                  node_F(j,x,y,z,n)  =  n_equ_e(j) * F_outlet(n)
                end do
              end if
            end do
          end do
        end do
      end do
      
      return
      end
