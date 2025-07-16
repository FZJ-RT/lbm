      subroutine  diff_error
      use para
      implicit none
      integer x,y,z,i

        m_sum  =0.d0
        y=y_m

        do x=1,lx
          do z=1,lz
            if ( obst(x,y_m,z).eq.0 .and. &
              m_sum < abs(Old_diff(x,y_m,z,1)-cenc_diff(x,y_m,z,1))/abs(cenc_diff(x,y_m,z,1))) then
              m_sum= abs(Old_diff(x,y_m,z,1)-cenc_diff(x,y_m,z,1))/abs(cenc_diff(x,y_m,z,1))
            end if
          end do
        enddo
        !write(*,*) '   Step error :', m_sum
        !write(*,*) '***************************************'
        steperror_diff= m_sum
        Old_diff=cenc_diff

      return
      end
