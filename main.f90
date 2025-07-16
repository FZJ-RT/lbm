! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is the MRT Lattice Boltzmann code for 3D convective diffusion by D3Q7 scheme
! used to simulate ionic transport processes with charge balance.
! This code package can be used for: 
! 1. Multi-species transport in heterogeneous porous media;
! 2. Anisotropic diffusion
! The code is developed by Dr. Yuankai Yang
! PNP model is used presnted in : 
! Tournassat C, Steefel C I, Gimmi T., Water Resources Research, 2020;
! MRT-LBM is presented in:          
! Yoshida H, Nagaoka M., Journal of Computational Physics, 2010.
!                         Edited @2020/3/17
!------------copyright-------------------
! TransLBM: a Lattice Boltzmann codes to simulate species transport in 
! porous media at pore or continuum scale
!
! Copyright (C) <2022>  <Yuankai Yang>
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!     
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! If you have any question about TransLBM, please feel free to send a email
! to yangyuankai@outlook.com

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      program main

      use para
      implicit none

      integer clock1,clock2,total_cost

      call system_clock(count=clock1)

      boundary_check=1 !Use Bounce-back Boundary

! ***Read Data & Initialization potential&concentration*****
      call dat
      call initial

! *** Flow field and other fields computing *********************************
      if ( F_N .gt. 0 ) then
        write(*,'(a80)') '********************************************************************************' 
        write(*,'(a40)') '*** Solving Flow Field! ***'
        call F_main 
        write(*,'(a80)') '********************************************************************************' 
        write(*,'(a40)') '*** Flow Field Finished ! ***'
        write(*,'(a80)') '********************************************************************************' 
        !write(*,'(a40,E11.4)') 'Averaged Velocity: ',u_avg
        !write(*,'(a40,E11.4)') 'Reynolds number: ',Re
        !write(*,'(a40,E11.4)') 'Peclet number: ',Pe
      end if
! **********************************************************
      write(*,"('**** Totle Time of LPM is' ,I7, '********')") tmaxcirc
      do  numbers=1,tmaxcirc
        
! *************Ioinic transport computing*******************
        call diff_main !solve concentration
        if ( mod(numbers,1000)==0 ) then
          !call diff_results
          call diff_error
          write(*,"( '                  ',I5, ' times has run')") numbers
          write(*,*) '   Step error :', steperror_diff
          write(*,*) '***************************************'
        end if      
!**************   Output unsteady results    ************
        if ( mod(numbers,tmaxcirc/Ncirc)==0 ) then
          call diff_results
        end if
!********************************************************
! Steady
        if ( 1 .eq. 1 ) then
          !call diff_error
          if (steperror_diff < 1e-6) then
            call diff_results
            write(*,*) '   Step error :', steperror_diff
            write(*,*) 'Steady state'
            stop
          end if
        end if
      end do
!********************************************************
      call system_clock(count=clock2)
      total_cost=(clock2-clock1)/1000.0
      write(*,*) 'Total Time Cost:',total_cost,'s'

      write(*,*) '*********  Congratulation ***********'
      write(*,*) '*********       end       ***********'

      stop
      end
