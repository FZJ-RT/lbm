      subroutine F_results
      use para
      implicit none

      integer  x,y,z,j,n,M
      character*80 filename,names
      real*8 M_porosity
      real*8,allocatable,dimension(:):: Flux_inlet,Flux_outlet
      real*8,allocatable,dimension(:,:,:,:)::  Flux_x

      filename=trim(path)//'/F_tecplot_3D.dat'
                        open(10,file=filename)

      write(10,*)  'VARIABLES = X, Y, Z, P,u_x, obst'
      write(10,*) 'ZONE I=', lx, ', J=', ly,',K=,', lz, ', F=POINT'

      do z = 1, lz
        do y = 1, ly
          do x = 1, lx

            write(10,*) x,y,z,F_M(x,y,z,1),u_x(x,y,z),obst(x,y,z)

          enddo
        enddo
      enddo
      close(10)
      if(F_N .gt. 1)then
        ! VTK FORMAT for Paraview
        filename=trim(path)//'/plot3d_field.vtk'
        open(23,file=filename,form = 'formatted',status='replace')
        write(23,'(a26)') '# vtk DataFile Version 2.0'
        write(23,'(a31)') 'Lattice_Boltzmann_flow'
        write(23,'(a5)') 'ASCII'
        write(23,'(a24)') 'DATASET RECTILINEAR_GRID'
        write(23,'(a10,i5,i5,i5)') 'DIMENSIONS', lx, ly, lz
        write(23,'(a13,i5,a10)') 'X_COORDINATES', lx, 'float'
        write(23,*)  x_phy
        write(23,'(a13,i5,a10)') 'Y_COORDINATES', ly, 'float'
        write(23,*)  y_phy
        write(23,'(a13,i5,a10)') 'Z_COORDINATES', lz, 'float'
        write(23,*)  z_phy
        write(23,'(a10,i12)') 'POINT_DATA', lx*ly*lz
        write(23,'(a7,a10,a10)') 'SCALARS ','Structure',' float 1'
        write(23,'(a20)') 'LOOKUP_TABLE default'
        write(23,*) (((obst(x,y,z),x=1,lx),y=1,ly),z=1,lz)
        do n=2,F_N
          write(names,'(I10)') n
          write(23,'(a7,a30,a10)') 'SCALARS ','  '//trim(adjustl(names))//'_F',' float 1'
          write(23,'(a20)') 'LOOKUP_TABLE default'
      
          write(23,*) (((F_M(x,y,z,n),x=1,lx),y=1,ly),z=1,lz)
        end do
        close(23)
      endif

      ! VTK FORMAT for Paraview
      filename=trim(path)//'/plot3d_flow.vtk'
      open(24,file=filename,form = 'formatted',status='replace')
      write(24,'(a26)') '# vtk DataFile Version 2.0'
      write(24,'(a32)') 'Lattice_Boltzmann_fluid_Density'
      write(24,'(a5)') 'ASCII'
      write(24,'(a24)') 'DATASET RECTILINEAR_GRID'
      write(24,'(a10,i5,i5,i5)') 'DIMENSIONS', lx, ly, lz
      write(24,'(a13,i5,a10)') 'X_COORDINATES', lx, 'float'
      write(24,*)  x_phy
      write(24,'(a13,i5,a10)') 'Y_COORDINATES', ly, 'float'
      write(24,*)  y_phy
      write(24,'(a13,i5,a10)') 'Z_COORDINATES', lz, 'float'
      write(24,*)  z_phy
      write(24,'(a10,i12)') 'POINT_DATA', lx*ly*lz

      write(24,'(a7,a10,a10)') 'SCALARS ',' Pression',' float 1'
      write(24,'(a20)') 'LOOKUP_TABLE default'
      write(24,*) (((F_M(x,y,z,1),x=1,lx),y=1,ly),z=1,lz)
      write(24,'(a7,a10,a7)') 'VECTORS',' Velocity',' float'
      !write(23,'(a20)') 'LOOKUP_TABLE default'
      do z=1,lz
        do y=1,ly
          do x=1,lx
            write(24,'(E13.5,E13.5,E13.5)') u_x(x,y,z),u_y(x,y,z),u_z(x,y,z)
          enddo
        enddo
      enddo
    close(24)

!>>>>>>>>>>>> flux along the x axis >>>>>>>>>>>>>>
      allocate(flux_x(lx,ly,lz,F_N))
      allocate(Flux_inlet(F_N),Flux_outlet(F_N))
      filename=trim(path)//'/Sum_flow_Flux.dat'
      open(400, file=filename, iostat=ios)
      if (ios /= 0 ) stop "Error opening file flux inlet and outlet"
      do x=1, lx
        do n=1,F_N
          do  y = 1, ly
             do z=  1, lz
                if(obst(x,y,z).eq.0) then
                     flux_x(x,y,z,n)=(tao_F_x(n,x,y,z)-0.5d0)/tao_F_x(n,x,y,z)*(n_hlp_F(1,x,y,z,n)-n_hlp_F(2,x,y,z,n))*c_sqrt_F(n)
                end if
             end do
          end do
        end do

        do n=1,F_N
          Flux_inlet(n)=0.d0
          do  y = 1, ly
            do z=  1, lz
                  if(obst(x,y,z).eq.0) then
                Flux_inlet(n)=Flux_inlet(n)+flux_x(x,y,z,n)
              end if
            end do
          end do
        end do
        Flux_inlet=Flux_inlet /ly /lz
        write(400,*) x," Inlet: ", Flux_inlet(1:F_N)
      end do
      close(unit=400, iostat=ios)
      if ( ios /= 0 ) stop "Error closing file unit flux inlet and outlet"

      return
      end