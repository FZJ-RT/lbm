      subroutine  diff_results
      use para
      implicit none

      integer  x,y,z,j,n,M
      character*80 filename,names
      real*8 M_porosity
      real*8,allocatable,dimension(:,:):: effe_c
      real*8,allocatable,dimension(:):: Flux_inlet,Flux_outlet
		  real*8,allocatable,dimension(:,:,:,:)::  Flux_x
      M=100+numbers/(tmaxcirc/Ncirc)
      write(names,'(I5)') M
      print*, 'Output No.', M

      ! VTK FORMAT for Paraview
      filename=trim(path)//'/plot3d_concentration'//trim(adjustl(names))//'.vtk'
      open(23,file=filename,form = 'formatted',status='replace')
      write(23,'(a26)') '# vtk DataFile Version 2.0'
      write(23,'(a31)') 'Lattice_Boltzmann_Concentration'
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
      do n=1,N_ion
        write(names,'(I10)') n
        write(23,'(a7,a30,a10)') 'SCALARS ','  '//trim(adjustl(names))//'_ion',' DOUBLE 1'
        write(23,'(a20)') 'LOOKUP_TABLE default'
      
        write(23,*) (((cenc_diff(x,y,z,n),x=1,lx),y=1,ly),z=1,lz)
      end do
      do n=0,6
        write(names,'(I10)') n
        write(23,'(a7,a40,a10)') 'SCALARS ','  '//trim(adjustl(names))//'_direction',' DOUBLE 1'
        write(23,'(a20)') 'LOOKUP_TABLE default'
        write(23,*) (((n_hlp_diff(n,x,y,z,1),x=1,lx),y=1,ly),z=1,lz)
      end do
      close(23)

      filename=trim(path)//'/dataseries_concentration.vtk.series'
      open(43,file=filename,form = 'formatted',status='replace')
      write(43,'(a1)') '{'
      write(43,'(a32)') '  "file-series-version" : "1.0",'
      write(43,'(a13)') '  "files" : ['
      if(M >= 102) then
        do n=101,M-1
          write(names,'(I5)') n
          write(43,'(a16,a27,a12,E10.3,a3)') '    { "name" : "', 'plot3d_concentration'//trim(adjustl(names))//'.vtk','", "time" : ', (n-100)*(tmaxcirc/Ncirc)*deitat_diff(1),' },'
        end do
      endif
      write(names,'(I5)') M
      write(43,'(a16,a27,a12,E10.3,a2)') '    { "name" : "', 'plot3d_concentration'//trim(adjustl(names))//'.vtk','", "time" : ', (n-100)*(tmaxcirc/Ncirc)*deitat_diff(1),' }'
      write(43,'(a3)') '  ]'
      write(43,'(a1)') '}'
      close(43)
!>>>>>>>>>>>> Average concentration along the x axis >>>>>>>>>>>>>>
      allocate(effe_c(lx,N_ion))
      write(names,'(I5)') M
      filename=trim(path)//'/Sum_'//trim(names)//'_diff.dat'
	  open(20,file=filename)
      do x = 1, lx
        do n=1,N_ion
          effe_c(x,n)=0.d0
          M_porosity=0
          do  y = 1, ly
            do z=  1, lz
              if(obst(x,y,z).eq.0 ) then
                effe_c(x,n)=effe_c(x,n)+cenc_diff(x,y,z,n)
                M_porosity=M_porosity+1
              end if
            end do
          end do
        end do
        write(20,*) x, effe_c(x,1:N_ion)/M_porosity
      end do
      close(20)


!>>>>>>>>>>>> Ionic flux along the x axis >>>>>>>>>>>>>>
      allocate(flux_x(lx,ly,lz,N_ion))
      allocate(Flux_inlet(N_ion),Flux_outlet(N_ion))
      filename=trim(path)//'/Sum_'//trim(names)//'_Flux.dat'
      open(400, file=filename, iostat=ios)
      if (ios /= 0 ) stop "Error opening file flux inlet and outlet"
      do x=1, lx
      	do n=1,N_ion
          do  y = 1, ly
          	do z=  1, lz
            	if(obst(x,y,z).eq.0) then
           			flux_x(x,y,z,n)=(tao_x(n,x,y,z)-0.5d0)/tao_x(n,x,y,z)* (n_hlp_diff(1,x,y,z,n) - n_hlp_diff(2,x,y,z,n))*c_sqrt_diff(n)
            	end if
          	end do
          end do
        end do

        do n=1,N_ion
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
        write(400,*) x," Inlet: ", Flux_inlet(1:N_ion) ,"mol /m2 /s"
      end do
      close(unit=400, iostat=ios)
      if ( ios /= 0 ) stop "Error closing file unit flux inlet and outlet"

      ! VTK FORMAT for Paraview
      !filename=trim(path)//'/plot3d_ionic_flux'//trim(adjustl(names))//'.vtk'
      !open(24,file=filename,form = 'formatted',status='replace')
      !write(24,'(a26)') '# vtk DataFile Version 2.0'
      !write(24,'(a31)') 'Lattice_Boltzmann_Ionic_Flux'
      !write(24,'(a5)') 'ASCII'
      !write(24,'(a24)') 'DATASET RECTILINEAR_GRID'
      !write(24,'(a10,i5,i5,i5)') 'DIMENSIONS', lx, ly, lz
      !write(24,'(a13,i5,a10)') 'X_COORDINATES', lx, 'float'
      !write(24,*)  x_phy
      !write(24,'(a13,i5,a10)') 'Y_COORDINATES', ly, 'float'
      !write(24,*)  y_phy
      !write(24,'(a13,i5,a10)') 'Z_COORDINATES', lz, 'float'
      !write(24,*)  z_phy
      !write(24,'(a10,i12)') 'POINT_DATA', lx*ly*lz
      !do n=1,N_ion
      !  write(names,'(I10)') n
      !  write(24,'(a7,a30,a10)') 'SCALARS ','  '//trim(adjustl(names))//'_ion',' DOUBLE 1'
      !  write(24,'(a20)') 'LOOKUP_TABLE default'
     ! 
     !   write(24,*) (((flux_x(x,y,z,n),x=1,lx),y=1,ly),z=1,lz)
     ! end do
     ! close(24)      
     !! filename=trim(path)//'/dataseries_flux.vtk.series'
     !!! open(44,file=filename,form = 'formatted',status='replace')
     !! write(44,'(a1)') '{'
     ! write(44,'(a32)') '  "file-series-version" : "1.0",'
     ! write(44,'(a13)') '  "files" : ['
     ! if(M >= 102) then
     !   do n=101,M-1
     !     write(names,'(I5)') n
     !     write(44,'(a16,a24,a12,E10.3,a3)') '    { "name" : "', 'plot3d_ionic_flux'//trim(adjustl(names))//'.vtk','", "time" : ', (n-100)*deitat_diff(1),' },'
     !   end do
     ! endif
     ! write(names,'(I5)') M
     ! write(44,'(a16,a24,a12,E10.3,a2)') '    { "name" : "', 'plot3d_ionic_flux'//trim(adjustl(names))//'.vtk','", "time" : ', (M-100)*deitat_diff(1),' }'
     ! write(44,'(a3)') '  ]'
     ! write(44,'(a1)') '}'
     ! close(44)

      return
      end
