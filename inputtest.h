!dx= Lattice size
9.0475e-06
512
!grid number along y
512
!grid number along z
1
!Number of ions <=20
1
!temperature
298.15
!Diffusivity in precipitation Ds/D0= 
0.1
!  IF OPENMP 1 yes 0 no and  Parallel computing cores number  (not yet) 
0
1
data
structure.txt
! Diffusivity m2/s
${base_diffusivity}$
! Inlet&outlet concentration mol/m3 （-1 mean no flux）
1.0
0.0
! Inital concentration mol/m3
0.0
! Valences
0
!F_N equation number:0 No other equation, 1 flow, (2 temperature not yet) 
0
!Inlet & Outlet pressure (-1 means -flux)
100.0
100.0

!  IF MRT (1 MRT, 0 SRT)
0
!  IF PNP (1 PNP, 0 FICK)
0
!Total steps for simulation
500000000
!Number of Outputs 
100000
