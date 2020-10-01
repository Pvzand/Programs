subroutine  potBapHeXZPADfit(Vpot, rdist,ndist)
! Ba+ - He potential in the ground state
! version CCSDT from Fausto,
  ! averaged over He ZPAD wave function by potconv-baphe.f90
  ! then fitted to an analytical form with gnuplot

!!! distances in Angstroms, energies in cm-1
  
!
! implicit real (kind=8) :: (a-h,o-z)
logical first, lordered
integer :: ndist, idist
integer, parameter :: ndistx=10000
double precision, dimension (ndist) :: Vpot, rdist,vlr
double precision, dimension (ndistx) :: vsr,rdistsqinv
double precision :: C4, C6, C8, a1, alp1, a2, alp2, a, b, D, r0, del

parameter (C4=-7856.28d0,C6=-778533d0,C8=1.6402d7) ! long range part
parameter (a1=448387d0,alp1=2.12257d0,a2=512378d0,alp2=0.862913d0) ! short range part
parameter (a=1.54219d0,b=4.53313d0) ! switching function
parameter (D=279.582d0,r0=0.246636d0,del=1.74349d0) ! additional gaussian

if (ndist.gt.ndistx) then
   print *,'  >>>>> ERROR STOP IN POTBAPHEXZPADfit  <<<<'
   print *,' ndist = ',ndist,' > ndistx = ',ndistx
   print *,' increase ndistx '
   STOP
endif
Vpot=-D*dexp(-((rdist-r0)/del)**2)  ! additional gaussian
vsr(1:ndist)=a1*dexp(-alp1*rdist)+a2*dexp(-alp2*rdist**2)  !  short range 
rdistsqinv(1:ndist)=1.d0/rdist**2
vlr(1:ndist)=((C8*rdistsqinv(1:ndist)+C6)*rdistsqinv(1:ndist)+C4)*rdistsqinv(1:ndist)**2
Vpot=Vpot+vsr(1:ndist)+0.5d0*(1.d0+dtanh(a*(rdist-b)))*(vlr(1:ndist)-vsr(1:ndist))

return

end subroutine potBapHeXZPADfit

