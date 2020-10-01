subroutine  deranpotBapHeXZPAD(dVpot, rdist,ndist)
! potential derivative for Ba+ - He in the ground state
! version CCSDT from Fausto,
  ! averaged over He ZPAD wave function by potconv-baphe.f90
  ! then fitted to an analytical form with gnuplot

!!! distances in Angstroms, energies in cm-1
  
!
  ! implicit real (kind=8) :: (a-h,o-z)
  implicit none
logical first, lordered
integer :: ndist, idist
integer, parameter :: ndistx=10000
double precision, dimension (ndist) :: dVpot, rdist
double precision, dimension (ndistx) :: vsr,dvsr,rdistsqinv,vlr,dvlr
double precision, dimension (ndistx) :: tanhr
double precision :: C4, C6, C8, a1, alp1, a2, alp2, a, b, D, r0, del
double precision :: fourC4m, sixC6m, eightC8m

parameter (C4=-7856.28d0,C6=-778533d0,C8=1.6402d7) ! long range part
parameter (fourC4m=-4.d0*C4,sixC6m=-6.d0*C6,eightC8m=-8.d0*C8)
parameter (a1=448387d0,alp1=2.12257d0,a2=512378d0,alp2=0.862913d0) ! short range part
parameter (a=1.54219d0,b=4.53313d0) ! switching function
parameter (D=279.582d0,r0=0.246636d0,del=1.74349d0) ! additional gaussian

if (ndist.gt.ndistx) then
   print *,'  >>>>> ERROR STOP IN deranPOTBAPHEXZPADfit  <<<<'
   print *,' ndist = ',ndist,' > ndistx = ',ndistx
   print *,' increase ndistx '
   STOP
endif
dVpot=-2.d0*(rdist-r0)/del**2*D*dexp(-((rdist-r0)/del)**2)  ! additional gaussian
vsr(1:ndist)=a1*dexp(-alp1*rdist) + a2*dexp(-alp2*rdist**2)  !  short range
dvsr(1:ndist)=-alp1*a1*dexp(-alp1*rdist)-2.d0*alp2*rdist*a2*dexp(-alp2*rdist**2)
rdistsqinv(1:ndist)=1.d0/rdist**2
vlr(1:ndist) = ((C8*rdistsqinv(1:ndist)+C6)*rdistsqinv(1:ndist)+C4)*rdistsqinv(1:ndist)**2
dvlr(1:ndist)=((eightC8m*rdistsqinv(1:ndist)+SixC6m)*rdistsqinv(1:ndist)+fourC4m)*rdistsqinv(1:ndist)**2/rdist
tanhr(1:ndist) = dtanh(a*(rdist-b))

dVpot=dVpot+dvsr(1:ndist)+0.5d0*(1.d0+tanhr(1:ndist))*(dvlr(1:ndist)-dvsr(1:ndist))  &
     +0.5d0*a*(1.d0-tanhr(1:ndist)**2)*(vlr(1:ndist)-vsr(1:ndist))

return

end subroutine deranpotBapHeXZPAD



