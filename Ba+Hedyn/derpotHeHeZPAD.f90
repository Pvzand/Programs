!!! Calculation of derivatives of the He-He interaction potential:  d/dx_k( Sum_(i<j) V_He-He(R_ij))
!!! Forme LM2M2 pour le potentiel de He-He
!!! R.A. Aziz and M. J. Slaman, J. Chem. Phys., 94 (12), 8047-8053 (1991).

!!! UNITS: Angstroms, cm-1
 
subroutine derpotHeHeZPAD(dV22,V22,dist,ndist)
  
  implicit none

!!! Parameters for He-He effective potential from ZPAD Ne4He100.
!!!  Rcut = 1.6 Angs.

  integer :: ndist
  real*8, dimension(ndist) :: dist, dV22,V22
  integer, parameter :: ndistx=1000
  real*8 :: asr,alphasr,atr,btr,C6,c8,c10,re,reinv,de,beta,xe
  real*8 :: sixC6m, eightC8m, tenC10m, twelveC12m,twobetade

  parameter(asr =  -430128.0d0, alphasr = 18.742d0, atr = 4.67503d0, btr = 4.71596d0)
  parameter(C6 = -1.1876d0, c8 = -0.728204d0, c10 = -0.411248d0)
  parameter(sixC6m=-6.d0*C6, eightC8m=-8.d0*C8, tenC10m=-10.d0*C10)
  parameter(re = 4.2d0)
  parameter(de=1.45437d0,beta=-8.30807d0,xe=1.00088d0) 

  integer :: i,ipos,ipos2,j,jpos,jpos2,k,ndim,nbashe,ndimtot
  real (kind=8), dimension(ndistx) ::  distred,distredinv,distredinvsq,vsr,dvsr,vlr,dvlr,vtrans,dvtrans,vsr1,vM,expM,dvM,dvsr1

  if (ndist.gt.ndistx) then
     print *,'  >>>>> ERROR STOP IN derpotHeHeZPAD  <<<<'
     print *,' ndist = ',ndist,' > ndistx = ',ndistx
     print *,' increase ndistx '
     STOP
  endif
!     cre = bohr/re   ! to be used for input distances in Bohrs
  reinv = 1.d0/re
  distred(1:ndist) = dist(1:ndist)*reinv
  distredinv(1:ndist) = 1.d0/distred(1:ndist)
  distredinvsq(1:ndist) = distredinv(1:ndist)**2
  
  twobetade = 2.d0*beta*de

! Short range potential
  vsr1(1:ndist) = asr*dexp(-alphasr*distred(1:ndist))
  
  !Morse potential
  
  vM(1:ndist) = de*(exp(beta*(distred(1:ndist)-xe)))*(exp(beta*(distred(1:ndist)-xe))-2.d0)

  !Total short range potential

  vsr(1:ndist)=vsr1(1:ndist)+vM(1:ndist)


! Long range potential 
  vlr(1:ndist) = ((C10*distredinvsq(1:ndist) + C8)*distredinvsq(1:ndist)+C6)*distredinvsq(1:ndist)**3
  
! Switching function
  vtrans(1:ndist) = 0.5d0*(1.0d0+dtanh(atr*(dist(1:ndist)-btr)))
 
! Derivative of the short-range potential
  dvsr1(1:ndist) = -asr*alphasr*reinv*dexp(-alphasr*distred(1:ndist))

!Derivative of Morse potential
  
  expM(1:ndist) = dexp(beta*(distred(1:ndist)-xe))
  dvM(1:ndist) = twobetade*expM(1:ndist)*(expM(1:ndist)-1.d0)

!Derivative of total short range potential
 dvsr(1:ndist) = dvsr1(1:ndist)+dvM(1:ndist)

! Derivative of the long-range potential
  dvlr(1:ndist) =((tenC10m*distredinvsq(1:ndist) + &
       eightC8m)*distredinvsq(1:ndist)+sixC6m)*reinv*distredinv(1:ndist)*distredinvsq(1:ndist)**3

! Derivative of the switching function
  dvtrans(1:ndist) = 0.5d0*atr*(1.0d0-dtanh(atr*(dist(1:ndist)-btr))**2)

! He-He potential
  v22(1:ndist) =  vsr(1:ndist) + vtrans(1:ndist)*(vlr(1:ndist) - vsr(1:ndist))  ! for result in cm-1

! Total derivative

  dV22 = ( dvsr(1:ndist) +dvtrans(1:ndist)*(vlr(1:ndist) - vsr(1:ndist)) + vtrans(1:ndist)*(dvlr(1:ndist) - dvsr(1:ndist)) )

end subroutine derpotHeHeZPAD
