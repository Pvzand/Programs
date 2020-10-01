!!! Calculation of derivatives of the He-He interaction potential:  d/dx_k( Sum_(i<j) V_He-He(R_ij))
!!! Forme LM2M2 pour le potentiel de He-He
!!! R.A. Aziz and M. J. Slaman, J. Chem. Phys., 94 (12), 8047-8053 (1991).

!!! UNITS: Angstroms, cm-1, use iunit to change units
!!! (so far only iunit=0 (unchanged) or 1 (atomic units) is coded)
 
subroutine derpotHeHeZPADu(dV22,V22,dist,ndist,iunit)
  
  implicit none

!!! Parameters for He-He effective potential from ZPAD Ne4He100.
!!!  Rcut = 1.6 Angs.

  integer :: ndist, iunit
  real*8, dimension(ndist) :: dist, dV22,V22
  integer, parameter :: ndistx=1000
  real*8 :: asrmu
  real*8 :: atrmu,btrmu,remu,demu,betamu,xemu
  real*8 :: C6mu,c8mu,c10mu
  real*8 :: asr,alphasr,atr,btr,C6,c8,c10,re,reinv,de,beta,xe,twobetade
  real*8 :: sixC6m, eightC8m, tenC10m, twelveC12m

  logical :: lfirst=.true.

  parameter(asrmu = -430128.d0, alphasr = 18.742d0)   !!! cm-1,  dimensionless
  parameter(atrmu = 4.67503d0, btrmu = 4.71596d0)     !!! Angstrom-1, Angstrom
  parameter(C6mu =  -1.1876d0, c8mu = -0.728204d0, c10mu =  -0.411248d0) !!! cm-1
  parameter(remu = 4.2d0)    !!! Angstroms
  parameter(demu=1.45437d0,betamu=-8.30807d0,xemu=1.00088d0)

  include 'units.f90'

  integer :: i,ipos,ipos2,j,jpos,jpos2,k,ndim,nbashe,ndimtot
  real (kind=8), dimension(ndistx) ::  distred,distredinv,distredinvsq,vsr,dvsr,vlr,dvlr,vtrans,dvtrans,vs1,vM,dvM,dvsr1,vsr1,expM

  if (lfirst) then
     print *,' in derpotHeHeZPADu, iunit = ',iunit
     print *,' He-He potential parameters '
        print *,'Original values in molecular units (Angstroms, cm-1)'
        print *,' asr,alphasr = ',asrmu,alphasr
        print *,' atr, btr, re = ',atrmu, btrmu, remu
        print *,' C6,C8,C10,C12 = ',C6mu,C8mu,C10mu
     select case (iunit)
     case(0)   !!! here stick to molecular units
        print *,' potential parameters are used in original units'
        asr = asrmu
        atr = atrmu
        btr = btrmu
        re = remu
        C6 = C6mu
        C8 = C8mu
        C10 = C10mu
        
        de = demu
        xe = xemu
        beta = betamu

     case(1)   !!! here convert to atomic units (Bohr, Hartree)
        print *,' potential parameters converted to atomic units (Bohr, Hartree): '
        asr = asrmu*cm2Hartree
        atr = atrmu*Bohr
        btr = btrmu*Angs2Bohr
        re = remu*Angs2Bohr
        C6 = C6mu * cm2Hartree
        C8 = C8mu * cm2Hartree
        C10 = C10mu * cm2Hartree
        de = demu* cm2Hartree
        xe = xe* cm2Hartree
        beta =betamu* cm2Hartree

    case default
      write(6,*)' >>>>> ERROR STOP in derpotHeHeZPADu <<<<<<'
      write(6,*)' unit case iunit = ',iunit,' invalid '
      write(6,*)' choose 0 for molecular units (cm-1, Angstroms) '
      write(6,*)'   or 1 for atomic units  (Hartree, Bohr)'
   end select
   print *,' asr,alphasr = ',asr,alphasr
   print *,' atr, btr, re = ',atr, btr, re
   print *,' C6,C8,C10 = ',C6,C8,C10
   reinv = 1.d0/re
   sixC6m=-6.d0*C6
   eightC8m=-8.d0*C8
   tenC10m=-10.d0*C10
   lfirst = .false.
endif


  if (ndist.gt.ndistx) then
     print *,'  >>>>> ERROR STOP IN derpotHeHeZPAD  <<<<'
     print *,' ndist = ',ndist,' > ndistx = ',ndistx
     print *,' increase ndistx '
     STOP
  endif
!     cre = bohr/re   ! to be used for input distances in Bohrs
  distred(1:ndist) = dist(1:ndist)*reinv
  distredinv(1:ndist) = 1.d0/distred(1:ndist)
  distredinvsq(1:ndist) = distredinv(1:ndist)**2
  
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
  dvlr(1:ndist) =(((twelveC12m*distredinvsq(1:ndist)+tenC10m)*distredinvsq(1:ndist) + &
       eightC8m)*distredinvsq(1:ndist)+sixC6m)*reinv*distredinv(1:ndist)*distredinvsq(1:ndist)**3

! Derivative of the switching function
  dvtrans(1:ndist) = 0.5d0*atr*(1.0d0-dtanh(atr*(dist(1:ndist)-btr))**2)

! He-He potential
  v22(1:ndist) =  vsr(1:ndist) + vtrans(1:ndist)*(vlr(1:ndist) - vsr(1:ndist))  ! for result in cm-1
! print *,'dist, v22 = ',dist(1:ndist), v22(1:ndist)

!!! Total derivative
  dV22 = ( dvsr(1:ndist) + dvtrans(1:ndist)*(vlr(1:ndist) - vsr(1:ndist)) + vtrans(1:ndist)*(dvlr(1:ndist) - dvsr(1:ndist)) )
! print *,' dV22 = ',dV22(1:ndist)
end subroutine derpotHeHeZPADu

