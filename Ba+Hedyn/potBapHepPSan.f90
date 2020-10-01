subroutine  potBapHepPSan(VpSig, VpPi, rdist,ndist,iunit)
! Ba+ - He potentials in the 4p excited states Sigma and Pi
! fitted from Fausto Cargnoni's ab initio potentials
  ! to analytical forms with gnuplot

!!! original potential: distances in Angstroms, energies in cm-1
!!! iunit = 0 to work in these units;
!!! iunit = 1 to work in atomic units
!
! implicit real (kind=8) :: (a-h,o-z)
logical first
integer :: ndist, idist
integer, parameter :: ndistx=10000
double precision, dimension (ndist) :: VpSig, VpPi, rdist
double precision, dimension (ndistx) :: vsrSig,vsrPi,VMSig,VMPi,VdWSig,VdWPi,rdistsqinv
double precision :: a1Sig, alp1Sig, a2Sig, alp2Sig
double precision :: DeSig, betSig, reSig
double precision :: aSig, bSig
double precision :: C4Sig, C6Sig
double precision :: a2Pi, alp2Pi
double precision :: DePi, betPi, rePi
double precision :: aPi, bPi
double precision :: C4Pi, C6Pi, C8Pi
double precision :: a1Sigmu, alp1Sigmu, a2Sigmu, alp2Sigmu
double precision :: DeSigmu, betSigmu, reSigmu
double precision :: aSigmu, bSigmu
double precision :: C4Sigmu, C6Sigmu
double precision :: a2Pimu, alp2Pimu
double precision :: DePimu, betPimu, rePimu
double precision :: aPimu, bPimu
double precision :: C4Pimu, C6Pimu, C8Pimu

real (kind=8) :: Vpasympcm, Vpasymp
real (kind=8) :: Ep3h, Ep1h, DeltaSOp

include 'units.f90'

!!! Sigma potential parameters (molecular units)
parameter (a1Sigmu=6.81044d6,alp1Sigmu=-3.05154d0,a2Sigmu=-143888.d0,alp2Sigmu=-0.565436d0) ! short range part
parameter (DeSigmu=1.8366d0,betSigmu=-0.643702d0,reSigmu=9.1086d0) ! Morse potential
parameter (aSigmu=0.48248d0,bSigmu=5.13892d0) ! switching function
parameter (C4Sigmu=-11128.4d0,C6igmu=-130333.d0) ! long range part
!!! Pi potential parameters (molecular units)
parameter (a2Pimu=-1.47083d6,alp2Pimu=-2.03809d0) ! short range part
parameter (DePimu=473.896d0,betPimu=-1.85676d0,rePimu=2.88162d0) ! Morse potential
parameter (aPimu=2.49651d0,bPimu=3.69071d0) ! switching function
parameter (C4Pimu=-11858.d0,C6Pimu=-225244.d0,C8Pimu=-2.74938d6) ! long range part

!  atomic (spin-orbit) NIST energy levels (cm-1)
parameter (Ep3h=21952.404d0,Ep1h=20261.561d0)

data first /.true./

save first

if (first) then
   ! asymptotic values (cm-1) from atomic (spin-orbit) NIST energy levels
   DeltaSOp = Ep3h-Ep1h
   print *,'E(P(3/2)) = ',Ep3h,' E(P(1/2)) = ',Ep1h,' DeltaSO(P) = ',DeltaSOp
   Vpasympcm = (2.d0*Ep3h+Ep1h)/3.d0
   print *,' potBapHepPSan: '
   print *,' asymptotic energies (no SO) in cm-1: (P)',Vpasympcm

   select case (iunit)

   case(0) ! molecular units
      a1Sig = a1Sigmu
      alp1Sig = alp1Sigmu
      a2Sig = a2Sigmu
      alp2Sig = alp2Sigmu
      DeSig = DeSigmu
      betSig = betSigmu
      reSig = reSigmu
      aSig = aSigmu
      bSig = bSigmu
      C4Sig = C4Sigmu
      C6Sig = C6igmu
      a2Pi = a2Pimu
      alp2Pi = alp2Pimu
      DePi = DePimu
      betPi = betPimu
      rePi = rePimu
      aPi = aPimu
      bPi = bPimu
      C4Pi = C4Pimu
      C6Pi = C6Pimu
      C8Pi = C8Pimu
      Vpasymp = Vpasympcm
      
   case(1) ! atomic units
      print *,'  potential coefficients converted to au'
      a1Sig = a1Sigmu*cm2Hartree
      alp1Sig = alp1Sigmu*Bohr
      a2Sig = a2Sigmu*cm2Hartree
      alp2Sig = alp2Sigmu*Bohr**2
      DeSig = DeSigmu*cm2Hartree
      betSig = betSigmu*Bohr
      reSig = reSigmu*Angs2Bohr
      aSig = aSigmu*Bohr
      bSig = bSigmu*Angs2Bohr
      C4Sig = C4Sigmu*cm2Hartree*Bohr**4
      C6Sig = C6igmu*cm2Hartree*Bohr**6
      a2Pi = a2Pimu*cm2Hartree
      alp2Pi = alp2Pimu*Bohr**2
      DePi = DePimu*cm2Hartree
      betPi = betPimu*Bohr
      rePi = rePimu*Angs2Bohr
      aPi = aPimu*Bohr
      bPi = bPimu*Angs2Bohr
      C4Pi = C4Pimu*cm2Hartree*Bohr**4
      C6Pi = C6Pimu*cm2Hartree*Bohr**6
      C8Pi = C8Pimu*cm2Hartree*Bohr**8
      Vpasymp = Vpasympcm*cm2Hartree
      print *,' asymptotic energies (no SO) in au: (P)',Vpasymp
   case default
      write(6,*)' >>>>> ERROR STOP in potBapHepPSan <<<<<<'
      write(6,*)' unit case iunit = ',iunit,' invalid '
      write(6,*)' choose 0 for molecular units (cm-1, Angstroms) '
      write(6,*)'   or 1 for atomic units  (Hartree, Bohr)'
   end select
   
   first=.false.
endif

if (ndist.gt.ndistx) then
   print *,'  >>>>> ERROR STOP IN PotBapHepPSan  <<<<'
   print *,' ndist = ',ndist,' > ndistx = ',ndistx
   print *,' increase ndistx '
   STOP
endif
VMSig(1:ndist)=DeSig*dexp(betSig*(rdist-reSig))*(dexp(betSig*(rdist-reSig))-2.d0)  ! Morse Sigma
VMPi(1:ndist)=DePi*dexp(betPi*(rdist-rePi))*(dexp(betPi*(rdist-rePi))-2.d0)  ! Morse Pi
vsrSig(1:ndist) = a1Sig*dexp(alp1Sig*rdist) + a2Sig*dexp(alp2Sig*rdist**2)  !  short range Sigma 
vsrPi(1:ndist) = a2Pi*dexp(alp2Pi*rdist**2)  !  short range Pi
rdistsqinv(1:ndist)=1.d0/rdist**2

VpSig=VPasymp+vsrSig(1:ndist)+VMSig(1:ndist)+0.5d0*(1.d0+dtanh(aSig*(rdist-bSig)))&
     *((C6Sig*rdistsqinv(1:ndist)+C4Sig)*rdistsqinv(1:ndist)**2-VMSig(1:ndist))
VpPi=VPasymp+vsrPi(1:ndist)+VMPi(1:ndist)+0.5d0*(1.d0+dtanh(aPi*(rdist-bPi)))*(((C8Pi*rdistsqinv(1:ndist)+C6Pi) &
     *rdistsqinv(1:ndist)+C4Pi)*rdistsqinv(1:ndist)**2-VMPi(1:ndist))
do idist=1, ndist
   write(25,'(9(1x,g15.8))') rdist(idist), VpSig(idist), VpPi(idist)
   write(26,'(4(1x,g15.8))') rdist(idist),vsrPi(idist),VMPi(idist),VpPi(idist)
enddo

!do idist=1, ndist
!   VdWSig(idist)=(C6Sig*rdistsqinv(idist)+C4Sig)*rdistsqinv(idist)**2
!   VpSig(idist)=VPasymp+vsrSig(idist)+VMSig(idist)+0.5d0*(1.d0+dtanh(aSig*(rdist(idist)-bSig)))&
!        *(VdWSig(idist)-VMSig(idist))
!   VdWPi(idist)=((C8Pi*rdistsqinv(idist)+C6Pi)*rdistsqinv(idist)+C4Pi)*rdistsqinv(idist)**2
!   VpPi(idist)=VPasymp+vsrPi(idist)+VMPi(idist)+0.5d0*(1.d0+dtanh(aPi*(rdist(idist)-bPi)))*(VdWPi(idist)-VMPi(idist))
!!   VpPi(idist)=VPasymp+vsrPi(idist)+VpPi(idist)+0.5d0*(1.d0+dtanh(aPi*(rdist(idist)-bPi)))*(((C8Pi*rdistsqinv(idist)&
!!        +C6Pi)*rdistsqinv(idist)+C4Pi)*rdistsqinv(idist)**2-VpPi(idist))
!   write(80,*) rdist(idist),VpSig(idist),VpPi(idist)
!enddo

return

end subroutine potBapHepPSan


