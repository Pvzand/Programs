subroutine  derpotpotBapHepPSan(VpSig, VpPi, dVpSig, dVpPi, rdist,ndist,iunit)
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
double precision, dimension (ndist) :: VpSig, VpPi, dVpsig, dVpPi, rdist
double precision, dimension (ndistx) :: vsr1Sig,vsr2Sig,vsrPi,VMSig,VMPi,expMSig,expMPi,VdWSig,VdWPi&
                                       ,tanhSig,tanhPi,fSig,fPi,rdistsqinv
double precision, dimension (ndistx) :: dvsrSig,dvsrPi,dVMSig,dVMPi,dVdWSig,dVdWPi,dfSig,dfPi,rdistinv
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

double precision :: twoalp2Sig
double precision :: twoalp2Pi
double precision :: twobetDeSig
double precision :: twobetDePi
double precision :: fourC4mSig, sixC6mSig
double precision :: fourC4mPi, sixC6mPi, eightC8mPi

real (kind=8) :: Vpasymp
real (kind=8) :: Ep3h, Ep1h, DeltaSOp

double precision :: convE, convL
double precision :: Rydberg,Hartree2cm, Bohr,cm2Hartree,Angs2Bohr


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

! Rydberg cst= 10 973 731.568 527(73) m-1; Hartree2cm=2*Ryd/100
! a0=0.529 177 208 59(36)×10−10m;
parameter (Rydbergcst=109737.31568527d0,Hartree2cm=2.d0*Rydbergcst,Bohr=0.52917720859d0)
parameter (cm2Hartree=1.d0/Hartree2cm,Angs2Bohr=1.d0/Bohr) 

data first /.true./

save first

if (first) then
   ! asymptotic values (cm-1) from atomic (spin-orbit) NIST energy levels
   DeltaSOp = Ep3h-Ep1h
   print *,'E(P(3/2)) = ',Ep3h,' E(P(1/2)) = ',Ep1h,' DeltaSO(P) = ',DeltaSOp
   Vpasymp = (2.d0*Ep3h+Ep1h)/3.d0
   print *,' asymptotic energies (no SO) in cm-1: (P)',Vpasymp

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
      
   case(1) ! atomic units
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
   case default
      write(6,*)' >>>>> ERROR STOP in potBapHepPSan <<<<<<'
      write(6,*)' unit case iunit = ',iunit,' invalid '
      write(6,*)' choose 0 for molecular units (cm-1, Angstroms) '
      write(6,*)'   or 1 for atomic units  (Hartree, Bohr)'
   end select
   twoalp2Sig = 2.d0*alp2Sig
   twoalp2Pi = 2.d0*alp2Pi
   twobetDeSig = 2.d0*betSig*DeSig
   twobetDePi = 2.d0*betPi*DePi
   fourC4mSig = -4.d0*C4Sig
   sixC6mSig = -6.d0*C6Sig
   fourC4mPi = -4.d0*C4Pi
   sixC6mPi = -6.d0*C6Pi
   eightC8mPi = -8.d0*C8Pi
   
   first=.false.
endif

if (ndist.gt.ndistx) then
   print *,'  >>>>> ERROR STOP IN PotBapHepPSan  <<<<'
   print *,' ndist = ',ndist,' > ndistx = ',ndistx
   print *,' increase ndistx '
   STOP
endif
!!! !  short range Sigma 
vsr1Sig(1:ndist) = a1Sig*dexp(alp1Sig*rdist)
vsr2Sig(1:ndist) = a2Sig*dexp(alp2Sig*rdist**2)
!!!!  short range Pi
vsrPi(1:ndist) = a2Pi*dexp(alp2Pi*rdist**2)
!!! ! Morse Sigma
expMSig(1:ndist) = dexp(betSig*(rdist-reSig))
VMSig(1:ndist)=DeSig*expMSig(1:ndist)*(expMSig(1:ndist)-2.d0)
!!! Morse Pi
expMPi(1:ndist) = dexp(betPi*(rdist-rePi))
VMPi(1:ndist)=DePi*expMPi(1:ndist)*(expMPi(1:ndist)-2.d0)
!!! switching function Sigma, Pi
tanhSig(1:ndist) = dtanh(aSig*(rdist-bSig))
fSig(1:ndist) = 0.5d0*(1.d0+tanhSig(1:ndist))
tanhPi(1:ndist) = dtanh(aPi*(rdist-bPi))
fPi(1:ndist) = 0.5d0*(1.d0+tanhPi(1:ndist))
!!! Van der Waals Sigma, Pi
rdistinv(1:ndist)=1.d0/rdist
rdistsqinv(1:ndist)=rdistinv(1:ndist)**2
VdWSig(1:ndist) = (C6Sig*rdistsqinv(1:ndist)+C4Sig)*rdistsqinv(1:ndist)**2
VdWPi(1:ndist) = ((C8Pi*rdistsqinv(1:ndist)+C6Pi)*rdistsqinv(1:ndist)+C4Pi)*rdistsqinv(1:ndist)**2
!!! Total potential Sigma, Pi
VpSig = VPasymp + vsr1Sig(1:ndist) + vsr2Sig(1:ndist) + VMSig(1:ndist) + fSig(1:ndist)&
     *(VdWSig(1:ndist)-VMSig(1:ndist))
VpPi = VPasymp + vsrPi(1:ndist) + VMPi(1:ndist) + fPi(1:ndist)*(VdWPi(1:ndist)-VMPi(1:ndist))
do idist = 1, ndist
   write(80,'(5(1x,g15.8))') rdist(idist), vsrPi(idist),VMPi(idist),VdWPi(idist),VpPi(idist)
enddo


!!! short range derivative Sigma, Pi
dvsrSig(1:ndist) = alp1Sig*vsr1Sig(1:ndist) + twoalp2Sig*rdist*vsr2Sig(1:ndist)
dvsrPi(1:ndist) = twoalp2Pi*rdist*vsrPi(1:ndist)

!!! Morse derivative Sigma, Pi
dVMSig(1:ndist) = twobetDeSig*expMSig(1:ndist)*(expMSig(1:ndist)-1.d0)
dVMPi(1:ndist) = twobetDePi*expMPi(1:ndist)*(expMPi(1:ndist)-1.d0)

!!! Van der Waals derivative Sigma, Pi
dVdWSig(1:ndist) =  (sixC6mSig*rdistsqinv(1:ndist)+fourC4mSig)*rdistinv(1:ndist)*rdistsqinv(1:ndist)**2
dVdWPi(1:ndist) = ((eightC8mPi*rdistsqinv(1:ndist)+sixC6mPi)*rdistsqinv(1:ndist)+fourC4mPi)*rdistinv(1:ndist)*rdistsqinv(1:ndist)**2

!!! switching function derivative, Sigma, Pi
dfSig(1:ndist) = 0.5d0*aSig*(1.d0-tanhSig(1:ndist)**2)
dfPi(1:ndist) = 0.5d0*aPi*(1.d0-tanhPi(1:ndist)**2)

!!! total potential derivative Sigma, Pi
dVpSig = dvsrSig(1:ndist) + dVMSig(1:ndist) + fSig(1:ndist)*(dVdWSig(1:ndist)-dVMSig(1:ndist)) &
     + dfSig(1:ndist)*(VdWSig(1:ndist)-VMSig(1:ndist))
dVpPi = dvsrPi(1:ndist) + dVMPi(1:ndist) + fPi(1:ndist)*(dVdWPi(1:ndist)-dVMPi(1:ndist)) &
     + dfPi(1:ndist)*(VdWPi(1:ndist)-VMPi(1:ndist))

!do idist=1, ndist
!   write(25,'(5(1x,g15.8))') rdist(idist), VpSig(idist), VpPi(idist), dVpSig(idist), dVpPi(idist)
!enddo

!do idist=1, ndist
!   VdWSig(idist)=(C6Sig*rdistsqinv(idist)+C4Sig)*rdistsqinv(idist)**2
!   VpSig(idist)=VPasymp+vsrSig(idist)+VpSig(idist)+0.5d0*(1.d0+dtanh(aSig*(rdist(idist)-bSig)))&
!        *(VdWSig(idist)-VpSig(idist))
!   VdWPi(idist)=((C8Pi*rdistsqinv(idist)+C6Pi)*rdistsqinv(idist)+C4Pi)*rdistsqinv(idist)**2
!   VpPi(idist)=VPasymp+vsrPi(idist)+VpPi(idist)+0.5d0*(1.d0+dtanh(aPi*(rdist(idist)-bPi)))*(VdWPi(idist)-VpPi(idist))
!!   VpPi(idist)=VPasymp+vsrPi(idist)+VpPi(idist)+0.5d0*(1.d0+dtanh(aPi*(rdist(idist)-bPi)))*(((C8Pi*rdistsqinv(idist)&
!!        +C6Pi)*rdistsqinv(idist)+C4Pi)*rdistsqinv(idist)**2-VpPi(idist))
!   write(80,*) rdist(idist),VpSig(idist),VpPi(idist)
!enddo

return

end subroutine derpotpotBapHepPSan



