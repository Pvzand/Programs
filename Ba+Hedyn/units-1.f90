real (kind=8) :: convgrau
real (kind=8) :: bohr, ang2bohr
real (kind=8) :: convt
real (kind=8) :: EK2cm1,convcmau,convaucm
real (kind=8) :: dVdrcmAng2au
real (kind=8) :: Rydbergcst,Hartree2cm !, Bohr,
real (kind=8) :: cm2Hartree,Angs2Bohr
real (kind=8) :: xm2Ang, Ang2m
real (kind=8) :: Bohr2m, xm2Bohr
real (kind=8) :: Hartree2EJ,EJ2Hartree
real (kind=8) :: hbarSI
real (kind=8) :: pau2pSI,pSI2pau
real (kind=8) :: xmsinv2au, au2xmsinv

parameter(convgrau = 1822.88852d0) ! mass in g/mol to au
! parameter(bohr = 0.529177249d0,ang2bohr=1.d0/bohr) ! length from Bohr (au) to Angstroms
! NIST version from non SI atomic units:
parameter(bohr = 0.52917720859d0,ang2bohr=1.d0/bohr) ! length from Bohr (au) to Angstroms
parameter(convt = 2.418884327d-5) ! time from au to picoseconds
!C Kelvins to cm-1 (from NIST Energy equivalents)
PARAMETER(EK2cm1 = 0.6950356d0)
!C FROM cm-1 TO au
PARAMETER(CONVCMAU=4.5563353D-6,CONVAUCM=1.d0/CONVCMAU)
parameter(dVdrcmAng2au = convcmau*bohr)
! Rydberg cst= 10 973 731.568 527(73) m-1; Hartree2cm=2*Ryd/100
! a0=0.529 177 208 59(36)×10−10m;
parameter (Rydbergcst=109737.31568527d0,Hartree2cm=2.d0*Rydbergcst)!,Bohr=0.52917720859d0)
parameter (cm2Hartree=1.d0/Hartree2cm,Angs2Bohr=1.d0/Bohr)
parameter (xm2Ang=1.d10, Ang2m=1.d-10) ! meter to Angstrom, Angstrom to meter
parameter (Bohr2m=bohr*Ang2m,xm2Bohr=xm2Ang*Ang2Bohr)  ! Bohr to meter, meter to Bohr
parameter (Hartree2EJ=4.35974394d-18,EJ2Hartree=1.d0/Hartree2EJ) ! Hartree to Joule, Joule to Hartree
!!! atomic unit of momentum = hbar/a0
parameter (hbarSI=1.054571628d-34) !!! in J s (SI units)
parameter (pau2pSI=hbarSI/Bohr2m,pSI2pau=1.d0/pau2pSI) ! momentum from au to SI, from SI to au
!!! atomic units of speed: m/s to au
parameter (xmsinv2au=xm2Bohr*(convt*1.d-12),au2xmsinv=1.d0/xmsinv2au)

    
