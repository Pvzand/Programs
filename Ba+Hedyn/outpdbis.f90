SUBROUTINE OUTPdbis(T,Y,DERY,ndimy,prmt,nprmt)
!!! output subroutine for HPCG integrator
!!! in Re(dk), Im(dk) for the coefficients
!!! rather than mod(ck), phase(ck)
!!! where dk = ck exp(i int_0^t Vk/hbar dt'), Vk being the k-th eigenvalue attached to eigenvector k
!!! This form should avoid describing many unnecessary oscillations compared to Re(ck), Im(ck),
!!! while avoiding numerical undeterminacy when a coefficient tends to zero with mod(ck), phase(ck)
!!! Added energy phases as variables: phik=int_0^t Vk/hbar dt'

!!! This version to be used before beginning a trajectory (no previous call to fct)
  
  implicit none

  integer ndimy, nprmt,iat
  real (kind=8) :: t
  real (kind=8), dimension (ndimy) :: Y, DERY
  real (kind=8), dimension (nprmt) :: prmt

  include 'nbrat12d.f90'
  include 'units.f90'

  real (kind=8) :: Etot, Ekin, Epot, EpotHe

  integer, parameter :: iunit = 1   !!! to change input potentials from (Angs,cm-1) to atomic units

!!! for common/masses/
  real (kind=8) :: xmass1, xmass2
  real (kind=8) :: xmass     !!!  dimension (nbrat12)   (dimensioned in the COMMON itself)
  real (kind=8) :: xmassinv  !!!  dimension (ncoord12)  (dimensioned in the COMMON itself)
!!! for common/dyn/
  integer :: ihop, ioldhop, ihopinic, nbhop
!!! for common nadcplng (non-adiabatic couplings)
  real (kind=8), dimension (ndimorb,ndimorb) :: rdkj_SOR, rdkj_SOI
!!! for common/DIM
  real (kind=8), dimension(ndimorb) ::  eigenEcomm
  real (kind=8), dimension(ndimorb) ::  eigenE   !!! added for the call to DIMpSOdynanZPAD
  real (kind=8), dimension(ndimorb,ndimorb) ::  vecr_SOR, vecr_SOI
  real (kind=8), dimension(ndimorb,ndimorb) ::  vecrprim_SOR, vecrprim_SOI
!!! for common/potHeHe
  real (kind=8), dimension(ndist22) :: V22, dV22
!!! for common/pos/
  real (kind=8), dimension(3,nbrat1) :: xyz1
  real (kind=8), dimension(3,nbrat2) :: xyz2
!!! for common/dist/
  real (kind=8), dimension(3,nbrat2,nbrat1) :: distvec12
  real (kind=8), dimension(nbrat2,nbrat1) :: dist12
!!! for common/dist2/
  real (kind=8), dimension(3,nbrat2,nbrat2) :: distvec22
  real (kind=8), dimension(nbrat2,nbrat2) :: dist22
  real (kind=8), dimension(3,ndist22) :: distvec221D
  real (kind=8), dimension(ndist22) :: dist221D
!!! for common/output/
  integer :: ittt
!!! for common/printinge
  real (kind=8) ::  deltatwe
!!! for common cevalhop
  logical :: levalhop ! set to .false. to skip hopping probability evaluation at the 1st step
  ! (passed in common cevalhop to outp subroutine)
!!! for COMMON randomc
    real (kind=8) :: yrand

!!Opening files
  logical lfirst
  data lfirst/.true./
  save lfirst
!  character (len=13) :: title 
!!! general use
  integer ivar, jvar, idist22
  logical   :: linit,linitoutp
  integer   :: itraj
  common/inittraj/linit,itraj,linitoutp

  common/masses/xmass1,xmass2,xmass(nbrat12),xmassinv(ncoord12)

  common/dyn/ihop, ioldhop, ihopinic, nbhop   ! index of the current, previous step, DIM eigenvalue (hence surface for classical coordinates)
  common/nadcplng/ rdkj_SOR, rdkj_SOI

  common/DIM/EigenEcomm,vecr_SOR,vecr_SOI,vecrprim_SOR,vecrprim_SOI
  common/potHeHe/V22
  common/Eclass/Etot,Ekin,Epot,EpotHe
  common/pos/xyz1,xyz2
  common/dist/distvec12,dist12
  common/dist2/distvec22,dist22,distvec221D,dist221D
  common/output/ittt
  common/printinge/deltaTwe
  common/cevalhop/levalhop
!  common/writefile/title
  !!! for random subroutine: TO BE INCLUDED IN ALL PROGRAMS OR SUBROUTINES CALLING RANDOM
  common/randomc/yrand


!  print *
!  print *,' OUTP: t = ',t,' au',' = ',t*convt,' ps',' ihop = ',ihop
!  call flush(6)
!  print *,' current PES : ihop = ',ihop
!  print *
!  print *,' current variables : Y'
!  print *,' positions '
!  print *,(y(ivar),ivar=1,ncoord12)
!  print *,' momenta '
!  print *,(y(ivar),ivar=ncoord12+1,ncoordmom)
!  print *,' coefficients modulus '
!  print *,(y(ivar),ivar=ncoordmom+1,ncoordmomcoef)
!  print *,' coefficients phases '
!  print *,(y(ivar),ivar=ncoordmomcoef+1,nvar)

!  if (lfirst) then
!    open(59,file='movie.'//title)   deleted
!   open(31,file=title//'.coords') opened in main
!   open(32,file=title//'.moments') opened in main
!   open(33,file=title//'.coefts') opened in main
!   open(30,file=title//'.energies') opened in main
!   open(41,file=title//'.hopprob_allow') ! already opened in MAIN; written in evalhop
!   lfirst=.false.
   

!   open(51,file=title//'.couplings')   ! opened in main
!   open(50,file=title//'.eigenvect')   ! opened in main
  !  end if
  

!!!
  call distN(y,ndimy)

  call  DIMpSOdynanZPAD(eigenE,ndimorb,iunit)
  EigenEcomm = EigenE
  
  EpotHe = 0.d0
  if (nbrat2.gt.1) then
     call derpotHeHeZPADu(dV22,V22,dist221D,ndist22,iunit)
     EpotHe = sum(V22(1:ndist22))
  endif
  
  Epot = EpotHe + EigenEcomm(ihop)
  ekin = 0.5d0*sum(xmassinv(1:ncoord12)*y(ncoord12+1:ncoordmom)**2)
  Etot = Epot + Ekin
  ! tests

!!! check for degeneracies:
!!! Assuming a real crossing only occurs punctually,
!!!  there is very little chance to stumble on one.
!!! Hence integration of the current trajectory is stopped
!!! in case 2 pairs have the same energy (case where Ba+ is isolated, no more He around)
  do ivar = 1, ndimorb-1,2
     do jvar = ivar+2, ndimorb-1,2
        if (dabs(EigenEcomm(ivar)-EigenEcomm(jvar)).lt.1.d-11) then
           prmt(5) = 2.d0
           print *,' itraj, t = ',itraj,t
           print *,'  >>>>> DEGENERACY between ',ivar,' and ',jvar,' eigenvectors detected'
           print *,'  ie, between Kramer-s pair ',(ivar+1)/2,' and ',(jvar+1)/2
           print *,' trajectory nb ',itraj,' ended at t = ',t/convt,' ps'
           print *,' Ba+ has no more He around'
        endif
     enddo
  enddo

  return
end SUBROUTINE OUTPDBIS
