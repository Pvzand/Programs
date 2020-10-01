SUBROUTINE OUTP(T,Y,DERY,ndimy,prmt,nprmt)

  implicit none

  integer ndimy, nprmt
  real (kind=8) :: t
  real (kind=8), dimension (ndimy) :: Y, DERY
  real (kind=8), dimension (nprmt) :: prmt

  include 'nbrat12.f90'
  include 'units.f90'

  real (kind=8) :: Etot, Ekin, Epot, EpotHe

!!! no need to re-calculate DIM, outp is always called right after a call to FCT
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
  real (kind=8), dimension(ndimorb,ndimorb) ::  vecr_SOR, vecr_SOI
  real (kind=8), dimension(ndimorb,ndimorb) ::  vecrprim_SOR, vecrprim_SOI
!!! for common/potHeHe
  real (kind=8), dimension(ndist22) :: V22
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
  real (kind=8) ::  deltatwe,deltatwg
!!! for common cevalhop
  logical :: levalhop ! set to .false. to skip hopping probability evaluation at the 1st step
  ! (passed in common cevalhop to outp subroutine)

!!! general use
  integer ivar, jvar, idist22
  logical   :: linit
  integer   :: itraj
  common/inittraj/linit,itraj

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


!  print *
!  print *,' OUTP: t = ',t,' au',' = ',t*convt,' ps'
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
  if (dabs(t/deltaTwe-anint(t/deltaTwe))< prmt(3)/(2.d0*deltaTwe) ) then
    write(60,*) t, (y(ivar),ivar=1,3+3*nbrat2)
  endif
  call flush(60)

!!! outp is always called after a call to fct: hence potential energies are known
  EpotHe = 0.d0
  if (nbrat2.gt.1) EpotHe = sum(V22(1:ndist22))
  Epot = EpotHe + EigenEcomm(ihop)
  ekin = 0.5d0*sum(xmassinv(1:ncoord12)*y(ncoord12+1:ncoordmom)**2)
  Etot = Epot + Ekin
  write(61,*) t, Etot, Ekin, Epot, EigenEcomm(ihop), EpotHe
  call flush(61)
!  write(62,*) t, EigenEcomm(1:nDIMorb)
  write(62,*) t, (EigenEcomm(ivar),ivar=1,ndimorb)
  call flush(62)

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


  prmt(5)=0.d0
!!! write components of the eigenvectors
  
!  do ivar = 1, nDIMorb
!     write(70+ivar,*) t, abs(dcmplx(vecr_SOR(1:nDIMorb,ivar),vecr_SOI(1:nDIMorb,ivar)))
!     write(90+ivar,*) t, abs(dcmplx(rdkj_SOR(1:nDIMorb,ivar),rdkj_SOI(1:nDIMorb,ivar)))
!  enddo

!!! Test for possible hops
!!! If a hop occurs and it is allowed, prmt(5) is set to non zero which makes HPCG return
!!! and prmt(6) is set to 1.0 to let the MAIN program know that it is because of a hop
  if (levalhop) then   !!! calculate hopping probability and prepares hops
     call evalhop(t,y,dery,ndimy,ittt,prmt,nprmt)
  else  !!! do not attempt hops at initialization step but resets levalhop for next step
     levalhop=.true.
  endif

  return
end SUBROUTINE OUTP
