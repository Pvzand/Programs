SUBROUTINE OUTP(T,Y,DERY,ndimy,prmt,nprmt)
!!! output subroutine for HPCG integrator
!!! in Re(dk), Im(dk) for the coefficients
!!! rather than mod(ck), phase(ck)
!!! where dk = ck exp(i int_0^t Vk/hbar dt'), Vk being the k-th eigenvalue attached to eigenvector k
!!! This form should avoid describing many unnecessary oscillations compared to Re(ck), Im(ck),
!!! while avoiding numerical undeterminacy when a coefficient tends to zero with mod(ck), phase(ck)
!!! Added energy phases as variables: phik=int_0^t Vk/hbar dt'

  implicit none

  integer ndimy, nprmt,iat
  real (kind=8) :: t
  real (kind=8), dimension (ndimy) :: Y, DERY
  real (kind=8), dimension (nprmt) :: prmt

  include 'nbrat12d.f90'
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
  
  if (linitoutp) then  !!! 2 empty lines to separate trajectories in output files
     write(30,*)
     write(30,*)
     write(31,*)
     write(31,*)
     write(32,*)
     write(32,*)
     write(33,*)
     write(33,*)
     write(50,*)
     write(50,*)
!     write(51,*)
!     write(51,*)
     write(52,*)
     write(52,*)
     
  linitoutp=.false.
  end if

!!! outp is always called after a call to fct: hence potential energies are known
  EpotHe = 0.d0
  if (nbrat2.gt.1) EpotHe = sum(V22(1:ndist22))
  Epot = EpotHe + EigenEcomm(ihop)
  ekin = 0.5d0*sum(xmassinv(1:ncoord12)*y(ncoord12+1:ncoordmom)**2)
  Etot = Epot + Ekin
  ! tests
  if (dabs(t/deltaTwe-anint(t/deltaTwe)) < prmt(3)/(2.d0*deltaTwe) ) then
     print *,' OUTP: t = ',t,' au',' = ',t*convt,' ps',' ihop = ',ihop
     call flush(6)
!
 !    write(59,*) nbrat12
 !    write(59,*) 'Ba  ', y(1:3)*bohr
 !    do iat=1,nbrat2
 !    write(59,*) 'He  ', y(3*iat+1:3*iat+3)*bohr
     !    enddo
     write(30,*) ittt, t, Etot, Ekin, Epot, EigenEcomm(1:6), EpotHe
     write(31,*) ittt, t, (y(ivar),ivar=1,ncoord12)
     write(32,*) ittt, t, (y(ivar),ivar=ncoord12+1,ncoordmom)
     write(33,*) ittt, t, (y(ivar),ivar=ncoordmom+1,nvar)
!     write(34,*) t, t*convt, cdabs(dcmplx(y(ncoordmom+1:ncoordmomcoef:2),y(ncoordmom+2:ncoordmomcoef:2))),&
!     datan2(y(ncoordmom+2:ncoordmomcoef:2),y(ncoordmom+1:ncoordmomcoef:2))
!     write(35,*) t, t*convt, cdabs(dcmplx(y(ncoordmom+1:ncoordmomcoef:2),y(ncoordmom+2:ncoordmomcoef:2))),&
!     datan2(y(ncoordmom+2:ncoordmomcoef:2),y(ncoordmom+1:ncoordmomcoef:2))+y(ncoordmomcoef+1:nvar)
     write(50,*) ittt, t,(&
           (vecr_SOR(ivar,jvar),vecr_SOI(ivar,jvar),ivar=1,ndimorb),jvar=1,nDIMorb)
!        abs(dcmplx(vecr_SOR(1:nDIMorb,ivar),vecr_SOI(1:nDIMorb,ivar))),ivar=1,ndimorb)
!        write(51,*) t,(&
!                     abs(dcmplx(rdkj_SOR(1:nDIMorb,ivar),rdkj_SOI(1:nDIMorb,ivar))),ivar=1,ndimorb)
     write(52,*) ittt, t,(&
                     (rdkj_SOR(ivar,jvar),rdkj_SOI(ivar,jvar),ivar=1,ndimorb),jvar=1,nDIMorb)
!call flush(59)
     call flush(30)  
     call flush(31)
     call flush(32)
     call flush(33)
     call flush(50)
! call flush(51)
     call flush(52)
! tests
  end if

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

!!! write components of the eigenvectors
  

!!! Test for possible hops
!!! If a hop occurs and it is allowed, prmt(5) is set to non zero which makes HPCG return
!!! and prmt(6) is set to 1.0 to let the MAIN program know that it is because of a hop
  if (levalhop) then   !!! calculate hopping probability and prepares hops
     call evalhopd(t,y,dery,ndimy,ittt,prmt,nprmt)
  else  !!! do not attempt hops at initialization step but resets levalhop for next step
     levalhop=.true.
  endif

  return
end SUBROUTINE OUTP
