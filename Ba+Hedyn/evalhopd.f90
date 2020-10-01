!!!  SUBROUTINE TO EVALUATE the PROBABILITY for SURFACE HOPPING
!!!  returns IHOP = surface index for the classical trajectory
!!!
!!! This version for propagation in Re(dk), Im(dk) for the coefficients
!!! rather than mod(ck), phase(ck)
!!! where dk = ck exp(i Ephase_k) where Ephase_k=int_0^t (Vk/hbar dt'),
!!!            Vk being the k-th eigenvalue attached to eigenvector k
!!! This form should avoid describing many unnecessary oscillations compared to Re(ck), Im(ck),
!!! while avoiding numerical undeterminacy when a coefficient tends to zero with mod(ck), phase(ck)
!!! Added energy phases as variables: phik=int_0^t Vk/hbar dt'


subroutine evalhopd(t,y,dery,ndimy,ittt,prmt,nprmt)

  implicit none

  include 'nbrat12d.f90'
  include 'units.f90'


  integer :: ndim, ndimy, ittt, nprmt
  real (kind=8) :: t
  real (kind=8), dimension(ndimy) :: y, dery
  real (kind=8), dimension(nprmt) :: prmt
!!! for COMMON randomc
  real (kind=8) :: yrand  
!! for common /DIM/: eigenvalue; eigenvectors (real, imaginary part); eigenvectors at previous time step
  real (kind=8), dimension(ndimorb) ::  eigenEcomm
  real (kind=8), dimension(ndimorb,ndimorb) ::  vecr_SOR, vecr_SOI
  real (kind=8), dimension(ndimorb,ndimorb) ::  vecrprim_SOR, vecrprim_SOI
  common/DIM/eigenEcomm,vecr_SOR,vecr_SOI,vecrprim_SOR,vecrprim_SOI

!!! for common /nadcplng/: real and imaginary part of the non-adiabatic couplings (rd_{k,j})
  real (kind=8), dimension (ndimorb,ndimorb) :: rdkj_SOR, rdkj_SOI
  common/nadcplng/ rdkj_SOR, rdkj_SOI

!!! for common/masses/ : Ba+ mass, He mass, mass for each atomic coordinate, 1/mass for each atomic coordinate
  real (kind=8) :: xmass1, xmass2, xmass, xmassinv
  common/masses/xmass1,xmass2,xmass(nbrat12),xmassinv(ncoord12)
!!! for common/dyn/: current, previous, and initial index for the surface on which classical dynamics is run
  integer :: ihop, ioldhop, ihopinic, nbhop             
  common/dyn/ihop, ioldhop, ihopinic, nbhop   

  logical :: lforbid
  integer :: i, ihoptest
!!! for random sampling 
  integer iout, iseed
  integer, parameter :: nrand=1
  real (kind=8), dimension(nDIMorb) :: probhop
  real (kind=8) :: probhopsum, sumhop
  real (kind=8), dimension(nrand) :: randhop
!!! for random subroutine: TO BE INCLUDED IN ALL PROGRAMS OR SUBROUTINES CALLING RANDOM
  common/randomc/yrand

  real (kind=8) :: ccc
  real (kind=8), dimension(ndimorb) :: Ephase
  double complex, dimension(ndimorb) :: zdcoef
  double complex, dimension(ndimorb,ndimorb) :: zrdkj, za
  
  
! work: initial       integer i,j,k,iout,nbas,ittt,ihpcg,ncoef,nbhop
! work: initial       integer inhopv,oldhop,hop,nphase,ndim,ndimtot
! work: initial      integer moment,momenthe,nbashe,iat,iadj,iadjhe
! work: initial       integer iat,iadj,iadjhe
! work: initial       real*8 g,randhop,dhdyvib,rdkjvib,sumg,suplim,vv
! work: initial       real*8 x2,valp,rdkj,vecr,dum

! work: initial      parameter(nbas = 3*nbrat, moment = 3*nbrat, ncoef = 6*nbrat, 
! work: initial      &            nphase = 9*nbrat, ndim = 12*nbrat, nbashe = 3*nbrathe,
! work: initial      &            momenthe = ndim + nbashe, ndimtot = ndim + 6*nbrathe)
! in units.f90      parameter(convt = 2.418884327d-5)

! work: initial       dimension g(nbas,nbas),randhop(1)
! work: initial       dimension dhdyvib(nbas,nbas),rdkjvib(nbas)
! work: initial       dimension vecr(nbas,nbas)
      
! work: initial       common/adjust/iadj,iadjhe
! work: initial       common/sauts/hop,oldhop 
! work: initial       common/time/x2
! work: initial       common/nbsauts/nbhop
! work: initial       common/couplages/rdkj(nbas,nbas)

!!! probability for a hop (from ihop to ihopnew)
!!! Prob(ihop -> ihopnew) = -2*delta_t*n_j/n_hop cos(q_j-q_hop) (dR/dt).d(j,hop)
  probhop = 0.d0   ! initialize array for hopping probabilities
  probhopsum = 0.d0
!!!  coef= 2.d0*prmt(3)/y(ncoordmom+ihop)

!!! define za = density matrix: za_{jk}=d_j d_k* exp(i(phi_k-phi_j))
  Ephase = y(ncoordmomcoef+1:nvar)
  zdcoef = dcmplx(y(ncoordmom+1:ncoordmomcoef:2),y(ncoordmom+2:ncoordmomcoef:2))
  zrdkj = dcmplx(rdkj_SOR,rdkj_SOI)
  do i=1, nDIMorb
     za(i,:) = zdcoef(i)*dconjg(zdcoef(:))*cdexp(dcmplx(0.d0,(ephase(:)-ephase(i))))
  enddo
  ccc= 2.d0*prmt(3)/cdabs(za(ihop,ihop))
  
  do i = 1, nDIMorb
     if (i.eq.ihop) cycle
!!!     probhop(i)=(y(ncoordmom+i)*(rdkj_SOI(i,ihop)*dsin(y(ncoordmomcoef+i)-y(ncoordmomcoef+ihop)) &
!!!          -rdkj_SOR(i,ihop)*dcos(y(ncoordmomcoef+i)-y(ncoordmomcoef+ihop))))*coef
!!!     probhop(i) = dimag(za(i,ihop)*dconjg(zrdkj(i,ihop)))*ccc
          probhop(i) = -dreal(dconjg(za(i,ihop))*zrdkj(i,ihop))*ccc
!     print *,' i, yi, rdkj_SOI,rdkj_SOR, probhop = ',i,y(ncoordmom+i),rdkj_SOI(i,ihop),rdkj_SOR(i,ihop),probhop(i)
     if (probhop(i).lt.0.d0) probhop(i)=0.d0
  enddo
!  probhop = -probhop*2.d0*prmt(3)/y(ncoordmom+ihop)
  probhopsum = sum(probhop)
  if (probhopsum.gt.1.d0) then
     probhop = probhop/probhopsum
     probhopsum = 1.d0
  endif
  probhop(ihop) = 1.d0-probhopsum
!  print *,' probhop: ',probhop(1:nDIMorb)
  write(40,*) ittt,t,probhop(1:nDIMorb)

!!! TEST
  do i = 1, nDIMorb
!     print *,' i, probhop(i) = ',i,probhop(i)
!       if (probhop(i).lt.0.d0) then
!        print *,' >>>>>> in EVALHOP: TEST ERROR <<<<< '
!        print *,' i = ',i,' probability for hop is negative: ',-probhop(i)*2.d0*prmt(3)/y(ncoordmom+ihop)
!     endif
  enddo
!  write(40,*) t, (probhop(i),i=1,nDIMorb)
  
!!! test the possibility to hop
  iout = 6
!  iseed = 0 !!! to continue sampling from previous call to vranf
!  call vranf(randhop,nrand,iseed,iout)
  call random(yrand,randhop(1))
  sumhop = 0.d0
  i=0
  do while (sumhop.lt.randhop(1))
     i=i+1
     if (i.gt.nDIMorb) then
        print *,' >>>>> ERROR STOP in EVALHOP <<<<<'
        print *,' i = ',i,' = nDIMorb = ',nDIMorb,' and no solution found '
        print *,' randhop = ',randhop(1),' probhopsum = ',probhopsum,' sumhop = ',sumhop
        STOP
     endif
     sumhop = sumhop + probhop(i)
  end do

!!! TEST
!!!  print *,' in evalhop: randhop, ihop, probhop = ',randhop(1),ihop,probhop(ihop)
!!! test if there is a possible hop
  if (i.ne.ihop) then
     ihoptest = i
     print *,' evalhop, t = ',t,'   test hop from ', ihop,' to ',ihoptest
!!! test if the hop is energetically allowed, adjust momenta to conserve energy if it is
     call globveladj(y,dery,ndimy,ihoptest,lforbid)
     if (.not.lforbid) then   !!! hop is classically allowed
        nbhop = nbhop + 1
        ioldhop = ihop
        ihop = ihoptest
        prmt(5) = 1.d0   !!! flag set to >0 to exit from HPCG in order to restart dynamics on new surface
        prmt(6) = 1.d0   !!! flag transmitted to MAIN via HPCG to indicate restart was due to a hop
!        open(41,file='result/probsaut_allow',status='unknown') !opened in main
        write(41,"(I5,G16.8,3(I5),19(G16.8))") &
             ittt,t,nbhop,ioldhop,ihop,randhop(1), (probhop(i),i=1,nDIMorb)
        call flush(41)
     else     !!! hop is classically forbidden, stay on the same surface
        print *, ' not enough energy to hop '
        call flush(6)
     endif
  endif

  return
end subroutine evalhopd

