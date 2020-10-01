!!! PROGRAM FOR ZPAD-MDQT (quasi-CLASSI!AL, with Tully's surface hopping) TRAJECTORIES for HeN-Ba+
!!! (MDQT in the adiabatic basis set)
   
!!! This version for propagating Re(dk), Im(dk) for the coefficients rather than mod(ck), phase(ck)
!!! where dk = ck exp(i int_0^t Vk/hbar dt'), Vk being the k-th eigenvalue attached to eigenvector k
!!! This form should avoid describing many unnecessary oscillations compared to Re(ck), Im(ck),
!!! while avoiding numerical undeterminacy when a coefficient tends to zero with mod(ck), phase(ck)
!!! Added energy phases as variables: phik=int_0^t Vk/hbar dt'

  implicit none

  include 'nbrat12d.f90'
  include 'units.f90'
!     Y contains all the variables that are going to be propagated
!     and DERY the corresponding equations of motion
!   such that dY(i)/dt = DERY(i)      
!     These equations of motion are evaluated in subroutine ffct (ground state) or fctd (excited state)
!     which is called by the propagator.
!     Here we use fhpcg (ground) and hpcg (excited state) Hamming's predictor-corrector propagator
!!! Order of the variables: 
!!!  - all atomic (atom1 then atom2) cartesian coordinates;
!!!  - then all atomic (atom1 then atom2) momenta;
!!!  - then coefficients: Re(d1), Im(d1); Re(d2), Im(d2); .... for the electronic wave packet coefficients
!!!    in the adiabatic basis set
!!!  - then the energy phases phik

  logical :: lterminate, lfail
  integer :: ntraje
  integer :: i, iat, ipos
  integer :: itraj, itrajoutp
  integer :: limpa, initpcg, ihlf
  integer, parameter :: nprmte=7 !!! minimum is 5 for HPCG, but it is set to more for output handling

  real (kind=8) :: rhardsphBa,rhardsphHe

!!! for common /dyn/
  integer :: ihop, ioldhop, ihopinic, nbhop
  integer :: iopt
  real (kind=8), dimension(nvar) :: ye, derye
  real (kind=8), dimension(ncoordmom) :: yg,deryg
  real (kind=8) :: Emin, temp, etemp, etot
  real (kind=8) :: t0,tg,te,tstep0g,tstep0e,accurg,accure,tfgps,tfeps,dtgps,dtg,deltatwe,deltatwg
  real (kind=8), dimension(5) :: prmtg
  real (kind=8), dimension(nprmte) :: prmte
  real (kind=8) :: xmass1, xmass2, xmass, xmassinv
  real (kind=8) :: hyperad0, hyperad

!!! for COMMON randomc0, randomc
  real(kind=8) :: yrand0, yrand
  real (kind=8), allocatable ::  yrand0traj(:)

!!! for common cevalhop
  logical :: levalhop ! set to .false. to skip hopping probability evaluation at the 1st step
  ! (passed in common cevalhop to outp subroutine)

!!! for writing ground state energies for the initially sampled conditions
  real (kind=8) :: etotg,eking,epotg

!!! for common Eclass
  real (kind=8) :: Etote,Ekine,Epote,EpotHe

  common/masses/xmass1,xmass2,xmass(nbrat12),xmassinv(ncoord12)
  common/printingg/deltatwg
  common/printinge/deltatwe
  common/output/itrajoutp
!!!  ! index of the current, previous, and initial DIM eigenvalue (hence surface for classical coordinates)
  common/dyn/ihop, ioldhop, ihopinic, nbhop
  common/cevalhop/levalhop
!!! index set to .T. when starting a new trajectory (for reinitializing adiabatic following)
  logical :: linit,linitoutp
  common/inittraj/linit, itraj, linitoutp
!!! for random subroutine: TO BE INCLUDED IN ALL PROGRAMS OR SUBROUTINES CALLING RANDOM
  common/randomc0/yrand0
  common/randomc/yrand

  common/Eclass/Etote,Ekine,Epote,EpotHe

  character*255 :: filenameresults
  common/writefile/filenameresults ! read from Bap.in (in subroutine config); all output files are then opened in MAIN

!!! create directory 'result' before running if it does not already exist
!!  !NH  read(5,'(a)') title
!!  title='result/0010traj100He'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!C


!!!   INITIALIZATION:  read input configuration
!!! and add internal energy corresponding to TEMP (read from Bap.in)
  call configdwf2(ntraje,iopt,yg,ncoordmom,Emin,Temp,ETemp,Etot &
              ,tstep0g,tstep0e,accurg,accure,tfgps,tfeps,dtgps,rhardsphBa,rhardsphHe)

!!! shuffles initial randomnumber values: 1 for each trajectory
  allocate (yrand0traj(ntraje))
  yrand = yrand0
  do itraj = 1, ntraje
     call RANDOM(yrand,yrand0traj(itraj))
  enddo
!*************************************************************************
! Open output files 
! (except Bap.in, input file read in config subroutine)
! and Bap.in qui est a part (lu une seule fois : config.f)
!  OPEN(60,FILE=title//'.BaHexyz.t0') ! file for xbs plot of configuration at t=0 for each traj
  OPEN(42,FILE=trim(filenameresults)//'.tf_info',STATUS='UNKNOWN')! fate of trajectory, number of hops, etc...
   write(42,*)' # itraj  lterminate   t(au)    ihopini   ihop   nbhop, yrand0traj'
!  OPEN(61,FILE=title//'.BaHexyz.tf',status='UNKNOWN') ! file for xbs plot of configuration at end of each traj

!!! files written in evalhop
  OPEN(40,FILE=trim(filenameresults)//'.evalhop') ! write in evalhop hopping probabilities at each time step
  OPEN(41,FILE=trim(filenameresults)//'.hopprob_allowed') ! write in evalhop each allowed hop
   write(40,*)' # itraj   t(au)    nbhop   (probhop(i),i=1,nDIMorb)'
   write(41,*)' # itraj   t(au)    nbhop    ioldhop   ihop   randhop(1) (probhop(i),i=1,nDIMorb)'

!!! files written in foutp:
  open(20,file=trim(filenameresults)//'.fenerg',status='unknown') ! t, pot, kin, and total ground state energies 
  OPEN(21,FILE=trim(filenameresults)//'.fcoords',STATUS='UNKNOWN')
  OPEN(22,FILE=trim(filenameresults)//'.fmoments',STATUS='UNKNOWN')
  write(20,*)'#   t(au)   etot    ekin   epot '
  write(21,*)'#   t(au)   xBa, yBa, zBa, xHe1, yHe1, zHe1, xHe2, yHe2, zHe2, ... (au)'
  write(22,*)'#   t(au)   pxBa, pyBa, pzBa, pxHe1, pyHe1, pzHe1, pxHe2, pyHe2, pzHe2, ... (moments, in au)'
  
!!! files written in outp:
   open(30,file=trim(filenameresults)//'.energies') ! ex-61 
   open(31,file=trim(filenameresults)//'.coords') ! ex 60 (contained also moments)
   open(32,file=trim(filenameresults)//'.moments') 
   open(33,file=trim(filenameresults)//'.dcoefts-eph')  ! coefficients (modulus and phase)
!   open(34,file=trim(filenameresults)//'.doefts-modph')  ! coefficients (modulus and phase)
!   open(35,file=trim(filenameresults)//'.coefts-modph')  ! coefficients (modulus and phase)
   open(50,file=trim(filenameresults)//'.eigenvect') ! ex 70
!   open(51,file=trim(filenameresults)//'.couplings')  ! ex 90
   open(52,file=trim(filenameresults)//'.couplingsReIm')  ! ex 90
  write(30,*)'#   itraj  t(au)   etot    ekin   epot   eigenecomm(1:6)   epotHe   (au)'
  write(31,*)'#   itraj  t(au)   xBa, yBa, zBa, xHe1, yHe1, zHe1, xHe2, yHe2, zHe2, ... (au)'
  write(32,*)'#   itraj  t(au)   pxBa, pyBa, pzBa, pxHe1, pyHe1, pzHe1, pxHe2, pyHe2, pzHe2, ... (moments, in au)'
  write(33,*)'#   itraj  t(au)   dcoefts (Re(1), Im(1), Re(2), Im(2),...,Ephase(1), Ephase(2)...)'
  write(50,*)'#   itraj  t(au)   eigvecRe(1,1),eigvecIm(1,1)  eigvecRe(2,1),eigvecIm(2,1) ...'
  write(52,*)'#   itraj  t(au)   rdkjRe(1,1),rdkjIm(1,1)  rdkjRe(2,1),rdkjIm(2,1) ...'

!!! ground state configuration after all trajectories
     open(92,file=trim(filenameresults)//'.finalconfig') ! 

!!! ground state energy for sampled configuration at the beginning of each excited state trajectory
     open(93,file=trim(filenameresults)//'.init-conds') ! itraj, yrand0traj(itraj) , Etotg, Epotg, Eking,Etote,Ekine,Epote,EpotHe, y


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dtg = dtgps/convt

!   calculation of initial hyperspherical radius 
!     as a benchmark to monitor dissociation in the ground state
!  Ba+ is the 1st atom, with coordinates (yg(1),yg(2),yg(3))
  hyperad0 = 0.d0
  do i = ncoord1+1, ncoord12-2, 3
     hyperad0 = hyperad0 &
         +  (yg(i)-yg(1))**2+(yg(i+1)-yg(2))**2+(yg(i+2)-yg(3))**2
  enddo
  print *
  print *,' initial value for hyperradius = ',dsqrt(hyperad0/nbrat2)*bohr, ' Angstrom'
  call flush(6)

!     initialization for propagation in the ground state, randomization step:
  t0 = 0.d0
  prmtg(1) = t0/convt        ! initial value of time (in au)
  prmtg(2) = tfgps/convt   ! final value of time (in au)
  prmtg(3) = tstep0g   ! initial value of the time step
  prmtg(4) = accurg   ! accuracy parameter for hpcg
  prmtg(5) = 0.d0   ! parameter to be set >0 for terminating propagation (in foutp)

  deryg(1:ncoordmom) = 1.d0/dfloat(ncoordmom) ! error weights, will be overwritten by time derivatives
  tg = t0
  call fhpcg(tg,yg,deryg,ncoordmom,prmtg)

      
!   calculation of hyperspherical radius 
!  as a benchmark to monitor dissociation after propagation in the ground state
  hyperad = 0.d0
  do i = ncoord1+1, ncoord12-2, 3
     hyperad = hyperad + (yg(i)-yg(1))**2+(yg(i+1)-yg(2))**2+(yg(i+2)-yg(3))**2
  enddo
  print *
  print *,' hyperradius after randomization = ',dsqrt(hyperad/nbrat2)*bohr,' Angstrom'
  print *
  print *,' values for the variables, t = ',tg,' au, = ',tg*convt,' ps'
  call flush(6)
!  print *,' positions:'
!  print *, yg(1:ncoord12)
!  print *,' momenta : '
!  print *, yg(ncoord12+1:ncoordmom)
 

!!!!!!!!!!
!!!   MAIN DO LOOP on EXCITED STATE TRAJECTORIES
!!!!!!!!!!
!!! choice of the potential energy surface: to be sampled later!!!
!!!  ihop is now read in config.f 
!!!ihop = 1
!!!!!!!!!!!
  
  do itraj = 1, ntraje
!!!   run ground state trajectory for dtgps ps (dtg au) to get new initial conditions
     !     initialization for propagation in the ground state, randomization step:
     print *
     print *,'  trajectory nb ', itraj
     print *,'  --------------'
     print *,' randomization step '
     call flush(6)
     t0 = tg
     prmtg(1) = t0          ! initial value of time  (in au)
     prmtg(2) = t0+dtg      ! final value of time (in au)
     prmtg(3) = tstep0g     ! initial value of the time step
     prmtg(4) = accurg      ! accuracy parameter for hpcg
     prmtg(5) = 0.d0        ! parameter to be set >0 for terminating propagation (in foutp)
     print *,' parameters for ground state randomization step : '
     print *, prmtg(1:5)
     deryg(1:ncoordmom) = 1.d0/dfloat(ncoordmom) ! error weights, will be overwritten by time derivatives
     call flush(6)
!     randomize by running dtg ps in the ground state before selecting new initial conditions
!     for next trajectory
     call fhpcg(tg,yg,deryg,ncoordmom,prmtg)  ! runs ground state dynamics for tg from prmtg(1) to prmtg(2)
     print *
     print *,' values for the variables, after randomization t = ',tg,' au, = ',tg*convt,' ps'
     call flush(6)
     hyperad = 0.d0
     do i = ncoord1+1, ncoord12-2, 3
        hyperad = hyperad + (yg(i)-yg(1))**2+(yg(i+1)-yg(2))**2+(yg(i+2)-yg(3))**2
     enddo
     print *
     print *,' hyperradius after randomization = ',dsqrt(hyperad/nbrat2)*bohr,' Angstrom'
     call flush(6)
!!  print *,' positions:'
!!  print *, yg(1:ncoord12)
!!  print *,' momenta : '
!!  print *, yg(ncoord12+1:ncoordmom)
!!  print *
!     write(60,*) 'atom Ba   ',yg(1:3)*bohr
!     do iat=1, nbrat2
!        write(60,*) 'atom He   ',yg(3*iat+1:3*iat+3)*bohr
!     enddo
!     write(60,*)
!     write(60,*)
!     write(60,*)' spec  Ba   1.0    0. 1.0 1.0'
!     write(60,*)' spec  He   0.5    1.0 1.0 1.0'
!     write(60,*)
!     write(60,*)' bonds Ba Ba  5.0 7. 0.15 0.0'
!     write(60,*)' bonds He He  5.0 7. 0.05 0.0'
!     write(60,*)' bonds Ba He  4.0 7. 0.1 0.0'
!     write(60,*)
!     write(60,*)' inc  5.00'
!     write(60,*)
!     write(60,*) 'tg = ',tg*convt
!     write(60,*)
!     call flush(60)

     write(40,*)
     write(40,*)
     
     print *,'  initiate trajectory nb ', itraj,' in the excited state'
     call flush(6)
!     use yg as initial conditions for a new trajectory in the excited states:
     te = 0.d0
     prmte(1) = te             ! initial value of time  (in au)
     prmte(2) = tfeps/convt    ! final value of time (in au)
     prmte(3) = tstep0e        ! initial value of the time step (in au)
     prmte(4) = accure         ! accuracy parameter for hpcg
     prmte(5) = 0.d0  ! return flag: HPCG returns if OUTP set it to non-zero
     prmte(6) = 0.d0          ! flag set to >0 in EVALHOP if there has been a hop (HPCG needs to be restarted)
     print *,' parameters for excited state trajectory propagation for hpcg:'
     print *, prmte(1:6)
!!!     yrand0traj = yrand  ! for bookkeeping, in case a specific trajectory needs to be rerun
     yrand = yrand0traj(itraj)  ! for reshuffling initial conditions: reinitialize sapling/each traj
     call flush(6)
     do i = 1, ncoordmom
        ye(i) = yg(i)
     enddo
     lterminate=.false.
!!! check angular momentum before sampling:
     print *
     print *,'check angular momentum before sampling:'
     call comJ(ye(1:ncoord12),ye(ncoord12+1:ncoordmom),ncoord12)
     print *
     print *,' sampling initial He position using ZPAD wf '
!!! sample He atom positions around their classical one using ZPAD wave function
     call samplewf(ye(1:ncoord12),ncoord12,rhardsphBa,rhardsphHe,lfail)
     if (lfail) then
        print *,' itraj = ',itraj,' SAMPLING for INITIAL CONDITIONS FAILED, skipping to next traj'
        write(42,*) itraj, lterminate, te, ihopinic,ihop,nbhop, yrand0traj(itraj)
     call flush(42)
        cycle
     endif
!!! the origin will no longer be exactly at the CoM
!!! and the angular momentum no longer exactly equal to zero
!!! check for this!!!
     call comJ(ye(1:ncoord12),ye(ncoord12+1:ncoordmom),ncoord12)
!     print *,' initial coordinates (au) : '
!     print *,ye(1:ncoord12)
!     call flush(6)
!     print *
!     print *,' initial momenta (au) : '
!     print *, ye(ncoord12+1:ncoordmom)
!     call flush(6)
!     ! select initial excited state for classical nuclear coordinates
!     ! (ihop, passed by COMMON/HOP/)
!     ! and initializes wave packet components accordingly
     call WPINI_SOd(ye,nvar,iopt)
     print *
     print *,' initial coefficients (real, imaginary part for d1, then for d2,...)'
     print *, ye(ncoordmom+1:ncoordmomcoef)
     call flush(6)
!     print *
     print *,' initial energy phases (should be zero)'
     print *, ye(ncoordmomcoef+1:nvar)
     call flush(6)
     
!     update itrajoutp for outp subroutine
     itrajoutp = itraj
     call foutpbis(te,ye,nvar,Etotg,Eking,Epotg)  !!! to record ground state energy corresponding to sampled initial conditions
     call outpdbis(te,ye,derye,nvar,prmte,nprmte) !!! to record excited state energies corresponding to initial conditions
     write(93,*) itraj,yrand0traj(itraj),Etotg,Eking,Epotg,Etote,Ekine,Epote,EpotHe,ye

     linit = .true.
     linitoutp = .true.
     nbhop = 0   ! reinitialize the number of hops for this trajectory
     initpcg = 0  ! for hpcg: number of times the current trajectory has been reinitialized
     limpa = 0  ! number of tentatives to half the step in addition to the 10 ones initially in hpcg
     ! ihlf = number of times the integration step has been halved (output): limit = 10+limpa
     do while (te.lt.prmte(2).and..not.lterminate.and.dabs(prmte(5)).lt.1.d-10)
        levalhop=.false.
!!! need to initialize DERY before calling HPCG!!! to error weights, with sum = 1

        derye(1:nvar) = 1.d0/dfloat(nvar)
        call hpcg(te,ye,derye,nvar,prmte,nprmte,limpa,initpcg,ihlf)
        if (ihlf.ge.10+limpa) then
           print *,' trajectory nb ',itraj,' stopped, step halved too many times: ',ihlf
           call flush(6)
           lterminate=.true.
           cycle
        endif
        if (prmte(6).gt.0.d0) then   !!! case of a hop
           print *,'te = ',te,'au, = ',te*convt,'ps:   new hop on ',ihop
           prmte(1) = te   ! reinitialize integrator but not time
           prmte(5) = 0.d0    ! reset flag for return 
           prmte(6) = 0.d0    ! reset flag signaling a hop
           call flush(6)
        endif
        !        if prmte(5) = 2., the trajectory is stopped
        !         (Degeneracy of 2 Kramer's pairs hence Ba+ has no more He)
           
     end do
     print *,' end of trajectory nb ',itraj,' t = ',te,' au, = ',te*convt,' ps'
     call flush(6)
     write(42,*) itraj, lterminate, te, ihopinic,ihop,nbhop, yrand0traj(itraj)
     call flush(42)
!     write(61,*) 'atom Ba   ',ye(1:3)*bohr
!     do iat=1, nbrat2
!        write(61,*) 'atom He   ',ye(3*iat+1:3*iat+3)*bohr
!     enddo
!     write(61,*)
!     write(61,*)
!     write(61,*)' spec  Ba   1.0    0. 1.0 1.0'
!     write(61,*)' spec  He   0.5    1.0 1.0 1.0'
!     write(61,*)
!     write(61,*)' bonds Ba Ba  5.0 7. 0.15 0.0'
!     write(61,*)' bonds He He  5.0 7. 0.05 0.0'
!     write(61,*)' bonds Ba He  4.0 9. 0.1 0.0'
!     write(61,*)
!     write(61,*)' inc  5.00'
!     write(61,*)
!     write(61,*) 'te = ',te*convt
!     write(61,*)
!     call flush(61)
         
  enddo
!!! output final ground state configuration as an input file for next runs
  write(92,*) epotg*Hartree2EJ, eking*Hartree2EJ
  write(92,*) (yg(i)*Ang2m,i=1,3), (yg(ncoord12+i)/(xmass1*xmsinv2au),i=1,3)
  do iat=1, nbrat2
     ipos = 3*(nbrat1+iat-1)
     write(92,*) (yg(ipos+i)*Ang2m,i=1,3), (yg(ncoord12+ipos+i)/(xmass1*xmsinv2au),i=1,3)
  enddo
  close(92)      
  STOP

END program






