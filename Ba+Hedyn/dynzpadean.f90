!!! PROGRAM FOR ZPAD-MDQT (quasi-CLASSICAL, with Tully's surface hopping) TRAJECTORIES for HeN-Ba+
!!! (MDQT in the adiabatic basis set)
  implicit none

!N      include 'PARAMETER.f'
!N      INCLUDE 'COMMON.f'

  include 'nbrat12.f90'
  include 'units.f90'
!     Y contains all the variables that are going to be propagated
!     and DERY the corresponding equations of motion
!   such that dY(i)/dt = DERY(i)      
!     These equations of motion are evaluated in subroutine ffct (ground state) or fct (excited state)
!     which is called by the propagator.
!     Here we use fhpcg (ground) and hpcg (excited state) Hamming's predictor-corrector propagator
!N: changed order of the variables compared to David's version
!     In David's version y contained atom1 coordinates in the first 3*nbrat positions,
!     then atom1 momenta, then coefficients and phases for the electronic wave packet
!     then after ndim all the atom2 coordinates and then their momenta
!N:  Changed to all atomic (atom1 then atom2) cartesian coordinates
!N:  then all atomic (atom1 then atom2) momenta
!N:  then coefficients modulus and then phases for the electronic wave packet     

  logical :: lterminate
  integer :: ntraje, irand
  integer :: i, iat
  integer :: itraj, itrajoutp
  integer :: limpa, initpcg,ihlf
  integer, parameter :: nprmte=7 !!! minimum is 5 for HPCG, but it is set to more for output handling
!!! for common /dyn/
  integer :: ihop, ioldhop, ihopinic, nbhop
  integer :: iopt
  real (kind=8), dimension(nvar) :: yg,deryg, ye, derye
  real (kind=8) :: Emin, temp, etemp, etot
  real (kind=8) ::t0,tg,te,tstep0g,tstep0e,accurg,accure,tfgps,tfeps,dtgps,dtg,deltaTwe,deltaTwg
  real (kind=8), dimension(5) :: prmtg
  real (kind=8), dimension(nprmte) :: prmte
  real (kind=8) :: xmass1, xmass2, xmass, xmassinv
  real (kind=8) :: hyperad0, hyperad
!!! for common cevalhop
  logical :: levalhop ! set to .false. to skip hopping probability evaluation at the 1st step
  ! (passed in common cevalhop to outp subroutine)
      
  common/masses/xmass1,xmass2,xmass(nbrat12),xmassinv(ncoord12)
  common/printingg/deltaTwg
  common/printinge/deltaTwe
  common/output/itrajoutp
!!!  ! index of the current, previous, and initial DIM eigenvalue (hence surface for classical coordinates)
  common/dyn/ihop, ioldhop, ihopinic, nbhop
  common/cevalhop/levalhop
  
  character (len=7) :: title
  common/writefile/title

  read(5,'(a)') title

!*************************************************************************
! Open all output files 
! except  eigvec, rd, traj???? et ftraj???? (depend de NBAS)
! and Bap.in qui est a part (lu une seule fois : config.f)
! = tests pour Xe3+ 40+IFILE, 50+IBAS, 60 et 61 + Xe.in -> 90
  OPEN(20,FILE=title//'.initraj',STATUS='UNKNOWN')
!      OPEN(21,FILE='result/ihlftest',STATUS='UNKNOWN')
  OPEN(22,FILE=title//'.traj',STATUS='UNKNOWN')
  OPEN(23,FILE=title//'.dtraj',STATUS='UNKNOWN')
!!! moved to foutp  OPEN(24,FILE='result/ftraj',STATUS='UNKNOWN')
  OPEN(25,FILE=title//'.hoptraj',STATUS='UNKNOWN')
  OPEN(26,FILE=title//'.hoptrajlife',STATUS='UNKNOWN')
!      OPEN(27,FILE='result/hoptrajlife_sigma',STATUS='UNKNOWN')
!      OPEN(28,FILE='result/hoptrajlife_pi',STATUS='UNKNOWN')
  OPEN(29,FILE=title//'.momang',STATUS='UNKNOWN')
!      OPEN(30,FILE='result/etot',STATUS='UNKNOWN')
      OPEN(31,FILE=title//'.probsaut',STATUS='UNKNOWN')
  OPEN(32,FILE=title//'.probsaut_allow',STATUS='UNKNOWN')
  OPEN(33,FILE=title//'.hoptrajion',STATUS='UNKNOWN')
  OPEN(34,FILE=title//'.hoptrajcan',STATUS='UNKNOWN')
  OPEN(35,FILE=title//'.etotfrag',STATUS='UNKNOWN')
!      OPEN(36,FILE='result/coords',STATUS='UNKNOWN')
!      OPEN(37,FILE='result/eigen',STATUS='UNKNOWN')
!      OPEN(38,FILE='result/coeff',STATUS='UNKNOWN')
!      OPEN(39,FILE='result/phases',STATUS='UNKNOWN')
  OPEN(62,FILE=title//'.nano',STATUS='UNKNOWN')
  OPEN(110,FILE=title//'.saut_fin',STATUS='UNKNOWN')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!C


!!!   INITIALIZATION:  read input configuration
!!! and add internal energy corresponding to TEMP (read from Bap.in)
  call config(ntraje,irand,iopt,yg,deryg,Emin,Temp,ETemp,Etot &
              ,tstep0g,tstep0e,accurg,accure,tfgps,tfeps,dtgps)
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
  print *,' initial value for hyperradius = ',dsqrt(hyperad0)*bohr, ' Angstrom'

!     initialization for propagation in the ground state, randomization step:
  t0 = 0.d0
  prmtg(1) = t0/convt        ! initial value of time (in au)
  prmtg(2) = tfgps/convt   ! final value of time (in au)
  prmtg(3) = tstep0g   ! initial value of the time step
  prmtg(4) = accurg   ! accuracy parameter for hpcg
  prmtg(5) = 0.d0   ! parameter to be set >0 for terminating propagation (in foutp)
      
  tg = t0
  call fhpcg(tg,yg,deryg,nvar,prmtg)

      
!   calculation of hyperspherical radius 
!  as a benchmark to monitor dissociation after propagation in the ground state
  hyperad = 0.d0
  do i = ncoord1+1, ncoord12-2, 3
     hyperad = hyperad + (yg(i)-yg(1))**2+(yg(i+1)-yg(2))**2+(yg(i+2)-yg(3))**2
  enddo
  print *
  print *,' hyperradius after randomization = ',dsqrt(hyperad)*bohr,' Angstrom'
  print *
  print *,' values for the variables, t = ',tg,tg*convt
  print *,' positions:'
  print *, yg(1:ncoord12)
  print *,' momenta : '
  print *, yg(ncoord12+1:ncoordmom)
  
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
     t0 = tg
     prmtg(1) = t0          ! initial value of time  (in au)
     prmtg(2) = t0+dtg      ! final value of time (in au)
     prmtg(3) = tstep0g     ! initial value of the time step
     prmtg(4) = accurg      ! accuracy parameter for hpcg
     prmtg(5) = 0.d0        ! parameter to be set >0 for terminating propagation (in foutp)
!     randomize by running dtg ps in the ground state before selecting new initial conditions
!     for next trajectory
     call fhpcg(tg,yg,deryg,nvar,prmtg)  ! runs ground state dynamics for tg from prmtg(1) to prmtg(2)
  print *
  print *,' values for the variables, after randomization t = ',tg,tg*convt
  print *,' positions:'
  print *, yg(1:ncoord12)
  print *,' momenta : '
  print *, yg(ncoord12+1:ncoordmom)
  print *
!  write(30,*) 'atom Ba   ',yg(1:3)*bohr
  do iat=1, nbrat2
!     write(30,*) 'atom He   ',yg(3*iat+1:3*iat+3)*bohr
  enddo
!  write(30,*)
!  write(30,*)
!  write(30,*)' spec  Ba   1.0    0. 1.0 1.0'
!  write(30,*)' spec  He   0.5    1.0 1.0 1.0'
!  write(30,*)
!  write(30,*)' bonds Ba Ba  5.0 7. 0.15 0.0'
!  write(30,*)' bonds He He  5.0 7. 0.05 0.0'
!  write(30,*)' bonds Ba He  4.0 7. 0.1 0.0'
!  write(30,*)
!  write(30,*)' inc  5.00'
!  write(30,*)
!  write(30,*) 'tg = ',tg
!  write(30,*)


         
     print *,'  jumping to excited state, trajectory nb ', itraj
!     use yg as initial conditions for a new trajectory in the excited states:
     te = 0.d0
     prmte(1) = te             ! initial value of time  (in au)
     prmte(2) = tfeps/convt    ! final value of time (in au)
     prmte(3) = tstep0e        ! initial value of the time step (in au)
     prmte(4) = accure         ! accuracy parameter for hpcg
!!!     prmte(5) = 0.d0  ! return flag: HPCG returns if OUTP set it to non-zero
     ! it is set to zero initially by HPCG
     prmte(6) = 0.d0          ! flag set to >0 in EVALHOP if there has been a hop (HPCG needs to be restarted)
     do i = 1, ncoordmom
        ye(i) = yg(i)
     enddo
     ! select initial excited state for classical nuclear coordinates
     ! (ihop, passed by COMMON/HOP/)
     ! and initializes wave packet components accordingly
     call WPINI_SO(ye,nvar,iopt)
!     update itrajoutp for outp subroutine
     itrajoutp = itraj

!
     
!     ye = 0.d0

!     ye(4)=5.6220*ang2bohr
     
!     ihop=3
     
!     ye(12+ihop)=1.d0

     
     initpcg = 0  ! for hpcg: number of times the current trajectory has been reinitialized
     limpa = 0  ! number of tentatives to half the step in addition to the 10 ones initially in hpcg
     ! ihlf = number of times the integration step has been halved (output):limit = 10+limpa
     lterminate=.false.
     do while (te.lt.prmte(2).and..not.lterminate)
        levalhop=.false.       
        derye(1:nvar) = 0.d0
        derye(1:ncoordmom) = 1.d0/dfloat(ncoordmom)   !!! error can be too large on phase
        call hpcg(te,ye,derye,nvar,prmte,nprmte,limpa,initpcg,ihlf)
        if (ihlf.ge.10+limpa) then
           print *,' trajectory nb ',itraj,' stopped, step halved too many times: ',ihlf
           lterminate=.true.
    !       prmte(1)=te
           cycle
        endif
        if (prmte(6).gt.0.d0) then
        print *,'te = ',te,'au, = ',te*convt,'ps:   new hop on ',ihop   
        prmte(6) = 0.d0
        prmte(1)=te
        endif
     end do
  print *,' end of trajectory nb ',itraj,' t = ',te,' au, = ',te*convt,' ps'         
  write(110,*) itraj, lterminate, te*convt,ihop,nbhop
  nbhop=0 ! On initialise nbhop à 0 pour chaque trajectoire
  ihop=ihopinic
  write(31,*) 'atom Ba   ',ye(1:3)*bohr
  do iat=1, nbrat2
     write(31,*) 'atom He   ',ye(3*iat+1:3*iat+3)*bohr
  enddo
  write(31,*)
  write(31,*)
  write(31,*)' spec  Ba   1.0    0. 1.0 1.0'
  write(31,*)' spec  He   0.5    1.0 1.0 1.0'
  write(31,*)
  write(31,*)' bonds Ba Ba  5.0 7. 0.15 0.0'
  write(31,*)' bonds He He  5.0 7. 0.05 0.0'
  write(31,*)' bonds Ba He  4.0 9. 0.1 0.0'
  write(31,*)
  write(31,*)' inc  5.00'
  write(31,*)
  write(31,*) 'te = ',te
  write(31,*)
         
  enddo


      
  STOP

END program






