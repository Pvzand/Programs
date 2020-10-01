 subroutine config(ntraje,irand,iopt,y,dery,Emin,Temp,ETemp,Etot,tstep0f,tstep0e,accurg,accure,tfgps,tfeps,dtgps)
!   read minimum energy configuration (with energy Emin) from input file
!   add internal energy corresponding to T=0.4K
!   by displacing Ba+ along a random direction
!   until the potential energy difference with Emin is equal to ETemp
!   where ETemp = (3(NHe+1)-6)*kT/2
!   Etot=Emin+ETemp

  implicit none

  include 'nbrat12.f90'
  include 'units.f90'

  character*4 :: atom1, atom2
  character*255 :: filenameinpcfg
  integer :: ntraje,i,j,moment,ndim,momenthe,itest,iadj,iadjhe
  integer :: ipos,irand,iopt,ndimtot,ixyz
  integer :: nhecheck, iat
!!! for common /dyn/
 integer :: ihop, ioldhop, ihopinic, nbhop


  real (kind=8) :: tinf,rcut,rpol,alfa,vit
  real (kind=8) :: xmass1gmol,xmass1,xmass2gmol,xmass2,xmass, xmassinv
  real (kind=8) :: energ0,energ0cm
  real (kind=8) :: rtemp1,rtemp2,rtemp3
  real (kind=8) :: Emin,Temp,ETemp,Etot,dx0
  real (kind=8) :: Epotread, Ekinread, Epot, Ekin
  real (kind=8) :: tstep0f,tstep0e,accurg,accure,tfgps,tfeps,dtgps,deltaTwg,deltaTwe
  real (kind=8) :: detotxg, detotxe, detotxgcm, detotxecm
  parameter (alfa = 2.6696d0)

! induced dipole-induced dipole interaction no longer calculated beyond rcut
! Use anadipol.f and uncomment the following line
!     parameter(tinf = 0.0d0, rcut = 1.0d3, rpol = 2.527d0)
!  induced dipole-induced dipole interaction calculated at all distances
! Use anadip.f and uncomment the following line
!      parameter(tinf = 0.0d0, rcut = 1.0d3, rpol = 0.0d0)

! array Dimensions
  real (kind=8), dimension (3) :: pos, vel, ycm, pcm
  real (kind=8), dimension (nvar) :: y, dery
  real (kind=8), dimension (ncoord12) :: ysv, psv, derysv, derpsv

  common/printingg/deltaTwg
  common/printinge/deltaTwe
  common/adjust/iadj,iadjhe
  common/masses/xmass1,xmass2,xmass(nbrat12),xmassinv(ncoord12)
  common/econservg/detotxg
  common/econserve/detotxe
!!!  ! index of the current, previous, and initial DIM eigenvalue (hence surface for classical coordinates)
  common/dyn/ihop, ioldhop, ihopinic, nbhop

! 1. read input data 

  open(90,file='Bap.in',status='unknown')
  read (90,*)
  read (90,*) atom1
  
  read (90,*)
  read (90,*) xmass1gmol
 
  xmass1 = xmass1gmol*convgrau
  do iat = 1, nbrat1
     xmass(iat) = xmass1
  enddo
  read (90,*)
  read (90,*) atom2

  read (90,*)
  read (90,*) xmass2gmol

  xmass2 = xmass2gmol*convgrau
  do iat = nbrat1+1, nbrat12
     xmass(iat) = xmass2
  enddo
  do i = 1, ncoord1
     xmassinv(i) = 1.d0/xmass1
  enddo
  do i = ncoord1+1, ncoord12
     xmassinv(i) = 1.d0/xmass2
  enddo
  
! Temperature (in Kelvin); if <0, read input positions AND velocities and check that T=|Temperature|;
!         if >=0, read input positions only and pull away atom 1 until reaching the required kinetic energy.
  read (90,*)
  read (90,*)
  read (90,*) Temp
! input configuration file name (file 91)
  read (90,*)
  read (90,*) filenameinpcfg

! convert temperature to energy
  if (nbrat12.GT.2) then
     etemp = 0.5d0 * dabs(temp) * EK2cm1 * DBLE(ncoord12-6)
  else
     etemp = 0.5d0 * dabs(temp) * EK2cm1
  endif
  etemp = etemp * convcmau
  
  read (90,*)
  read (90,*) irand
  read (90,*)
  read (90,*) iopt, ihop
  ihopinic = ihop  ! save initial value for iopt=3 
  read (90,*)
  read (90,*) ntraje

  
  print *,' DISSOCIATION DYNAMICS OF ',atom2,nbrat2, atom1,nbrat1
  print *,' ===================================================='
  print *,' nbrat1 = ',nbrat1,' (number of',atom1,' atoms)'
  print *,' m(',atom1,') = ',xmass1gmol,' g/mol,= ',xmass1,' au'
  print *,' m(',atom2,') = ',xmass2gmol,' g/mol,= ',xmass2,' au'
  print *
  print *,' randomization of initial conditions at T = ',temp,' K'
  print *,' Initialization for random sampling: irand = ',irand
  print *,' initialization for wave packet: iopt = ',iopt,' ihop = ',ihop
  print *
  print *,' running ntraje = ',ntraje,' trajectories'
  print *

  read (90,*)
  read (90,*) iadj
  print *,' iadj = ',iadj
  read (90,*)
  read (90,*) iadjhe
  print *,' iadjhe = ',iadjhe
  read (90,*)
  read (90,*)
  read (90,*) tstep0f
  print *,' tstep0f = ',tstep0f
  read (90,*)
  read (90,*) tstep0e
  print *,' tstep0e = ',tstep0e
  read (90,*)
  read (90,*) accurg
  print *,' accurg = ',accurg
  read (90,*)
  read (90,*) accure
  print *,' accure = ',accure
  read (90,*)
  read (90,*) tfgps
  print *,' tfgps = ',tfgps
  read (90,*)
  read (90,*) dtgps
  print *,' dtgps = ',dtgps
  read (90,*)
  read (90,*) tfeps
  print *,' tfeps = ',tfeps
  read (90,*)
  read (90,*) deltatwg, deltatwe
  print *,' deltatwg, deltatwe = ',deltatwg, deltatwe,' in ps'
  deltatwg=deltatwg/convt
  deltatwe=deltatwe/convt
  print *,' deltaTwg, deltaTwe = ',deltaTwg, deltaTwe,' in au'
  read(90,*)
  read(90,*) detotxgcm, detotxecm
  print *,' detotxgcm, detotxecm = ', detotxgcm, detotxecm
  detotxg = detotxgcm*convcmau
  detotxe = detotxecm*convcmau

  print *,' DIM + induced dipole-induced dipole P.E.S. ',' for the ionic cluster '
  print *,'--------------------------------------------','----------------------'
  print *,' tinf = ',tinf,'  (zero for energies in the DIM ','Hamiltonian calculation)'
  print *,' rcut = ',rcut,'  (distance from which calculation of',' the diatomic potential is assumed to be zero)'
  print *,' (commented out for the time being)'
  print *,' rpol = ',rpol,'  (distance below which polarizability',' is screened in the induced dipole-induced dipole interaction)'
  print *,' alfa = ',alfa,'  (polarizability in Bohr**3)'
  print *

  print *,"Momentum adjustment vector for Ba+ and He's upon hop: "
  print *,'Delta P_alpha prop. to  g_jk = grad_alpha (V_j V_k)'

  close(90)


!  initialize the arrays of propagated variables
  do i = 1, nvar
     y(i) = 0.d0
     dery(i) = 0.d0
  enddo
!  initialize position and momenta arrays used for backups
  do i=1,ncoord12
     ysv(i) = 0.0d0
     psv(i) = 0.0d0
     derysv(i) = 0.0d0
     derpsv(i) = 0.0d0
  enddo

! 2. Lecture de la configuration Ba+@HeN (en USI)

  print *
  print *,' reading input positions from file 91:',trim(filenameinpcfg)
!  print *,'./data/baphe100_pow4.cfg'
!  open(91,file='./data/baphe100_pow4.cfg',status='old')
  open(91,file=trim(filenameinpcfg),status='old')
  if (Temp.ge.0.d0) then
     read(91,'(15x,I5,32x,f19.13)') nhecheck,energ0cm
     energ0=energ0cm*convcmau
     print *,' nhe = ',nhecheck,' Energ0 = ',energ0,' au, = ',energ0cm,' cm-1'
     if (nhecheck.ne.nbrat2) then
        if (nhecheck.gt.nbrat2) then
           print *
           print *,'   >>> WARNING <<<  '
           print *,' input configuration file has nhe = ',nhecheck
           print *,' nbrat2 = ',nbrat2,' expected'
           print *,' using only the first ',nbrat2,' ones'
        else
           print *
           print *,'   >>> ERROR <<<  '
           print *,' input configuration file has nhe = ',nhecheck
           print *,' nbrat2 = ',nbrat2,' expected'
           STOP
        endif
     endif

!     read input coordinates for atom(s) 1 (in Angstroms):
     print *
     print *,' Initial atom coordinates in Angstroms :'
     do i = 1, nbrat1
        ipos = 3*(i-1)
        read(91,'(A4,I3,3(1x,g15.8))')atom1,iat,(pos(j),j=1,3)
        print *,atom1,iat,(pos(j),j=1,3)
        call flush(6)
        !     convert to atomic units and store in y
        do j = 1, 3
           y(ipos+j) = pos(j)/bohr
        enddo
     enddo
!         print *,(pos(j),j=1,3),(vit(j),j=1,3)
!         print *,ipos,(y(ipos+j),j=1,3),(dery(ipos+j),j=1,3),
!     &           (y(moment+ipos+j),j=1,3)
!     call flush(6)

!  read input coordinates for atoms 2 (in Angstroms):
     do i=1,nbrat2
        ipos = 3*(nbrat1+i-1)
        read(91,'(A4,I3,3(1x,g15.8))') atom2,iat,(pos(j),j=1,3)
        print *,atom2,iat,(pos(j),j=1,3)
        call flush(6)
!  convert to atomic units and store in y
        do j=1,3
           y(ipos+j) = pos(j)/bohr
        enddo
     enddo
     close(91)

!     calculate initial potential energy:     
     print *
     call fpottot2(y,nvar,energ0)
     energ0cm = energ0*convaucm
     print *,' Initial potential energy: ',energ0,' au, = ',energ0cm,'cm-1'

!  pull the Ba atom away from equilibrium in any random direction
!     until reaching internal energy corresponding to Temp
     dx0=2.d0 ! unit of distance to pull Ba initially (in a.u.)
     call tire(y,nvar,dx0,irand,etemp,energ0,etot)


  else
!!! here read positions AND velocities from input configuration file
!!! then checks kinetic energy (temperature)
     print *,' reading input positions AND velocities from file 91: ',trim(filenameinpcfg)
     read (91,*) Epotread, Ekinread
     print *,' potential, kinetic energies as recorded in file 91: ',Epotread, Ekinread,' J'
     Epotread = Epotread * EJ2Hartree
     Ekinread = Ekinread * EJ2Hartree
     print *,' potential, kinetic energies as recorded in file 91: ',Epotread, Ekinread,' Hartree'
     print *,' potential, kinetic energies as recorded in file 91: ',Epotread*convaucm, Ekinread*convaucm,' cm-1'
     print *,' potential, kinetic energies as recorded in file 91: ',Epotread*convaucm/EK2cm1, Ekinread*convaucm/EK2CM1,' K'
!     read input coordinates for atom(s) 1 (in meters):
     print *
     print *,' Initial atom coordinates in meters, velocities in m/s :'
     do iat = 1, nbrat1
        ipos = 3*(i-1)
        read(91,*) (pos(j),j=1,3), (vel(j),j=1,3)
        print *,atom1,iat,(pos(j),j=1,3), (vel(j),j=1,3)
        call flush(6)
        !     convert to atomic units and store in y
        do j = 1, 3
           y(ipos+j) = pos(j)*xm2bohr
           y(ncoord12+ipos+j) = vel(j)*xmass1*xmsinv2au
        enddo
     enddo

!  read input coordinates for atoms 2 (in Angstroms):
     do iat=1,nbrat2
        ipos = 3*(nbrat1+i-1)
        read(91,*) (pos(j),j=1,3), (vel(j),j=1,3)
        print *,atom2,iat,(pos(j),j=1,3), (vel(j),j=1,3)
        call flush(6)
!  convert to atomic units and store in y
        do j=1,3
           y(ipos+j) = pos(j)*xm2bohr
           y(ncoord12+ipos+j) = vel(j)*xmass2*xmsinv2au
        enddo
     enddo
     close(91)
  endif

!     3. put all atoms in the center of mass frame
  print *
  print *, '  determine center of mass position,'
  print *,' and reset all atoms in the center of mass frame'

!     do i=1, ncoord12
!        ysv(i) = y(i)
!     enddo
     ysv(1:ncoord12) = y(1:ncoord12)
     rtemp1 = dfloat(nbrat1)*xmass1 + dfloat(nbrat2)*xmass2 ! total mass

! calculate center of mass position ycm(3) and momentum pcm(3)
     do ixyz=1,3
        rtemp2 = 0.0d0
        do i=1,nbrat1
           ipos = 3*(i-1)
           rtemp2 = rtemp2 + xmass1*ysv(ipos+ixyz)
        enddo
        do i=1, nbrat2
           ipos = ncoord1 + 3*(i-1)
           rtemp2 = rtemp2 + xmass2*ysv(ipos+ixyz)
        enddo
        ycm(ixyz) = rtemp2/rtemp1
     enddo
     print *,' center of mass coordinates : '
     print *,(ycm(ixyz),ixyz=1, 3),' bohr'
     print *,(ycm(ixyz)*bohr,ixyz=1,3),' Angstrom'
                                                   
!! subtract center of mass position from all atom positions;
     do ixyz=1,3
        do i = 1, nbrat12
           ipos = 3*(i-1)
           y(ipos+ixyz) = ysv(ipos+ixyz) - ycm(ixyz)
        enddo
        
     enddo

     print *,'The origin is now at the center of mass:'
     print *,'POSITIONS in Bohr = '
     print *,atom1,' : ',(y(i),i=1,ncoord1)
     do iat = 1, nbrat2
        print *,atom2,iat,' : ',(y(ncoord1+i),i=3*(iat-1)+1,3*iat)
     enddo
     print *,'POSITIONS in Angstrom = '
     print *,atom1,' : ',(y(i)*bohr,i=1,ncoord1)
     do iat = 1, nbrat2
        print *,atom2,iat,' : ',(y(ncoord1+i)*bohr,i=3*(iat-1)+1,3*iat)
     enddo
     call flush(6)
     print *
     call fpottot2(y,nvar,energ0)
     energ0cm = energ0*convaucm
     print *,' Check Initial total energy after pulling and COM: ', energ0,' au, = ',energ0cm,'cm-1'

!!! Case where velocities were read in addition to positions (Temp <0)
!!! calculate momenta in the COM system
     if (Temp.lt.0.d0) then
        do ixyz = 1, 3
           pcm(ixyz) = sum(y(ncoord12+ixyz:ncoordmom:3))
        enddo
        print *,'   Center of Mass Momentum before setting the origin at the COM'
        print *, pcm(1:3), 'au (hbar/a0)'
        Ekin = 0.5d0*sum(xmassinv(1:ncoord12)*y(ncoord12+1:ncoordmom)**2)
        print *,'   Total kinetic energy = ',Ekin,' Hartree, = ',Ekin*Hartree2cm, ' cm-1, ',Ekin*Hartree2cm/EK2cm1,' K'
        do ixyz = 1, 3
           y(ncoord12+ixyz:ncoordmom) = y(ncoord12+ixyz:ncoordmom) - pcm(ixyz)
        enddo
        print *,'   Center of Mass Momentum before setting the origin at the COM'
        print *,'   (should be zero) '
        print *, pcm(1:3), 'au (hbar/a0)'
!  check total kinetic energy
        Ekin = 0.5d0*sum(xmassinv(1:ncoord12)*y(ncoord12+1:ncoordmom)**2)
        print *,'   Total kinetic energy = ',Ekin,' Hartree, = ',Ekin*Hartree2cm, ' cm-1, ',Ekin*Hartree2cm/EK2cm1,' K'
        print *,'  COM kinetic energy = ',0.5d0*sum(pcm(1:3)**2)/rtemp1,' Hartree'
        print *,'  check by comparing COM K.E. + K.E. in the COM system to K.E. before '
        print *, Ekin + 0.5d0*sum(pcm(1:3))**2/rtemp1,' Hartree'
        
!!! If total kinetic energy
!!!!!!!!!!! TO BE CONTINUED !!!!!!!!!!!!!!
        !!!! reflechir s'il faut chercher le minimum ici ou s'il faut dans foutp !!!!

     endif
     
        
  return
  end subroutine config
