!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                !!!!
!!!!              PROGRAM KINETICS  (version 1.0)                 !!!!
!!!!                                                                !!!!
!!!! Program to follow the formation of fragments by generating 2D  !!!!
!!!! and 3D mass spectra from movie files (XYZ format) coming from  !!!!
!!!! DYNIMAT runs. This program is freely inspired from             !!!!
!!!! MC2-SPECTRA-dv2.2 (21/11/2016).                                !!!!
!!!!                                                                !!!!
!!!!                                                                !!!!
!!!! List of routines called by the main program:                   !!!!
!!!! 1) SUBROUTINE INCSORT generates configurations with particles  !!!!
!!!!    (coordinates, mass, symbol) sorted in such a way            !!!!
!!!!    that their distance to the molecular COM increases (l. 893).!!!!
!!!!                                                                !!!!
!!!!                                                                !!!!
!!!! Author: Dr David A. Bonhommeau(1)                              !!!!
!!!! modified by Nadine
!!!! to get other observables as a function of time
!!!! (adapted to other trajectory output files;
!!!! mass spectra treated separately
!!!!                                                                !!!!
!!!! (1)  University of Reims Champagne-Ardenne (France).           !!!!
!!!!                                                                !!!!
!!!!                                                                !!!!
!!!! Coding start date: 25/04/2018                                  !!!!
!!!! Last modification: 22/11/2018                                  !!!!   
!!!!                                                                !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! General approach to generate the mass spectra
!
! 1.   Reading of input configurations
! 1.1. Reading of the pieces of information provided by the user. 
! 1.2. Reading of configurations and move into the molecular frame attached to
!!! either Ba+ or COM.
! 1.3. Sorting of particles based on their distance to the molecular frame. 
!
! 2.   Identification of fragments 
! 2.1. Distances between particles are computed.
! 2.2. If the distance dist(1,i) between particle 1 and particles i 
!      is below a threshold d_frag particle 1 and its affiliated particles i 
!      are assumed to form fragment 1 with mass Sum_i m_i (mass of particles 
!      composing the fragment). Fragment 2 is obtained by testing the first 
!      particle that do not belong to fragment 1, etc. At the end of this 
!      process we obtain N_frag<=N primary trial fragments.  
! 2.3. To ensure that these fragments have no particle in common, we generate
!      a vector frag(i,j) for each fragment where i is the fragment number
!      and j<=N is its component (=1 if particle j is included in the cluster, 
!      0 otherwise). Independent fragments should have zero overlap. In the 
!      contrary case (ie, non-zero overlap), the two fragments are assumed 
!      to form the same fragment. This process should lead to the real 
!      fragments after several iterations. The resulting frag vectors are 
!      reordered and the coordinates of the largest fragment are saved.
!
! 3.  Observables:
!      NHe(t) attached to Ba+
!      Number of He's (t) attached to Ba+ as an exciplex
!      Kinetic energy of dissociated He's
!      They are averaged over the N_traj 
!      trajectories at the end of the program.

program spectradyn

 implicit none

! Some important parameters (hard-coded or defined by the user):
! N_tot: Total number of atoms.
! N_at = Number of atoms of each type.
! N_traj = Number of trajectories per set of MD simulations.
! N_sim = Number of sets of MD simulations.
! N_step = Number of time steps printed out for each trajectory of each set of MD simulations
!          (= "configurations in movie files" for David).
! d_frag = distance above which fragmentation is assumed. Beware: a too large 
!          value of d_frag would lead to the detection of only one fragment. 
!          Selecting a smaller distance is therefore more appropriate than a large one, 
!          provided that the distance is significantly larger than typical equilibrium distances.
!          Here the value of d_frag should properly represent the fragmentation of 
!          Rgn+@HeN clusters. 
!          For Ar, R_eq(He-Ar) = 6.6 a0, R_eq(He-Ar+) = 5.4 a0, 
!                  R_eq(He-He,eff) = 4.3-4.4 A = 8.1-8.3 a0,
!                  R_eq(He-Ar,eff) = 4.1-4.2 A = 7.8-7.9 a0, 
!                  R_eq(Ar-Ar+) = 2.4-5.8 A = 5.5-11.0 a0 (depending on the molecular state).
!          In general it is better to use a d_frag about twice as large as R_eq, that is about 
!          12-16 a0 in the case of argon clusters embedded in helium nanodroplets.
!          For Ba+, R_eq(He-Ba+) = 5 A ~ 10 a0
!                   Req(He-He, eff) = 4.1-4.4 A ~ 8-9 a0
!          A value of d_frag ~ 20 aO is probably suitable for Ba+.
 INTEGER :: i, j, k, l, m, n, N_ttraj, itraj, iconv, icpt, lat
 INTEGER :: jdum   !!!NH added
 INTEGER :: imass, imass_max, iprint, ifrag, ifrag_real, prodscal, fragsize_max
 INTEGER :: N_tot, N_tot2, N_print, N_sim, N_traj, N_type, NHe_tot, NBa_tot, N_rmv
 INTEGER, ALLOCATABLE, DIMENSION(:) :: N_at, NHe, NBa, fragsize, fragsize_temp, ideloc
 INTEGER, ALLOCATABLE, DIMENSION(:,:) :: N_step, frag, frag_temp
 DOUBLE PRECISION :: mass_bar, perc, mass_tot, time, time_f, time_step, time_lim, dt_frame
 DOUBLE PRECISION :: massfrag_tot
 double precision :: NHeexcplx   !!! nb of He in Ba+ exciplex f(t) for each traj
 double precision :: timeunit, distunit
 double precision :: timeau
 DOUBLE PRECISION, PARAMETER :: hundred=100.0d0, zero=0.0d0, one=1.0d0, tenth=0.1d0
!!! DOUBLE PRECISION, PARAMETER :: d_frag=15.0d0, d_excplx=4.d0
!!! DOUBLE PRECISION, PARAMETER :: d_frag=17.5d0, d_excplx=4.5d0
!!! DOUBLE PRECISION, PARAMETER :: d_frag=20.d0, d_excplx=5.d0
 DOUBLE PRECISION, PARAMETER :: d_frag=20.d0, d_excplx=4.5d0

 DOUBLE PRECISION, DIMENSION(3) :: rcom
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: mass, mass_at, mass_com, time_ref
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: coords, coords_com
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: dist, mass_frag
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: mspec
 CHARACTER(len=1000) :: name_frame, name_path, name_in, name_out, name_comp
 CHARACTER(len=3), ALLOCATABLE, DIMENSION(:) :: atsymb, atsymb_com, typesymb
! CHARACTER(len=3), ALLOCATABLE, DIMENSION(:,:) :: atsymb_frag
 LOGICAL, ALLOCATABLE, DIMENSION(:) :: lfrag, lprod, Bapres
 LOGICAL, ALLOCATABLE, DIMENSION(:) :: sflag
 LOGICAL :: first

 double precision, ALLOCATABLE, dimension(:) :: NHehist, NHe_excplx
 
 write(6,*)
 write(6,*) '                 WELCOME to Kinetics-Ba'
 write(6,*) 'A program that generates time-dependent observables'
 write(6,*) ' upon photoexcitation of Ba+HeN clusters.'
 write(6,*)
 write(6,*) 'NB: A folder called RESULTS must be created for storing'
 write(6,*) 'the output files before running the program. Please, check it'
 write(6,*) 'before proceeding.'
 write(6,*)
 print *,'   read input in atomic units, outputs in ps and Angstroms '
 print *,'   atom-atom distance above which dissociation is assumed: d_frag = ',d_frag,' Angs'
 print *,'   Ba+-He distance below which exciplex is assumed: d_excplx = ',d_excplx,' Angs'
 call flush(6)

 timeunit =  2.418884327d-5 ! time from au to picoseconds
 distunit = 0.52917720859d0 ! distance from Bohr to Angstrom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                !!!!
!!!! 1. READING OF THE SETUP FILE   !!!! 
!!!!                                !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1.1. General information on the simulation
 write(6,*) 'GENERAL INFORMATION ON THE SIMULATION'
 write(6,*) 'Number of MD runs?'
 read(5,*) N_sim
 print *, N_sim
 write(6,*) 'Number of trajectories by MD run?'
 read(5,*) N_traj
 print *, N_traj
! name_path should be a list of paths 
! MD simulation 1: as many rows as trajectories
! One blank
! MD simulation 2: as many rows as trajectories
! One blank
! ... 
! MD simulation N: as many rows as trajectories 
 write(6,*) 'Name of the file that collects the paths to configuration files?'
 read(5,*) name_path
 print *, trim(name_path)
 open(110, file=trim(name_path), status='old')
 write(6,*) 'Name of the file that collects the paths to input files for frame number evaluation?'
 read(5,*) name_frame
 print *, trim(name_frame)
 write(6,*) 'Number of lines to discard when counting frames (eg, related to comments)?'
 read(5,*) N_rmv
 print *, N_rmv
 call flush(6)
 
 open(111, file=trim(name_frame), status='old')

! Typical number of frames ~ 3 000 but this number may slightly vary for different trajectories.
! N_step is initialized to -N_rmv since N_rmv lines must be removed from the total number of lines in the input
! file that allows the calculation of the number of frames.
 ALLOCATE(N_step(N_sim,N_traj))
!!!NH removed N_step = -N_rmv
  N_step = 0   !!!NH added
 do i=1,N_sim
    read(111,'(a)') name_in   !NH added: only one file per simulation, N_traj each
    print *, trim(name_in)    !NH added
    call flush(6)
    open(112,file=trim(name_in),status='old')   !NH added
    do j=1,N_rmv    !NH added
       read(112,*)   ! NH added
    enddo           ! NH added
    do j=1,N_traj
!!!NH removed      read(111,'(a)') name_in
!!!NH removed       print *, trim(name_in)
!!!NH removed       open(112,file=trim(name_in),status='old')
       read(112,*) !!!NH added  (2 blank lines before each new trajectory)
       read(112,*) !!!NH added

!!!NH removed       do  
!!!NH removed         read(112,*,end=998)
       k=j   !!!NH added
       do while (k == j)   !!!NH added
          read(112,*,end=998) k, timeau   !!!NH added
          N_step(i,j) = N_step(i,j) + 1
       enddo
       N_step(i,j) = N_step(i,j) - 1   !!!NH added
       backspace(112)
       backspace(112)
       backspace(112)
998    continue   !!!NH added
       print *,' trajectory nb ',j,': ',N_step(i,j),' time steps'
       call flush(6)
!!!NH removed 998    close(112)
    enddo
    close(112)   !!!NH added
    read(111,*,end=999)
 enddo
999 close(111)

! XXX files are closed after evaluating the number of frames.
! They are open once again for later data treatment.
!!! NH removed  open(111, file=trim(name_frame), status='old')  (not needed with my files where t is printed in the trajectory files)

! Time interval between frames: (simulation time step)*(printing period) ~ 1e-4*100 = 1e-2
! in typical MD simulations.
 write(6,*) 'Time interval between two frames (ps) in input trajectories?'
 call flush(6)
 read(5,*) dt_frame
 print *, dt_frame   !!NH added
 call flush(6)
! The time limit is fixed to dt_frame
  time_lim = dt_frame
! time_lim = tenth*dt_frame  
 write(6,*) 'Elapsed time of the dynamics (ps)?'
call flush(6)
  read(5,*) time_f
 print *, time_f   !!NH added
call flush(6)
  
 write(6,*) 'Name of the file that details the cluster composition '
 write(6,*) '(particle symbols, numbers, and masses)?' 
call flush(6)
  read(5,*) name_comp
 print *, trim(name_comp)   !!NH added
call flush(6)
  open(100,file=trim(name_comp), status='old')
 read(100,*) N_type
 ALLOCATE(typesymb(N_type),N_at(N_type),mass_at(N_type)) 
 N_tot = 0
 mass_tot = zero
! WARNING: We assume that only He and Ba atoms are present.
 do i=1,N_type
    read(100,*) typesymb(i),N_at(i),mass_at(i)
    print *,N_at(i),typesymb(i),' atoms with mass ',mass_at(i)   !!NH added
    N_tot = N_tot + N_at(i)
    if ( typesymb(i) == 'He' ) then 
       NHe_tot = N_at(i)
    else 
       NBa_tot = N_at(i)
    endif
    mass_tot = mass_tot + mass_at(i)*dble(N_at(i))
 enddo
 print *,' total mass = ',mass_tot   !!NH added
call flush(6)
 !!!NH added: here it is assumed that Ba atoms are first, then He atoms
 ALLOCATE(coords(N_tot,3),mass(N_tot),atsymb(N_tot))   !!!NH added
 mass(1:NBa_tot) = mass_at(1)   !!!NH added
 mass(NBa_tot+1:N_tot) = mass_at(2)   !!!NH added
 atsymb(1:NBa_tot) = typesymb(1)
 atsymb(NBa_tot+1:N_tot) = typesymb(2)
 print *,'>>>>> WARNING: Make sure that file 120 has Ba coordinates first, then Hes <<<<<'!!!NH added
 print *,' and that mass 1 is Ba mass and mass 2 He mass in ',trim(name_comp)   !!!NH added
 print *,' total number of atoms = ',N_tot,' total mass = ',mass_tot   !!!NH added
call flush(6)
 
! 1.2. Information on mass spectra

 write(6,*)
 write(6,*) 'GENERAL INFORMATION ON MASS SPECTRA'

 write(6,*) 'Number of time steps to be generated for observables from t=0 to ',time_f,' ?' 
call flush(6)
  read(5,*) N_print
 print *, N_print   !!NH added
call flush(6)
 ! N_print => N_print - 1 values excluding the first step.
 ALLOCATE(time_ref(N_print))
 time_ref(1) = zero
 time_step = time_f/dble(N_print-1)
 do i=2,N_print
    time_ref(i) = time_ref(1) + time_step*dble(i-1)
 enddo
 write(6,*) 'Size bin for mass spectra (in amu)?'
! write(6,*) 'NB: Typical values are 0.01, 0.1 or 1 amu.'
call flush(6)
  read(5,*) mass_bar
 print *, mass_bar   !!NH added
call flush(6)
 
! Opening output files
!!! do i=1,N_print
! Output file that collects information on charged fragments: number of He and Ba atoms and mass.
!    write(name_in,'("RESULTS/Info-ch-",i2.2,".dat")')i
 !    open(200+i,file=trim(name_in),status='unknown')
 open(201,file='RESULTS/kin_NHe.dat',status='unknown')
 write(201,"('#   t(ps)      NHe(t)   NHe-excplx(t) ')")
 call flush(201)
!    write(200+i,'(a6,a9,2a5,a16,a8,a14)') '# Traj','Frag nb','He','Ba','Mass (amu)' 
!! Output file that collects information on neutral fragments: number of He and Ba atoms and mass. 
!    write(name_in,'("RESULTS/Info-neut-",i2.2,".dat")')i
!    open(300+i,file=trim(name_in),status='unknown') 
!    write(300+i,'(a6,a9,2a5,a16)') '# Traj','Frag nb','He','Ba','Mass (amu)'
! enddo 
! Output file that collects information on charged fragments: number of He and Ba atoms and mass.
!    write(name_in,'("RESULTS/Info-ch-",i2.2,".dat")')i
 !    open(200+i,file=trim(name_in),status='unknown')
 open(401,file='RESULTS/kin_NHetraj.dat',status='unknown')
 write(401,"('#  itraj   t(ps)      NHe(t)   NHe-excplx(t) ')")
 call flush(401)
 ! Output file for fragments sizes
 open(501,file='RESULTS/kin_fragsize-traj.dat',status='unknown')
 write(501,"('#   itraj   t(ps)   Nb-fragments(t)   Frag-sizes(t)  ')")
 call flush(501)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                         !!!!
!!!! 2. READING OF THE CONFIGURATION FILES   !!!! 
!!!!                                         !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

 write(6,*)
 write(6,*) 'START OF THE JOB'

! Total number of trajectories
 N_ttraj = N_sim*N_traj
 print *,' total number of trajectories to be treated : ',N_ttraj   !!!NH added
 call flush(6)
 

! imass_max = int( mass_tot/mass_bar ) + 1
!! mspec(imass_max,N_print,1) = mass spectrum for charged fragments (experimentally accessible)
!! mspec(imass_max,N_print,2) = mass spectrum for neutral fragments (not experimentally accessible)
! ALLOCATE(mspec(imass_max,N_print,2))
 ! mspec = zero

!!! allocate number of time steps for histogram
 ALLOCATE(NHehist(N_print),NHe_excplx(N_print))
!!! initialize observables
 NHehist=0.d0
 NHe_excplx = zero

 ! Averages and array containing averages are initialized
!!! NH version: everything (time AND coords) is read from the same files,
!!! with one common file for each set of trajectories (1 to N_sim),
!!! in which 2 blank lines precede each trajectory
!!!
 do i=1,N_sim     !!! do loop on the number of MD runs
   
    write(6,*) 'Treatment of data set ',i,' underway ...'
!!!NHadd: moved here!! Opening of the MOVIE file for trajectory itraj = N_traj*(i-1)+j
    read(110,'(a1000)') name_in  !!!NH added
    open(120,file=trim(name_in),status='old') !!!NH added
    print *,' reading input trajectories from file :',trim(name_in)   !!!NH added
    call flush(6)
 
    do j=1,N_rmv    !NH added
       read(120,*)   ! NH added
    enddo           ! NH added

    do j=1,N_traj

       iprint = 1
       itraj = N_traj*(i-1)+j
       print *,' j, itraj = ',j,itraj
       print *,' -------------'
       call flush(6)
 
!!!! Opening of the MOVIE file for trajectory itraj = N_traj*(i-1)+j
!!!NH removed       read(110,'(a1000)') name_in
!!!NH removed       open(120,file=trim(name_in),status='old') 
!!!NH removed ! Opening of the input file for evaluating frames, for trajectory itraj = N_traj*(i-1)+j
!!!NH removed       read(111,'(a1000)') name_in
!!!NH removed       open(121,file=trim(name_in),status='old')              
!!!NH removed       read(121,*)
       read(120,*)   !!!NH added
       read(120,*)   !!!NH added
       print *,' number of time steps: ',N_step(i,j)
       print *,' time_lim = ',time_lim
       call flush(6)
 
       do k=1,N_step(i,j)

! 2.1 Reading of XXX.dat files to check whether the current time step
!     has to be considered to plot mass spectra.    
!
!!!NH removed            read(121,*) time
           read(120,*) jdum, timeau   !!!NH added
           time=timeau*timeunit   !!!NH added
          if ( dabs(time-time_ref(iprint)) < time_lim )  then
             backspace(120)   !!! back up 1 line to read coordinates as well
! 2.2. Reading of configurations (+ move into the COM frame).
!
!!!NH removed ! Reading of the MOVIE file for trajectory itraj
!!!NH removed             read(120,*) N_tot2
!!!NH removed               read(120,*)
!!!NH removed               if ( N_tot /= N_tot2 ) then 
!!!NH removed                  write(6,'(a,I6,a,I6,a)') 'Check the number of atoms since N_tot =',&
!!!NH removed                                           N_tot,' and N_tot2 = ',N_tot2,'.'
!!!NH removed                  stop
!!!NH removed               endif
!!!NH removed             ALLOCATE(coords(N_tot,3),mass(N_tot),atsymb(N_tot))
             ALLOCATE(coords_com(N_tot,3),mass_com(N_tot),atsymb_com(N_tot))
!!! reset atom order in the list (previous trajectory may have inverted Ba and He in INCSORT)
             mass(1:NBa_tot) = mass_at(1)   !!!NH added
             mass(NBa_tot+1:N_tot) = mass_at(2)   !!!NH added
             atsymb(1:NBa_tot) = typesymb(1)   !!!NH added
             atsymb(NBa_tot+1:N_tot) = typesymb(2)   !!!NH added

             read(120,*) jdum, timeau, (coords(l,1:3),l=1,N_tot)   !!!NH added
!             print *,jdum,timeau,coords(1,1:3), coords(2,1:3),' ... ',coords(N_tot,1:3)
             coords = coords*distunit   !!!NH added
             rcom = zero
             do m=1, 3   !!!NH added
                rcom(m)=sum(mass(1:N_tot)*coords(1:N_tot,m))   !!!NH added
             enddo   !!!NH added
!             print *,' COM coordinates : ',rcom(1:3)
!!!NH removed              do l=1,N_tot
!!!NH removed                 read(120,*) atsymb(l),(coords(l,m),m=1,3)
!!!NH removed                 do m=1,N_type
!!!NH removed                    if ( atsymb(l) == typesymb(m)) mass(l) = mass_at(m)
!!!NH removed                 enddo
!!!NH removed                 do m=1,3
!!!NH removed                    rcom(m) = rcom(m) + mass(l)*coords(l,m)
!!!NH removed                 enddo 
!!!NH removed              enddo
        
             do l=1,3
                rcom(l) = rcom(l)/dble(mass_tot)
             enddo
!!! change here reference to Ba+
             rcom(1:3) = coords(1,1:3)
!             print *,'   Ba+   coordinates : ',rcom(1:3)
             do l=1,N_tot
                do m=1,3
                   coords_com(l,m) = coords(l,m) - rcom(m)
                enddo
                atsymb_com(l) = atsymb(l)
                mass_com(l) = mass(l) 
             enddo

! 2.3. Sorting of particle positions
!
! For the subsequent identification of fragments to be effective it is more convenient 
! to sort clusters with respect to their distance to the molecular COM (increasing order).
! In the contrary case there is a small (but non negligible) probability for finding two
! trial fragments side by side without any particle in common.
! The new configurations are copied into coords(N_tot,3).
! In input files the first atoms are always Ba atoms.

             call incsort(N_tot,coords_com,mass_com,atsymb_com,coords,mass,atsymb)
             DEALLOCATE(coords_com,mass_com,atsymb_com)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                !!!!
!!!! 3. IDENTIFICATION OF FRAGMENTS !!!!
!!!!                                !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 3.1. Interatomic distances
             ALLOCATE(dist(N_tot,N_tot))
             dist = zero
             do l=1,N_tot-1
                do m=l+1,N_tot
                   dist(l,m) = dsqrt( ( coords(l,1) - coords(m,1) )**2 + &
                                      ( coords(l,2) - coords(m,2) )**2 + & 
                                      ( coords(l,3) - coords(m,3) )**2 ) 
                   dist(m,l) =  dist(l,m)
                enddo
             enddo
!!! 3.2.0. determination of exciplexes
             l=2
             NHeexcplx=0
             do while (dist(1,l).lt.d_excplx)
!                NHe_excplx(iprint)=NHe_excplx(iprint)+1
                NHeexcplx=NHeexcplx+1
                l=l+1
             end do
                NHe_excplx(iprint)=NHe_excplx(iprint)+NHeexcplx
             
! 3.2. Determination of trial fragments
!
! lfrag(i) = .false. -> atom i was not assigned to a fragment yet.
!          = .true. -> atom i has already been assigned to a fragment.
! lprod(i) = .false. -> fragment i was not involved in a scalar product during the 
!                       current iteration.
!           = .true. -> fragment i was involved in a scalar product during the 
!                       current iteration.

             ALLOCATE(lfrag(N_tot),frag(N_tot,N_tot),fragsize(N_tot))
             
             ifrag = 0
             frag = 0
             fragsize = 0
             lfrag = .false.

             do l=1,N_tot
                if ( .not.lfrag(l) ) then
  
                   ifrag = ifrag + 1
                   frag(ifrag,l) = 1
                   lfrag(l) = .true.
                   fragsize(ifrag) = fragsize(ifrag) + 1 

                   do m=1,N_tot
                      if ( ( m /= l ).and.( dist(l,m) < d_frag ) ) then 
                         frag(ifrag,m) = 1
                         lfrag(m) = .true.
                         fragsize(ifrag) = fragsize(ifrag) + 1 
                      endif
                   enddo

                endif
             enddo

! Partial deallocation of memory
             DEALLOCATE(lfrag,dist)

! 3.3. Determination of fragments.
!
!      Warning: Once fragment l has been updated, other fragments must be 
!               compared with this new fragment. That is the reason why we introduced 
!               an iterative convergence check characterized by "iconv". While iconv 
!               is not zero overlaps are detected between several fragments and the 
!               loop proceeds. ifrag_real is initialized to ifrag that is an upper 
!               bound of ifrag_real.

             iconv = 1
             ifrag_real = ifrag
             do while ( iconv /= 0 )

! 3.3.1. Check of the trial fragment by computing the scalar product of frag vectors.

                iconv = 0
                ALLOCATE(lprod(ifrag),frag_temp(ifrag,N_tot),fragsize_temp(ifrag))
                lprod = .false.
                do l=1,ifrag
                   do m=1,N_tot 
                      frag_temp(l,m) = frag(l,m)
                   enddo
                   fragsize_temp(l) = fragsize(l)
                enddo

                do l=1,ifrag
                   if ( .not.lprod(l) ) then
                      do m=1,ifrag
                         if ( m /= l ) then
                            prodscal = 0
                            do n = 1,N_tot
                               prodscal = prodscal + frag(l,n)*frag(m,n)
                            enddo
                            if ( prodscal /= 0 ) then 
! The two fragments share some atoms and consequently belong to the same fragment.
                               do n=1,N_tot
                                  if ( ( frag(l,n) == 0 ).and.( frag(m,n) == 1 ) ) then 
! The second IF loop is required in case of multiple fragments that overlap (at least
! three), it is not needed when overlaps only concern two primary fragments.
                                     if ( frag_temp(l,n) == 0 ) fragsize_temp(l) = fragsize_temp(l) + 1
                                     frag_temp(l,n) = 1
                                  endif
                                  frag_temp(m,n) = 0               
                               enddo
                               fragsize_temp(m) = 0
                               iconv = iconv + 1
                               if ( .not.lprod(m) ) ifrag_real = ifrag_real - 1
                               lprod(m) = .true.
                            endif
                         endif
                      enddo
                   endif
                enddo

                DEALLOCATE(frag,fragsize)
                ALLOCATE(frag(ifrag_real,N_tot),fragsize(ifrag_real))

! 3.3.2. Reordering of the frag vectors.
!
!        At this point, the frag vectors contain the "real" fragments but some 
!        fragments may be separated by blank vectors, ie frag(1,i) can contain the 
!        indices of the atoms that compose the first fragment, frag(2,i) may be blank, 
!        frag(3,i) may contain indices of the second fragment, etc. 
!        The frag vectors need to be reordered to fit the dimension ifrag_real.
!        NB: frag(icpt,m) is not converted into frag(icpt,icpt2) (double reordering) 
!        since we would lose the information on the atomic index by doing this. 
!        The information on the fragment size is stored in fragsize(icpt).

                icpt = 0  
                do l=1,ifrag
                   if ( fragsize_temp(l) > 0 ) then
                      icpt = icpt + 1
                      do m=1,N_tot
                         frag(icpt,m) = frag_temp(l,m)
                      enddo
                      fragsize(icpt) = fragsize_temp(l) 
                   endif
                enddo

                DEALLOCATE(lprod,frag_temp,fragsize_temp)

                ifrag = ifrag_real
             enddo
!             print *,' nb of fragments: ',ifrag_real

! 3.3.3. Storage of the fragment masses
! We calculate fragsize_max in advance to minimize the size of the mass_frag array.
             fragsize_max = 0
             do l=1,ifrag_real
                if ( fragsize(l) > fragsize_max ) fragsize_max = fragsize(l)
             enddo
!             print *,' fragment sizes: ',fragsize(1:ifrag_real)
!             print *,'    size of the fragment containing Ba+: ',fragsize(1)

! Storage of the mass of all the fragments in mass_frag
             ALLOCATE(mass_frag(ifrag_real,fragsize_max),NHe(ifrag_real),NBa(ifrag_real),Bapres(ifrag_real))
             NHe = 0
             NBa = 0
             Bapres=.false.
             do l=1,ifrag_real
                lat = 0 
             
                do m=1,N_tot
                   if ( frag(l,m) == 1 ) then 
                      lat = lat + 1
                      if ( atsymb(m) == 'He' ) then
                         NHe(l) = NHe(l) + 1
                      elseif (atsymb(m) == 'Ba') then
                         NBa(l) = NBa(l) + 1
                         Bapres(l) = .true.
                      else
                         write(6,*) 'At least one atomic symbol is not valid!'
                         stop
                      endif 
                      mass_frag(l,lat) = mass(m)
                   endif
                enddo  

             enddo

             write(401,*) itraj, time, NHe(1), NHeexcplx
             write(501,*) itraj, time, ifrag_real, (NHe(l)+NBa(l),l=1,ifrag_real)

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                 !!!!
!!!! 4. MASS SPECTRA !!!!
!!!!                 !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NHe(t) = number of He's in the 1st fragment (which is the one containing Ba+)
             NHehist(iprint) = NHehist(iprint) + NHe(1)

! Remarks:  
! mspec(i,j,k) => i = index of the fragment mass bin 
!                 j = mass spectrum number 
!                 k = 1 or 2 for charged or neutral fragments

! Mass spectra are plotted for neutral and charged fragments.
! At this point, frag(i,j) is an array of ifrag_real vectors with the real fragments.

!             do l=1,ifrag_real

! Total mass of the fragments
! NB: massfrag_tot is no longer the total mass of the parent cluster.
!                massfrag_tot = 0.d0   !!!NH added: .d0
!                do m=1,fragsize(l) 
!                   massfrag_tot  = massfrag_tot + mass_frag(l,m)
!                enddo 
!                imass = idnint( massfrag_tot/mass_bar ) + 1

! Mass spectra of charged fragments
!                if (Bapres(l)) then 
 
!                   mspec(imass,iprint,1) = mspec(imass,iprint,1) + one

!                   write(200+iprint,'(i6,2i7,i5,e16.8,i6,e16.8)') itraj,l,NHe(l),NBa(l),massfrag_tot

!                   print *,' itraj, fragment nb, NHe, NBa, mass = ',itraj,l,NHe(l),NBa(l),massfrag_tot
                   
!                   if ( itraj == N_ttraj ) then
!                      write(200+iprint,*)
!                      write(200+iprint,'(a27,e16.8,a1)') '# Fragment data at time t =',time,'.' 
!                   endif

! Mass spectra of neutral fragments
!                else

!                   mspec(imass,iprint,2) = mspec(imass,iprint,2) + one

!                   write(300+iprint,'(i6,2i7,i5,e16.8)') itraj,l,NHe(l),NBa(l),massfrag_tot

!                   if ( itraj == N_ttraj ) then
!                      write(300+iprint,*) 
!                      write(300+iprint,'(a27,e16.8,a1)') '# Fragment data at time t =',time,'.' 
!                   endif

!                endif

!             enddo

! iprint is incremented for later tests
             iprint = iprint + 1

! Partial deallocation of memory
!!!NH removed             DEALLOCATE(coords,atsymb,mass,NHe,NBa,Bapres,frag,fragsize,mass_frag)
             DEALLOCATE(NHe,NBa,Bapres,frag,fragsize,mass_frag)   !!!NH added

!!!NH removed          else

! We must jump N_tot+2 lines in XYZ movie files and 
!!!NH removed            do l=1,N_tot+2
!!!NH removed               read(120,*)
!!!NH removed            enddo

          endif

       enddo 
!!!NH removed       close(120)
       write(6,*) 'Trajectory ',j,'of data set ',i,'treated.'    
       print *   !!!NH added
       write(401,*)
       write(401,*)

    enddo
    
    close(120)   !!!NH added
    read(110,*)
    write(6,*) 'Data set ',i,'treated.'
 enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                       !!!!
!!!! 5. AVERAGE HISTOGRAMS !!!!
!!!!                       !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Creation of the file names in which the results are stored.

 print *
 print *,' final values for observables, total number of trajectories = ',N_ttraj
 
 do i=1,N_print
    print *,i, time_ref(i), NHehist(i), NHe_excplx(i)
    write(201,*) time_ref(i), NHehist(i)/N_ttraj, NHe_excplx(i)/N_ttraj

!    write(name_out,'("RESULTS/mass-spectrum-ch-",i3.3,".dat")')i
!    open(130,file=trim(name_out),status='unknown') 
!    write(130,*) '# Spectrum ',i,' at time t ~ ',time_ref(i)
!    write(130,'(2a18,a40)') '#        Ion mass','Intensity (tot)'
!    do j=1,imass_max
!       write(130,'(100f18.5)') dble((j-1)*mass_bar),mspec(j,i,1)/dble(N_ttraj)
!    enddo  
!    close(130)

!    write(name_out,'("RESULTS/mass-spectrum-neut-",i3.3,".dat")')i
!    open(130,file=trim(name_out),status='unknown') 
!    write(130,*) '# Spectrum ',i,' at time t ~ ',time_ref(i)
!    write(130,'(2a18,a40)') '#    Neutral mass','Intensity (tot)'
!    do j=1,imass_max
!       write(130,'(100f18.5)') dble((j-1)*mass_bar),mspec(j,i,2)/dble(N_ttraj)
!    enddo  
!    close(130)
 
 enddo

! Final deallocation of memory
!  DEALLOCATE(typesymb,mass_at,N_at,N_step,mspec,time_ref)  
 DEALLOCATE(typesymb,mass_at,N_at,N_step,time_ref)  


 write(6,*)
 write(6,*) 'Job completed. GOOD BYE!'

end program spectradyn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE 1: INCSORT
!
! DESCRIPTION: Returns the configurations, atomic masses and symbols 
!              sorted in such a way that the atomic distance to the COM 
!              increases.
!
! STRUCTURE OF THE SUBROUTINE (INPUTS/OUPUTS):
! nat: number of atoms (input).
! coords1(nat,3): initial coordinates (input).
! mass1(nat): initial mass (input).
! atsymb1(nat): initial atomic symbol (input). 
! coords2(nat,3): final coordinates (output).
! mass2(nat): final mass (output).
! atsymb2(nat): final atomic symbol (output).

 subroutine incsort(nat,coords1,mass1,atsymb1,coords2,mass2,atsymb2)

 implicit none

 INTEGER :: i, j, iconv, itemp
 INTEGER, INTENT(IN) :: nat
 INTEGER, DIMENSION(nat) :: ind 
 DOUBLE PRECISION :: temp
 DOUBLE PRECISION, DIMENSION(nat), INTENT(IN) :: mass1
 DOUBLE PRECISION, DIMENSION(nat), INTENT(OUT) :: mass2
 DOUBLE PRECISION, DIMENSION(nat,3), INTENT(IN) :: coords1
 DOUBLE PRECISION, DIMENSION(nat,3), INTENT(OUT) :: coords2   
 DOUBLE PRECISION, DIMENSION(nat) :: dist
 DOUBLE PRECISION, PARAMETER :: zero=0.0d0
 LOGICAL, DIMENSION(nat) :: lswap  
 CHARACTER(len=3), DIMENSION(nat), INTENT(IN) :: atsymb1
 CHARACTER(len=3), DIMENSION(nat), INTENT(OUT) :: atsymb2

 dist = zero 
 do i=1,nat
    dist(i) = dsqrt( coords1(i,1)**2 + coords1(i,2)**2 + coords1(i,3)**2 ) 
 enddo

 iconv = 1
 do i=1,nat
    ind(i) = i
 enddo
 do while ( iconv /= 0 ) 
    iconv = 0
    lswap = .false.      
    do i=1,nat-1
       if ( .not.lswap(i) ) then
          do j=i+1,nat
             if ( (.not.lswap(j)).and.( dist(i) > dist(j) ) ) then
                temp = dist(i)
                dist(i) = dist(j)
                dist(j) = temp
                itemp = ind(i) 
                ind(i) = ind(j)
                ind(j) = itemp
                lswap(j) = .true.    
                iconv = iconv + 1
             endif
          enddo
       endif
    enddo
 enddo
 do i=1,nat
    coords2(i,1) = coords1(ind(i),1)
    coords2(i,2) = coords1(ind(i),2)
    coords2(i,3) = coords1(ind(i),3)
    mass2(i) = mass1(ind(i))
    atsymb2(i) = atsymb1(ind(i))
 enddo

 end subroutine incsort
