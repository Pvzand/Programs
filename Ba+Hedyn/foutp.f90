!##########################################################################
!     subroutine called each time fhpcg has completed one step
!   for output purposes: check total energy conservation, dissociation, 
!   write various variables to output files...
!   set prmt(5)=1 for exiting propagation      
 subroutine foutp(t,y,dery,ndimtot,prmt)

 implicit none

 include 'nbrat12.f90'
 include 'units.f90'

 integer :: i
 integer :: ndimtot
 integer :: nstep, nstw, istw
 real (kind=8) :: t
 real (kind=8), dimension(ndimtot) :: y,dery
 real (kind=8), dimension(5) :: prmt
 real (kind=8) :: etot,epot,ekin,etotini,detotx
 real (kind=8) ::  xmass1,xmass2,xmass,xmassinv
 logical :: lfirst, lwrite
 character (len=7) :: title

 common/masses/xmass1,xmass2,xmass(nbrat12),xmassinv(ncoord12)
 common/printingg/nstep, nstw
 common/econservg/detotx
 common/writefile/title
 data etotini/0.0d0/

 data lfirst/.true./

 save etotini

 
 if (lfirst) then
    OPEN(24,FILE=title//'.ftraj',STATUS='UNKNOWN')
    lfirst = .false.
 endif
 
!     initialization
 lwrite=.false.
 if (dabs(t-prmt(1)).lt.1.d-8) then !  called for starting values, no step
    nstep = 0
    istw = 0
    lwrite = .true.
 else                 ! count steps for checking and printing
    istw = istw+1
    nstep = nstep + 1
    if (istw.eq.nstw) then
       lwrite = .true.
       istw = 0
    endif
 endif
 prmt(5) = 0.d0
      
!  Potential energy:     
!!      call fpottot(y,ndimtot,epot)
 call fpottot2(y,ndimtot,epot)

! Kinetic energy
 ekin = 0.0d0
 do i=1, ncoord12
    ekin = ekin + 0.5d0*xmassinv(i)*y(ncoord12+i)**2
 enddo

 etot = epot + ekin
!!      print *,' nstep, t, epot, ekin, etot : ',nstep,t,epot,ekin,etot
!!      print *,' ps, cm-1 ',t*convt,epot*convaucm,ekin*convaucm,
!!     & etot*convaucm
!!     write(6,10000) nstep, t, t*convt, etot, ekin, epot,
!!    &            etot*convaucm, ekin*convaucm, epot*convaucm

 if (nstep.eq.0) etotini=etot

!  conditions for setting up flag for error return in fhpcg
  if (dabs(etot-etotini).gt.detotx) prmt(5)=1.d0
!      if (nstep.gt.1) prmt(5) = 1.d0

   if (lwrite) then
     open(19,file=title//'.etotneut',status='unknown')
!     &          , access='append') 
     write(19,*) t*convt,epot*convaucm,ekin*convaucm,etot*convaucm
  call flush(19)
  write(24,'(19(1x,g15.8))') t*convt,(y(i)*Bohr,i=1,ncoord12),(y(i),i=ncoord12+1,2*ncoord12)
  endif

  return  

 
10000 format(I6,' t=',g15.8,' au, = ',g15.8,' ps,',3x,'etot,ekin,epot=',3(1x,g15.8),' au, =', 3(1x,g15.8),' cm-1')


end subroutine foutp
