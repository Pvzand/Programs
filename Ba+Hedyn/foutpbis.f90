!##########################################################################
!     subroutine called to get ground state energies for a given configuration
!    (positions, momenta)

 subroutine foutpbis(t,y,ndimtot,Etot,Ekin,Epot)

 implicit none

 include 'nbrat12d.f90'
 include 'units.f90'

 integer :: i
 integer :: ndimtot
 integer :: nstep, nstw, istw
 real (kind=8) :: t
 real (kind=8), dimension(ndimtot) :: y
 real (kind=8) :: etot,epot,ekin
 real (kind=8) ::  xmass1,xmass2,xmass,xmassinv

 logical :: lfirst, lwrite
! character (len=13) :: title 
 
! common/writefile/title

 common/masses/xmass1,xmass2,xmass(nbrat12),xmassinv(ncoord12)

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


  return  

 
10000 format(I6,' t=',g15.8,' au, = ',g15.8,' ps,',3x,'etot,ekin,epot=',3(1x,g15.8),' au, =', 3(1x,g15.8),' cm-1')


end subroutine foutpbis

