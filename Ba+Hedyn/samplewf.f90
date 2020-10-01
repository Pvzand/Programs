subroutine samplewf(yloc, nyloc,rhardsphBa,rhardsphHe,lfail)

  implicit none

  include 'nbrat12d.f90'

  integer :: nyloc
  real (kind=8), dimension (nyloc) :: yloc
  real (kind=8) :: rhardsphHe, rhardsphBa
  logical :: lfail

  integer  :: iat, jat, ixyz
  integer :: icount, icountat
  integer, parameter :: icountmax=1000, icountatmax=10000
  real (kind=8) :: rmax, pi, twopi, rshift, theta, phi, proba, probr, rdist
  real (kind=8) :: yrand, rand
  real (kind=8), dimension(ncoord2) :: ytemp
  real (kind=8), dimension(3) :: rdistvec

  logical  :: lfirst
  data lfirst /.true./
  
!!! for random subroutine: TO BE INCLUDED IN ALL PROGRAMS OR SUBROUTINES CALLING RANDOM
  common/randomc/yrand

  if (lfirst) then
     rmax = 2.d0
     pi = dacos(-1.d0)
     twopi = 2.d0*pi
     rshift = 0.d0
     lfirst = .false.
  endif
  icount = 0
  lfail = .false.
  
2000 continue
  icount = icount+1
!  write(6,*) ' icount = ',icount
!  call flush(6)
  if (icount.gt.icountmax) then
     print *,'   in samplewf, already ',icount,' attempts, try decreasing hard sphere radiuses'
     lfail = .true.
     return
  endif
  do iat = 1, nbrat2
     icountat=0
!     print *,' sampling for atom ',iat
1000 call RANDOM(yrand,rand)
     theta = rand*pi
     call RANDOM(yrand,rand)
     phi = rand*twopi
     proba = 1.d0
     probr = 0.d0
     do while (proba.gt.probr)
        call RANDOM(yrand,rand)
        rshift = rand*rmax
        call RANDOM(yrand,rand)
        proba = rand
        call wfhersq(probr,rshift,1,91)
!        write(99,*) rshift, probr
     enddo
     icountat=icountat+1
     if (icountat.gt.icountatmax) go to 2000
     ytemp(3*iat-2) = yloc(ncoord1+3*iat-2) + rshift*dsin(theta)*dcos(phi)
     ytemp(3*iat-1) = yloc(ncoord1+3*iat-1) + rshift*dsin(theta)*dsin(phi)
     ytemp(3*iat) = yloc(ncoord1+3*iat) + rshift*dcos(theta)
     do ixyz = 1, 3
        rdistvec(ixyz) = ytemp(3*(iat-1)+ixyz)-yloc(ixyz)   !!! distance vector to Ba+
     enddo
     rdist = sqrt(dot_product(rdistvec(1:3),rdistvec(1:3)))
     if (rdist.lt.rhardsphBa) go to 1000
     do jat = 1, iat-1
!!! calculate distance between new proposed atom position and previous ones
        do ixyz = 1, 3
           rdistvec(ixyz) = ytemp(3*(iat-1)+ixyz)-ytemp(3*(jat-1)+ixyz)
        enddo
        rdist = sqrt(dot_product(rdistvec(1:3),rdistvec(1:3)))
        if (rdist.lt.rhardsphHe) go to 1000
     enddo
!     write(6,*) 'iat, icountat = ',iat,icountat
!     flush(6)
  enddo
  yloc(ncoord1+1:ncoord12) = ytemp(1:ncoord2)        

!!! check CoM position and angular momentum: done in calling program
!!!
  return
end subroutine samplewf
