!!! calculate -gradient_r(Vtot)
!!! from the y array of variables used in the propagation program
!!! by calling the appropriate subroutines
!!! input: y = atom cartesian coordinates in Bohr (atomic units)
!!!            then cartesian momenta (au)
!!! output: y unchanged;
!!!       dery(ncoord12+i)=-grad_i(Vtot) for i=1, ncoord12
!!!       with Vtot = total potential energy
!!!       and -grad_i = d/(dy(i))
!!! potential energy in Hartree (atomic units)

SUBROUTINE derfpottot(y,dery,ndimtot)
      
  IMPLICIT NONE

  include 'nbrat12.f90'
  include 'units.f90'

  INTEGER :: iat, jat, ndist, ndimtot, idist, ip, jp, ip0
  integer   ireadpot12, i, ixyz
  real (kind=8), dimension (ndimtot) :: y, dery
  real (kind=8), dimension (ncoord1) :: xyzBap
  real (kind=8), dimension (3,nbrat2) :: xyzHe,r12
  real (kind=8), dimension (nbrat2) :: dist12, dv12
  real (kind=8), dimension ((nbrat2*(nbrat2-1))/2) :: dist22, dv22, V22
  real (kind=8), dimension (3,(nbrat2*(nbrat2-1))/2) :: r22
  REAL*8 :: dvint


  if (nbrat1.gt.1) then
     print *,' ***** ERROR STOP in fpottot ***** '
     print *,' nbrat1 = ',nbrat1
     print *,'   so far only nbrat1=1 is provided '
     stop
  endif

  xyzBap = y(1:ncoord1)*bohr   ! atom 1 (Ba+) cartesian coordinates in Angstroms
  xyzHe=reshape(y(ncoord1+1:ncoord12)*bohr, (/ 3, nbrat2 /))  ! atoms 2 (He's) cartesian coordinates in Angstroms
  
  dery(ncoord12+1:2*ncoord12) = 0.d0

  ! atom1-atom2's distances (Ba+ - He's)  (in Angstroms)
 
  do iat = 1, nbrat2
     r12(1:3,iat) = xyzHe(1:3,iat)-xyzBap
     dist12(iat)=dot_product(r12(1:3,iat),r12(1:3,iat))
  enddo
  dist12(1:nbrat2) = dsqrt(dist12(1:nbrat2))

  call deranpotBapHeXZPAD(dv12,dist12,nbrat2)   ! dV12/dr12 in cm-1/Angs

  dv12(1:nbrat2)=dv12(1:nbrat2)*dVdrcmAng2au   ! dv12 is now in atomic units

  do ixyz=1,3
     ip0 = ncoord12+ixyz
     dery(ip0) = 0.d0
     do iat = 1, nbrat2
        ip = ip0+ncoord1+(iat-1)*3
        dery(ip) = -dv12(iat)*r12(ixyz,iat)/dist12(iat)
        dery(ip0) = dery(ip0) - dery(ip)
     enddo
  enddo
  
! atom2-atom2 distances  
  ndist = 0
  if (nbrat2.gt.1) then
     do iat = 1, nbrat2-1
        do jat = iat+1, nbrat2
           ndist = ndist + 1
           r22(1:3,ndist) = xyzHe(1:3,iat)-xyzHe(1:3,jat)
           dist22(ndist) = dot_product(r22(1:3,ndist),r22(1:3,ndist))
        enddo
     enddo
     dist22(1:ndist) = dsqrt(dist22(1:ndist)) ! dist22 is in Angstrom

     call derpotHeHeZPAD(dv22,V22,dist22,ndist) ! works in Angstrom and cm-1
     dv22(1:ndist)=dv22(1:ndist)*dVdrcmAng2au   ! dv22 is now in atomic units

     do ixyz = 1, 3
        idist = 0
        ip0 = ncoord12 + ncoord1 + ixyz
        do iat = 1, nbrat2-1
           ip = ip0 + (iat-1)*3 
           do jat = iat+1, nbrat2
              idist = idist+1
              jp = ip0 + (jat-1)*3
              dvint = dv22(idist)*r22(ixyz,idist)/dist22(idist)
              dery(ip) = dery(ip) - dvint
              dery(jp) = dery(jp) + dvint
           enddo
        enddo
     enddo

  endif
  
  RETURN
END SUBROUTINE derfpottot

