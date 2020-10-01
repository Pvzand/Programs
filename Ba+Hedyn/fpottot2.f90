!!! calculate total potential energy for atom1_n1-atom2_n2
!!! from the y array of variables used in the propagation program
!!! by calling the appropriate subroutines
!!! input: y = atom cartesian coordinates in Bohr (atomic units)
!!! output: epot = potential energy in Hartree (atomic units)

!!! this version uses the (fitted) analytic version of the ZPAD He-Ba+ potential

SUBROUTINE fpottot2(y,ndimy,epot)
      
  IMPLICIT NONE

  include 'nbrat12.f90'
  include 'units.f90'

  INTEGER :: ndimy
  integer, parameter :: ndim=3
  real (kind=8), dimension (ndimy) :: y
  real (kind=8), dimension (ncoord1) :: xyzBap
  real (kind=8), dimension (ndim,nbrat2) :: xyzHe
  integer    iat,i
  REAL*8 :: epot, epot11, epot22, epot12


  if (nbrat1.gt.1) then
     print *,' ***** ERROR STOP in fpottot2 ***** '
     print *,' nbrat1 = ',nbrat1
     print *,'   so far only nbrat1=1 is provided '
     stop
  endif

  xyzBap = y(1:ndim)*bohr
  xyzHe=reshape(y(ncoord1+1:ncoord12)*bohr, (/ 3, nbrat2 /))
 
!  print *,' xyzBap(Angs):',(xyzBap(i),i=1,ndim)
!  do iat = 1, nbrat2
!     print *,'xyzHe',iat,' (Angs):',(xyzHe(i,iat),i=1,3)
!  enddo
  call VBapHeNXtot2(xyzBap,ndim,xyzHe,nbrat2,epot) ! works in Angstrom and cm-1
!  print *,' epot = ',epot
  epot = epot*convcmau

  RETURN
END SUBROUTINE fpottot2


