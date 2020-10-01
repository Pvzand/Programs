!!!      SUBROUTINE WPINI 
!!!  to select at random the initial excited state surface
!!!  for propagation of the nuclear coordinates
!!! to give initial values to wavepacket coefficients (modulus and phase)
!!!
!!! iopt = option for initial wavepacket coefficients:
!!! iopt = 1: "adiabatic": initial state ihop selected at random,
!!!                        the coefficient for ihop eigenstate is set to 1
!!!                        (modulus=1, phase=0) and all the others are set to 0
!!! iopt = 2: "white light": all the coefficients modulus are set to 1/sqrt{3}
!!!                          all phases are set to zero.
!!! iopt = 3: "test": initial state ihop read from input file (same line as iopt)
!!!                   the coefficient for ihop eigenstate is set to 1
!!!                  (modulus=1, phase=0) and all the others are set to 0
!!! ADD OTHER OPTIONS FOR TIME-DEPENDENT LASER PULSES !!!
!!!
!!!
!!!  returns ihop in COMMON/dyn/ 
!!!       IF ihop=1, then the nuclear trajectory will start on the lowest excited state surface,
!!!       IF ihop=2, then it will start on the 2nd lowest excited state surface,
!!!	     etc...

!!! returns wave packet coefficients at their place in dynamics variable Y
!!!

SUBROUTINE WPINI_SO(Y0,ndimy,iopt)

  IMPLICIT NONE

  integer :: ndimy, iopt
  real (kind=8), dimension(ndimy) :: Y0

  INCLUDE 'nbrat12.f90'

  integer ihop, ioldhop, ihopinic, nbhop
  
  integer, parameter :: nrand=1
  integer :: iseed
!!! He-He and He-Ba+ potentials are given in molecular units:
!!! convert them to atomic units
  
  integer, parameter :: iunit=1

  real (kind=8) :: coef
  

  INTEGER J,IOUT
!!!  REAL (kind=8), dimension(nvar) :: Y

  real (kind=8), dimension(nrand) :: RANDHOP
!!! real, imaginary part of eigenvectors (for common/DIM/)
double precision, dimension (ndimorb,ndimorb):: eigvcRe, eigvcIm, eigvcoldRe, eigvcoldIm
!!! eigenenergies (for common/DIM/)
double precision, dimension (ndimorb):: eigenE

!!! index of the current DIM eigenvalue (hence surface for classical coordinates)
  common/dyn/ihop, ioldhop, ihopinic, nbhop  

  common/DIM/eigenE,eigvcRe,eigvcIm,eigvcoldRe,eigvcoldIm


  print *,'in wpini_so: iopt, ihop = ',iopt, ihop
  
!!! initialize all phases:
  y0(ncoordmomcoef+1:nvar) = 0.d0
  
  select case (iopt)

  case(1) !!! select one initial state ihop at random
!!!           randhop(1)=random number between 0 and 1
     IOUT = 6
     iseed = 0 !!! to continue sampling from previous call to vranf
     CALL VRANF(RANDHOP,nrand,iseed,IOUT)

!     call RANDOM(iseed,RANDHOP)

     iHOP = IDINT((RANDHOP(1)*NDIMorb)+1)
     y0(ncoordmom+1:ncoordmomcoef) = 1.d-8   ! instead of 0, to avoid numerical problems
     y0(ncoordmom+ihop) = 1.d0

  case(2) !!! all coefficients modulus are set equal to 1/sqrt(nDIMorb)
     !!! initial surface for classical dynamics of nuclear coordinates chosen at random
     coef=1.d0/dsqrt(dfloat(nDIMorb))
     y0(ncoordmom+1:ncoordmomcoef) = coef
     IOUT = 6
     iseed = 0 !!! to continue sampling from previous call to vranf
     CALL VRANF(RANDHOP,nrand,iseed,IOUT)
    
!     call RANDOM(iseed,RANDHOP)
  
     iHOP = IDINT((RANDHOP(1)*NDIMorb)+1)

  case(3) !!! impose one initial state (from input file, ihop, same line as iopt)
     ihop = ihopinic ! reset ihop in case previous trajectory has hopped to another surface
     y0(ncoordmom+1:ncoordmomcoef) = 1.d-8    ! instead of 0, to avoid numerical problems
     y0(ncoordmom+ihop) = 1.d0

  case default
     print *
     print *,' >>>>> ERROR STOP in  wpini_so : <<<<< '
     print *,' iopt = ',iopt,' is not possible yet '
     STOP
  end select

  WRITE (*,*) 'wpini_so : initial state = ',ihop
   
  print *, 'in wpini_so: initial coefficient modulus: '
  print *,y0(ncoordmom+1:ncoordmomcoef)
  print *, 'in wpini_so: initial coefficient phases: '
  print *,y0(ncoordmomcoef+1:nvar)
  
  RETURN
END SUBROUTINE WPINI_SO



       
