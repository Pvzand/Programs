!!! Pull atom 1 away from equilibrium configuration (with energy E0)
!!! along a random direction
!!!   until reaching the desired potential energy Egoal=E0+Etemp
!
!!!   1st try: x, y, z of atom 1 are incremented by +/- r * DX
!!!   with a random sign and r a random number between 0 and 1



 SUBROUTINE TIRE(Y,ndimy,DX0,IALEA,etemp,E0,egoal)      
 
 IMPLICIT NONE

 include 'nbrat12.f90'
 include 'units.f90'
!INCLUDE 'PARAMETERPOT_so.f'
!INCLUDE 'COMMON_so.f'

 INTEGER :: ndimy,IALEA,IEN,ISORT,I,J,IOUT,ICONV,NCONV
 REAL (KIND=8) :: Eneut,Emin,E0,PRECIS,DIFF,RAND
 REAL (KIND=8) :: etemp, egoal, epot, dx0, dx
 REAL (KIND=8),DIMENSION(NDIMy) :: Y
 REAL (KIND=8),DIMENSION(ncoord12) :: YSHIFT,YINI,YSAVE
 REAL (KIND=8),DIMENSION(2*ncoord1) :: XRAND
 REAL (KIND=8),DIMENSION(ncoord1) :: DELTAR

!DIMENSION Eneut(20),E0(20),Emin(20)

!NH CALL MINE0(Eneut,Emin,E0)

! Critere de precision
  PRECIS = 1.0D-15

! save initial values for y

 DO I=1, ncoord1
    YINI(I) = Y(I)
    YSAVE(I) = 0.0D0
 ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 STEP 1                               !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! counter for the number of attempts to get epot>Egoal
 IEN = 0
! index set to 1 after 5 trials: reset and start from a new sampling
 ISORT = 0
! Parameter of the convergency of the stretching + test integer
 NCONV = 10
 ICONV = 0
! energy to be obtained:
 print *
 egoal = E0 + etemp
 print *,'Initial energy for the neutral cluster : '
 print *, 'E0 = ',E0,' etemp = ',etemp,' egoal = ',egoal,' a.u.'
 print *, 'E0 = ',E0/convcmau,' etemp = ',etemp/convcmau,' egoal = ',egoal/convcmau,' cm-1'
 print *
!      print *,' initial positions:'
!      do i = 1, nbrat12
!         print *,i,(y(3*(i-1)+j),j=1,3)
!      enddo
!  initial value for the potential energy
 print *
 call fpottot2(y,ndimy,epot)
 print *,' check potential energy: ',epot,' au, = ',epot*convaucm,'cm-1'
 if (dabs(epot-E0).gt.1.d-10) then
    print *
    print *,' >>>>> ERROR STOP in tire <<<<< '
    print *,' potential energy of the input configuration = ',epot
    print *,' should be equal to E0 = ',E0
    STOP
 endif

 DO I=1,Ncoord1
    YSHIFT(I) = 0.0D0
    XRAND(I) = 0.0D0
    xrand(ncoord1+i) = 0.d0
    DELTAR(I) = 0.0D0 
 ENDDO 

! PULL RANDOMLY ATOM 1 away FROM its position at THE MINIMUM ENERGY CONFIGURATION
! TO GET Epot REQUIRED

!INITIALIZATION OF THE RANDOM NUMBERS GENERATOR (Marius)
 IOUT = 6
CALL VRANF(RAND,0,IALEA,IOUT)
!   CALL RANDOM(IALEA,RAND)

 DO I=1,NCONV      
!     at each new attempt in dowhile, a new random direction is sampled and DX is multiplied by 2
!     after 5 such attempts, ISORT is set to one and I is incremented, DX is reset  to DX0
  DX = DX0
  DO WHILE (epot.lt.Egoal.AND.ISORT.eq.0)

!C Sample 6 random numbers: 3 for direction, 3 for signs in each direction
  CALL VRANF(XRAND,2*ncoord1,0,IOUT)
!print *,' xrand: ',(xrand(j),j=1,2*ncoord1)

! RANDOM CHOICE OF A VECTOR YSHIFT
!     pull the ATOM until getting above the desired energy
! sign(ar1,arg2) Performs sign transfer:􏰒 if arg2<0 the result is -arg1; if arg2 >= 0 the result is arg1􏰅 􏰂
   DO J=1, ncoord1
     YSHIFT(J) = XRAND(J)*sign(DX,xrand(J+ncoord1)-1.d0)
     Y(J) = Y(J) + YSHIFT(J)
!               print *,j,'yshift, new y:',yshift(j),y(j)
   ENDDO

!total potential energy for the new configuration
 call fpottot2(y,ndimy,epot)
  
!PRINT *,' epot = ',epot,'(E0 = ',E0,' egoal = ',egoal,')'
! CALL FLUSH(6)

  IEN = IEN + 1
  IF (IEN.gt.5) THEN
    PRINT *,'more than 5 attempts -> ',' double the step and restart samplings'
    ISORT = 1   ! forces end of do while and next step in i
    ICONV = ICONV + 1
    DX = DX*2.d0
  ENDIF
 
 ENDDO  ! end do while

   isort = 0
ENDDO   ! end do on i

 WRITE(*,*) 'Total number of infructuous trials before reaching ',' Egoal = ',ICONV ,' final DX = ',DX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                 ETAPE 2                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reinitialisation de IEN
 IEN = 0

 DIFF = epot - Egoal

! Determination of steps along x, y and z used to obtain epot > E0(NBRATDEF)
 DO I=1,Ncoord1
    DELTAR(I) = Y(I) - YINI(I)
 ENDDO

! reach desired energy E0 DICHOTOMIcally

 DO WHILE (DIFF.gt.PRECIS)

    IF (epot.lt.Egoal) THEN   !  pull the atom
       DO I=1, ncoord1
          Y(I) = Y(I) + DELTAR(I)
       ENDDO
    ELSE                       ! push the atom back
       DO I=1, ncoord1
          Y(I) = Y(I) - DELTAR(I)
       ENDDO
    ENDIF

    DO I=1, ncoord1
       DELTAR(I) = DELTAR(I)*0.5D0
    ENDDO

! CALCULates the new POTENTIAL and the mismatch with the desired energy
   CALL FPOTtot2(Y,ndimy,epot)
   DIFF= DABS(epot-Egoal)
! WRITE (*,*) ien, epot,DIFF

    IEN = IEN + 1 
    IF (IEN.gt.100) THEN 
       PRINT *,'Pb : more than 100 dichotomies, goal not reached'
       CALL FLUSH(6) 
       STOP
    ENDIF

 ENDDO

!  print *,' end of tire_dicho_gen4 '
!  print *,'y = ',(y(i),i=1,nbasdef)
!  print *,dsqrt((y(4)-y(1))**2+(y(5)-y(2))**2+(y(6)-y(3))**2)
!  call flush(6)

!  STOP

 RETURN
 END SUBROUTINE TIRE 

