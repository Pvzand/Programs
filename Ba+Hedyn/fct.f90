!! ******************************* FCT *********************************

SUBROUTINE FCT(T,Y,DERY,ndimy)
  
!!! equations of motion for HPCG integrator

!!! in this version calls first for diagonalization of the HDIM+HSO matrix
!!! and then separately for adiabatic following
!!! No treatment for degeneracy between Kramer's pairs of eigenvectors
!!! but adiabatic followingfor the treatment of the degeneracy in each pair AND
!!! for reordering of the eigenvectors in case of a crossing between pairs
!!! Assuming a real crossing only occurs punctually
!!! hence there is very little chance to stumble on one,
!!! integration of the current trajectory will be stopped in outp
!!! in case 2 pairs have the same energy (case where Ba+ is isolated, no more He around)

  IMPLICIT NONE
  integer :: ndimy
  real (kind=8) :: T
  real (kind=8), dimension (NDIMy) :: Y,DERY

  real (kind=8), parameter :: splim=0.8d0
      
  include 'nbrat12.f90'
  include 'units.f90'

  INTEGER :: I,J,K,IPOS,IL
  INTEGER :: iat, jat, idist, ip, jp, ip0, ixyz
  integer, parameter :: iunit = 1   !!! to change input potentials from (Angs,cm-1) to atomic units
  integer :: ifail
  real (kind=8), dimension(3,nbrat1) :: xyz1
  real (kind=8), dimension(3,nbrat2) :: xyz2

  real (kind=8) :: DeltaV, dvint
  real (kind=8) :: derytemp
  real (kind=8), dimension(ndimorb,ndimorb)  :: D,DR,DI
!!! for common/dist/
  real (kind=8), dimension (nbrat2,nbrat1) :: dist12
  real (kind=8), dimension (3,nbrat2,nbrat1) ::  distvec12
!!! for common /dist2/
  real (kind=8), dimension (3,nbrat2,nbrat2) ::  distvec22
  real (kind=8), dimension (nbrat2,nbrat2) :: dist22
  real (kind=8), dimension(3,ndist22) :: distvec221D
  real (kind=8), dimension(ndist22) :: dist221D
!!!  for common/dyn/
  integer :: ihop, ioldhop, ihopinic, nbhop            

  real (kind=8), dimension(ndist22) :: V22,dV22

!!! for common potder
  real (kind=8), dimension(nbrat2):: VpSig,VpPi,VpSigder,VpPider
  
  real (kind=8), dimension (ndimorb0,ndimorb0,ncoord12) :: dhnoSOdy
  real (kind=8), dimension (ndimorb,ndimorb,ncoord12) :: dhdy, dhdy_sor, dhdy_soi
  real (kind=8), dimension (ndimorb,ndimorb) :: rdkj_SOR, rdkj_SOI

  real (kind=8), dimension(ndimorb) ::  eigenE, eigenEcomm
  real (kind=8), dimension(ndimorb,ndimorb) ::  vecr_SOR, vecr_SOI
  real (kind=8), dimension(ndimorb,ndimorb) ::  vecrprim_SOR, vecrprim_SOI

!!! for energy conservation test:
  real (kind=8), dimension(ndimorb,ndimorb) :: URe, UIm
  
!!! for common/masses/
  real (kind=8) :: xmass1, xmass2, xmass, xmassinv

  integer :: ittt

  integer :: nstep, nstw
  
  common/masses/xmass1,xmass2,xmass(nbrat12),xmassinv(ncoord12)
!!! index of the current DIM eigenvalue (hence surface for classical coordinates)
  common/dyn/ihop, ioldhop, ihopinic, nbhop
!!! index set to .T. when starting a new trajectory (for reinitializing adiabatic following)
  logical :: linit
  integer :: itraj
  common/inittraj/linit,itraj
  common/nadcplng/ rdkj_SOR, rdkj_SOI
  common/DIM/eigenEcomm,vecr_SOR,vecr_SOI,vecrprim_SOR,vecrprim_SOI
  common/HDIMder/dHdy_SOr, dHdy_SOi
  common/potHeHe/V22
  common/potder/VpSig,VpPi,VpSigder,VpPider
  common/pos/xyz1,xyz2
  common/dist/distvec12,dist12
  common/dist2/distvec22,dist22,distvec221D,dist221D
  common/output/ittt
  common/printinge/nstep, nstw

!  print *
!  print *,' in FCT,  t = ',t
!  print *,' positions ',y(1:ncoord12)
!  print *,' momenta ',y(ncoord12+1:ncoordmom)
!  print *,' coefs mod ',y(ncoordmom+1:ncoordmomcoef)
!  print *,' coefs phases ',y(ncoordmomcoef+1:nvar)
!  call flush(6)
  
!!! calculate all distance vectors and squared modulus for common /dist/      
  call distN(y,nvar)

!       get the eigenvectors vecr_SOR,vecr_SOI
!     CALL POT(Y,nvar)
!!!  call DIMpSOadiabfdynanZPAD(eigenE,ndimorb,iunit,linit)
!!! build Htot = HDIM+HSO and diagonalizes it to get eigenvalues and eigenvectors
  call DIMpSOdynanZPAD(eigenE,ndimorb,iunit)
  eigenEcomm = eigenE   !!! to pass to outp
!!! reordering of eigenvectors in case of a crossing, and adiabatic following for Kramer's pairs
  call adiabfK(eigenE,ndimorb,splim, ifail, linit)
!!!  if (ifail.ne.0) then ???

! Validity of the propagation 
  DO I = ncoordmom+1,ncoordmom+NDIMorb
     IF (Y(I).lt.1.0D-10.AND.Y(I).gt.0.0D0) THEN
        WRITE(*,*)  ' traj ',ITTT,' t = ',T,T*CONVT &
     ,    'Warning Y(',I,') = ',Y(I),' <10**(-10) (hpcg)'
        PRINT *,'******** Coefficients tend to zero ********'
  call flush(6)
     ENDIF
  ENDDO


!!! Method1 : Hellmann-Feynman 

!!!   NONADIABATIC COUPLING VECTOR D(k,j)
!!!  calculated as (Gam**(-1) (dH/dR) Gam)/(E(j)-E(k))
!!!  where dH/dR is calculated in the diabatic basis
!!!  and Gam is the eigenvector matrix (eigenvectors in columns)
!!! (Hellman-Feynman formula)

!!!     COMPUTE ALL dH/dY in the diabatic basis set
!!!  (= basis set in which the DIM matrix is written before diagonalization,
!!!  (i.e., p orbitals (or s,d,p))

  D = 0.d0
  DR = 0.d0
  DI = 0.d0
  DHDY = 0.d0
  DHDY_SOR = 0.d0
  DHDY_SOI = 0.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! OPTIMIZE: the DIM matrix has already been calculated
!!! this is not yet taken into account in DIMpderan
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call DIMpderanZPAD(distvec12,dist12,dhnoSOdy,ndimorb0,nbrat2,ncoord12,iunit)

!!! CHANGE BASIS to the adiabatic basis set in which the DIM matrix has been diagonalized
!!!  D = initial matrix (unchanged), 
!!!  DR, DI = real part, imaginary part of the transformed matrix 
!!!           (complex because spin-orbit is included)
  DO K=1,Ncoord12
     D(1:ndimorb0,1:ndimorb0) = dHnoSOdy(:,:,K)     ! dH_{ij}/dr_k  , 1st diagonal block
     D(ndimorb0+1:ndimorb,ndimorb0+1:ndimorb) = dHnoSOdy(:,:,K)     ! dH_{ij}/dr_k, 2nd diagonal block  
!!!! TEST ENERGY CONSERVATION:
!!!! WRITE A FORTRAN 90 EQUIVALENT to CHBSIS!!!
     
     URe = 0.d0
     do i = 1, ndimorb
        URe(i,i) = 1.d0
     enddo     
     UIm = 0.d0
!     print *
!     print *,' before CHBSIS_SO: dH/dR_k,k = ', K
!     do i = 1, ndimorb
!        print *,(D(i,j),j=1,ndimorb)
!     enddo

     CALL CHBSIS_SO(D,DR,DI,NDIMORB,VECR_SOR,VECR_SOI,NDIMORB,NDIMORB)             
!!!     CALL CHBSIS_SO(D,DR,DI,NDIMORB,URe,UIm,NDIMORB,NDIMORB)
     
!     print *
!     print *,' after CHBSIS_SO: dH/dR'
!     do i = 1, ndimorb
!        print *,' Re:',(DR(i,j),j=1,ndimorb)
!        print *,' Im:',(DI(i,j),j=1,ndimorb)
!     enddo

     
!!! Symmetrization (is it really needed???)
     dHdY_SOr(:,:,K) = 0.5d0*(DR+TRANSPOSE(DR))
     dHdY_SOi(:,:,K) = 0.5d0*(DI-TRANSPOSE(DI))
     
!     print *
!     print *,' after CHBSIS_SO and symmetrization : dH/dR'
!     do i = 1, ndimorb
!        print *,' Re:',(DHdY_SOr(i,j,K),j=1,ndimorb)
!        print *,' Im:',(DHdY_SOi(i,j,K),j=1,ndimorb)
!     enddo

  ENDDO
 

!!!   Equations of motion for positions and momenta
!!! DERY[1..3N] = positions -> p/m
  DERY(1:ncoord12) = Y(ncoord12+1:ncoordmom)*xmassinv
!  print *,' dery for positions '
!  print *, dery(1:ncoord12)
!    call flush(6)

  
!!! DERY[3N+1..6N] = momenta -> -<IHOP|dH/dY(i)|IHOP> for Ba+;
  !!!  -<IHOP|dH/dY(i)|IHOP>  - dV{He-He}/dY(i)  for any He
  DERY(ncoord12+1:ncoordmom) = -DHDY_SOR(IHOP,IHOP,:)
!  print *,' dery for momenta without He-He pot '
!  print *, dery(ncoord12+1:ncoordmom)
!    call flush(6)


  if (nbrat2.gt.1) then
     call derpotHeHeZPADu(dV22,V22,dist221D,ndist22,iunit)
!     print *,' ndist22, V22, dV22 = ',ndist22, V22, dV22
     do ixyz = 1, 3
        idist = 0
!        print *,' ixyz = ',ixyz
        ip0 = ncoord12 + ncoord1 + ixyz
        do iat = 1, nbrat2-1
           ip = ip0 + (iat-1)*3
!           print *,' iat, ip = ',iat,ip
           do jat = iat+1, nbrat2
              idist = idist+1
              jp = ip0 + (jat-1)*3
              dvint = dv22(idist)*distvec221D(ixyz,idist)/dist221D(idist)
!              print *,' idist, dv22,distvec221D, dist221D, dvint = ',idist,dv22(idist), &
!                   distvec221D(ixyz,idist), dist221D(idist), dvint
!              print *,' dery(ip), (jp), before dvint = ',dery(ip), dery(jp)
              dery(ip) = dery(ip) + dvint
              dery(jp) = dery(jp) - dvint
!              print *,' dery(ip), (jp), after dvint = ',dery(ip), dery(jp)
!              call flush(6)
           enddo
        enddo
     enddo
!     print *,' dery for positions after He-He pot (check)'
!     print *, dery(1:ncoord12)
!     print *,' dery for momenta including He-He pot '
!     print *,dery(ncoord12+1:ncoordmom)
!       call flush(6)

   endif

!!! ***  RDkj = (dR/dt).D(k,j)
!!!  If the adiabatic surfaces are degenerate, the coupling
!!!  between them is set to 0 (Kramer's pairs in particular)

!!! For optimization:calculate RDKJ(J,I) by antisymmetrization since DHDY_SO has been symmetrized
   
  RDKJ_SOR = 0.0D0
  RDKJ_SOI = 0.0D0
  DO J=1,NDIMORB-1
     DO I=J+1,NDIMORB
        DeltaV = eigenE(J)-eigenE(I)
        IF (DABS(DeltaV).gt.1.0D-11) THEN  
           RDKJ_SOR(I,J) = dot_product(DERY(1:ncoord12),DHDY_SOR(I,J,1:ncoord12))/DeltaV
           RDKJ_SOI(I,J) = dot_product(DERY(1:ncoord12),DHDY_SOI(I,J,1:ncoord12))/DeltaV
           RDKJ_SOR(J,I) = -dot_product(DERY(1:ncoord12),DHDY_SOR(J,I,1:ncoord12))/DeltaV
           RDKJ_SOI(J,I) = -dot_product(DERY(1:ncoord12),DHDY_SOI(J,I,1:ncoord12))/DeltaV
        ENDIF
     ENDDO
  ENDDO

      
!!!  EQUATIONS OF MOTION for wave packet coefficients (modulus and phase)
!!!   C_i=sqrt(ksi_k)*exp(-i*q_k)

!!!!!!!!!!!!!!! TO BE OPTIMIZED LATER !!!!!!!!!!!!!!!!!
  
  dery(ncoordmom+1:nvar) = 0.d0

!!! DERY[ncoordmom+1...(ncoordmom+ndimorb=ncoordmomcoef)] -> Derivatives of the modulus ksi_k
  DO I=1,NDIMORB
!!! Modulus n_k         
!     if (y(ncoordmom+i).gt.1.d-10) then  
     DO J=1,NDIMORB
        IF (J.ne.I) then
           DERY(I+Ncoordmom) = DERY(I+ncoordmom)  &
                  +(-RDKJ_SOR(I,J)*DCOS(Y(Ncoordmomcoef+I)-Y(Ncoordmomcoef+J)) &
                    +RDKJ_SOI(I,J)*DSIN(Y(Ncoordmomcoef+I)-Y(Ncoordmomcoef+J)))*Y(Ncoordmom+J)
        endif
     enddo
  enddo
!!! DERY[ncoordmomcoef+1... nvar] -> Derivatives of the phases q_k
!!! if the jth-modulus is zero, skip this part and only use V_j/hbar for the time derivative of the jth-phase
  do i=1, ndimorb
     DERY(I+NCOORDMOMCOEF) = eigenE(I) 
     if (dabs(y(ncoordmom+i)).lt.1.d-10) cycle
     derytemp = 0.d0
     do j=1, ndimorb
!        IF (dabs(y(ncoordmom+j)).gt.1.d-10) then
           DERYtemp = DERYtemp &
                + Y(NCOORDMOM+J)*(RDKJ_SOR(I,J)*DSIN(Y(NCOORDMOMCOEF+I)-Y(NCOORDMOMCOEF+J)) &
                                 +RDKJ_SOI(I,J)*DCOS(Y(NCOORDMOMCOEF+I)-Y(NCOORDMOMCOEF+J)))
!        endif
     ENDDO
     DERY(I+NCOORDMOMCOEF) = DERYtemp/Y(NCOORDMOM+I) + eigenE(I) 
  ENDDO

!  if (itraj.eq.2.and.t.lt.20.d0) then
!     write(6,*) 'IN FCT : ', itraj, t
!     write(6,*) 'positions : ', y(1:ncoord12)
!     write(6,*) ' moments : ',y(ncoord12+1:ncoordmom)
!     write(6,*) ' modulus : ',y(ncoordmom+1:ncoordmomcoef)
!     write(6,*) ' phases :',y(ncoordmomcoef+1:nvar)
!     write(6,*) ' dery-positions: ',dery(1:ncoord12)
!     write(6,*) ' dery-moments : ',dery(ncoord12+1:ncoordmom)
!     write(6,*) ' dery-modulus :',dery(ncoordmom+1:ncoordmomcoef)
!     write(6,*) ' dery-phases : ',dery(ncoordmomcoef+1:nvar)
!  endif
!  call flush(6)
  
  RETURN
END SUBROUTINE FCT

