subroutine DIMpSOdynanZPAD(eigenSO,ndimeigenSO,iunit)
!!! set up DIM matrix p states + spin-orbit (HeN-Ba+) and diagonalizes it
!!! adiabatic following is made in another subroutine
!!! potentials for DIM p states: ZPAD version fitted to analytical form

!!! iunit=1 to convert potentials AND spin-orbit splitting from (Angs.,cm-1) to au

!!! ifail flag set to 0 if no problem, set to 1 if some vectors are not assigned

  
  implicit none

integer               :: ndimeigenSO, iunit
real (kind=8), dimension(ndimeigenSO) :: eigenSO   ! eigenenergies

logical :: lfirst
integer               :: il, ic

include 'nbrat12.f90'
include 'units.f90'

!!!parameter(ndimp=3, ndimpSO=2*ndimp)
integer               :: iHe, ixyz, idimp, idimpSO, iwrte, iwrteSO, ifail
real (kind=8), dimension(nbrat2) :: rdist
real (kind=8), dimension(3,nbrat2) :: rdistvec
real (kind=8), dimension(ndimorb0,ndimorb0) :: HDIMp
double complex, dimension(ndimorb,ndimorb) :: zHDIMpSO,zHSO
real (kind=8), dimension(ndimorb,ndimorb) :: HSORe, HSOIm, hdimpsoRe, hdimpsoIm

! for subroutine CH (diagonalization of a complex Hermitian matrix)
integer   :: moption, ierr
! real, imaginary part of eigenvectors
double precision, dimension (ndimorb,ndimorb):: eigvcRe, eigvcIm, eigvcoldRe, eigvcoldIm
double precision, dimension (ndimorb):: eigenE
double precision, dimension (ndimorb) :: work1, work2
double precision, dimension (2,ndimorb) :: work3
! (complex) eigenvector (output),
double complex, dimension(ndimorb,ndimorb) :: zeigvc

integer :: iat

common/DIM/eigenE,eigvcRe,eigvcIm,eigvcoldRe,eigvcoldIm

!  atomic (spin-orbit) NIST energy levels (cm-1)
double precision :: Ep3h, Ep1h, DeltaSOcm, DeltaSO    ! atomic spin-orbit (J=3/2, 1/2) energies and splitting
double precision :: gSOh    ! half atomic spin-orbit coupling constant (deduced from splitting)
parameter (Ep3h=21952.404d0,Ep1h=20261.561d0)

double complex :: zzero, zi, z1
parameter(zzero=(0.d0,0.d0),zi=(0.d0,1.d0),z1=(1.d0,0.d0))
  real (kind=8), dimension(3,nbrat2,nbrat1) :: distvec12
  real (kind=8), dimension(nbrat2,nbrat1) :: dist12

common/dist/distvec12,dist12

data lfirst /.true./

save lfirst, HSORe, HSOIm

if (ndimeigenSO.lt.ndimorb) then
   print *,' >>>>> DIMENSION ERROR IN DIMpSOadiabfdyn  <<<<<'
   print *,' ndimeigenSO = ',ndimeigenSO,' should be at least : ndimorb = ',ndimorb
   STOP
endif

ifail = 0

if (lfirst) then
   
   DeltaSOcm = Ep3h - Ep1h
   select case (iunit)
      case(0)  !!! here stick to molecular units
         DeltaSO = DeltaSOcm
      case(1) !!!  here convert to atomic units (Bohr, Hartree)
         DeltaSO = DeltaSOcm*cm2Hartree
    case default
      write(6,*)' >>>>> ERROR STOP in DIMpSOadiabfdynan <<<<<<'
      write(6,*)' unit case iunit = ',iunit,' invalid '
      write(6,*)' choose 0 for molecular units (cm-1, Angstroms) '
      write(6,*)'   or 1 for atomic units  (Hartree, Bohr)'
   end select
         
   gSOh = DeltaSO/3.d0
   zHSO = zzero

! call SOmatpxyz(zHSO,ndimorb,DeltaSO)
  ! zHSO: (px,-1/2) (py,-1/2) (pz,-1/2) (px,1/2) (py,1/2) (pz,1/2)

   zHSO(1,2) = zi
   zHSO(2,1) = -zi
   zHSO(1,6) = -z1
   zHSO(6,1) = -z1
   zHSO(2,6) = -zi
   zHSO(6,2) = zi
   zHSO(3,4) = z1
   zHSO(4,3) = z1
   zHSO(3,5) = zi
   zHSO(5,3) = -zi
   zHSO(4,5) = -zi
   zHSO(5,4) = zi

   zHSO = zHSO*gSOh
   HSORe = dreal(zHSO)
   HSOIm = dimag(zHSO)

   lfirst = .false.
endif


!!! for compatibility
rdist=dist12(:,1)
do iat=1,nbrat2
   rdistvec(:,iat)=distvec12(:,iat,1)
enddo
call DIMpdynanZPAD(rdistvec,rdist,HDIMp,ndimorb0,nbrat2,iunit)

!!! copy HDIMp into the whole (double-sized) H matrix

hdimpsoRe = 0.d0
hdimpsoIm = 0.d0

do ic=1, ndimorb0
   do il=1, ndimorb0
      hdimpsoRe(il,ic)=HDIMp(il,ic)
      hdimpsoRe(ndimorb0+il,ndimorb0+ic)=HDIMp(il,ic)
   enddo
enddo


hdimpsoRe = hdimpsoRe + HSORe
hdimpsoIm = HSOIm

! set Moption to 0 for only eigenvalues, to any other value for eigenvectors and eigenvalues
Moption = 1
ierr = 0
call CH(ndimorb,ndimorb,HDIMpSORe,HDIMpSOIm,eigenSO,Moption,eigvcRe,eigvcIm,work1,work2,work3,ierr)

return
end subroutine DIMpSOdynanZPAD


