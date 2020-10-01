subroutine  DIMpdynanZPAD(xyzdist,rdist,HDIMp,ndimp,nHe,iunit)
!!! build DIM matrix for Ba+ - He excited state potentials: 6p sigma and pi
!!! using analytical form for the 6p Sigma and Pi potentials
!!! iunit = 0 for molecular units (Angstrom, cm-1); 1 for atomic units
!
! implicit real (kind=8) :: (a-h,o-z)
double precision, dimension (3,nHe) :: xyzdist
double precision, dimension (nHe) :: rdist
double precision, dimension (ndimp,ndimp) :: HDIMp
integer :: iunit
logical first
integer :: ndimp, nHe, iHe, iread, il, ic
real (kind=8) :: delV

include 'nbrat12.f90'

! double precision, dimension (nbrat2) :: rdistsqinv, Vpsig, Vppi
double precision, dimension (nbrat2) :: rdistsqinv, Vpsig, Vppi, VpSigder, VpPider

! do iHe = 1, nHe
!   rdistsq=xyzdist(1,iHe)**2+xyzdist(2,iHe)**2+xyzdist(3,iHe)**2
!   print *,' rdistsq = ',rdistsq
!   rdist(iHe) = dsqrt(rdistsq)
!    rdistsqinv(iHe) = 1.d0/rdistsq
! enddo
rdistsqinv(1:nHe) = 1.d0/(rdist*rdist)

!call potBapHepPSan(Vpsig,Vppi,rdist,nHe,iunit)
call derpotpotBapHepPSanZPAD(Vpsig,Vppi,VpSigder,VpPider,rdist,nHe,iunit)
! print *,' in DIMpdynan, potentials from potBapHepPSan '
! print *,'rdist, Vpsig, Vppi = ',rdist(1),Vpsig(1), Vppi(1),VpSigder(1), VpPider(1)

HDIMp = 0.d0
do iHe = 1, nHe
   delV = Vpsig(iHe)-Vppi(iHe)
   do il = 1, ndimp
      do ic = il, ndimp
         HDIMp(il,ic)=HDIMp(il,ic)+xyzdist(il,iHe)*xyzdist(ic,iHe)*delV*rdistsqinv(iHe)
      enddo
      HDIMp(il,il) = HDIMp(il,il) + Vppi(iHe)
   enddo
enddo

! complete by symmetry

do il = 1, ndimp-1
   do ic = il+1, ndimp
      HDIMp(ic,il) = HDIMp(il,ic)
   enddo
enddo

!do  il = 1, ndimp
!   print *,(HDIMp(il,ic),ic=1,3)
!enddo

return
end subroutine DIMpdynanZPAD

