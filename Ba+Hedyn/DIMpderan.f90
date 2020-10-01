subroutine  DIMpderan(xyzdist,rdist,HDIMpder,ndimp,nHe,nxyz,iunit)
!!! build DIM matrix derivatives of Ba+ - He excited state with: 6p sigma and pi potentials
!!! with respect to all the Ba+ and He's (cartesian) coordinates
!!!
!!! This version uses analytical form for the 6p Sigma and Pi potentials
  
!!! iunit = 0 for molecular units (Angstrom, cm-1); 1 for atomic units

! implicit real (kind=8) :: (a-h,o-z)
integer :: ndimp, nHe, nxyz, iunit
double precision, dimension (3,nHe) :: xyzdist
double precision, dimension (nHe) :: rdist
double precision, dimension (ndimp,ndimp,nxyz) :: HDIMpder

logical lfirst
integer :: iHe, il, ic, iat2, ixyz
real (kind=8) :: delV

include 'nbrat12.f90'
double precision, dimension (nDIMorb0,nDIMorb0) :: Unitmtx, HDIMpred, tempMtx, tempMtxp

double precision, dimension (nbrat2) :: rdistsqinv, Vpsig, Vppi, Vpsigder, Vppider

! for checks:
double precision, dimension (nbrat2) :: Vpsigp, Vpsigm, Vppip, Vppim, rdistm, rdistp
!double precision, dimension (3,3) :: test
double precision :: eps

data lfirst /.true./

if (lfirst) then
   ! Unit matrix
   Unitmtx = 0.d0
   do il = 1, ndimp
      Unitmtx(il,il) = 1.d0
   enddo
   lfirst = .false.
endif

! do iHe = 1, nHe
!   rdistsq=xyzdist(1,iHe)**2+xyzdist(2,iHe)**2+xyzdist(3,iHe)**2
!   print *,' rdistsq = ',rdistsq
!   rdist(iHe) = dsqrt(rdistsq)
!    rdistsqinv(iHe) = 1.d0/rdistsq
! enddo
rdistsqinv(1:nHe) = 1.d0/(rdist*rdist)

!!! Calculate Sigma and Pi Ba+-He potentials and derivatives for each Helium
call derpotpotBapHepPSan(Vpsig,Vppi,Vpsigder,Vppider,rdist,nHe,iunit)
! print *
! print *,' in DIMpderan '
! print *,' Vpsig, Vppi = ',Vpsig(1), Vppi(1)
! print *,' Vpsigder, Vppider = ',Vpsigder(1), Vppider(1)
! print *,' check potential derivatives'
! eps = 1.d-6
! rdistp(1) = rdist(1) + eps
! call  potBapHepPSan(Vpsigp,Vppip, rdistp,nHe,iunit)
! print *,' rdistp, Vpsigp, Vppip = ',rdistp(1),Vpsigp(1),Vppip(1)
! rdistm(1) = rdist(1) - eps
! call  potBapHepPSan(Vpsigm,Vppim, rdistm,nHe,iunit)
! print *,' rdistm, Vpsigm, Vppim = ',rdistm(1),Vpsigm(1),Vppim(1)
! print *,' Vpsigder, Vppider, NUMERICAL = ',(Vpsigp(1)-Vpsigm(1))/(2.d0*eps), (Vppip(1)-Vppim(1))/(2.d0*eps)
! print *,' eps = ',eps
! print *

HDIMpder = 0.d0

do iHe = 1, nHe
!!! Calculate the "reduced" contribution of helium atom iHe to the DIM matrix
   HDIMpred = 0.d0
   delV = Vpsig(iHe)-Vppi(iHe)
   do il = 1, 3
      do ic = 1, 3
         HDIMpred(il,ic)=xyzdist(il,iHe)*xyzdist(ic,iHe)
      enddo
   enddo
   do ixyz = 1, 3
      indx = ncoord1+3*(iHe-1)+ixyz
      rho = xyzdist(ixyz,iHe)/rdist(iHe)
      Tempmtx = 0.d0
      Tempmtxp = 0.d0
      HDIMpder(:,:,indx) = Vppider(iHe)*rho*Unitmtx + &
         ((Vpsigder(iHe)-Vppider(iHe))*rho-2.d0*delV*xyzdist(ixyz,iHe)*rdistsqinv(iHe))*rdistsqinv(iHe)*HDIMpred
      TempMtx(ixyz,:) = xyzdist(:,iHe) !!! ex ixyz=1: (x,y,z//0,0,0//0,0,0)
      TempMtxp(:,ixyz) = xyzdist(:,iHe) !!! ex ixyz=1: (x,0,0//y,0,0//z,0,0)
      HDIMpder(:,:,indx) = HDIMpder(:,:,indx) + delV*(TempMtx+TempMtxp)*rdistsqinv(iHe)
   enddo

enddo
!!! Derivatives with respect to Ba+ coordinates
!print *
!print *,' in DIMpder, check derivative / Ba+ coordinates:'

do ixyz = 1, 3
!   test = 0.d0
!   do il = 1, 3
!      do ic = 1, 3
!         do iHe = 1, nHe
!            test(il,ic) = test(il,ic) - HDIMpder(il,ic,ncoord1+3*(iHe-1)+ixyz)
!         enddo
!      enddo
!   enddo
   
   HDIMpder(:,:,ixyz) = -sum(HDIMpder(:,:,ncoord1+ixyz:ncoord12:3),3)
!   do il = 1, 3
!      print *, ' test :',test(il,:)
!      print *,' HDIMpder ',HDIMpder(il,:,ixyz)
!   enddo
   
enddo

return
end subroutine DIMpderan

  
  


