! test derivatives of the potentials for DIM p states + spin-orbit (He-Ba+)
!!! using analytical potentials and derivatives for the p Sigma and Pi potentials
!!! iunit = 0 for molecular units (Angstrom, cm-1); 1 for atomic units
implicit none

logical, parameter  :: lderiv=.false.   !!! set to .true. for checking derivatives (choose npts small, lots of prints!)
integer, parameter    :: npts=101
integer               :: nat, nHe, nxyz, ndimp, ndimpSO, il, ic
parameter(nat=2, nHe=nat-1, nxyz=3*nat, ndimp=3, ndimpSO=2*ndimp)
integer               :: iunit
integer               :: iHe, ixyz, ipt, idimp, idimpSO, ifail
integer               :: ireadpotexc, iwrteDIM, iwrteDIMder, iwrteDIMeigE
integer               :: iwrteDIMSO, iwrteDIMSOder, iwrteDIMSOeigE

real (kind=8), dimension(3,nat) :: rat, ratp, ratm  ! cartesian coordinates Ba, He
real (kind=8) :: deg2rad
real (kind=8), dimension(3,nHe) :: xyzdist, xyzdistp,xyzdistm
real (kind=8), dimension(nHe) :: rdist, rdistm, rdistp
real (kind=8), dimension(ndimp,ndimp) :: HDIMp, eigvc
real (kind=8), dimension(ndimp,ndimp,3*nat) :: HDIMpder
real (kind=8), dimension(ndimp,ndimp) :: HDIMpp,HDIMpm
double complex, dimension(ndimpSO,ndimpSO) :: zHDIMpSO,zHSO,zeigvec
real (kind=8), dimension(ndimpSO,ndimpSO) :: zHSORe, zHSOIm, zHDIMpSORe, zHDIMpSOIm
real (kind=8), dimension(ndimp) :: eigen, work
! for subroutine CH (diagonalization of a complex Hermitian matrix)
double precision, dimension(ndimpSO) :: eigenSO
integer   :: moption, ierr
double precision, dimension (ndimpSO,ndimpSO):: zeigvcRe, zeigvcIm
double precision, dimension (ndimpSO) :: work1, work2
double precision, dimension (2,ndimpSO) :: work3
!  atomic (spin-orbit) NIST energy levels (cm-1)
double precision :: Ep3h, Ep1h, DeltaSO
parameter (Ep3h=21952.404d0,Ep1h=20261.561d0)
double complex :: zzero

real (kind=8) :: delr, rmin, rmax, theta, phi, pi,rho
real (kind=8) :: eps

pi = dacos(-1.d0)
deg2rad=pi/180.d0
zzero = dcmplx(0.d0,0.d0)
DeltaSO = Ep3h - Ep1h

!!! small number for numerical derivatives
eps = 1.d-6

!!! iunit = 0 for molecular units (Angstrom, cm-1); 1 for atomic units
iunit = 0

theta = 42.d0*deg2rad
phi = 138.d0*deg2rad
!theta = 0.d0*deg2rad
!phi = 0.d0*deg2rad
rat(1,1) = 0.d0
rat(2,1) = 0.d0
rat(3,1) = 0.d0
! rat(1,1) = -1.d0
! rat(2,1) = 0.5d0
! rat(3,1) = 0.3d0
print *,' Ba+ xyz coordinates :',(rat(ixyz,1),ixyz=1,3)
print *,' nat = ',nat,' nHe = ',nhe

ireadpotexc=11
iwrteDIM = 20
iwrteDIMSO = 21
iwrteDIMeigE = 22
iwrteDIMSOeigE = 23
iwrteDIMder = 30
iwrteDIMSOder = 31

write(iwrteDIM,*) '# He-Ba+ DIM p mtx (no SO), analytic fit of Fausto"s inputs, Angstoms, cm-1'
write(iwrteDIM,*) '#      r      6p xx      6p xy      6p xz      6pyx     6pyy    6pyz    6pzx    6pzy    6pzz'

write(iwrteDIMSO,*) '# He-Ba+ DIM+SO p mtx, analytic fit of Fausto"s inputs, Angstoms, cm-1'
write(iwrteDIMSO,*) '#   r      6p x-x-      6p x-y-      6p x-z-      6p x-x+      6p x-y+      6p x-z+'&
     ,'     6p y-x-      6p y-y-      6p y-z-      6p y-x+      6p y-y+      6p y-z+'&
     ,'     6p z-x-      6p z-y-      6p z-z-      6p z-x+      6p z-y+      6p z-z+'&
     ,'     6p x+x-      6p x+y-      6p x+z-      6p x+x+      6p x+y+      6p x+z+'&
     ,'     6p y+x-      6p y+y-      6p y+z-      6p y+x+      6p y+y+      6p y+z+'&
     ,'     6p z+x-      6p z+y-      6p z+z-      6p z+x+      6p z+y+      6p z+z+'


write(iwrteDIMeigE,*) '# He-Ba+ DIM p eigenvalues (no SO), analytic fit of Fausto"s inputs, Angstoms, cm-1'
write(iwrteDIMeigE,*) '#        r            6p 1          6p 2          6p 3  '

write(iwrteDIMSOeigE,*) '# He-Ba+ DIM+SO p eigenvalues, analytic fit of Fausto"s inputs, Angstoms, cm-1'
write(iwrteDIMSOeigE,*) '#    r          6p 1          6p 2          6p 3     6p4      6p5    6p6'

! if (lderiv) then
! write(iwrteDIMder,*) '# He-Ba+ DIM p mtx deriv (no SO), analytic fit of Fausto"s inputs, Angstoms, cm-1'
! write(iwrteDIMder,*) '#          '

! write(iwrteDIMSOder,*) '# He-Ba+ DIM+SO p mtx deriv, analytic fit of Fausto"s inputs, Angstoms, cm-1'
! write(iwrteDIMSOder,*) '#        '
! endif
   
rmin = 2.d0
! rmin = 2.7d0
rmax = 12.d0
delr = 0.d0
if (npts.gt.1) delr = (rmax-rmin)/(npts-1.d0)

do ipt = 1, npts
   print *
   print *,' ipt = ',ipt
   rho = rmin + (ipt-1.d0)*delr
   rat(1,2) = rho*dsin(theta)*dcos(phi)
   rat(2,2) = rho*dsin(theta)*dsin(phi)
   rat(3,2) = rho*dcos(theta)
   do iHe = 1, nHe
      do ixyz=1, 3
         xyzdist(ixyz,iHe) = rat(ixyz,1+iHe)-rat(ixyz,1)
      enddo
      rdist(iHe) = sqrt(dot_product(xyzdist(:,iHe),xyzdist(:,iHe)))
      print *,' xyzdist, rdist = ',xyzdist,rdist
   enddo
   call DIMpdynanZPAD(xyzdist,rdist,HDIMp,ndimp,nHe,iunit)

   if (lderiv) then
      call DIMpderanZPAD(xyzdist,rdist,HDIMpder,ndimp,nHe,nxyz,iunit)

      do iHe = 1, nHe
         do ixyz = 1, 3
            print *
            print *,'Derivative/He', iHe, ' coord  ixyz  = ',ixyz
            print *,' xyzdist, rdist = ', xyzdist, rdist
            xyzdistp = xyzdist
            rdistp = rdist
            xyzdistp(ixyz,iHe) = xyzdist(ixyz,iHe) + eps
            rdistp(iHe) = sqrt(dot_product(xyzdistp(:,iHe),xyzdistp(:,iHe)))
            print *,' xyzdistp, rdistp = ', xyzdistp, rdistp
            call DIMpdynanZPAD(xyzdistp,rdistp,HDIMpp,ndimp,nHe,ireadpotexc)
            xyzdistm = xyzdist
            rdistm = rdist
            xyzdistm(ixyz,iHe) = xyzdist(ixyz,iHe) - eps
            rdistm(iHe) = sqrt(dot_product(xyzdistm(:,iHe),xyzdistm(:,iHe)))
            print *,' xyzdistm, rdistm = ', xyzdistm, rdistm
            call DIMpdynanZPAD(xyzdistm,rdistm,HDIMpm,ndimp,nHe,ireadpotexc)
            do il = 1, 3
               print *,' num-deriv:',(HDIMpp(il,:)-HDIMpm(il,:))/(2.d0*eps)
               print *,' DIMpder  :',HDIMpder(il,:,3+ixyz)
            enddo
         enddo
      enddo

      print *
      print *,' check derivatives for Ba+ '
      do ixyz = 1, 3
         ratp = rat
         ratm = rat
         print *,' Ba+ ixyz  = ',ixyz
         print *,' xyzdist, rdist = ', xyzdist, rdist
         ratp(ixyz,1) = rat(ixyz,1) + eps
         do iHe = 1, nHe
            xyzdistp(:,iHe) = ratp(:,iHe+1)-ratp(:,1)
            rdistp(iHe) = sqrt(dot_product(xyzdistp(:,iHe),xyzdistp(:,iHe)))
         enddo
         print *,' xyzdistp, rdistp = ', xyzdistp, rdistp
         call DIMpdynanZPAD(xyzdistp,rdistp,HDIMpp,ndimp,nHe,ireadpotexc)
         ratm(ixyz,1) = rat(ixyz,1) - eps
         do iHe = 1, nHe
            xyzdistm(:,iHe) = ratm(:,iHe+1)-ratm(:,1)
            rdistm(iHe) = sqrt(dot_product(xyzdistm(:,iHe),xyzdistm(:,iHe)))
            print *,' xyzdistm, rdistm = ', xyzdistm, rdistm
         enddo
         call DIMpdynanZPAD(xyzdistm,rdistm,HDIMpm,ndimp,nHe,ireadpotexc)
         do il = 1, 3
            print *,' num-deriv:',(HDIMpp(il,:)-HDIMpm(il,:))/(2.d0*eps)
            print *,' DIMpder  :',HDIMpder(il,:,ixyz)
         enddo
      enddo
   endif

! diadonalization of HDIMp
   ifail=0
   call f02abf(HDIMp,ndimp,ndimp,eigen,eigvc,ndimp,work,ifail)
   print *,' diagonalization of HDIMp: ifail = ',ifail
   print *,(eigen(idimp),idimp=1,ndimp)
   write(iwrteDIMeigE,'(4(1x,g15.8))') rdist, (eigen(idimp),idimp=1,ndimp)


! copy HDIMp into the whole (double-sized) H matrix

   zHDIMpSORe = 0.d0
   zHDIMpSOIm = 0.d0

   do ic=1, ndimp
      do il=1, ndimp
         zHDIMpSORe(il,ic)=HDIMp(il,ic)
         zHDIMpSORe(ndimp+il,ndimp+ic)=HDIMp(il,ic)
      enddo
   enddo

   call SOmatpxyz(zHSO,ndimpSO,DeltaSO)      ! Matrice Spin-Orbite 
   zHSORe = dreal(zHSO)
   zHSOIm = dimag(zHSO)

   zHDIMpSORe = zHDIMpSORe + zHSORe
   zHDIMpSOIm = zHSOIm

! set Moption to 0 for only eigenvalues, to any other value for eigenvectors and eigenvalues
   Moption = 1
   ierr = 0
   call CH(ndimpSO,ndimpSO,zHDIMpSORe,zHDIMpSOIm,eigenSO,Moption,zeigvcRe,zeigvcIm,work1,work2,work3,ierr)    ! diagonalisation de la matrice (DIM+SO)
   print *,'  CH-diagonalization of zHDIMpSO: ierr = ', ierr
   print *,(eigenSO(idimp),idimp=1,ndimpSO)
   write(iwrteDIMSOeigE,'(7(1x,g15.8))') rdist, (eigenSO(idimp),idimp=1,ndimpSO)
enddo

stop
end program

