! test potentials and derivatives

implicit none

integer, parameter    :: npts=136

logical               :: lordered=.TRUE.
integer               :: ipt, ireadpot, ireadpotder, iwrte
integer               :: iunit  ! 0 for molecular, 1 for atomic units

real (kind=8), dimension (npts) :: rdist, rdistau, Vpot
real (kind=8), dimension (npts) :: Vpsigan, Vpsiganp,Vpsiganm,dpsigan, dpsignum
real (kind=8), dimension (npts) :: Vpsigspl, Vpsigsplp,Vpsigsplm,dpsigspl, dpsigsplnum
real (kind=8), dimension (npts) :: Vppian, Vppianp,Vppianm,dppian, dppinum
!real (kind=8), dimension (npts) :: Vppispl, Vppisplp,Vppisplm,dppispl, dppisplnum
real (kind=8), dimension (npts) :: worksig, workpi

real (kind=8) :: delr, rmin, rmax, eps, epsau

include 'units.f90'

ireadpot = 10
ireadpotder = 12
iwrte = 20
eps=1.d-10
epsau = eps*ang2bohr

rmin = 1.5d0
rmax = 15.d0
delr = 0.d0
if (npts.gt.1) delr = (rmax-rmin)/(npts-1.d0)

do ipt = 1, npts
   rdist(ipt) = rmin + (ipt-1.d0)*delr
enddo

!!! fitted (hence analytical) potential
iunit = 0
print *,' calling potBapHeXZPADfit'
call potBapHeXZPADfit(Vpot,rdist,npts)
write(iwrte,*)'# He-Ba+ ground potential ZPAD fitted to analytical form '
write(iwrte,*)'#    r(A)         Vground(ZPADan)(cm-1)'
do ipt = 1, npts
   print *,rdist(ipt), Vpot(ipt)
   write(iwrte,'(3(1x,g22.15))') rdist(ipt), Vpot(ipt)
enddo

end program

