! test potentials and derivatives

implicit none

integer, parameter    :: npts=136

logical               :: lordered=.TRUE.
integer               :: ipt, ireadpot, ireadpotder, iwrte
integer               :: iunit  ! 0 for molecular, 1 for atomic units

real (kind=8), dimension (npts) :: rdist, rdistau
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
print *,' calling derpotpotBapHepPSZPADan'
call  derpotpotBapHepPSanZPAD(VpSigan, VpPian, dpSigan,dpPian, rdist,npts,iunit)
write(iwrte,*)'# He-Ba+ p potentials ZPAD fitted to analytical form '
write(iwrte,*)'#    r(A)         VpSig(ZPADan)(cm-1)            VpPi(ZPADan)(cm-1)'
do ipt = 1, npts
   print *,rdist(ipt), Vpsigan(ipt), Vppian(ipt)
   write(iwrte,'(3(1x,g22.15))') rdist(ipt), Vpsigan(ipt), Vppian(ipt)
enddo

!!! numerical derivative:
print *
print *,'  compare analytical vs numerical derivatives '
call  derpotpotBapHepPSanZPAD(Vpsiganp, Vppianp, worksig,workpi,rdist+eps,npts, iunit)
call  derpotpotBapHepPSanZPAD(Vpsiganm, Vppianm, worksig,workpi,rdist-eps,npts, iunit)
dpsignum = (Vpsiganp-Vpsiganm)/(2.d0*eps)
dppinum = (Vppianp-Vppianm)/(2.d0*eps)

!!! analytical derivative
!!!call  derpotpotBapHepPSan(VpSigan, VpPian, dpSigan, dpPian, rdist,npts,iunit)
write(iwrte+1,*)'# He-Ba+ p potentials ZPAD fitted to analytical form and derivatives '
write(iwrte+1,*)'#    r    VpSigZPADan    VpPiZPADan   dVpSigZPADan   dVpSigZPADnum   dVpPiZPADan   dVpPiZPADnum '
write(iwrte+1,*)'#    units: length = Angstroms, energy = cm-1 '

print *,'  r   VSig.an   VPi.an   dVSig.an  dVSig.num  dVPi.an  dVPi.num '
do ipt = 1, npts
   print *, rdist(ipt), Vpsigan(ipt),  Vppian(ipt), dpsigan(ipt),dpsignum(ipt), dppian(ipt), dppinum(ipt)  
   write(iwrte+1,'(7(1x,g22.15))') rdist(ipt), Vpsigan(ipt),  Vppian(ipt), dpsigan(ipt),dpsignum(ipt), dppian(ipt), dppinum(ipt) 
enddo

stop
end program

