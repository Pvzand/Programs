subroutine globveladj(y,dery,ndimy,ihoptest,lforbid)

!!! Test if a hop to new adiabatic surface ihoptest is energetically allowed:
!!!  - if not, set lforbid to .true. and returns
!!!  - if yes, set lforbid to .false. and adjust Ba+ and He-s momenta to conserve energy.
!!! The momentum adjustment for coordinate R_j is made along the d[V(ihoptest)-V(ihop)]/dR_j vector:
!!! Pnew(j)-P(j) = a d[V(ihoptest)-V(ihop)]/dR_j
!!! The a coefficient is determined by energy conservation:
!!! V(ihoptest) + Sum_j Pnew(i)**2/(2M_j) = V(ihop) + Sum_j P(i)**22/(2 M_j)
!!! This results in a 2nd degree equation in a
!!! If its discriminant is <0, there is no solution (hop "classically forbidden")
!!! If it is >0 there are 2 solutions: the one with the smallest absolute value is selected
  
  implicit none

  include 'nbrat12.f90'
  include 'units.f90'

  logical :: lforbid
  integer :: ndimy, ihoptest
  real (kind=8), dimension(ndimy) :: y, dery

!!! for common /DIM/: eigenvalue; eigenvectors (real, imaginary part); eigenvectors at previous time step
  real (kind=8), dimension(ndimorb) ::  eigenEcomm
  real (kind=8), dimension(ndimorb,ndimorb) ::  vecr_SOR, vecr_SOI
  real (kind=8), dimension(ndimorb,ndimorb) ::  vecrprim_SOR, vecrprim_SOI
  common/DIM/eigenEcomm,vecr_SOR,vecr_SOI,vecrprim_SOR,vecrprim_SOI
!!! for common/masses/ : Ba+ mass, He mass, mass for each atomic coordinate, 1/mass for each atomic coordinate
  real (kind=8) :: xmass1, xmass2, xmass, xmassinv
  common/masses/xmass1,xmass2,xmass(nbrat12),xmassinv(ncoord12)
!!! for common/dyn/: current, previous, and initial index for the surface on which classical dynamics is run
  integer :: ihop, ioldhop, ihopinic, nbhop             
  common/dyn/ihop, ioldhop, ihopinic, nbhop
!!! for common/Eclass/: total energy, kinetic energy, potential energy, sum of He-He potential energies
!!! (energies calculated in subroutine outp)
  real (kind=8) :: Etot,Ekin,Epot,EpotHe
  common/Eclass/Etot,Ekin,Epot,EpotHe
!!! for common/HDIMder/: derivative of HDIM in adiabatic basis (Re,Im)
  real (kind=8), dimension (ndimorb,ndimorb,ncoord12) :: dHdy_SOr, dHdy_SOi
  common/HDIMder/dHdy_SOr, dHdy_SOi


  real (kind=8) :: deltav, discr, acoef, bcoef, coef, coefp, coefm
  real (kind=8), dimension(ncoord12) :: deltap

!!! Calculate the potential energy difference between current (ihop) adiabatic surface
!!! and the one to which it is proposed to hop (ihoptest)
  deltav = eigenEcomm(ihoptest) - eigenEcomm(ihop)
  write(*,*) 'Energy gap for hop: deltav = ',deltav
  call flush(6)
  if (Ekin.lt.deltav) then
     print *,' not enough energy to hop '
     print *,' Ekin = ',Ekin
     lforbid=.true.
     return
  endif

  deltap = dHdY_SOr(ihoptest,ihoptest,:)-dHdY_SOr(ihop,ihop,:)
  bcoef = 0.5d0*sum(y(ncoord12+1:ncoordmom)*deltap*xmassinv)
  acoef = 0.5d0*sum(deltap**2*xmassinv)
  discr=bcoef**2-acoef*deltav
  if (discr.lt.0.d0) then   !!! here there is no solution
     print *,' not enough energy to hop with this velocity adjustment '
     print *,' discr = ',discr
     lforbid=.true.
     return
  else !!! here there are 1 or 2 solutions: choose the one with the smaller absolute value
     lforbid = .false.
     coefp = (-bcoef+dsqrt(discr))/acoef
     coefm = (-bcoef-dsqrt(discr))/acoef
     if (dabs(coefp).le.dabs(coefm)) then
        coef=coefp
     else
        coef=coefm
     endif
!!!   calculate new moments and their time-derivatives
     y(ncoord12+1:ncoordmom) = y(ncoord12+1:ncoordmom) + coef*deltap
     dery(ncoord12+1:ncoordmom) = -dHdY_SOR(ihoptest,ihoptest,:)
  endif
      
  return
end subroutine globveladj

