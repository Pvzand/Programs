! Calculation of the energy associated with the He-He interaction potential:
! Sum_(i<j) V_He-He(R_ij).
! Forme LM2M2 pour le potentiel de He-He
! R.A. Aziz and M. J. Slaman, J. Chem. Phys., 94 (12), 8047-8053 (1991).

subroutine potheheZPADtot(disthehe,vhehetot,ndist,iunit)

  implicit none
 
! Parameters for the He-He effective (ZPAD) potential for Ne4He100.
!     New parameters for ZPAD potential for Ba+4He from David 
!     Distance de coupure Rcut = 1.6 Angs.
!!! distances in Angstrom, potential in cm-1
!!! iunit = 0 to work in these units;
!!! iunit = 1 to work in atomic units
!

  include './nbrat12.f90'
  include './units.f90'

  integer :: ndist, iunit
  real (kind=8), dimension(ndist) :: disthehe
  real (kind=8) :: vhehetot

  real (kind=8) :: ashmu,alphash,atrmu,btrmu,c6lgmu,c8lgmu,c10lgmu,remu,xemu,demu,betamu
  real (kind=8) :: ash,atr,btr,c6lg,c8lg,c10lg,c12lg,re,xe,de,beta

  parameter(ashmu = -430128.0d0, alphash = 18.742d0, &
      atrmu = 4.67503d0, btrmu = 4.71596d0)
  parameter(c6lgmu = -1.18767d0, c8lgmu = -0.728204d0, &
      c10lgmu = -0.411248d0)
   
  parameter(remu = 4.2d0) 
  parameter(demu=1.45437d0,betamu=-8.30807d0,xemu=1.00088d0)

  integer, parameter :: ndistx=(nbrat2*(nbrat2-1))/2
  integer :: idist
  real (kind=8) :: cre
  real (kind=8), dimension(ndistx) :: xx, vshort,vlong,vtrans,vhehe,vM,vshort1

  logical :: first

!!!  parameter(bohr = 0.529177249d0, convcmua = 4.5563353d-6)


data first /.true./

save first

if (first) then

   select case (iunit)

   case(0) ! molecular units
      re = remu
      cre = 1.d0/re
      ash = ashmu
      atr = atrmu
      btr = btrmu
      c6lg = c6lgmu
      c8lg = c8lgmu
      c10lg = c10lgmu
      de = demu
      beta = betamu
      xe = xemu

   case(1) ! atomic units (input distances in Bohrs, output potential in Hartree)
      print *,'  potential coefficients converted to au'
      re = remu/bohr
      cre = 1.d0/re  
      ash = ashmu*cm2Hartree
      atr = atrmu*bohr
      btr = btrmu/bohr
      c6lg = c6lgmu*cm2Hartree
      c8lg = c8lgmu*cm2Hartree
      c10lg = c10lgmu*cm2Hartree
      xe = xemu*cm2Hartree
      de = demu*cm2Hartree
      beta = betamu*cm2Hartree


   case default
      write(6,*)' >>>>> ERROR STOP in potBapHepPSan <<<<<<'
      write(6,*)' unit case iunit = ',iunit,' invalid '
      write(6,*)' choose 0 for molecular units (cm-1, Angstroms) '
      write(6,*)'   or 1 for atomic units  (Hartree, Bohr)'
   end select
   
   first=.false.
endif

xx(1:ndist) = disthehe(1:ndist)*cre

! Short range potential
vshort1(1:ndist) = ash*dexp(-alphash*xx(1:ndist)*xx(1:ndist))

!Morse potential

vM(1:ndist) = de*(exp(beta*(xx(1:ndist)-xe)))*(exp(beta*(xx(1:ndist)-xe))-2.d0)

!Total short range potential

vshort(1:ndist)=vshort1(1:ndist)+vM(1:ndist)

! Long range potential 
vlong(1:ndist) = c6lg/xx(1:ndist)**6 + c8lg/xx(1:ndist)**8 + c10lg/xx(1:ndist)**10
  
! Switching function
vtrans(1:ndist) = 0.5d0*(1.0d0+dtanh(atr*(disthehe(1:ndist)-btr)))
 
! He-He potential
vhehe(1:ndist) = vshort(1:ndist) + vtrans(1:ndist)*(vlong(1:ndist) - vshort(1:ndist))

vhehetot=sum(vhehe(1:ndist))
return
end subroutine potheheZPADtot

