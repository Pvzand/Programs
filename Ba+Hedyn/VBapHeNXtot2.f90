subroutine VBapHeNXtot2(xyzBap,nxyz,xyzHe,natHe,epot)

!!! Analytic version of the ZPAD potential (resulting from a fit)
!!! xyzBap = Ba+ cartesian coordinates (in Angstrom)
!!! xyzHe(ixyz,iHe) = iHe-th He cartesian coordinates ixyz (1: x; 2: y; 3: z) in Angstrom
!!! output potential energy epot in cm-1
!!! = sum of He-Ba+(X) and He-He interactions
  
  
  integer nxyz, natHe
  real (kind=8), dimension(nxyz) :: xyzBap
  real (kind=8), dimension(nxyz,natHe) :: xyzHe
  real (kind=8) :: epot, vhehetot

  integer :: iat, jat, ndist

  integer, parameter :: nathex=100, ndistx=(nathex*(nathex-1))/2
  real (kind=8), dimension(natHex) :: pot12
  real (kind=8), dimension(nathex) :: dist12
  real (kind=8), dimension(ndistx) :: dist22

  integer :: iunit

  do iat = 1, natHe
     dist12(iat)=dot_product(xyzHe(1:3,iat)-xyzBap,xyzHe(1:3,iat)-xyzBap)
  enddo
  dist12(1:natHe) = dsqrt(dist12(1:natHe))

  call potBapHeXZPADfit(pot12,dist12,natHe)

  
  epot = sum(pot12(1:natHe))
!  write(52,*) dist12(1), epot

  ndist = 0
  do iat = 1, natHe-1
     do jat = iat+1, natHe
        ndist = ndist + 1
        dist22(ndist) = dot_product(xyzHe(1:3,iat)-xyzHe(1:3,jat),xyzHe(1:3,iat)-xyzHe(1:3,jat))
     enddo
  enddo

  dist22(1:ndist) = dsqrt(dist22(1:ndist))

  iunit=0
  call potheheZPADtot(dist22,vhehetot,ndist,iunit)

  epot = epot + vhehetot


  return

end subroutine VBapHeNXtot2

        
