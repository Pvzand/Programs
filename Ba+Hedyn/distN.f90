!!!     compute all atom1-atom2 distance vectors and modulus squared
!!!    Results dist12,distsq12,dist22,distsq22 passed in common
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!subroutine DISTN(ydyn,distvec12,dist12,distvec22,dist22,nvar,ncoord1,ncoord2,nbrat1,nbrat2)
subroutine DISTN(ydyn,ndimydyn)
     
  implicit none

  include 'nbrat12.f90'

  real (kind=8), dimension(ndimydyn) :: ydyn
  integer                            :: ndimydyn
  
  real (kind=8), dimension(3,nbrat1) :: xyz1
  real (kind=8), dimension(3,nbrat2) :: xyz2
  real (kind=8), dimension(3,nbrat2,nbrat1) :: distvec12
  real (kind=8), dimension(nbrat2,nbrat1) :: dist12
  real (kind=8), dimension(3,nbrat2,nbrat2) :: distvec22
  real (kind=8), dimension(nbrat2,nbrat2) :: dist22
  real (kind=8), dimension(3,ndist22) :: distvec221D
  real (kind=8), dimension(ndist22) :: dist221D

  integer iat1,iat2,jat2,ixyz,idist22

  common/pos/xyz1,xyz2
  common/dist/distvec12,dist12
  common/dist2/distvec22,dist22,distvec221D,dist221D

  ! xyz1(ixyz,iat1)=cartesian coordinate ixyz of (Ba+) atom 1 = iat1
  xyz1 = reshape(ydyn(1:ncoord1),(/ 3, nbrat1 /))
  ! xyz2(ixyz,iat2)=cartesian coordinate ixyz of (He) atom 2 = iat2
  xyz2 = reshape(ydyn(ncoord1+1:ncoord1+ncoord2),(/ 3, nbrat2 /))  

!!! calculate all distance vectors and modulus squared between atoms 2 and atoms 1
  do iat1 = 1, nbrat1
     do ixyz = 1, 3
        distvec12(ixyz,1:nbrat2,iat1) = xyz2(ixyz,1:nbrat2)-xyz1(ixyz,iat1)
     enddo

     do iat2 = 1, nbrat2
        dist12(iat2,iat1) = sqrt(dot_product(distvec12(1:3,iat2,iat1),distvec12(1:3,iat2,iat1)))
     enddo
  enddo
  
!!! calculate all distance vectors and modulus squared between atoms 2
  do iat2 = 1, nbrat2
     do ixyz = 1, 3
        distvec22(ixyz,1:nbrat2,iat2) = xyz2(ixyz,1:nbrat2)-xyz2(ixyz,iat2)
     enddo

     do jat2 = 1, nbrat2
        dist22(jat2,iat2) = sqrt(dot_product(distvec22(1:3,jat2,iat2),distvec22(1:3,jat2,iat2)))
     enddo
  enddo
! atom2-atom2 distances arranged in a 1-D vector
  idist22 = 0
  do iat2 = 1, nbrat2-1
     do jat2 = iat2+1, nbrat2
        idist22 = idist22 + 1
        distvec221D(:,idist22) = distvec22(:,jat2,iat2)
        dist221D(idist22) = dist22(jat2,iat2)
     enddo
  enddo

  return
end subroutine DISTN



