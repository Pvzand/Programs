subroutine SOmatpxyz(zHSO,ndimzHSO,DeltaSO)
  implicit none
  logical :: lfirst
  integer :: ndimzHSO
  double complex, dimension (ndimzHSO,ndimzHSO) :: zHSO
  double complex :: zzero, zi, z1
  real (kind=8) :: DeltaSO, gSOh
  ! useful data
  parameter(zzero=(0.d0,0.d0),zi=(0.d0,1.d0),z1=(1.d0,0.d0))

  gSOh = DeltaSO/3.d0
  zHSO = zzero

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

  return
end subroutine SOmatpxyz

  
  
