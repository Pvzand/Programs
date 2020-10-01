!!! version for propagating Re(dk), Im(dk) for the coefficients
!!! rather than mod(ck), phase(ck)
!!! where dk = ck exp(i int_0^t Vk/hbar dt'), Vk being the k-th eigenvalue attached to eigenvector k
!!! This form should avoid describing many unnecessary oscillations compared to Re(ck), Im(ck),
!!! while avoiding numerical undeterminacy when a coefficient tends to zero with mod(ck), phase(ck)
!!! Added energy phases as variables: phik=int_0^t Vk/hbar dt'
  integer, parameter :: nbrat1=1, nbrat2=1000, nDIMorb0=3,nDIMorb=2*nDIMorb0,nDIMorb2 = 2*nDIMorb
  integer :: nbrat12
  integer :: ncoord1, ncoord2, ncoord12
  integer :: ncoordmom1, ncoordmom, ncoordmomcoefr, ncoordmomcoef, nvar
  integer :: ndist22
!ccc  nbrat1 = number of atom1 (Ba+) atoms
!ccc  nbrat2 = number of atom2 (He) atoms
!ccc  nDIMorb = number of basis functions for the electronic wave packet (including SO)
!ccc  = size of the DIM matrix
!ccc  = number of coefficients modulus = number of coefficients phases
  !      parameter(nbrat1=1,nbrat2=2,nDIMorb=6)
!!! nDIMorb2 = number of coefficients modulus + phases
  parameter (nbrat12=nbrat1+nbrat2)   ! total number of atoms
  parameter (ndist22=max(((nbrat2*(nbrat2-1))/2),1))  ! number of atom2-atom2 distances

!ccc   variables:
!ccc 1 to 3*nbrat1 : cartesian coordinates for type 1 atoms (atom1=Ba)
!ccc (3*nbrat1+1) to 3*(nbrat1+nbrat2) : cartesian coordinates for type 2 atoms (atom2=He)
!ccc 3*nbrat12+1 to 3*nbrat12+3*nbrat1 : momenta for type 1 atoms
!ccc 3*nbrat12+3*nbrat1+1 to 6*nbrat12 : momenta for type 2 atoms
!ccc 6*nbrat12+1 to 6*nbrat12+nDIMorb  : modulus of the wave packet coefficients
!ccc 6*nbrat12+nDIMorb+1 to 6*nbrat12+2*nDIMorb : phase of the wave packet coefficients
  parameter(ncoord1=3*nbrat1,ncoord2=3*nbrat2)
  parameter(ncoord12=ncoord1+ncoord2) ! number of cartesian coordinates
  parameter (ncoordmom1=ncoord12+ncoord1,ncoordmom=2*ncoord12) ! number of coordinates and momenta
  parameter (ncoordmomcoef=ncoordmom+nDIMorb2) ! number of coords, mom., and coeff.modulus
  parameter (nvar=ncoordmomcoef+nDIMorb) ! total number of variables
