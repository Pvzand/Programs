subroutine comJ(pos,pmom,ncoord)
!!! sets the origin at the center of mass
!!! and angular momentum to zero using the same idea as Jellinek and Li, PRL 62, 241 (1989).
!!!
  implicit none
  integer ncoord
  real (kind=8), dimension(ncoord) :: pos, pmom

  include 'nbrat12d.f90'
  include 'units.f90'

  logical :: lfirst
  integer :: ixyz, jxyz, iat, ipos
  real (kind=8) :: xmasstot, angmomx, angmomy, angmomz, Jtot, Erot, Ekin, TCM
  real (kind=8) :: temp
  real (kind=8), dimension(3) :: rcm, pcm, angmom, omega
  real (kind=8), dimension(3,3) :: xinertM
  
  !!! for CoM position and total angular momentum
  real (kind=8) :: xmass1,xmass2,xmass, xmassinv
  common/masses/xmass1,xmass2,xmass(nbrat12),xmassinv(ncoord12)
  data lfirst /.true./

  if (lfirst) then
     xmasstot = dfloat(nbrat1)*xmass1 + dfloat(nbrat2)*xmass2 ! total mass
     print *,' Total mass = ',xmasstot,' amu, = ',xmasstot/convgrau,' g/mol'
     call flush(6)
     lfirst = .false.
  endif

!!! center of mass position:
  do ixyz = 1, 3
     rcm(ixyz) = sum(pos(ixyz:ncoord1:3))*xmass1 + sum(pos(ncoord1+ixyz:ncoord12:3))*xmass2
     pcm(ixyz) = sum(pmom(ixyz:ncoord12:3))
  enddo
  rcm = rcm/xmasstot
  print *
  print *,' Center of mass position : ',rcm(1:3),' (coordinates in Bohr)'
  print *,' Center of mass (total) momentum : ', pcm(1:3), ' in a.u.'
  Ekin = 0.5d0*(sum(pmom(1:ncoord1)**2)/xmass1+sum(pmom(ncoord1+1:ncoord12)**2)/xmass2)
  print *,' Total kinetic energy = ',Ekin,' Hartree, = ',Ekin*Hartree2cm,' cm-1'

  do ixyz = 1, 3
     pos(ixyz:ncoord12:3) = pos(ixyz:ncoord12:3) - rcm(ixyz)
     pmom(ixyz:ncoord1:3) = pmom(ixyz:ncoord1:3) - pcm(ixyz)*xmass1/xmasstot
     pmom(ncoord1+ixyz:ncoord12:3) = pmom(ncoord1+ixyz:ncoord12:3) - pcm(ixyz)*xmass2/xmasstot
  enddo
  print *,' new origin set at the center of mass'
  call flush(6)
  Ekin =  0.5d0*(sum(pmom(1:ncoord1)**2)/xmass1+sum(pmom(ncoord1+1:ncoord12)**2)/xmass2)
  print *,' Total kinetic energy in the CoM frame = ',Ekin,' Hartree, = ',Ekin*Hartree2cm,' cm-1'
  TCM = 0.5d0*sum(pcm(1:3)**2)/xmasstot
  print *,' CoM kinetic energy = ', TCM,' Hartree, = ',TCM*Hartree2cm,' cm-1'
  print *,' total kinetic energy (for a check) = ',Ekin+TCM,' Hartree, = ',(Ekin+TCM)*Hartree2cm,' cm-1'
  

!!! total angular momentum
!!! Jtot = Sum_alpha R_alpha x P_alpha
  angmom = 0.d0
  angmomx = 0.d0
  do ixyz = 1, ncoord12, 3
     angmomx = angmomx + pos(ixyz+1)*pmom(ixyz+2) - pos(ixyz+2)*pmom(ixyz+1)
  enddo
  angmom(1) = angmomx
  angmomy = 0.d0
  do ixyz = 2, ncoord12, 3
     angmomy = angmomy + pos(ixyz+1)*pmom(ixyz-1) - pos(ixyz-1)*pmom(ixyz+1)
  enddo
  angmom(2) = angmomy
  angmomZ = 0.d0
  do ixyz = 3, ncoord12, 3
     angmomz = angmomz + pos(ixyz-2)*pmom(ixyz-1) - pos(ixyz-1)*pmom(ixyz-2)
  enddo
  angmom(3) = angmomZ
  print *
  print *,' Jtot(x,y,z) = ', angmom,' Jtot = ',dsqrt(angmom(1)**2+angmom(2)**2+angmom(3)**2)
  call flush(6)

!!! calculate instantaneous inertia matrix
!!!  I_{i,j} =  Sum_alp M_alp [- R_{i,alp) R_(j,alp) + delta_{i,j} R_alp**2]
  xinertM = 0.d0
  temp = sum(pos(1:ncoord1)*pos(1:ncoord1))*xmass1 + &
             sum(pos(ncoord1+1:ncoord12)*pos(ncoord1+1:ncoord12))*xmass2 ! Sum_alp  M_alp R_alp**2]
  do ixyz = 1, 3
     do jxyz = 1, 3
        xinertM(ixyz,jxyz) = sum(pos(ixyz:ncoord1:3)*pos(jxyz:ncoord1:3))*xmass1 + &
             sum(pos(ncoord1+ixyz:ncoord12:3)*pos(ncoord1+jxyz:ncoord12:3))*xmass2
     enddo
     xinertM(ixyz,ixyz) = xinertM(ixyz,ixyz) - temp
  enddo
  xinertM = -xinertM
  

!!! invert inertia matrix to get instantaneous angular velocity omega, from Jtot = I omega
  call matinv(xinertM,3)
  omega = matmul(xinertM,angmom)
  Erot = DOT_PRODUCT(angmom,omega)
  print *,' Erot = ',Erot,' Hartree, = ',Erot*Hartree2cm,' cm-1'
  call flush(6)

!!! each atom contributes M_alp omega x R_alp to the angular momentum:
!!! hence if this quantity is subtracted from its momentum the total J will be zero.
!!! the total momentum will also remain zero
!!! only the total kinetic energy will be slightly affected (of the order of several milli-cm-1 for 1000 atoms)

!!! decide if I really want to do it!!!
  do iat = 1, nbrat1
     ipos = 3*(iat-1)   !!! so that ipos+1 = x(iat), ipos+2 = y(iat), ipos+3 = z(iat)
     pmom(ipos+1) = pmom(ipos+1) - xmass1*(omega(2)*pos(ipos+3)-omega(3)*pos(ipos+2))
     pmom(ipos+2) = pmom(ipos+2) - xmass1*(omega(3)*pos(ipos+1)-omega(1)*pos(ipos+3))
     pmom(ipos+3) = pmom(ipos+3) - xmass1*(omega(1)*pos(ipos+2)-omega(2)*pos(ipos+1))
  enddo
  do iat = 1, nbrat2
     ipos = ncoord1 + 3*(iat-1)
     pmom(ipos+1) = pmom(ipos+1) - xmass2*(omega(2)*pos(ipos+3)-omega(3)*pos(ipos+2))
     pmom(ipos+2) = pmom(ipos+2) - xmass2*(omega(3)*pos(ipos+1)-omega(1)*pos(ipos+3))
     pmom(ipos+3) = pmom(ipos+3) - xmass2*(omega(1)*pos(ipos+2)-omega(2)*pos(ipos+1))
  enddo
     

  return
end subroutine comJ



  
  
