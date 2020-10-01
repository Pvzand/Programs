!Equations of motion for the ground state

subroutine ffct(t,y,dery,ndimtot) 
                                        
  implicit none

  include 'nbrat12.f90'

  real (kind=8) :: t
  real (kind=8) :: xmass1,xmass2,xmass,xmassinv
  
  integer :: ndimtot
  integer :: i

  real (kind=8), dimension (ndimtot):: y,dery
  
  common/masses/xmass1,xmass2,xmass(nbrat12),xmassinv(ncoord12)

!     dpi/dt = - dV/dri: calculates dery(i), i=ncoord12+1, 2*ncoord12
!   (equations of motion for the momenta)
  call derfpottot(y,dery,ndimtot)

! Equations of motion for positions: dxi/dt = pi/mi

  do i = 1, ncoord12
     dery(i) = y(ncoord12+i)*xmassinv(i)
  enddo
 
  return
  end subroutine ffct
