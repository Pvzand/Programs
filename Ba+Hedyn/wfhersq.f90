subroutine wfhersq (prob, rdist,ndist,iread)
!!! He ZPAD wave function
!!! read from file iread, as  r psi(r) from  0  to +rl 
!!! and send back |r * psi(r)|**2 for sampling (prob),
!!! normalized to 1 so that the probability gives a nomber between 0 and 1

! implicit real (kind=8) :: (a-h,o-z)
logical :: lfirst
character (len=80) :: aline
integer, parameter :: ngridx=1000
integer :: ndist, idist, iread, iline, ngrid, igrid, indxin, indxout
real (kind=8) :: Hartree2cm1, Hartreemicro2cm1, splinr
double precision, dimension (ndist) :: prob, rdist
double precision, dimension (ngridx) :: rgrid
double precision, dimension (2*ngridx) :: fgrid
real (kind=8) :: probmax



data lfirst /.true./
parameter (Hartree2cm1=219474.63,Hartreemicro2cm1=1.0d-6*Hartree2cm1)

save rgrid, fgrid, ngrid

if (lfirst) then   
! read input points and set spline interpolation coefficients
!   read(iread,'(A80)') aline
   !   print *, aline
   ngrid = 0
   do iline = 1, ngridx
      read(iread,*,END=10000) rgrid(iline), fgrid(iline)
      ngrid = ngrid + 1
!      print *, rgrid(iline), fgrid(iline)
   enddo
   print *,' >>> WARNING from subroutine wfheb <<<'
   print *,' on file ',iread,' CAUTION: there may be more points to read '
10000 print *, ngrid,' wave function grid points read from file ',iread
   fgrid(1:ngrid) = fgrid(1:ngrid)**2 !!! in order to get probability
   probmax = MAXVAL(fgrid(1:ngrid))
   print *,' maximum value of fgrid = ',probmax
   fgrid(1:ngrid) = fgrid(1:ngrid)/probmax

   if (rgrid(1).lt.0.d0) then
      print *, 'the input wave function is read as:  psi(r)'
      print *,' from rmin = ',rgrid(1),' to rmax = ',rgrid(ngrid)
      print *,' >>> ERROR STOP, this is not the correct input file! <<<'
      STOP
   endif

   fgrid(ngrid+1:2*ngrid) = 0.d0
   
   call splset(fgrid, rgrid, ngrid)

   lfirst=.FALSE.
endif

indxin = 2
indxout = 2
   do idist = 1, ndist
      prob(idist) = SPLINR(fgrid,rgrid,indxin, indxout, ngrid, rdist(idist))
   enddo

return

end subroutine wfhersq


