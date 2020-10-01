program histo1d

  implicit none
  
  integer, parameter :: nxhisto = 20, nyhisto=20, nzhisto=20, nrhisto=20   ! number of intervals for histogram

  logical   :: lauto=.true.   !!! for automatic determination of the histogram bounds
  logical :: lout
  
  real     (kind=8) :: xhisto(0:nxhisto), yhisto(0:nyhisto), zhisto(0:nzhisto), rhisto(0:nrhisto)
  real     (kind=8) :: fxhisto(nxhisto), fyhisto(0:nyhisto), fzhisto(0:nzhisto),frhisto(0:nrhisto)
  
!!! variables
  real     (kind=8), allocatable :: xdat(:), ydat(:), zdat(:), rdat(:)
  real     (kind=8), allocatable :: rdat2(:) 
!!! min, max values of the data to be represented
  real     (kind=8) :: xdatmin, xdatmax
  real     (kind=8) :: ydatmin, ydatmax
  real     (kind=8) :: zdatmin, zdatmax
  real     (kind=8) :: rdatmin, rdatmax
  
!!! min, max values for the histogram   (overwritten if lauto=.true.)
  real     (kind=8) :: xhmin=-10.d0, xhmax=10.d0
  real     (kind=8) :: yhmin=-10.d0, yhmax=10.d0
  real     (kind=8) :: zhmin=-10.d0, zhmax=100.d0
  real     (kind=8) :: rhmin=-491.d-5, rhmax=-334.d-5
  real     (kind=8) :: delxh, delyh, delzh, delrh,energie

  integer           :: ndat, nxout, nyout, nzout, nxlow, nxhigh, nylow, nyhigh, nzlow, nzhigh,nrout,nrhigh,nrlow
  integer           :: idat, ix, iy, iz, ir,traj,nHe,nHe2,ndat2,i
  open(101,file='exploitation.in',status='unknown')
  open(100,file='fort.100',status='old')
!!! read input data from file 20
  read(101,*) nHe
  read(101,*) 
  read(101,*) 
  read(101,*) nHe2

!!! skip first lines (header)  
  read(100,*)
  !read(20,*)
  !read(20,*)
  !read(20,*)
  !read(20,*)

!!! read number of data points
  read(100,*) ndat2
  read(100,*)
  !ndat2=99
  !nHe=99
  print *,'  reading from file 20: ',ndat,' points'
  !allocate(xdat(ndat)); allocate(ydat(ndat)); allocate(zdat(ndat))
  
ndat=ndat2*nHe
  allocate(rdat(ndat))   ! a modifer pour revenir au trjectoire
  allocate(rdat2(nHe))
 do idat=1, ndat2
     read(100,*) traj,rdat2(1:nHe)
     do i=1,nHe
  !      read(100,*)traj, energie
        rdat(i+(idat-1)*nHe)=rdat2(i)
  !      rdat(i)=energie
     end do
 end do
print *,rdat
!!     print *,idat,xdat(idat)

  !read(20,*) rdat
  !xdat(:)=rdat(:,1)
  !ydat(:) = rdat(:,2)
  !zdat(:) = rdat(:,3)
  
  rdatmin = minval(rdat)
  rdatmax = maxval(rdat)
  print *,' xdatmin, xdatmax = ',rdatmin, rdatmax
  !ydatmin = minval(ydat)
  !ydatmax = maxval(ydat)
  !print *,' ydatmin, ydatmax = ',ydatmin, ydatmax
  !zdatmin = minval(zdat)
  !zdatmax = maxval(zdat)
  !print *,' zdatmin, zdatmax = ',zdatmin, zdatmax

   if (lauto) then
     rhmin = rdatmin
     rhmax = rdatmax
     !yhmin = ydatmin
     !yhmax = ydatmax
     !zhmin = zdatmin
     !zhmax = zdatmax
  endif
  delrh = (rhmax-rhmin)/dfloat(nrhisto)      
  write(*,*) delrh
  !delxh = (xhmax-xhmin)/dfloat(nxhisto)
  !delyh = (yhmax-yhmin)/dfloat(nyhisto)
  !delzh = (zhmax-zhmin)/dfloat(nzhisto)

  print *
  print *,' Histogram: '
  print *, ' in x: ',nrhisto,' intervals from ',rhmin,' to ',rhmax,' interval = ',delrh
  !print *, ' in x: ',nyhisto,' intervals from ',yhmin,' to ',yhmax,' interval = ',delyh
  !print *, ' in x: ',nzhisto,' intervals from ',zhmin,' to ',zhmax,' interval = ',delzh

  do ir = 0, nrhisto
     rhisto(ir) = rhmin + ir*delrh
  enddo
  !do iy = 0, nyhisto
  !   yhisto(iy) = yhmin + iy*delyh
  !enddo
  !do iz = 0, nzhisto
  !   zhisto(iz) = zhmin + iz*delzh
  !enddo
  
!!! determine grids for the 2D histogram
  write(80,*)'#   histogram initial for distance He-Ba '
 !write(80,*) '# built from',ndat,' He atoms'

  frhisto = 0.d0; fyhisto = 0.d0; fzhisto = 0.d0
  
  nrlow = 0; nrhigh =0 
  nylow = 0; nyhigh = 0
  nzlow = 0; nzhigh = 0
  nrout = 0; nyout = 0; nzout = 0
  do idat = 1, ndat2
     lout = .false.
     ir = nint((rdat(idat)-rhmin)/delrh)
     if (ir.lt.0) then
        nrlow=nrlow+1
        lout = .true.
     elseif (ir.gt.nrhisto) then
        nrhigh = nrhigh+1
        lout = .true.
     endif
     if (lout) then
        nrout = nrout + 1
      else
        frhisto(ir) = frhisto(ir) + 1.d0
     endif     
     lout = .false.
     !iy = nint((ydat(idat)-yhmin)/delyh)
     !if (iy.lt.0) then
     !   nylow=nylow+1
     !   lout = .true.
     !elseif (iy.gt.nyhisto) then
     !   nyhigh = nyhigh+1
     !   lout = .true.
     !endif
     !if (lout) then
     !   nyout = nyout + 1
     !else
     !   fyhisto(iy) = fyhisto(iy) + 1.d0
     !endif     
     !lout = .false.
     !iz = nint((zdat(idat)-zhmin)/delzh)
     !if (iz.lt.0) then
     !   nzlow=nzlow+1
     !   lout = .true.
     !elseif (iz.gt.nzhisto) then
     !   nzhigh = nzhigh+1
     !   lout = .true.
     !endif
     !if (lout) then
     !   nzout = nzout + 1
     !else
     !   fzhisto(iz) = fzhisto(iz) + 1.d0
     !endif     
  enddo

!!!  fhisto = fhisto/ndat

  print *
  print *,'   Histogram completed '
  print *,' in R:'
  print *,' total number of data points out of range: nrout = ',nrout
  print *,'    number of points below minimum r value: nrlow = ',nrlow
  print *,'    number of points above maximum r value: nrhigh = ',nrhigh
  !print *,' in Y:'
  !print *,' total number of data points out of range: nyout = ',nyout
  !print *,'    number of points below minimum y value: nylow = ',nylow
  !print *,'    number of points above maximum y value: nyhigh = ',nyhigh
  !print *,' in Z:'
  !print *,' total number of data points out of range: nzout = ',nzout
  !print *,'    number of points below minimum z value: nzlow = ',nzlow
  !print *,'    number of points above maximum z value: nzhigh = ',nzhigh!

  print *
  print *,'   writing output histogram to file 80'

  do ir = 0, nrhisto
     write(80,*) rhisto(ir), frhisto(ir)
  enddo
  !do iy = 0, nyhisto
  !      write(81,*) yhisto(iy), fyhisto(iy)
  !enddo
  !do iz = 0, nzhisto
  !      write(82,*) zhisto(iz), fzhisto(iz)
  !enddo
  
  stop
end program histo1d



