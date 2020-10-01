      program ppp
      use SHTOOLS

      implicit none

      real*8     pi
      parameter (pi=3.141592635d0)

      integer    k,m

      integer    i,j,n,ndum
      parameter (n=48)

      real*8     r
      real*8     th,ph
      real*8     potCH3I
      real*8     pott
      real*8     grid(n,n)

      integer    sampling,csphase,lmax,norm,lmax_calc

      parameter (sampling=1)       ! use nxn grid
      parameter (csphase=1)        ! use nxn grid
      parameter (norm=4)           ! use orthonormal sh
      parameter (lmax_calc=n/2-1)  ! use nxn grid

      real*8 cilm(2,lmax_calc+1,lmax_calc+1)


      r=9.0d0



      do i=1,n
      do j=1,n

      write(*,*) i,j


      th=      3.1415927d0*dble(i-1)/dble(n)
      ph= 2.d0*3.1415927d0*dble(j-1)/dble(n)
    
      grid(i,j)=pott(th,ph)
c     grid(i,j)=potCH3I(0,r,th,ph)

      write(22,*) th,ph,grid(i,j)
c     write(23,*) i,j,grid(i,j)
 
      enddo
      write(22,*) 
      enddo

      call SHExpandDH(grid,n,cilm,lmax,norm,sampling,csphase,lmax_calc)

      do k=0,lmax_calc
      do m=-k,k

       if (m.lt.0) write(33,*) k,m,cilm(2,k+1,-m+1)
c     .                           ,dimag(cilm(2,k+1,-m+1))
      if (m.ge.0) write(33,*) k,m,cilm(1,k+1, m+1)
c     .                           ,dimag(cilm(1,k+1, m+1))


c      if (m.lt.0) write(33,*) k,m,dreal(cilm(2,k+1,-m+1))
c     .                           ,dimag(cilm(2,k+1,-m+1))
c      if (m.ge.0) write(33,*) k,m,dreal(cilm(1,k+1, m+1))
c     .                           ,dimag(cilm(1,k+1, m+1))

      enddo
      enddo

c     write(*,*) potCH3I(11,r,0.1d0,0.d0)
c     write(*,*) potCH3I(12,r,0.9*pi,0.d0)

      end


      


      function pott(th,ph)
      implicit none

      real*8 fpi
      parameter (fpi=1.d0/dsqrt(4.d0*3.1415927d0))

      real*8 th,ph

      real*8  pott
c      complex*16 pott
      complex*16 y00,y10,y20,y30,y3m3,y3p3
   

      y00=fpi

      y10=fpi*dsqrt(3.d0)*cos(th)

      y20=fpi*dsqrt(5.d0/4.d0)*(3.d0*cos(th)**2-1.d0)

      y3m3=fpi*dsqrt(35.d0/16.d0)*sin(th)**3
     .    *cdexp(dcmplx(0.d0,-3.d0*ph))

      y30 = fpi*dsqrt(7.d0/4.d0)*(5.d0*cos(th)**3-3.d0*cos(th))

      y3p3=-fpi*dsqrt(35.d0/16.d0)*sin(th)**3
     .    *cdexp(dcmplx(0.d0, 3.d0*ph))


      pott=(0.1*y00+0.012*y10+0.2*y20+0.1*y30)+0.4*(y3m3-y3p3)


      return

      end



