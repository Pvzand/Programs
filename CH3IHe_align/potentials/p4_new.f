      program ppp
      use SHTOOLS

      implicit none

      integer    k,l,m

      integer    i,j,n,ndum,nr
      parameter (n=256)     ! number of points in th,ph, needs to be even 

      real*8     r,rmax,rmin

      real*8     th,ph,pi,fact
      real*8     potCH3I
      real*8     rvec(3)
      complex*16 pott
      complex*16 grid(n,n)

      integer    sampling,csphase,lmax,norm,lmax_calc
   
      parameter (sampling=1)       ! use nxn grid
      parameter (csphase=1)        ! use nxn grid
      parameter (norm=4)           ! use orthonormal sh
      parameter (lmax_calc=n/2-1)  ! use nxn grid

      complex*16 cilm(2,lmax_calc+1,lmax_calc+1)

 
      pi=dacos(-1.d0)
      open(1,file='grid')
      read(1,*)nr
      read(1,*)rmin,rmax

      open(45,file='potdev.dat')
      open(46,file='pdev.dat')
      open(50,file='cilm_neg.dat')
      open(51,file='cilm_pos.dat')
      open(52,file='cilm_0.dat')


c_____loop over radius

      do k=1,nr
c       r=rmin
      r=rmin+(rmax-rmin)*dble(k-1)/dble(nr-1)
!      write(*,*) r
c       r=5.669177966
c       r=8.503766949
c      r=11.338355931
c       r=15.117807909
c        r=22.676711863

      do i=1,n
      do j=1,n

c      write(*,*) i,j

       th= pi*dble(i-1)/dble(n-1)
       ph= 2.d0*pi*dble(j-1)/dble(n-1)
!       ph=0.0

       rvec(1)=r*sin(th)*cos(ph)
       rvec(2)=r*sin(th)*sin(ph)
       rvec(3)=r*cos(th)

    
c       grid(i,j)=pott(th,ph)        ! for test-model function
      grid(i,j)=dcmplx(potCH3I(0,r,th,ph),0.d0)

      write(120,*) th,ph,dreal(grid(i,j))  ! for gnu-plot
      enddo
      write(120,*)   ! for gnu-plot
      enddo
      
      call SHExpandDHC(grid,n,cilm,lmax,norm,sampling,csphase,lmax_calc)

       do m=-10,10
!       m=0
       fact = (-1)**m       
       if (m.lt.0) write(50,*) (fact*dreal(cilm(2,j+1,-m+1)),j=-m,40)
       if (m.eq.0) write(52,*) r,(fact*dreal(cilm(1,j+1,m+1)),j=m,40)
     
       if (m.gt.0) write(51,*) (fact*dreal(cilm(1,j+1,m+1)),j=m,40)
       enddo

       do j=0,39
       write(300,*) r, dreal(cilm(1,j+1,m+1))
       enddo

      write(45,88) r
     .  ,dreal(cilm(1,1,1))  ! 0,0
c     .  ,dreal(cilm(1,2,1))  ! 1,0
c     .  ,dreal(cilm(1,2,2))  ! 1,1
c     .  ,dreal(cilm(2,2,2))  ! 1,-1
c     .  ,dreal(cilm(1,3,1))  ! 2,0
c     .  ,dreal(cilm(1,3,2))  ! 2,1
c     .  ,dreal(cilm(1,3,3))  ! 2,2
c     .  ,dreal(cilm(1,4,1))  ! 3,0
c     .  ,dreal(cilm(1,4,4))  ! 3,3
c     .  ,dreal(cilm(1,5,1))  ! 4,0
c     .  ,dreal(cilm(1,5,4))  ! 4,3
c     .  ,dreal(cilm(1,6,1))  ! 5,0
c     .  ,dreal(cilm(1,6,4))  ! 5,3
c     .  ,dreal(cilm(1,7,1))  ! 6,0
c     .  ,dreal(cilm(1,7,4))  ! 6,3
c     .  ,dreal(cilm(1,7,7))  ! 6,6
c     .  ,dreal(cilm(1,8,1))  ! 7,0
c     .  ,dreal(cilm(1,8,4))  ! 7,3
c     .  ,dreal(cilm(1,8,7))  ! 7,6
c     .  ,dreal(cilm(1,9,1))  ! 8,0
c     .  ,dreal(cilm(1,9,4))  ! 8,3
c     .  ,dreal(cilm(1,10,1))   ! 9,0
c     .  ,dreal(cilm(1,10,4))   ! 9,3
c     .  ,dreal(cilm(1,11,1))   ! 10,0
c     .  ,dreal(cilm(1,12,1))   ! 11,0
c     .  ,dreal(cilm(1,13,1))   ! 12,0
c     .  ,dreal(cilm(1,14,1))   ! 13,0
c     .  ,dreal(cilm(1,15,1))   ! 13,0
c     .  ,dreal(cilm(1,16,1))   ! 13,0
c     .  ,dreal(cilm(1,17,1))   ! 13,0
c     .  ,dreal(cilm(1,18,1))   ! 13,0
c     .  ,dreal(cilm(1,19,1))   ! 13,0
c     .  ,dreal(cilm(1,20,1))   ! 13,0

c      write(47,88) r
c     .  ,dreal(cilm(2,1,1)),dimag(cilm(2,1,1))   ! 0,0
c     .  ,dreal(cilm(2,2,1)),dimag(cilm(2,2,1))   ! 1,0
c     .  ,dreal(cilm(2,3,1)),dimag(cilm(2,3,1))   ! 2,0
c     .  ,dreal(cilm(2,4,1)),dimag(cilm(2,4,1))   ! 3,0
c     .  ,dreal(cilm(2,4,4)),dimag(cilm(2,4,4))   ! 3,-3
c     .  ,dreal(cilm(2,5,1)),dimag(cilm(2,5,1))   ! 4,0
c     .  ,dreal(cilm(2,5,4)),dimag(cilm(2,5,4))   ! 4,-3
c     .  ,dreal(cilm(2,6,1)),dimag(cilm(2,1,1))   ! 5,0
c     .  ,dreal(cilm(2,6,4)),dimag(cilm(2,1,1))   ! 5,-3
c     .  ,dreal(cilm(2,7,1)),dimag(cilm(2,7,1))   ! 6,0
c     .  ,dreal(cilm(2,7,4)),dimag(cilm(2,7,4))   ! 6,-3
c     .  ,dreal(cilm(2,7,7)),dimag(cilm(2,7,7))   ! 6,-6
c     .  ,dreal(cilm(2,8,1)),dimag(cilm(2,8,1))   ! 7,0
c     .  ,dreal(cilm(2,8,4)),dimag(cilm(2,8,4))   ! 7,-3
c     .  ,dreal(cilm(2,8,8)),dimag(cilm(2,8,8))   ! 7,-7
c     .  ,dreal(cilm(2,9,1)),dimag(cilm(2,9,1))   ! 8,0
c     .  ,dreal(cilm(2,9,4)),dimag(cilm(2,9,4))   ! 8.-3



c      write(46,88) r
c     .  ,dreal(cilm(1,1,1))   ! 0,0
c     .  ,dreal(cilm(1,2,1))   ! 1,0
c     .  ,dreal(cilm(1,3,1))   ! 2,0
c     .  ,dreal(cilm(1,4,1))   ! 3,0
c     .  ,dreal(cilm(1,4,4))   ! 3,3
c     .  ,dreal(cilm(1,5,1))   ! 4,0
c     .  ,dreal(cilm(1,5,4))   ! 4,3
c     .  ,dreal(cilm(1,6,1))   ! 5,0
c     .  ,dreal(cilm(1,6,4))   ! 5,3
c     .  ,dreal(cilm(1,7,1))   ! 6,0
c     .  ,dreal(cilm(1,7,4))   ! 6,3
c     .  ,dreal(cilm(1,7,7))   ! 6,6
c  

      
      enddo   
c end of loop over radius

      
      
88    format(26F14.8)

      end


      


      function pott(th,ph)
      implicit none

      real*8 fpi,pi
      parameter (pi=dacos(-1.d0))
      parameter (fpi=1.d0/dsqrt(4.d0*pi))

      real*8 th,ph

      complex*16 pott
      complex*16 y00,y10,y20,y30,y3m3,y3p3
      complex*16 y2m2,y2p2  

      y00=fpi
 
      y10=fpi*dsqrt(3.d0)*cos(th)

      y20=fpi*dsqrt(5.d0/4.d0)*(3.d0*cos(th)**2-1.d0)

      y2m2=fpi*dsqrt(15.d0/8.d0)*sin(th)**2
     .    *cdexp(dcmplx(0.d0,-2.d0*ph))

      y2p2=fpi*dsqrt(15.d0/8.d0)*sin(th)**2
     .    *cdexp(dcmplx(0.d0,2.d0*ph))

      y3m3=fpi*dsqrt(35.d0/16.d0)*sin(th)**3
     .    *cdexp(dcmplx(0.d0,-3.d0*ph))

      y30 = fpi*dsqrt(7.d0/4.d0)*(5.d0*cos(th)**3-3.d0*cos(th))

      y3p3=-fpi*dsqrt(35.d0/16.d0)*sin(th)**3
     .    *cdexp(dcmplx(0.d0, 3.d0*ph))


      pott=(0.1*y00+0.012*y10+0.2*y20)+0.2*(y2m2-y2p2)+0.4*(y3m3-y3p3)
!      pott=0.4*(y3m3-y3p3)

      return

      end



