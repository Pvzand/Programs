      program pot_ylmb
      implicit none
    
      integer i,n
      real*8 fpi,pi
      parameter (pi=dacos(-1.d0))
      parameter (n=256)
      parameter (fpi=1.d0/dsqrt(4.d0*pi))

      real*8 th,ph

      complex*16 pott
      complex*16 y00,y10,y20,y30,y3m3,y3p3
      complex*16 y2m2,y2p2
   
      open(1,file='pot_ylm_comp')


      do i=1,n
     

      th=pi*dble(i-1)/dble(n)
      ph=0.d0
c     ph=pi 

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


c      pott=(0.1*y00+0.012*y10+0.2*y20+0.1*y30)+0.4*(y3m3-y3p3)
      pott=(0.1*y00+0.012*y10+0.2*y20)+0.2*(y2m2-y2p2)+0.4*(y3m3-y3p3)

      write(1,*) th,dreal(pott),dreal(y00),dreal(y10),dreal(y20),
     . dreal(y3m3),dreal(y30),dreal(y3p3)

      enddo
      end



