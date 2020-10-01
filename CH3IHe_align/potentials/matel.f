      function matel(j,m,k,jp,mp,kp,rvec)
      implicit none

      real*8 pi
      parameter (pi=3.141592635d0)

      integer j,m,k,jp,mp,kp
      real*8  rvec(3)
      real*8  cutoff
      parameter (cutoff=0.001d0)

      integer count, nu
      real*8  rr,radpot,threeJ
      complex*16 Y(5*5),matel,csum

      integer lambda

      call  Ylmb(rvec,4,Y)

      csum=(0.d0,0.d0)

      rr=dsqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)

      count=0
      do lambda=0,4
      do nu=-lambda,lambda
      count=count+1

      csum=csum+dsqrt(4.d0*pi/(2*lambda+1))
     . *threej(j,lambda,jp,k, 0,-kp)
     . *threej(j,lambda,jp,m,nu,-mp)
     . *radpot(lambda,rr)
     . *Y(count)

      enddo
      enddo

      matel=dreal(csum)

      if (dreal(matel).gt.cutoff) matel=cutoff

      return
      end
    
      
