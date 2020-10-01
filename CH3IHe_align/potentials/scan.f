      program scan
      implicit none

      integer i,j,n
      parameter (n=256)

      real*8 rvec(3),dx,dy,dz
      real*8 matel


      open(22,file='pp.dat')
     
      dx=0.1
      dy=0.1
      dz=0.1


      do i=1,n
      do j=1,n

      rvec(1)=dble(i-n/2)*dx
      rvec(2)=dble(j-n/2)*dy
      rvec(3)=dble(j-n/2)*dz
 
c     rvec(3)=0.d0
c     rvec(3)=0.d0
      rvec(3)=0.d0

      write(22,*) rvec(1),rvec(2),matel(0,0,0,0,0,0,rvec)

      enddo
      write(22,*) 
      enddo

      close(22)


      return 

      end
