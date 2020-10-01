      program pot_cYLM
     
      implicit none
      integer np,npo
      parameter (npo=128)
      integer i,j,m,k,jp,mp,kp
      integer count, nu, lambda
      real*8  pi
      real*8  rold,coef,coeff
      real*8  rr,cylm,xx
      real*8  x,xmin,xmax
      real*8  y,ymin,ymax
      real*8  z,zmin,zmax
      real*8  rvec(3)
      real*8  c(npo*2)
      real*8  splint
      complex*16 Ylm(5*5),csum
     
 
      dimension rold(npo),coef(npo,13,13)
      dimension coeff(2400)

      pi=dacos(-1.d0)

      open(1,file="cilm.dat")
      open(10,file="cilm_test.dat")
      
       do i=1,npo
         m=0
          
          read(1,*) rold(i),(coef(i,j,m),j=m,10)
          write(10,*) rold(i),(coef(i,j,m),j=m,10)
          
           
         enddo
      

      open(2,file="grid3D")
      read(2,*)np
      read(2,*)xmin,xmax
      read(2,*)ymin,ymax
      read(2,*)zmin,zmax


 
      open(3,file="c_Ylmb.dat")


      do i=1,np
        do j=1,np
         
         rvec(1)=xmin+(xmax-xmin)*dble(i-1)/dble(np-1)
         rvec(2)=ymin+(ymax-ymin)*dble(j-1)/dble(np-1)
         rvec(3)=zmin+(zmax-zmin)*dble(j-1)/dble(np-1)
c         rvec(3)=0.0
        
         call  Ylmb(rvec,4,Ylm)
  
         csum=(0.d0,0.d0)

         rr=dsqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)
     
         count=0
          do lambda=1,7
          nu=1
          count=count+1
c         
c          c(i)=coef(i,lambda,nu)       
          
          write(150,*)rr,lambda,nu,c(i)
          cylm=splint(c,rold,np,rr) 

          write(151,*)rr,cylm
          csum=csum+cylm*Ylm(count)
         
          
                 
c        enddo
      enddo

      write(3,*)rvec(1),rvec(2),dreal(csum),dimag(csum)

      enddo
      write(3,*)
      enddo

      end

c      enddo
C***********************************************************************
      FUNCTION SPLINT(F,X,NX,R)
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION U(4)
      DIMENSION F(2*NX),X(NX)
      DATA UN/1D0/,TWO/2D0/,THREE/3D0/
C
      IF(R.GE.X(NX)) GO TO 30
c      IF (R.EQ.RMEM) GO TO 40
c      RMEM = R
      DO 10 IDOL=2,NX
      IF(R.LT.X(IDOL)) GOTO 20
   10 CONTINUE
   20 HI=X(IDOL)-X(IDOL-1)
      XR=(R-X(IDOL-1))/HI
      U(1)=XR*XR*(-TWO*XR+THREE)
      U(3)=HI*XR*XR*(XR-UN)
      U(2)=UN-U(1)
      U(4)=HI*XR*((XR-TWO)*XR+UN)
   40 SPLINT=U(1)*F(IDOL)+U(2)*F(IDOL-1)+U(3)*F(NX+IDOL)
     &+U(4)*F(NX+IDOL-1)
      RETURN
30    RO=X(NX)
      YO=F(NX)
      N2X=2*NX
      YP=F(N2X)
      AIN=YO+YP*RO/6.D0
      C8=-AIN*3.D0*RO**8
      C6=YO*RO**6-C8/RO/RO
      SPLINT=C6/R**6+C8/R**8
         RETURN
      END

