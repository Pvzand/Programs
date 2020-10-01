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
