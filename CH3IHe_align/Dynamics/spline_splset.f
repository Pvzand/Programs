ccccccccccc

      SUBROUTINE SPLSET(F,X,NX)
C
C         THIS ROUTINE SETS THE SPLINE INTERPOLATION ON THE GRID
C     (X(I),I=1,NX) FOR THE FUNCTION (F(I),I=1,NX).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NXMAX=8000, NXMX2=NXMAX-2)
      DIMENSION F(1),X(1)
CCC   DIMENSION HX(99),RLX(98),RMUX(98),XI(98),AB(4),YZ(4),A(4),B(98)
      DIMENSION HX(NXMAX-1)
      DIMENSION RLX(NXMX2), RMUX(NXMX2), XI(NXMX2), B(NXMX2)
      DIMENSION AB(4),YZ(4),A(4)
      DATA UN/1.0D0/,THREE/3.0D0/
C
      IF (NX.GT.NXMAX) GO TO 9999
C
      NX2=NX-2
C     COMPUTE THE HX'S
      DO 10 I=2,NX
   10 HX(I-1)=X(I)-X(I-1)
C     COMPUTE LAMBDA'S & MU'S
      DO 40 I=1,NX2
      RLX(I)=HX(I+1)/(HX(I)+HX(I+1))
   40 RMUX(I)=UN-RLX(I)
      MAN=NX-3
      DO 60 I=1,4
      A(I)=X(MAN)
   60 MAN=MAN+1
C
C     SPLINE-FIT DE P(X)
      DO 110 I=1,4
      AB(I)=F(I)
  110 YZ(I)=F(NX+I-4)
      P0=DLAGRA(X,AB,4,1)
      F(NX+1)=P0
      PN=DLAGRA(A,YZ,4,4)
      F(2*NX)=PN
C     CALCUL SECOND MEMBRE
      DO 120 I=1,NX2
  120 B(I)=THREE*RLX(I)/HX(I)*(F(I+1)-F(I))
     & +THREE*RMUX(I)/HX(I+1)*(F(I+2)-F(I+1))
      B(1)=B(1)-RLX(1)*P0
      B(NX2)=B(NX2)-RMUX(NX2)*PN
      CALL JORDAN(RMUX,RLX,XI,NX2,B)
      DO 100 I=1,NX2
  100 F(NX+I+1)=XI(I)
      RETURN
C
 9999 CONTINUE
      PRINT 90000
      PRINT 99999, NX, NXMAX
      STOP
C
90000 FORMAT(//,5X,15('*'),' STOP IN SPLSET ',15('*'),/)
99999 FORMAT(5X,' NX = ',I5,
     &      ' GREATER THAN DIMENSION PARAMETER NXMAX = ',I5)
C
      END
C
      SUBROUTINE JORDAN(MU,LAMBDA,X,N,B)
      IMPLICIT DOUBLE PRECISION (A-H,L-M,O-Z)
C
      PARAMETER (NXMAX=8000)
C
CCC   DIMENSION MU(N),LAMBDA(N),X(N),PIV(98),B(N)
      DIMENSION MU(N),LAMBDA(N),X(N),B(N)
      DIMENSION PIV(NXMAX)
C
      IF (NX.GT.NXMAX) GO TO 9999
C
C     CALCUL DES PIVOTS
      PIV(1)=2.D0
      DO 10 I=2,N
      PIV(I)=2.D0-LAMBDA(I)*MU(I-1)/PIV(I-1)
   10 B(I)=B(I)-LAMBDA(I)/PIV(I-1)*B(I-1)
C
      X(N)=B(N)/PIV(N)
      I=N-1
   20 X(I)=(B(I)-X(I+1)*MU(I))/PIV(I)
      I=I-1
      IF(I.GT.0) GOTO 20
      RETURN
C
 9999 CONTINUE
      PRINT 90000
      PRINT 99999, NX, NXMAX
      STOP
C
90000 FORMAT(//,5X,15('*'),' STOP IN JORDAN ',15('*'),/)
99999 FORMAT(5X,' NX = ',I5,
     &      ' GREATER THAN DIMENSION PARAMETER NXMAX = ',I5)
C
      END
C
      DOUBLE PRECISION FUNCTION DLAGRA (X,Y,MIN,IP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MIN),Y(MIN)
      DLAGRA=0.D0
      DO 10 I=1,MIN
      IF(I.EQ.IP) GOTO 10
      YP=Y(I)
      DO 20 J=1,MIN
      IF(J.EQ.IP) GOTO 20
      IF(J.EQ.I) GOTO 20
      YP=YP*(X(IP)-X(J))
   20 CONTINUE
      DO 30 J=1,MIN
      IF(J.EQ.I) GOTO 30
      YP=YP/(X(I)-X(J))
   30 CONTINUE
      DLAGRA=DLAGRA+YP
   10 CONTINUE
      DO 40 I=1,MIN
      IF(I.EQ.IP) GOTO 40
      DLAGRA=DLAGRA+Y(IP)/(X(IP)-X(I))
   40 CONTINUE
      RETURN
      END
C
C


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
cccccc

c**********************************************************************
c      FUNCTION SPLINT(F,X,IOLD,NX,R)
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C **********************************************************************
C  USE IOLD TO SAVE TIME IF SPLINT IS GOING TO BE CALLED FOR INCREASING
C  VALUES OF R:
C    - FOR 1ST CALL, SET IOLD TO 2
C    - FOR THE FOLLOWING CALLS, USE THE OUT VALUE OF IOLD
C  IF THE VALUES OF R ARE NOT MONOTONICALLY INCREASING, SET IOLD TO 2
C  FOR ALL THE CALLS
C **********************************************************************
c      DIMENSION U(4)
c      DIMENSION F(1),X(1)
c      DATA UN/1D0/,TWO/2D0/,THREE/3D0/
cC
c      IF(R.GE.X(NX)) GO TO 30
c      DO 10 IDOL = IOLD, NX
c      IF(R.LT.X(IDOL)) GOTO 20
c   10 CONTINUE
c   20 IOLD = IDOL
c      HI=X(IDOL)-X(IDOL-1)
c      XR=(R-X(IDOL-1))/HI
c      U(1)=XR*XR*(-TWO*XR+THREE)
c      U(3)=HI*XR*XR*(XR-UN)
c      U(2)=UN-U(1)
c      U(4)=HI*XR*((XR-TWO)*XR+UN)
c      SPLINT=U(1)*F(IDOL)+U(2)*F(IDOL-1)
c     &      +U(3)*F(NX+IDOL)+U(4)*F(NX+IDOL-1)
c      RETURN
c30    CONTINUE
c         SPLINT = 0.D0
c         RETURN
c      END
