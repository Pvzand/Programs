
      function threeJ(l1,l2,l3,m1,m2,m3)
      implicit none

      integer i

      integer  l1,l2,l3
      integer  m1,m2,m3
      double precision  al1,al2,al3
      double precision  am1,am2,am3
      double precision  threeJ
      double precision  FG01BD

      al1=l1
      al2=l2
      al3=l3
      am1=m1
      am2=m2
      am3=m3


c   relation between 3j and CG
c   / l1  l2  l3 \
c   \ m1  m2  m3 /  =  [(-1)**(l1-l2-m3)/sqrt(2*l3+1) ] *  < l1 m1, l2 m2 | l3, -m3 >
c

      threeJ=(-1.d0)**(l1-l2-m3)*1.d0/dsqrt(dble(2*l3)+1.d0)
     .         *FG01BD(al1,al2,al3,am1,am2,-am3)

      return
      end


      DOUBLE PRECISION FUNCTION FG01BD(A,B,C,XX,YY,ZZ)
C***********************************************************************
C*    CLEBSH-GORDAN COEFFICIENT  < A,XX; B,YY I C,ZZ > *
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z),INTEGER*2(I-N)
C
C
C     DATA JJJ/0/
      DIMENSION H(501),J(501)
      DIMENSION AY(4),IAY(4)
C
      INTPTF(Q)=Q+Q+DSIGN(.1D0,Q)
      IPARF(I)=4*(I/4)-I+1
C
    1 K1=INTPTF(A)
      K2=INTPTF(B)
      K3=INTPTF(C)
      K4=INTPTF(XX)
      K5=INTPTF(YY)
      K6=INTPTF(ZZ)
C
      H(1)=1.0D0
      J(1)=0
      X=0.D0
      DO 400 I=2,501
      X=X+1.0D0
      H(I)=H(I-1)*X
      J(I)=J(I-1)
  200 IF(H(I).LT.10.0D0) GO TO 400
      H(I)=0.01D0*H(I)
      J(I)=J(I)+2
      GO TO 200
  400 CONTINUE
C
      IF((K4+K5-K6).NE.0) GO TO 710
      M1=K1+K2-K3
      M2=K2+K3-K1
      M3=K3+K1-K2
      M4=K1+K4
      M5=K1-K4
      M6=K2+K5
      M7=K2-K5
      M8=K3+K6
      M9=K3-K6
      M10=K1+K2+K3+2
C
      IF(M1.LT.0) GO TO 710
      IF(M2.LT.0) GO TO 710
      IF(M3.LT.0) GO TO 710
      IF(M4.LT.0) GO TO 710
      IF(M5.LT.0) GO TO 710
      IF(M6.LT.0) GO TO 710
      IF(M7.LT.0) GO TO 710
      IF(M8.LT.0) GO TO 710
      IF(M9.LT.0) GO TO 710
      IF((M4-(M4/2)-(M4/2)).NE.0) GO TO 710
      IF((M4-(M4/2)-(M4/2)).NE.0) GO TO 710
      IF((M4-(M4/2)-(M4/2)).NE.0) GO TO 710
      IF((M10-(M10/2)-(M10/2)).NE.0) GO TO 710
C
      Y=K3+1
      M1=M1/2+1
      M2=M2/2+1
      M3=M3/2+1
      M4=M4/2+1
      M5=M5/2+1
      M6=M6/2+1
      M7=M7/2+1
      M8=M8/2+1
      M9=M9/2+1
      M10=M10/2+1
C
      Y=DSQRT(Y*H(M1)*H(M2)*H(M3)*H(M4)*H(M5)*
     X H(M6)*H(M7)*H(M8)*H(M9)/H(M10))
      IY=(J(M1)+J(M2)+J(M3)+J(M4)+J(M5)+
     X J(M6)+J(M7)+J(M8)+J(M9)-J(M10))/2
C
      N4=M1
      IF(N4.GT.M5)N4=M5
      IF(N4.GT.M6)N4=M6
      N4=N4-1
      M2=K2-K3-K4
      M3=K1+K5-K3
      N5=0
      IF(N5.LT.M2) N5=M2
      IF(N5.LT.M3) N5=M3
      N5PAR=IPARF(N5)
      N5=N5/2
      Z=0.0D0
      GO TO 610
C
  700 MM1=M1-N5
      MM2=M5-N5
      MM3=M6-N5
      MM4=N5-(M2/2)+1
      MM5=N5-(M3/2)+1
C
      X=1.D0/(H(MM1)*H(MM2)*H(MM3)*H(MM4)*H(MM5)*H(N5+1))
      IX=-J(MM1)-J(MM2)-J(MM3)-J(MM4)-J(MM5)-J(N5+1)
C
  800 IF(IX+IY)900,210,110
  900 X=0.1D0*X
      IX=IX+1
      GO TO 800
  110 X=10.0D0*X
      IX=IX-1
      GO TO 800
C
  210 IF(N5PAR.LT.0) X=-X
      Z=Z+X
  510 N5PAR=-N5PAR
      N5=N5+1
C
  610 IF(N5-N4)700,700,810
C
 710  CLEBSH=0.0D0
      GO TO 910
C
 810  CLEBSH=Z*Y
  910 CONTINUE
C
  220 FG01BD=CLEBSH
      RETURN
      END
