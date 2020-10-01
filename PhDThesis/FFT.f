c Subrutina FFT valida para potencias de 2
      subroutine tfft(a,nmax,n,tipo)
      implicit none
      integer n,nmax,m,mtwo,tipo,i
      complex*16 a(nmax)

c     Comprobar q n es potencia de 2 y guardar la potencia en m
      mtwo=2
      m=1
      do while (mtwo.lt.n)
       mtwo=mtwo+mtwo
       m=m+1
c      write(6,*) m,mtwo
      enddo

      if (mtwo.eq.n) then
       if (tipo.eq.1) call creaIWK(a,m,nmax)
       if (tipo.eq.-1) then
                     DO 10 I=1,N                                       
                        A(I) = CONJG(A(I))                            
10                   CONTINUE                                        
                     CALL creaiwk (A,M,nmax)                           
                     DO 20 I=1,N                                   
                        A(I) = CONJG(A(I))/N                      
20                   CONTINUE                                    
        endif
      else
       write(6,*) n,' no es una potencia de 2'
       stop
      endif
      return
      end
       
      
      subroutine creaIWK(a,m,nmax)
      implicit none
      integer m,nmax
      integer iwk(m+1)
      complex*16 a(nmax-1)
      call FFT2C(a,m,iwk)
      return
      end

C   IMSL ROUTINE NAME   - FFT2C                                         FFT20010
C                                                                       FFT20020
C-----------------------------------------------------------------------FFT20030
C                                                                       FFT20040
C   COMPUTER            - PRIME/DOUBLE                                  FFT20050
C                                                                       FFT20060
C   LATEST REVISION     - JANUARY 1, 1978                               FFT20070
C                                                                       FFT20080
C   PURPOSE             - COMPUTE THE FAST FOURIER TRANSFORM OF A       FFT20090
C                           COMPLEX VALUED SEQUENCE OF LENGTH EQUAL TO  FFT20100
C                           A POWER TWO                                 FFT20110
C                                                                       FFT20120
C   USAGE               - CALL FFT2C (A,M,IWK)                          FFT20130
C                                                                       FFT20140
C   ARGUMENTS    A      - COMPLEX VECTOR OF LENGTH N, WHERE N=2**M.     FFT20150
C                           ON INPUT A CONTAINS THE COMPLEX VALUED      FFT20160
C                           SEQUENCE TO BE TRANSFORMED.                 FFT20170
C                           ON OUTPUT A IS REPLACED BY THE              FFT20180
C                           FOURIER TRANSFORM.                          FFT20190
C                M      - INPUT EXPONENT TO WHICH 2 IS RAISED TO        FFT20200
C                           PRODUCE THE NUMBER OF DATA POINTS, N        FFT20210
C                           (I.E. N = 2**M).                            FFT20220
C                IWK    - WORK VECTOR OF LENGTH M+1.                    FFT20230
C                                                                       FFT20240
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         FFT20250
C                       - SINGLE/H36,H48,H60                            FFT20260
C                                                                       FFT20270
C   REQD. IMSL ROUTINES - NONE REQUIRED                                 FFT20280
C                                                                       FFT20290
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           FFT20300
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      FFT20310
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  FFT20320
C                                                                       FFT20330
C   REMARKS  1.  FFT2C COMPUTES THE FOURIER TRANSFORM, X, ACCORDING     FFT20340
C                TO THE FOLLOWING FORMULA;                              FFT20350
C                                                                       FFT20360
C                  X(K+1) = SUM FROM J = 0 TO N-1 OF                    FFT20370
C                           A(J+1)*CEXP((0.0,(2.0*PI*J*K)/N))           FFT20380
C                  FOR K=0,1,...,N-1 AND PI=3.1415...                   FFT20390
C                                                                       FFT20400
C                NOTE THAT X OVERWRITES A ON OUTPUT.                    FFT20410
C            2.  FFT2C CAN BE USED TO COMPUTE                           FFT20420
C                                                                       FFT20430
C                  X(K+1) = (1/N)*SUM FROM J = 0 TO N-1 OF              FFT20440
C                           A(J+1)*CEXP((0.0,(-2.0*PI*J*K)/N))          FFT20450
C                  FOR K=0,1,...,N-1 AND PI=3.1415...                   FFT20460
C                                                                       FFT20470
C                BY PERFORMING THE FOLLOWING STEPS;                     FFT20480
C                                                                       FFT20490
C                     DO 10 I=1,N                                       FFT20500
C                        A(I) = CONJG(A(I))                             FFT20510
C                  10 CONTINUE                                          FFT20520
C                     CALL FFT2C (A,M,IWK)                              FFT20530
C                     DO 20 I=1,N                                       FFT20540
C                        A(I) = CONJG(A(I))/N                           FFT20550
C                  20 CONTINUE                                          FFT20560
C                                                                       FFT20570
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       FFT20580
C                                                                       FFT20590
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN FFT20600
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    FFT20610
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        FFT20620
C                                                                       FFT20630
C-----------------------------------------------------------------------FFT20640
C                                                                       FFT20650
      SUBROUTINE  FFT2C (A,M,IWK)                                       FFT20660
C                                  SPECIFICATIONS FOR ARGUMENTS         FFT20670
      INTEGER            M,IWK(1)                                       FFT20680
      COMPLEX*16         A(1)                                           FFT20690
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   FFT20700
      INTEGER            I,ISP,J,JJ,JSP,K,K0,K1,K2,K3,KB,KN,MK,MM,MP,N, FFT20710
     1                   N4,N8,N2,LM,NN,JK                              FFT20720
      DOUBLE PRECISION   RAD,C1,C2,C3,S1,S2,S3,CK,SK,SQ,A0,A1,A2,A3,    FFT20730
     1                   B0,B1,B2,B3,TWOPI,TEMP,                        FFT20740
     2                   ZERO,ONE,Z0(2),Z1(2),Z2(2),Z3(2)               FFT20750
      COMPLEX*16         ZA0,ZA1,ZA2,ZA3,AK2                            FFT20760
      EQUIVALENCE        (ZA0,Z0(1)),(ZA1,Z1(1)),(ZA2,Z2(1)),           FFT20770
     1                   (ZA3,Z3(1)),(A0,Z0(1)),(B0,Z0(2)),(A1,Z1(1)),  FFT20780
     2                   (B1,Z1(2)),(A2,Z2(1)),(B2,Z2(2)),(A3,Z3(1)),   FFT20790
     3                   (B3,Z3(2))                                     FFT20800
      DATA               SQ/.70710678118655D0/,                         FFT20810
     1                   SK/.38268343236509D0/,                         FFT20820
     2                   CK/.92387953251129D0/,                         FFT20830
     3                   TWOPI/6.2831853071796D0/                       FFT20840
      DATA               ZERO/0.0D0/,ONE/1.0D0/                         FFT20850
C                                  SQ=SQRT2/2,SK=SIN(PI/8),CK=COS(PI/8) FFT20860
C                                  TWOPI=2*PI                           FFT20870
C                                  FIRST EXECUTABLE STATEMENT           FFT20880
      MP = M+1                                                          FFT20890
      N = 2**M                                                          FFT20900
      IWK(1) = 1                                                        FFT20910
      MM = (M/2)*2                                                      FFT20920
      KN = N+1                                                          FFT20930
C                                  INITIALIZE WORK VECTOR               FFT20940
      DO 5  I=2,MP                                                      FFT20950
         IWK(I) = IWK(I-1)+IWK(I-1)                                     FFT20960
    5 CONTINUE                                                          FFT20970
      RAD = TWOPI/N                                                     FFT20980
      MK = M - 4                                                        FFT20990
      KB = 1                                                            FFT21000
      IF (MM .EQ. M) GO TO 15                                           FFT21010
      K2 = KN                                                           FFT21020
      K0 = IWK(MM+1) + KB                                               FFT21030
   10 K2 = K2 - 1                                                       FFT21040
      K0 = K0 - 1                                                       FFT21050
      AK2 = A(K2)                                                       FFT21060
      A(K2) = A(K0) - AK2                                               FFT21070
      A(K0) = A(K0) + AK2                                               FFT21080
      IF (K0 .GT. KB) GO TO 10                                          FFT21090
   15 C1 = ONE                                                          FFT21100
      S1 = ZERO                                                         FFT21110
      JJ = 0                                                            FFT21120
      K = MM - 1                                                        FFT21130
      J = 4                                                             FFT21140
      IF (K .GE. 1) GO TO 30                                            FFT21150
      GO TO 70                                                          FFT21160
   20 IF (IWK(J) .GT. JJ) GO TO 25                                      FFT21170
      JJ = JJ - IWK(J)                                                  FFT21180
      J = J-1                                                           FFT21190
      IF (IWK(J) .GT. JJ) GO TO 25                                      FFT21200
      JJ = JJ - IWK(J)                                                  FFT21210
      J = J - 1                                                         FFT21220
      K = K + 2                                                         FFT21230
      GO TO 20                                                          FFT21240
   25 JJ = IWK(J) + JJ                                                  FFT21250
      J = 4                                                             FFT21260
   30 ISP = IWK(K)                                                      FFT21270
      IF (JJ .EQ. 0) GO TO 40                                           FFT21280
C                                  RESET TRIGONOMETRIC PARAMETERS       FFT21290
      C2 = JJ * ISP * RAD                                               FFT21300
      C1 = DCOS(C2)                                                     FFT21310
      S1 = DSIN(C2)                                                     FFT21320
   35 C2 = C1 * C1 - S1 * S1                                            FFT21330
      S2 = C1 * (S1 + S1)                                               FFT21340
      C3 = C2 * C1 - S2 * S1                                            FFT21350
      S3 = C2 * S1 + S2 * C1                                            FFT21360
   40 JSP = ISP + KB                                                    FFT21370
C                                  DETERMINE FOURIER COEFFICIENTS       FFT21380
C                                    IN GROUPS OF 4                     FFT21390
      DO 50 I=1,ISP                                                     FFT21400
         K0 = JSP - I                                                   FFT21410
         K1 = K0 + ISP                                                  FFT21420
         K2 = K1 + ISP                                                  FFT21430
         K3 = K2 + ISP                                                  FFT21440
         ZA0 = A(K0)                                                    FFT21450
         ZA1 = A(K1)                                                    FFT21460
         ZA2 = A(K2)                                                    FFT21470
         ZA3 = A(K3)                                                    FFT21480
         IF (S1 .EQ. ZERO) GO TO 45                                     FFT21490
         TEMP = A1                                                      FFT21500
         A1 = A1 * C1 - B1 * S1                                         FFT21510
         B1 = TEMP * S1 + B1 * C1                                       FFT21520
         TEMP = A2                                                      FFT21530
         A2 = A2 * C2 - B2 * S2                                         FFT21540
         B2 = TEMP * S2 + B2 * C2                                       FFT21550
         TEMP = A3                                                      FFT21560
         A3 = A3 * C3 - B3 * S3                                         FFT21570
         B3 = TEMP * S3 + B3 * C3                                       FFT21580
   45    TEMP = A0 + A2                                                 FFT21590
         A2 = A0 - A2                                                   FFT21600
         A0 = TEMP                                                      FFT21610
         TEMP = A1 + A3                                                 FFT21620
         A3 = A1 - A3                                                   FFT21630
         A1 = TEMP                                                      FFT21640
         TEMP = B0 + B2                                                 FFT21650
         B2 = B0 - B2                                                   FFT21660
         B0 = TEMP                                                      FFT21670
         TEMP = B1 + B3                                                 FFT21680
         B3 = B1 - B3                                                   FFT21690
         B1 = TEMP                                                      FFT21700
         A(K0) = DCMPLX(A0+A1,B0+B1)                                    FFT21710
         A(K1) = DCMPLX(A0-A1,B0-B1)                                    FFT21720
         A(K2) = DCMPLX(A2-B3,B2+A3)                                    FFT21730
         A(K3) = DCMPLX(A2+B3,B2-A3)                                    FFT21740
   50 CONTINUE                                                          FFT21750
      IF (K .LE. 1) GO TO 55                                            FFT21760
      K = K - 2                                                         FFT21770
      GO TO 30                                                          FFT21780
   55 KB = K3 + ISP                                                     FFT21790
C                                  CHECK FOR COMPLETION OF FINAL        FFT21800
C                                    ITERATION                          FFT21810
      IF (KN .LE. KB) GO TO 70                                          FFT21820
      IF (J .NE. 1) GO TO 60                                            FFT21830
      K = 3                                                             FFT21840
      J = MK                                                            FFT21850
      GO TO 20                                                          FFT21860
   60 J = J - 1                                                         FFT21870
      C2 = C1                                                           FFT21880
      IF (J .NE. 2) GO TO 65                                            FFT21890
      C1 = C1 * CK + S1 * SK                                            FFT21900
      S1 = S1 * CK - C2 * SK                                            FFT21910
      GO TO 35                                                          FFT21920
   65 C1 = (C1 - S1) * SQ                                               FFT21930
      S1 = (C2 + S1) * SQ                                               FFT21940
      GO TO 35                                                          FFT21950
   70 CONTINUE                                                          FFT21960
C                                  PERMUTE THE COMPLEX VECTOR IN        FFT21970
C                                    REVERSE BINARY ORDER TO NORMAL     FFT21980
C                                    ORDER                              FFT21990
      IF(M .LE. 1) GO TO 9005                                           FFT22000
      MP = M+1                                                          FFT22010
      JJ = 1                                                            FFT22020
C                                  INITIALIZE WORK VECTOR               FFT22030
      IWK(1) = 1                                                        FFT22040
      DO 75  I = 2,MP                                                   FFT22050
         IWK(I) = IWK(I-1) * 2                                          FFT22060
   75 CONTINUE                                                          FFT22070
      N4 = IWK(MP-2)                                                    FFT22080
      IF (M .GT. 2) N8 = IWK(MP-3)                                      FFT22090
      N2 = IWK(MP-1)                                                    FFT22100
      LM = N2                                                           FFT22110
      NN = IWK(MP)+1                                                    FFT22120
      MP = MP-4                                                         FFT22130
C                                  DETERMINE INDICES AND SWITCH A       FFT22140
      J = 2                                                             FFT22150
   80 JK = JJ + N2                                                      FFT22160
      AK2 = A(J)                                                        FFT22170
      A(J) = A(JK)                                                      FFT22180
      A(JK) = AK2                                                       FFT22190
      J = J+1                                                           FFT22200
      IF (JJ .GT. N4) GO TO 85                                          FFT22210
      JJ = JJ + N4                                                      FFT22220
      GO TO 105                                                         FFT22230
   85 JJ = JJ - N4                                                      FFT22240
      IF (JJ .GT. N8) GO TO 90                                          FFT22250
      JJ = JJ + N8                                                      FFT22260
      GO TO 105                                                         FFT22270
   90 JJ = JJ - N8                                                      FFT22280
      K = MP                                                            FFT22290
   95 IF (IWK(K) .GE. JJ) GO TO 100                                     FFT22300
      JJ = JJ - IWK(K)                                                  FFT22310
      K = K - 1                                                         FFT22320
      GO TO 95                                                          FFT22330
  100 JJ = IWK(K) + JJ                                                  FFT22340
  105 IF (JJ .LE. J) GO TO 110                                          FFT22350
      K = NN - J                                                        FFT22360
      JK = NN - JJ                                                      FFT22370
      AK2 = A(J)                                                        FFT22380
      A(J) = A(JJ)                                                      FFT22390
      A(JJ) = AK2                                                       FFT22400
      AK2 = A(K)                                                        FFT22410
      A(K) = A(JK)                                                      FFT22420
      A(JK) = AK2                                                       FFT22430
  110 J = J + 1                                                         FFT22440
C                                  CYCLE REPEATED UNTIL LIMITING NUMBER FFT22450
C                                    OF CHANGES IS ACHIEVED             FFT22460
      IF (J .LE. LM) GO TO 80                                           FFT22470
C                                                                       FFT22480
 9005 RETURN                                                            FFT22490
      END                                                               FFT22500
C   IMSL ROUTINE NAME   - FFTCC                                         FFTP0010
C                                                                       FFTP0020
C-----------------------------------------------------------------------FFTP0030
C                                                                       FFTP0040
C   COMPUTER            - PRIME/DOUBLE                                  FFTP0050
C                                                                       FFTP0060
C   LATEST REVISION     - JANUARY 1, 1978                               FFTP0070
C                                                                       FFTP0080
C   PURPOSE             - COMPUTE THE FAST FOURIER TRANSFORM OF A       FFTP0090
C                           COMPLEX VALUED SEQUENCE                     FFTP0100
C                                                                       FFTP0110
C   USAGE               - CALL FFTCC (A,N,IWK,WK)                       FFTP0120
C                                                                       FFTP0130
C   ARGUMENTS    A      - COMPLEX VECTOR OF LENGTH N. ON INPUT A        FFTP0140
C                           CONTAINS THE COMPLEX VALUED SEQUENCE TO BE  FFTP0150
C                           TRANSFORMED. ON OUTPUT A IS REPLACED BY THE FFTP0160
C                           FOURIER TRANSFORM.                          FFTP0170
C                N      - INPUT NUMBER OF DATA POINTS TO BE             FFTP0180
C                           TRANSFORMED. N MAY BE ANY POSITIVE          FFTP0190
C                           INTEGER.                                    FFTP0200
C                IWK    - INTEGER WORK VECTOR OF LENGTH 6*N+150.        FFTP0210
C                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS) FFTP0220
C                WK     - REAL WORK VECTOR OF LENGTH 6*N+150.           FFTP0230
C                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS) FFTP0240
C                                                                       FFTP0250
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         FFTP0260
C                       - SINGLE/H36,H48,H60                            FFTP0270
C                                                                       FFTP0280
C   REQD. IMSL ROUTINES - NONE REQUIRED                                 FFTP0290
C                                                                       FFTP0300
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           FFTP0310
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      FFTP0320
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  FFTP0330
C                                                                       FFTP0340
C   REMARKS  1.  FFTCC COMPUTES THE FOURIER TRANSFORM, X, ACCORDING     FFTP0350
C                TO THE FOLLOWING FORMULA;                              FFTP0360
C                                                                       FFTP0370
C                  X(K+1) = SUM FROM J = 0 TO N-1 OF                    FFTP0380
C                           A(J+1)*CEXP((0.0,(2.0*PI*J*K)/N))           FFTP0390
C                  FOR K=0,1,...,N-1 AND PI=3.1415...                   FFTP0400
C                                                                       FFTP0410
C                NOTE THAT X OVERWRITES A ON OUTPUT.                    FFTP0420
C            2.  FFTCC CAN BE USED TO COMPUTE                           FFTP0430
C                                                                       FFTP0440
C                  X(K+1) = (1/N)*SUM FROM J = 0 TO N-1 OF              FFTP0450
C                           A(J+1)*CEXP((0.0,(-2.0*PI*J*K)/N))          FFTP0460
C                  FOR K=0,1,...,N-1 AND PI=3.1415...                   FFTP0470
C                                                                       FFTP0480
C                BY PERFORMING THE FOLLOWING STEPS;                     FFTP0490
C                                                                       FFTP0500
C                     DO 10 I=1,N                                       FFTP0510
C                        A(I) = CONJG(A(I))                             FFTP0520
C                  10 CONTINUE                                          FFTP0530
C                     CALL FFTCC (A,N,IWK,WK)                           FFTP0540
C                     DO 20 I=1,N                                       FFTP0550
C                        A(I) = CONJG(A(I))/N                           FFTP0560
C                  20 CONTINUE                                          FFTP0570
C                                                                       FFTP0580
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       FFTP0590
C                                                                       FFTP0600
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN FFTP0610
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    FFTP0620
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        FFTP0630
C                                                                       FFTP0640
C-----------------------------------------------------------------------FFTP0650
C                                                                       FFTP0660
      SUBROUTINE FFTCC (A,N,IWK,WK)                                     FFTP0670
C                                  SPECIFICATIONS FOR ARGUMENTS         FFTP0680
      INTEGER            N,IWK(1)                                       FFTP0690
      DOUBLE PRECISION   WK(1)                                          FFTP0700
      COMPLEX*16         A(N)                                           FFTP0710
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   FFTP0720
      INTEGER            I,IAM,IAP,IBM,IBP,IC,ICC,ICF,ICK,ID,IDM1,II,   FFTP0730
     1                   IJA,IKB,IKT,ILL,IM,IRD,ISF,ISK,ISP,ISS,ITA,ITB,FFTP0740
     2                   J,JA,JF,JJ,JK,K,K0,K1,K2,K3,KA,KB,KD2,KF,KH,KN,FFTP0750
     3                   KT,KTP,L,L1,M,MM,MM1,MP                        FFTP0760
      DOUBLE PRECISION   CM,SM,C1,C2,C3,S1,S2,S3,C30,RAD,A0,A1,A4,B4,   FFTP0770
     1                   A2,A3,B0,B1,B2,B3,ZERO,HALF,ONE,TWO,Z0(2),     FFTP0780
     2                   Z1(2),Z2(2),Z3(2),Z4(2)                        FFTP0790
      COMPLEX*16         ZA0,ZA1,ZA2,ZA3,ZA4,AK2                        FFTP0800
      EQUIVALENCE        (ZA0,Z0(1)),(ZA1,Z1(1)),(ZA2,Z2(1)),           FFTP0810
     1                   (ZA3,Z3(1)),(A0,Z0(1)),(B0,Z0(2)),(A1,Z1(1)),  FFTP0820
     2                   (B1,Z1(2)),(A2,Z2(1)),(B2,Z2(2)),(A3,Z3(1)),   FFTP0830
     3                   (B3,Z3(2)),(ZA4,Z4(1)),(Z4(1),A4),(Z4(2),B4)   FFTP0840
      DATA               RAD/6.2831853071796D0/,                        FFTP0850
     1                   C30/.86602540378444D0/                         FFTP0860
      DATA               ZERO,HALF,ONE,TWO/0.0D0,0.5D0,1.0D0,2.0D0/     FFTP0870
C                                  FIRST EXECUTABLE STATEMENT           FFTP0880
      IF (N .EQ. 1) GO TO 9005                                          FFTP0890
      K = N                                                             FFTP0900
      M = 0                                                             FFTP0910
      J = 2                                                             FFTP0920
      JJ = 4                                                            FFTP0930
      JF = 0                                                            FFTP0940
C                                  DETERMINE THE SQUARE FACTORS OF N    FFTP0950
      IWK(1) = 1                                                        FFTP0960
    5 I = K/JJ                                                          FFTP0970
      IF (I*JJ .NE. K) GO TO 10                                         FFTP0980
      M = M+1                                                           FFTP0990
      IWK(M+1) = J                                                      FFTP1000
      K = I                                                             FFTP1010
      GO TO 5                                                           FFTP1020
   10 J = J + 2                                                         FFTP1030
      IF (J .EQ. 4) J = 3                                               FFTP1040
      JJ = J * J                                                        FFTP1050
      IF (JJ .LE. K) GO TO 5                                            FFTP1060
      KT = M                                                            FFTP1070
C                                  DETERMINE THE REMAINING FACTORS OF N FFTP1080
      J = 2                                                             FFTP1090
   15 I = K / J                                                         FFTP1100
      IF (I*J .NE. K) GO TO 20                                          FFTP1110
      M = M + 1                                                         FFTP1120
      IWK(M+1) = J                                                      FFTP1130
      K = I                                                             FFTP1140
      GO TO 15                                                          FFTP1150
   20 J = J + 1                                                         FFTP1160
      IF (J .EQ. 3) GO TO 15                                            FFTP1170
      J = J + 1                                                         FFTP1180
      IF(J.LE.K) GO TO 15                                               FFTP1190
      K = IWK(M+1)                                                      FFTP1200
      IF (IWK(KT+1) .GT. IWK(M+1)) K = IWK(KT+1)                        FFTP1210
      IF(KT.LE.0) GO TO 30                                              FFTP1220
      KTP = KT + 2                                                      FFTP1230
      DO 25  I = 1,KT                                                   FFTP1240
         J = KTP - I                                                    FFTP1250
         M = M+1                                                        FFTP1260
         IWK(M+1) = IWK(J)                                              FFTP1270
   25 CONTINUE                                                          FFTP1280
   30 MP = M+1                                                          FFTP1290
      IC = MP+1                                                         FFTP1300
      ID = IC+MP                                                        FFTP1310
      ILL = ID+MP                                                       FFTP1320
      IRD = ILL+MP+1                                                    FFTP1330
      ICC = IRD+MP                                                      FFTP1340
      ISS = ICC+MP                                                      FFTP1350
      ICK = ISS+MP                                                      FFTP1360
      ISK = ICK+K                                                       FFTP1370
      ICF = ISK+K                                                       FFTP1380
      ISF = ICF+K                                                       FFTP1390
      IAP = ISF+K                                                       FFTP1400
      KD2 = (K-1) / 2 + 1                                               FFTP1410
      IBP = IAP + KD2                                                   FFTP1420
      IAM = IBP + KD2                                                   FFTP1430
      IBM = IAM + KD2                                                   FFTP1440
      MM1 = M-1                                                         FFTP1450
      I = 1                                                             FFTP1460
   35 L = MP - I                                                        FFTP1470
      J = IC - I                                                        FFTP1480
      IWK(ILL+L) = 0                                                    FFTP1490
      IF ((IWK(J-1) + IWK(J)) .EQ. 4) IWK(ILL+L) = 1                    FFTP1500
      IF (IWK(ILL+L) .EQ. 0) GO TO 40                                   FFTP1510
      I = I + 1                                                         FFTP1520
      L = L - 1                                                         FFTP1530
      IWK(ILL+L) = 0                                                    FFTP1540
   40 I = I + 1                                                         FFTP1550
      IF(I.LE.MM1) GO TO 35                                             FFTP1560
      IWK(ILL+1) = 0                                                    FFTP1570
      IWK(ILL+MP) = 0                                                   FFTP1580
      IWK(IC) = 1                                                       FFTP1590
      IWK(ID) = N                                                       FFTP1600
      DO 45  J = 1,M                                                    FFTP1610
         K = IWK(J+1)                                                   FFTP1620
         IWK(IC+J) = IWK(IC+J-1) * K                                    FFTP1630
         IWK(ID+J) = IWK(ID+J-1) / K                                    FFTP1640
         WK(IRD+J) = RAD/IWK(IC+J)                                      FFTP1650
         C1 = RAD/K                                                     FFTP1660
         IF (K .LE. 2) GO TO 45                                         FFTP1670
         WK(ICC+J) = DCOS(C1)                                           FFTP1680
         WK(ISS+J) = DSIN(C1)                                           FFTP1690
   45 CONTINUE                                                          FFTP1700
      MM = M                                                            FFTP1710
      IF (IWK(ILL+M) .EQ. 1) MM = M - 1                                 FFTP1720
      IF (MM .LE. 1) GO TO 50                                           FFTP1730
      SM = IWK(IC+MM-2) * WK(IRD+M)                                     FFTP1740
      CM = DCOS(SM)                                                     FFTP1750
      SM = DSIN(SM)                                                     FFTP1760
   50 KB = 0                                                            FFTP1770
      KN = N                                                            FFTP1780
      JJ = 0                                                            FFTP1790
      I = 1                                                             FFTP1800
      C1 = ONE                                                          FFTP1810
      S1 = ZERO                                                         FFTP1820
      L1 = 1                                                            FFTP1830
   55 IF (IWK(ILL+I+1) .EQ. 1) GO TO 60                                 FFTP1840
      KF = IWK(I+1)                                                     FFTP1850
      GO TO 65                                                          FFTP1860
   60 KF = 4                                                            FFTP1870
      I = I+1                                                           FFTP1880
   65 ISP = IWK(ID+I)                                                   FFTP1890
      IF (L1 .EQ. 1) GO TO 70                                           FFTP1900
      S1 = JJ * WK(IRD+I)                                               FFTP1910
      C1 = DCOS(S1)                                                     FFTP1920
      S1 = DSIN(S1)                                                     FFTP1930
C                                  FACTORS OF 2, 3, AND 4 ARE           FFTP1940
C                                  HANDLED SEPARATELY.                  FFTP1950
   70 IF (KF .GT. 4) GO TO 140                                          FFTP1960
      GO TO (75,75,90,115), KF                                          FFTP1970
   75 K0 = KB + ISP                                                     FFTP1980
      K2 = K0 + ISP                                                     FFTP1990
      IF (L1 .EQ. 1) GO TO 85                                           FFTP2000
   80 K0 = K0 - 1                                                       FFTP2010
      IF (K0 .LT. KB) GO TO 190                                         FFTP2020
      K2 = K2 - 1                                                       FFTP2030
      ZA4 = A(K2+1)                                                     FFTP2040
      A0 = A4*C1-B4*S1                                                  FFTP2050
      B0 = A4*S1+B4*C1                                                  FFTP2060
      A(K2+1) = A(K0+1)-ZA0                                             FFTP2070
      A(K0+1) = A(K0+1)+ZA0                                             FFTP2080
      GO TO 80                                                          FFTP2090
   85 K0 = K0 - 1                                                       FFTP2100
      IF (K0 .LT. KB) GO TO 190                                         FFTP2110
      K2 = K2 - 1                                                       FFTP2120
      AK2 = A(K2+1)                                                     FFTP2130
      A(K2+1) = A(K0+1)-AK2                                             FFTP2140
      A(K0+1) = A(K0+1)+AK2                                             FFTP2150
      GO TO 85                                                          FFTP2160
   90 IF (L1 .EQ. 1) GO TO 95                                           FFTP2170
      C2 = C1 * C1 - S1 * S1                                            FFTP2180
      S2 = TWO * C1 * S1                                                FFTP2190
   95 JA = KB + ISP - 1                                                 FFTP2200
      KA = JA + KB                                                      FFTP2210
      IKB = KB+1                                                        FFTP2220
      IJA = JA+1                                                        FFTP2230
      DO 110 II = IKB,IJA                                               FFTP2240
         K0 = KA - II + 1                                               FFTP2250
         K1 = K0 + ISP                                                  FFTP2260
         K2 = K1 + ISP                                                  FFTP2270
         ZA0 = A(K0+1)                                                  FFTP2280
         IF (L1 .EQ. 1) GO TO 100                                       FFTP2290
         ZA4 = A(K1+1)                                                  FFTP2300
         A1 = A4*C1-B4*S1                                               FFTP2310
         B1 = A4*S1+B4*C1                                               FFTP2320
         ZA4 = A(K2+1)                                                  FFTP2330
         A2 = A4*C2-B4*S2                                               FFTP2340
         B2 = A4*S2+B4*C2                                               FFTP2350
         GO TO 105                                                      FFTP2360
  100    ZA1 = A(K1+1)                                                  FFTP2370
         ZA2 = A(K2+1)                                                  FFTP2380
  105    A(K0+1) = DCMPLX(A0+A1+A2,B0+B1+B2)                            FFTP2390
         A0 = -HALF * (A1+A2) + A0                                      FFTP2400
         A1 = (A1-A2) * C30                                             FFTP2410
         B0 = -HALF * (B1+B2) + B0                                      FFTP2420
         B1 = (B1-B2) * C30                                             FFTP2430
         A(K1+1) = DCMPLX(A0-B1,B0+A1)                                  FFTP2440
         A(K2+1) = DCMPLX(A0+B1,B0-A1)                                  FFTP2450
  110 CONTINUE                                                          FFTP2460
      GO TO 190                                                         FFTP2470
  115 IF (L1 .EQ. 1) GO TO 120                                          FFTP2480
      C2 = C1 * C1 - S1 * S1                                            FFTP2490
      S2 = TWO * C1 * S1                                                FFTP2500
      C3 = C1 * C2 - S1 * S2                                            FFTP2510
      S3 = S1 * C2 + C1 * S2                                            FFTP2520
  120 JA = KB + ISP - 1                                                 FFTP2530
      KA = JA + KB                                                      FFTP2540
      IKB = KB+1                                                        FFTP2550
      IJA = JA+1                                                        FFTP2560
      DO 135 II = IKB,IJA                                               FFTP2570
         K0 = KA - II + 1                                               FFTP2580
         K1 = K0 + ISP                                                  FFTP2590
         K2 = K1 + ISP                                                  FFTP2600
         K3 = K2 + ISP                                                  FFTP2610
         ZA0 = A(K0+1)                                                  FFTP2620
         IF (L1 .EQ. 1) GO TO 125                                       FFTP2630
         ZA4 = A(K1+1)                                                  FFTP2640
         A1 = A4*C1-B4*S1                                               FFTP2650
         B1 = A4*S1+B4*C1                                               FFTP2660
         ZA4 = A(K2+1)                                                  FFTP2670
         A2 = A4*C2-B4*S2                                               FFTP2680
         B2 = A4*S2+B4*C2                                               FFTP2690
         ZA4 = A(K3+1)                                                  FFTP2700
         A3 = A4*C3-B4*S3                                               FFTP2710
         B3 = A4*S3+B4*C3                                               FFTP2720
         GO TO 130                                                      FFTP2730
  125    ZA1 = A(K1+1)                                                  FFTP2740
         ZA2 = A(K2+1)                                                  FFTP2750
         ZA3 = A(K3+1)                                                  FFTP2760
  130    A(K0+1) = DCMPLX(A0+A2+A1+A3,B0+B2+B1+B3)                      FFTP2770
         A(K1+1) = DCMPLX(A0+A2-A1-A3,B0+B2-B1-B3)                      FFTP2780
         A(K2+1) = DCMPLX(A0-A2-B1+B3,B0-B2+A1-A3)                      FFTP2790
         A(K3+1) = DCMPLX(A0-A2+B1-B3,B0-B2-A1+A3)                      FFTP2800
  135 CONTINUE                                                          FFTP2810
      GO TO 190                                                         FFTP2820
  140 JK = KF - 1                                                       FFTP2830
      KH = JK/2                                                         FFTP2840
      K3 = IWK(ID+I-1)                                                  FFTP2850
      K0 = KB + ISP                                                     FFTP2860
      IF (L1 .EQ. 1) GO TO 150                                          FFTP2870
      K = JK - 1                                                        FFTP2880
      WK(ICF+1) = C1                                                    FFTP2890
      WK(ISF+1) = S1                                                    FFTP2900
      DO 145 J = 1,K                                                    FFTP2910
         WK(ICF+J+1) = WK(ICF+J) * C1 - WK(ISF+J) * S1                  FFTP2920
         WK(ISF+J+1) = WK(ICF+J) * S1 + WK(ISF+J) * C1                  FFTP2930
  145 CONTINUE                                                          FFTP2940
  150 IF (KF .EQ. JF) GO TO 160                                         FFTP2950
      C2 = WK(ICC+I)                                                    FFTP2960
      WK(ICK+1) = C2                                                    FFTP2970
      WK(ICK+JK) = C2                                                   FFTP2980
      S2 = WK(ISS+I)                                                    FFTP2990
      WK(ISK+1) = S2                                                    FFTP3000
      WK(ISK+JK) = -S2                                                  FFTP3010
      DO 155 J = 1,KH                                                   FFTP3020
         K = JK - J                                                     FFTP3030
         WK(ICK+K) = WK(ICK+J) * C2 - WK(ISK+J) * S2                    FFTP3040
         WK(ICK+J+1) = WK(ICK+K)                                        FFTP3050
         WK(ISK+J+1) = WK(ICK+J) * S2 + WK(ISK+J) * C2                  FFTP3060
         WK(ISK+K) = -WK(ISK+J+1)                                       FFTP3070
  155 CONTINUE                                                          FFTP3080
  160 K0 = K0 - 1                                                       FFTP3090
      K1 = K0                                                           FFTP3100
      K2 = K0 + K3                                                      FFTP3110
      ZA0 = A(K0+1)                                                     FFTP3120
      A3 = A0                                                           FFTP3130
      B3 = B0                                                           FFTP3140
      DO 175 J = 1,KH                                                   FFTP3150
         K1 = K1 + ISP                                                  FFTP3160
         K2 = K2 - ISP                                                  FFTP3170
         IF (L1 .EQ. 1) GO TO 165                                       FFTP3180
         K = KF - J                                                     FFTP3190
         ZA4 = A(K1+1)                                                  FFTP3200
         A1 = A4*WK(ICF+J)-B4*WK(ISF+J)                                 FFTP3210
         B1 = A4*WK(ISF+J)+B4*WK(ICF+J)                                 FFTP3220
         ZA4 = A(K2+1)                                                  FFTP3230
         A2 = A4*WK(ICF+K)-B4*WK(ISF+K)                                 FFTP3240
         B2 = A4*WK(ISF+K)+B4*WK(ICF+K)                                 FFTP3250
         GO TO 170                                                      FFTP3260
  165    ZA1 = A(K1+1)                                                  FFTP3270
         ZA2 = A(K2+1)                                                  FFTP3280
  170    WK(IAP+J) = A1 + A2                                            FFTP3290
         WK(IAM+J) = A1 - A2                                            FFTP3300
         WK(IBP+J) = B1 + B2                                            FFTP3310
         WK(IBM+J) = B1 - B2                                            FFTP3320
         A3 = A1 + A2 + A3                                              FFTP3330
         B3 = B1 + B2 + B3                                              FFTP3340
  175 CONTINUE                                                          FFTP3350
      A(K0+1) = DCMPLX(A3,B3)                                           FFTP3360
      K1 = K0                                                           FFTP3370
      K2 = K0 + K3                                                      FFTP3380
      DO 185 J = 1,KH                                                   FFTP3390
         K1 = K1 + ISP                                                  FFTP3400
         K2 = K2 - ISP                                                  FFTP3410
         JK = J                                                         FFTP3420
         A1 = A0                                                        FFTP3430
         B1 = B0                                                        FFTP3440
         A2 = ZERO                                                      FFTP3450
         B2 = ZERO                                                      FFTP3460
         DO 180  K = 1,KH                                               FFTP3470
            A1 = A1 + WK(IAP+K) * WK(ICK+JK)                            FFTP3480
            A2 = A2 + WK(IAM+K) * WK(ISK+JK)                            FFTP3490
            B1 = B1 + WK(IBP+K) * WK(ICK+JK)                            FFTP3500
            B2 = B2 + WK(IBM+K) * WK(ISK+JK)                            FFTP3510
            JK = JK + J                                                 FFTP3520
            IF (JK .GE. KF) JK = JK - KF                                FFTP3530
  180    CONTINUE                                                       FFTP3540
         A(K1+1) = DCMPLX(A1-B2,B1+A2)                                  FFTP3550
         A(K2+1) = DCMPLX(A1+B2,B1-A2)                                  FFTP3560
  185 CONTINUE                                                          FFTP3570
      IF (K0 .GT. KB) GO TO 160                                         FFTP3580
      JF = KF                                                           FFTP3590
  190 IF ( I .GE. MM ) GO TO 195                                        FFTP3600
      I = I + 1                                                         FFTP3610
      GO TO 55                                                          FFTP3620
  195 I = MM                                                            FFTP3630
      L1 = 0                                                            FFTP3640
      KB = IWK(ID+I-1) + KB                                             FFTP3650
      IF (KB .GE. KN) GO TO 215                                         FFTP3660
  200 JJ = IWK(IC+I-2) + JJ                                             FFTP3670
      IF (JJ .LT. IWK(IC+I-1)) GO TO 205                                FFTP3680
      I = I - 1                                                         FFTP3690
      JJ = JJ - IWK(IC+I)                                               FFTP3700
      GO TO 200                                                         FFTP3710
  205 IF (I .NE. MM) GO TO 210                                          FFTP3720
      C2 = C1                                                           FFTP3730
      C1 = CM * C1 - SM * S1                                            FFTP3740
      S1 = SM * C2 + CM * S1                                            FFTP3750
      GO TO 70                                                          FFTP3760
  210 IF (IWK(ILL+I) .EQ. 1) I = I + 1                                  FFTP3770
      GO TO 55                                                          FFTP3780
  215 I = 1                                                             FFTP3790
      JA = KT - 1                                                       FFTP3800
      KA = JA + 1                                                       FFTP3810
      IF(JA.LT.1) GO TO 225                                             FFTP3820
      DO 220  II = 1,JA                                                 FFTP3830
         J = KA - II                                                    FFTP3840
         IWK(J+1) = IWK(J+1) - 1                                        FFTP3850
         I = IWK(J+1) + I                                               FFTP3860
  220 CONTINUE                                                          FFTP3870
C                                  THE RESULT IS NOW PERMUTED TO        FFTP3880
C                                  NORMAL ORDER.                        FFTP3890
  225 IF (KT .LE. 0) GO TO 270                                          FFTP3900
      J = 1                                                             FFTP3910
      I = 0                                                             FFTP3920
      KB = 0                                                            FFTP3930
  230 K2 = IWK(ID+J) + KB                                               FFTP3940
      K3 = K2                                                           FFTP3950
      JJ = IWK(IC+J-1)                                                  FFTP3960
      JK = JJ                                                           FFTP3970
      K0 = KB + JJ                                                      FFTP3980
      ISP = IWK(IC+J) - JJ                                              FFTP3990
  235 K = K0 + JJ                                                       FFTP4000
  240 ZA4 = A(K0+1)                                                     FFTP4010
      A(K0+1) = A(K2+1)                                                 FFTP4020
      A(K2+1) = ZA4                                                     FFTP4030
      K0 = K0 + 1                                                       FFTP4040
      K2 = K2 + 1                                                       FFTP4050
      IF (K0 .LT. K) GO TO 240                                          FFTP4060
      K0 = K0 + ISP                                                     FFTP4070
      K2 = K2 + ISP                                                     FFTP4080
      IF (K0 .LT. K3) GO TO 235                                         FFTP4090
      IF (K0 .GE. K3 + ISP) GO TO 245                                   FFTP4100
      K0 = K0 - IWK(ID+J) + JJ                                          FFTP4110
      GO TO 235                                                         FFTP4120
  245 K3 = IWK(ID+J) + K3                                               FFTP4130
      IF (K3 - KB .GE. IWK(ID+J-1)) GO TO 250                           FFTP4140
      K2 = K3 + JK                                                      FFTP4150
      JK = JK + JJ                                                      FFTP4160
      K0 = K3 - IWK(ID+J) + JK                                          FFTP4170
      GO TO 235                                                         FFTP4180
  250 IF (J .GE. KT) GO TO 260                                          FFTP4190
      K = IWK(J+1) + I                                                  FFTP4200
      J = J + 1                                                         FFTP4210
  255 I = I + 1                                                         FFTP4220
      IWK(ILL+I) = J                                                    FFTP4230
      IF (I .LT. K) GO TO 255                                           FFTP4240
      GO TO 230                                                         FFTP4250
  260 KB = K3                                                           FFTP4260
      IF (I .LE. 0) GO TO 265                                           FFTP4270
      J = IWK(ILL+I)                                                    FFTP4280
      I = I - 1                                                         FFTP4290
      GO TO 230                                                         FFTP4300
  265 IF (KB .GE. N) GO TO 270                                          FFTP4310
      J = 1                                                             FFTP4320
      GO TO 230                                                         FFTP4330
  270 JK = IWK(IC+KT)                                                   FFTP4340
      ISP = IWK(ID+KT)                                                  FFTP4350
      M = M - KT                                                        FFTP4360
      KB = ISP/JK-2                                                     FFTP4370
      IF (KT .GE. M-1 ) GO TO 9005                                      FFTP4380
      ITA = ILL+KB+1                                                    FFTP4390
      ITB = ITA+JK                                                      FFTP4400
      IDM1 = ID-1                                                       FFTP4410
      IKT = KT+1                                                        FFTP4420
      IM = M+1                                                          FFTP4430
      DO 275 J = IKT,IM                                                 FFTP4440
         IWK(IDM1+J) = IWK(IDM1+J)/JK                                   FFTP4450
  275 CONTINUE                                                          FFTP4460
      JJ = 0                                                            FFTP4470
      DO 290 J = 1,KB                                                   FFTP4480
         K = KT                                                         FFTP4490
  280    JJ = IWK(ID+K+1) + JJ                                          FFTP4500
         IF (JJ .LT. IWK(ID+K)) GO TO 285                               FFTP4510
         JJ = JJ - IWK(ID+K)                                            FFTP4520
         K = K + 1                                                      FFTP4530
         GO TO 280                                                      FFTP4540
  285    IWK(ILL+J) = JJ                                                FFTP4550
         IF (JJ .EQ. J) IWK(ILL+J) = -J                                 FFTP4560
  290 CONTINUE                                                          FFTP4570
C                                  DETERMINE THE PERMUTATION CYCLES     FFTP4580
C                                  OF LENGTH GREATER THAN OR EQUAL      FFTP4590
C                                  TO TWO.                              FFTP4600
      DO 300  J = 1,KB                                                  FFTP4610
         IF (IWK(ILL+J) .LE. 0) GO TO 300                               FFTP4620
         K2 = J                                                         FFTP4630
  295    K2 = IABS(IWK(ILL+K2))                                         FFTP4640
         IF (K2 .EQ. J) GO TO 300                                       FFTP4650
         IWK(ILL+K2) = -IWK(ILL+K2)                                     FFTP4660
         GO TO 295                                                      FFTP4670
  300 CONTINUE                                                          FFTP4680
C                                  REORDER A FOLLOWING THE              FFTP4690
C                                  PERMUTATION CYCLES                   FFTP4700
      I = 0                                                             FFTP4710
      J = 0                                                             FFTP4720
      KB = 0                                                            FFTP4730
      KN = N                                                            FFTP4740
  305 J = J + 1                                                         FFTP4750
      IF (IWK(ILL+J) .LT. 0) GO TO 305                                  FFTP4760
      K = IWK(ILL+J)                                                    FFTP4770
      K0 = JK * K + KB                                                  FFTP4780
  310 ZA4 = A(K0+I+1)                                                   FFTP4790
      WK(ITA+I) = A4                                                    FFTP4800
      WK(ITB+I) = B4                                                    FFTP4810
      I = I + 1                                                         FFTP4820
      IF (I .LT. JK) GO TO 310                                          FFTP4830
      I = 0                                                             FFTP4840
  315 K = -IWK(ILL+K)                                                   FFTP4850
      JJ = K0                                                           FFTP4860
      K0 = JK * K + KB                                                  FFTP4870
  320 A(JJ+I+1) = A(K0+I+1)                                             FFTP4880
      I = I + 1                                                         FFTP4890
      IF (I .LT. JK) GO TO 320                                          FFTP4900
      I = 0                                                             FFTP4910
      IF (K .NE. J) GO TO 315                                           FFTP4920
  325 A(K0+I+1) = DCMPLX(WK(ITA+I),WK(ITB+I))                           FFTP4930
      I = I + 1                                                         FFTP4940
      IF (I .LT. JK) GO TO 325                                          FFTP4950
      I = 0                                                             FFTP4960
      IF (J .LT. K2) GO TO 305                                          FFTP4970
      J = 0                                                             FFTP4980
      KB = KB + ISP                                                     FFTP4990
      IF (KB .LT. KN) GO TO 305                                         FFTP5000
 9005 RETURN                                                            FFTP5010
      END                                                               FFTP5020
