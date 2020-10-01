      implicit double precision (a-h,o-z)
c lm = maximum value of l for expansion (all m's included)
      parameter (ntheta=3,nphi=3,ngrid=ntheta*nphi,lm=2)
      dimension thetar(ngrid),phir(ngrid)
     &    ,tes((lm+1)*(lm+1)), tesseral((lm+1)*(lm+1),ngrid)
     &    ,vv(ngrid),c((lm+1)*(lm+1)),tesseral0((lm+1)*(lm+1),ngrid)
     &    ,unity((lm+1)*(lm+1),ngrid)
      
      open(unit=2, file='expansion.out',form='formatted')

C GRID FOR AR-CH4 THETA, PHI
      
      pi=acos(-1.d0)
      
      thetamin = 10.d0
      thetamax = 170.d0
      dtheta = (thetamax-thetamin)/dble(ntheta-1)
      phimin = 0.5d0
      phimax = 359.5d0
      dphi = (phimax-phimin)/dble(nphi-1)

      write(2,*)'  Initial grid in theta, phi (radians)'
      j=0
      do itheta= 1, ntheta
         tthetar = (thetamin + dble(itheta-1)*dtheta)*pi/180.d0
         do iphi= 1, nphi
            j=j+1
            thetar(j) = tthetar
            phir(j) = (phimin + dble(iphi-1)*dphi)*pi/180.d0
            write (2,*)'j, theta,phi', j, thetar(j), phir(j)
         end do
      end do


C  TESSERAL HARMONICS FOR THETA Y PHI
CC  INCLUDE RACAH NORMALIZE 

      do i=1,(lm+1)*(lm+1)
         do j=1,ngrid
            tesseral(i,j)=0.d0
         end do
      end do

      pi4=4.d0*acos(-1.d0)
      do j=1,ngrid
         call tesser (tes,lm,thetar(j),phir(j))
CC  RACAH NORMALIZE (elements and restorage of 
CC  matrix tesseral have been checked after and
CC  before Racah normalize 
         ii=0
         do l=0,lm
            fac=sqrt(pi4/dble((2*l+1)))
            do m=-l,l
               ii=ii+1
               tes(ii)=tes(ii)*fac
cN             write (2,*) 'theta,phi,l,m,tesseral', 
cN     &                  thetar(j),phir(j),l,m,tes(ii)
            end do
         end do
         do i=1,(lm+1)*(lm+1)
cN            tesseral(i,j)=tes(i)
            tesseral(j,i)=tes(i)
         end do
      end do

      write(2,*)
      write(2,*) 'theta, phi, tesseral'
      do i=1,(lm+1)*(lm+1)
         do j=1, ngrid
            write (2,*) thetar(j),phir(j), tesseral(j,i)
         end do
      end do

      write(2,*)
      write(2,*) 'matrix tesseral before inversion'
      do igrid=1,ngrid
         write(2,*)
         write(2,*) (tesseral (igrid,lam),lam=1,(lm+1)*(lm+1))
      end do

      do i=1,(lm+1)*(lm+1)
         do j=1,ngrid
            tesseral0(i,j)=tesseral(i,j)
         end do
      end do

C  DIMENSION AND INVERSION TESSERAL
     
      ndimr=(lm+1)*(lm+1)
      ndimf=ngrid
      write(2,*)
      write (2,*) 'number of rows, l-m',ndimr
      write (2,*) 'number of columns, grid', ndimf

      call inv (tesseral,ndimr,ndimr)
   
      write(2,*) 'inverted Tesseral matrix'

      do i=1,(lm+1)*(lm+1)
         write (2,*) (tesseral(i,j), j=1,ngrid)
      end do  
    
      do i=1,(lm+1)*(lm+1)
         do j=1,ngrid
            unity(i,j)=0.d0
         end do
      end do

      do i=1,(lm+1)*(lm+1)
         do j=1,ngrid
            do k=1,ngrid 
               unity(i,j)=unity(i,j)+tesseral0(i,k)*tesseral(k,j)
            end do
         end do
      end do

      write(2,*)
      write (2,*) 'check: T * T**(-1) = '

      do i=1,(lm+1)*(lm+1)
         write (2,*) (unity(i,j), j=1,ngrid)
      end do  

C    CALCULATION OF MY POTENTIAL


      write(2,*) 'calculation of my potential'

      r=5.d0
         
         do i=1,ngrid
            call mypot (r,thetar(i),phir(i),vv(i))
         end do
 
         write(2,*)
         write(2,*) 'potential at r = ',r
         write(2,*)'vv', (vv(i),i=1,ngrid)

CC   COEFFICIENTS C ARE CALCULATED AS VV*INV MATRIX TESSERAL HARMONICS 

         do l=1,(lm+1)*(lm+1)
            c(l)=0.d0
            do m=1,ngrid
               c(l)=vv(m)*tesseral(l,m)+c(l)
            end do
         end do
         write(2,*)
         write (2,*) 'resulting coefficients at r=', r
         write (2,"(10(1x,g22.15))") (c(l), l=1,(lm+1)*(lm+1))

 
      stop
      end 


      SUBROUTINE TESSER(TES,N,THETA,PHI)

C     Subroutine computes tesseral harmonics as defined by J.L. Prather
C     [N.B.S. monograph 19 (1961)]. The tesserals are computed from
C     l=m=0 upwards to l=m=n for the polar angles theta and phi. The
C     results are stored in the linear array tes in increasing order of
C     l and m. Negative m-values refer to Prather's S(l,m)
C     [= P(l,m)*sin(m*phi)]. Positive m-values to Prather's C(l,m)
C     [= P(l,m)*cos(m*phi)].
C
C     Notes:
C        The tesseral T(l,m,theta,phi) is accessed via tes(l*(l+1)+m+1).
C        Angles in radians.
C
C     Author: P.E.S. Wormer (1982)

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TES(*),P(351),SN(26),CS(26)
C     DIMENSION TES( (N+1)*(N+1) ),P( (N+1)*(N+2)/2 ),SN(N+1),
C    1               CS(N+1)
      REAL*8 TPIH/3.989422804014327D-1/
C     (   TPIH = 1/DSQRT(2*PI)   )
      REAL*8 TWMH/7.071067811865475D-1/
C     (   TWMH = 1./DSQRT(2)    )
C
C
C     COMPUTE ASSOCIATED LEGENDRE FUNCTIONS
C
      CALL ASSLEG(P,N,DCOS(THETA))
C
      TES(1) = TPIH * TWMH * P(1)
c      write (2,*) 'tes(1)=',tes(1)
      IF ( N .EQ. 0) RETURN
C
C     COMPUTE NECESSARY SINES AND COSINES
C
      CS(1) = 1.D0
      SN(1) = 0.D0
      COSP = DCOS(PHI)
      SINP = DSIN(PHI)
      DO 10 M=1,N
      CS(M+1) = CS(M)*COSP - SN(M)*SINP
      SN(M+1) = SN(M)*COSP + CS(M)*SINP
   10 CONTINUE
C
      KP = 1
      LL = 1
      DL = 0.D0
C
      DO 40 L=1,N
      LL = LL + L + L
C     ( LL = L(L+1) + 1   )
      KP = KP + 1
      DL = DL + 1.D0
      FAC = DSQRT(DL+DL+1.D0) * TPIH
c      write(2,*)'ll=',ll,'kp=',kp,'fac=',fac
      TES(LL) = FAC*TWMH * P(KP)
c      write(2,*)'tes(ll)',tes(ll)
      KT1 = LL
      KT2 = LL
      DLM = DL
      DLM1 = DL + 1.D0
C
      DO 30 M=1,L
      KT1 = KT1 + 1
      KT2 = KT2 - 1
C     ( KT1 = L(L+1) +1 + M  )
C     ( KT2 = L(L+1) +1 - M )
      KP = KP + 1
      DLM  = DLM  + 1.D0
      DLM1 = DLM1 - 1.D0
C     ( DLM = L+M      )
C     ( DLM1= L+1-M    )
c      write(2,*)'dlm=',dlm,'dlm1=',dlm1
      FAC = FAC / DSQRT(DLM*DLM1)
c      write(2,*)'fac/dsqrt(dlm*dlm1)',fac
c      write(2,*)'kt1=',kt1,'kt2=',kt2,'kp=',kp,'m+1=',m+1
      TES(KT1) = FAC * P(KP) * CS(M+1)
      TES(KT2) = FAC * P(KP) * SN(M+1)
c      write(2,*)'tes(kt1)=',tes(kt1)
c      write(2,*)'tes(kt2)=',tes(kt2)
C     ( T(L,M) = FAC * P(L,M) * COS(M*PHI)  )
C     ( T(L,-M)= FAC * P(L,M) * SIN(M*PHI)  )
C
   30 CONTINUE
   40 CONTINUE
C
      RETURN
      END

      SUBROUTINE ASSLEG(P,N,X)
C
C     Author: P.E.S. Wormer (1982)
C
C
C     Subroutine computes associated Legendre polynomials as defined
C     by A.R. Edmonds (angular momentum in quantum mechanics).
C     x is the usual coordinate (-1 < x < +1 ), n is the maximum
C     quantum number. The associated Legendre functions are computed
C     in the order (0,0),(1,0),(1,1),(2,0),(2,1),(2,2),.... ,(n,n)
C     and returned in the array P.
C     The associated legendre function p(l,m,x) may be  accessed via
C     P( l*(l+1)/2 + m+1 )
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P((N+1)*(N+2)/2 )
C
      P(1) = 1.D0
      IF (N .EQ. 0) RETURN
      SINT = DSQRT(1.D0 - X*X)
      P(2) = X
      P(3) = SINT
      IF (N .EQ. 1) RETURN
C
      LM1 = 1
      LM  = 3
      DL  = 1.D0
C
      DO 20 L=2,N
      DL = DL + 1.D0
      LM1 = LM1 + 1
C     (   LM1 = L*(L-1)/2 + 1 )
      LM  = LM + 1
C     (   LM = L*(L+1)/2 + 1 )
C
      P(LM) = X*P(LM1) - SINT*P(LM1+1)/DL
C     ( P(L,0) = X*P(L-1,0) - DSQRT(1-X*X)*P(L-1,1)/L   )
      MMAX = L-1
      DLM = DL
      DO 10 M=1,MMAX
      LM1 = LM1 + 1
      LM  = LM + 1
C
      P(LM) = DLM*SINT*P(LM1-1) + X*P(LM1)
C     (  P(L,M) = (L+M-1)*DSQRT(1-X*X)*P(L-1,M-1) + X*P(L-1,M)   )
      DLM = DLM + 1.D0
   10 CONTINUE
      LM = LM + 1
      P(LM) = (DL+DL-1.D0)*SINT*P(LM1)
C     (  P(L,L) = (2*L-1)*DSQRT(1-X*X)*P(L-1,L-1)    )
C
   20 CONTINUE
C
      RETURN
      END


      SUBROUTINE INV(A,N,NDIMA)
CCC      SUBROUTINE INV(A,N)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NMAX=181 )
CCC       DIMENSION A(N,N)
       DIMENSION A(NDIMA,NDIMA)
      DIMENSION PIVOT(NMAX), IPIVOT(NMAX), INDEX(NMAX,2)
       ABS(X)=DABS(X)
      IF (N.GT.NMAX) GO TO 999
      IF(N.EQ.1) GO TO 741
   15 DO 20 J=1,N
   20 IPIVOT(J) = 0
   30 DO 550 I= 1,N
   40 AMAX =0.0
   45 DO 105 J=1,N
   50 IF(IPIVOT(J)-1) 60,105,60
   60 DO 100 K=1,N
   70 IF(IPIVOT(K)-1) 80,100,740
   80 IF( ABS(AMAX)- ABS(A(J,K))) 85,100,100
   85 IROW=J
   90 ICOLUM=K
   95 AMAX=A(J,K)
  100 CONTINUE
  105 CONTINUE
  110 IPIVOT(ICOLUM)=IPIVOT(ICOLUM)+1
  130 IF(IROW-ICOLUM) 150,260,150
  150 DO 200 L=1,N
  160 SWAP= A(IROW,L)
  170 A(IROW,L) = A(ICOLUM,L)
  200 A(ICOLUM,L) = SWAP
  260 INDEX (I,1) =IROW
  270 INDEX(I,2)=ICOLUM
  310 PIVOT(I) =A(ICOLUM,ICOLUM)
  330 A(ICOLUM,ICOLUM) =1.0 D0
  340 DO 350  L=1,N
  350  A(ICOLUM,L)=A(ICOLUM,L)/PIVOT(I)
  380 DO 550 L1 =1,N
  390 IF(L1-ICOLUM) 400,550,400
  400 T =A(L1,ICOLUM)
  420 A(L1,ICOLUM) = 0.0  D0
  430 DO 450  L=1,N
  450 A(L1,L) = A(L1,L)-A(ICOLUM,L)*T
  550 CONTINUE
  600 DO 710 I=1,N
  610 L = N+1-I
  620 IF(INDEX(L,1)-INDEX(L,2)) 630,710,630
  630 JROW = INDEX (L,1)
  640 JCOLUM = INDEX (L,2)
  650 DO  705 K= 1,N
  660 SWAP = A(K,JROW)
  670 A(K,JROW) = A(K,JCOLUM)
  700 A(K,JCOLUM) = SWAP
  705 CONTINUE
  710 CONTINUE
  740 RETURN
  741 A(1,1) = 1./A(1,1)
      RETURN
C
  999 CONTINUE
      PRINT 9000
      PRINT 9999, NMAX, N
      STOP
C
 9000 FORMAT(//,5X,10('*'),' STOP IN INV ',10('*'),/)
 9999 FORMAT(5X,'DIMENSION PARAMETER NMAX = ',I3,' TOO SMALL',/,
     &       5X,I3,' REQUIRED')
       END


      SUBROUTINE MYPOT(R, thetar, phir,vv)

      implicit double precision (a-h,o-z)
      logical lfirst

      parameter (lm=2)
      dimension tetra((lm+1)*(lm+1)), c((lm+1)*(lm+1))

      data lfirst/.true./

      if (lfirst) then

         write(2,*)
         write(2,*)' in mypot, test potential given by initial coeffts'

         do i=1,(lm+1)*(lm+1)
            c(i)= 0.d0
         end do

         c(1)=1.33d0
         c(2)=2.34d-1
         c(3)=2.78d-1
         c(4)=3.98d-1
         c(5)=1.22d-2
         c(6)=0.57d-2
         c(7)=0.75d-2
         c(8)=2.45d-2
         c(9)=9.10d-2
cN      c(1)=1.33d0
cN      c(2)=2.34d0
cN      c(3)=2.78d0
cN      c(4)=3.98d0
cN      c(5)=1.22d0
cN      c(6)=0.57d0
cN      c(7)=0.75d0
cN      c(8)=2.45d0
cN      c(9)=9.10d0
ccc      c(10)=6.777d0
ccc      c(11)=8.666d0
ccc      c(12)=9.875d0
ccc      c(13)=10.229d0
ccc      c(14)=11.777d0
ccc      c(15)=8.444d0
ccc      c(16)=8.999d0


         write (2,*) (c(i),i=1,(lm+1)*(lm+1))

         lfirst = .false.
      endif

      call tesser (tetra,lm,thetar,phir)

c     racah normalize

      pi4=4.d0*acos(-1.d0)
      ii=0
      do l=0,lm
        fac=sqrt(pi4/dble((2*l+1)))
        do m=-l,l
           ii=ii+1
           tetra(ii)=tetra(ii)*fac
        end do
      end do

ccc      write(2,*) 'tesseral in my pot', 
ccc     &                 (tetra(i),i=1,(lm+1)*(lm+1))

      vv=0.d0
      do i=1,(lm+1)*(lm+1)
cN         vv=vv+c(i)*tetra(i)*r**(i-1)
         vv=vv+c(i)*tetra(i)
      end do
   
c      write(2,*) 'vv in my pot at r,theta,phi = ',R,thetar,phir
c      write(2,*) vv

      return
      end





