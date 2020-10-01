c**************************************************
        subroutine prepr1d(n1,trigs,ifax,nfac)

        complex*16 trigs(0:*)
        integer ifax(0:*),nfac
	
        call ffttrig(n1,trigs)
	
        call myfactor(ifax,nfac,n1)
	
        return
        end
c---------------------------------------------------------------------------
	subroutine initak (ak,n,delta,iorder)

c	********************************************************************
c	*  The subroutine initak initializes an array ak , which can be    *
c	*  used for multiplication in the frequency domain of an FFT.      *
c	*  The array will contain the values  ((0,1)*k) ** iorder , where  *
c	*  the real variable k is the variable in the frequency domain.    *
c	*                                                                  *
c	*  ak	  is a complex one dimensional array of length n.          *
c	*  n	  is the length of the ak-array. n is a power of 2.        *
c	*  delta  is the distance between grid points in the time domain.  *
c	*  iorder is the power of (0,1)*k (equivalent to the order of the  *
c	*	  derivative when the FFT is used for differentiating).    *
c	*                                                                  *
c	********************************************************************

	IMPLICIT REAL*8(A-H,O-y)
	complex*16 ak(0:n-1)
	complex*16 c,ctemp

	pi=4*atan(1.0)
	anorm=2*pi/(n*delta)

c	***   Compute c=(0,1)**iorder     ***
	iord4=mod(iorder,4)
	if      (iord4.eq.0) then
		c=1.
	else if (iord4.eq.1) then
		c=(0.,1.)
	else if (iord4.eq.2) then
		c=-1.
	else if (iord4.eq.3) then
		c=(0.,-1.)
	end if
c
c	*** Compute sign=(-1)**iorder     ***
	iord2=mod(iorder,2)
	if (iord2.eq.0) then
		sign=1
	else
		sign=-1
	end if
	 
c	***   Initialize array            ***
	ak(0)=0
	do 10 j=1,n/2
		temp=(j*anorm)**iorder
		ctemp=c*temp
		ak(n-j)=ctemp
		ak(j)=sign*ctemp
10	continue
	 
	return
	end
c---------------------------------------------------------------------------
        subroutine myfactor(ifax,nfac,n)
        dimension ifax(*)
        nfac=0
        nn=n
100     continue
        if(mod(nn,6).eq.0) then
        nn=nn/6
        nfac=nfac+1
        ifax(nfac)=6
        go to 100
        else if(mod(nn,5).eq.0) then
        nn=nn/5
        nfac=nfac+1
        ifax(nfac)=5
        go to 100
        else if(mod(nn,4).eq.0) then
        nn=nn/4
        nfac=nfac+1
        ifax(nfac)=4
        go to 100
        else if(mod(nn,3).eq.0) then
        nn=nn/3
        nfac=nfac+1
        ifax(nfac)=3
        go to 100
        else if(mod(nn,2).eq.0) then
        nn=nn/2
        nfac=nfac+1
        ifax(nfac)=2
        go to 100
        else 
        if(nn.eq.1) return
        print *,' factorization failed'
        stop 999
        endif
        end
c------------------------------------------------------------------------------
        subroutine ffttrig(nfft,trig)
        complex*16 trig(*),cc

          dinc=2.*3.14159265/nfft
          do 1 n=1,nfft
           cc=dinc*(n-1)
           trig(n)=cdexp((0.,1.)*cc)
1       continue
          return
        end

C---------------------------------------------------------------------
C       MIXED RADIX FFT ROUTINE BASED ON TEMPERTON ,JOUR.COMP.
C       PHYS. VOL 52 NUMBER 1 OCT 1983 PAGE 1.
C
C       VARIABLES:
C                 N -     TRANSFORM LENGTH
C                 A -     COMPLEX INPUT ARRAY OF LENGTH N
C                 C -     ADDITIONAL COMPLEX WORK ARRAY OF LENGTH N
C                 TRIGS - COMPLEX ARRAY OF DFT EXPONENTIAL FACTORS
C                         OF LENGTH N (PREVIOUSLY GENERATED)
c                 IFAX  - INTEGER ARRAY CONTAINING FACTORIZATIONS OF N
C                         ACCORDING TO N=n1*n2....*nk
C                         WHERE NK=2,OR 3,4,5,6
C                 NFAC  - NUMBER OF FACTORS IN N (BOTH NFAC AND IFAX
C                         PREVIOUSLY CALCULATED IN FACTOR)
C                 ISKIP - STRIDE OF FFT (e.g IF=2  SKIP EVERY SECOND SAMPLE)
C                 ISIGN - FFT SIGN
C---------------------------------------------------------------------
         subroutine tfft(a,c,n,trigs,ifax,nfac,iskip,isign)
         implicit real*8(a-h,o-z)
         dimension a(1),c(1)
         iskip2=iskip+iskip
         call fftb(a(1),a(2),c(1),c(n+1),n,trigs,ifax,nfac,iskip2,isign)
         return
        end
        subroutine fftb(a,b,c,d,n,trigs,ifax,nfac,iskip,isign)
        implicit real*8(a-h,o-z)
        complex*16 trigs(1)
        integer ifax(1)
        dimension a(1),c(1),b(1),d(1)
        la=1
         ii=1
        do 1 i=1,nfac
         ifac=ifax(i)
         if(ii.gt.0) call pass(a,b,c,d,
     .     trigs,ifac,la,n,iskip,1,isign)
         if(ii.lt.0) call pass(c,d,a,b,
     .     trigs,ifac,la,n,1,iskip,isign)
         ii=-ii
         la=la*ifac
1       continue        
        if(ii.lt.0) then
        n2=n*iskip
         jj=0
         do 2 j=1,n2,iskip
          jj=jj+1
          a(j)=c(jj)
          b(j)=d(jj)
2       continue 
        end if
        return
        end
        subroutine pass(a,b,c,d,trigs,ifac,la,n,iskip,
     .     iskip1,isign)
        implicit real*8(a-h,o-z)
        dimension a(1),c(1),trigs(1),b(1),d(1)
        save sin60,sin72,sin36,sq54,sin3672
        data iflag/0/
        if(iflag.eq.0) then
         iflag=1
        sin60=sin(3.14159265/3.)
        sin72=sin(3.14159265*72./180.)
        sin36=sin(3.14159265*36./180.)
        sq54=0.25*sqrt(5.)
        sin3672=sin36/sin72
        end if
        m=n/ifac
        mskip=m*iskip
        ia=0
        ib=mskip
        ja=0
        laskip=la*iskip1
        jb=laskip
        i=1
        j=1
        jump=(ifac-1)*laskip
        go to(100,200,300,400,500,600),ifac
100     print *,'radix one encountered'
        stop 999
200     continue
c----------------- factor two ----------------------
         do 1 l=1,la
          aiQR=a(i)
          aiQI=b(i)
          iib=i+ib
          aibQR=a(iib)
          aibQI=b(iib)
          c(j)=aiQR+aibQR
          d(j)=aiQI+aibQI
          jjb=j+jb
          c(jjb)=aiQR-aibQR
          d(jjb)=aiQI-aibQI
         i=i+iskip
         j=j+iskip1
1       continue
         j=j+jump
        do 2 k=la,m-la,la
         kk=2*k+1
         tQR=trigs(kk)
         tQI=trigs(kk+1)
         if(isign.eq.-1) tQI=-tQI
         do 3 l=1,la
          aiQR=a(i)
          aiQI=b(i)
          iib=i+ib
          aibQR=a(iib)
          aibQI=b(iib)
          c(j)=aiQR+aibQR
          d(j)=aiQI+aibQI
          f1=aiQR-aibQR
          f2=aiQI-aibQI
          jjb=j+jb
          c(jjb)=tQR*f1-tQI*f2
          d(jjb)=tQI*f1+tQR*f2
         i=i+iskip
         j=j+iskip1
3       continue        
         j=j+jump
2       continue        
        return
300     continue
c       ------ factor three ---------
        ic=ib+mskip
        jc=jb+laskip
         if(isign.gt.0) then
          do 4 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           aibQR=a(iib)
           aibQI=b(iib)
           iic=i+ic
           aicQR=a(iic)
           aicQI=b(iic)
           t1QR=aibQR+aicQR
           t1QI=aibQI+aicQI
           t2QR=aiQR-0.5*t1QR
           t2QI=aiQI-0.5*t1QI
           t3QR=sin60*(aibQR-aicQR)
           t3QI=sin60*(aibQI-aicQI)
           c(j)=aiQR+t1QR
           d(j)=aiQI+t1QI
           jjb=j+jb
           c(jjb)=t2QR-t3QI
           d(jjb)=t2QI+t3QR
           jjc=j+jc
           c(jjc)=t2QR+t3QI
           d(jjc)=t2QI-t3QR
           i=i+iskip
           j=j+iskip1
4       continue 
          j=j+jump
         do 5 k=la,m-la,la
          k2=k+k
          kk=k2+1
          tr1QR=trigs(kk)
          tr1QI=trigs(kk+1)
          kk=kk+k2
          tr2QR=trigs(kk)
          tr2QI=trigs(kk+1)
          do 6 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           aibQR=a(iib)
           aibQI=b(iib)
           iic=i+ic
           aicQR=a(iic)
           aicQI=b(iic)
           t1QR=aibQR+aicQR
           t1QI=aibQI+aicQI
           t2QR=aiQR-0.5*t1QR
           t2QI=aiQI-0.5*t1QI
           t3QR=sin60*(aibQR-aicQR)
           t3QI=sin60*(aibQI-aicQI)
           c(j)=aiQR+t1QR
           d(j)=aiQI+t1QI
           x1QR=t2QR-t3QI
           x1QI=t2QI+t3QR
           x2QR=t2QR+t3QI
           x2QI=t2QI-t3QR
           jjb=j+jb
           c(jjb)=tr1QR*x1QR-tr1QI*x1QI
           d(jjb)=tr1QI*x1QR+tr1QR*x1QI
           jjc=j+jc
           c(jjc)=tr2QR*x2QR-tr2QI*x2QI
           d(jjc)=tr2QI*x2QR+tr2QR*x2QI
           i=i+iskip
           j=j+iskip1
6       continue
          j=j+jump
5       continue
        else
          do 7 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           aibQR=a(iib)
           aibQI=b(iib)
           iic=i+ic
           aicQR=a(iic)
           aicQI=b(iic)
           t1QR=aibQR+aicQR
           t1QI=aibQI+aicQI
           t2QR=aiQR-0.5*t1QR
           t2QI=aiQI-0.5*t1QI
           t3QR=sin60*(aibQR-aicQR)
           t3QI=sin60*(aibQI-aicQI)
           c(j)=aiQR+t1QR
           d(j)=aiQI+t1QI
           jjb=j+jb
           c(jjb)=t2QR+t3QI
           d(jjb)=t2QI-t3QR
           jjc=j+jc
           c(jjc)=t2QR-t3QI
           d(jjc)=t2QI+t3QR
           i=i+iskip
           j=j+iskip1
7        continue
          j=j+jump
         do 8 k=la,m-la,la
          k2=k+k
          kk=k2+1
          tr1QR=trigs(kk)
          tr1QI=-trigs(kk+1)
          kk=kk+k2
          tr2QR=trigs(kk)
          tr2QI=-trigs(kk+1)
          do 9 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           aibQR=a(iib)
           aibQI=b(iib)
           iic=i+ic
           aicQR=a(iic)
           aicQI=b(iic)
           t1QR=aibQR+aicQR
           t1QI=aibQI+aicQI
           t2QR=aiQR-0.5*t1QR
           t2QI=aiQI-0.5*t1QI
           t3QR=sin60*(aibQR-aicQR)
           t3QI=sin60*(aibQI-aicQI)
           c(j)=aiQR+t1QR
           d(j)=aiQI+t1QI
           x1QR=t2QR+t3QI
           x1QI=t2QI-t3QR
           x2QR=t2QR-t3QI
           x2QI=t2QI+t3QR
           jjb=j+jb
           c(jjb)=tr1QR*x1QR-tr1QI*x1QI
           d(jjb)=tr1QI*x1QR+tr1QR*x1QI
           jjc=j+jc
           c(jjc)=tr2QR*x2QR-tr2QI*x2QI
           d(jjc)=tr2QI*x2QR+tr2QR*x2QI
           i=i+iskip
           j=j+iskip1
9       continue
          j=j+jump
8       continue
        endif
        return
400     continue
c ------------  factor 4 ---------------------------
        ic=ib+mskip
        id=ic+mskip
        jc=jb+laskip
        jd=jc+laskip
         if(isign.gt.0) then
          do 10 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           aibQR=a(iib)
           aibQI=b(iib)
           iic=i+ic
           aicQR=a(iic)
           aicQI=b(iic)
           iid=i+id
           aidQR=a(iid)
           aidQI=b(iid)
           t1QR=aiQR+aicQR
           t1QI=aiQI+aicQI
           t2QR=aibQR+aidQR
           t2QI=aibQI+aidQI
           t3QR=aiQR-aicQR
           t3QI=aiQI-aicQI
           t4QR=aibQR-aidQR
           t4QI=aibQI-aidQI
           c(j)=t1QR+t2QR
           d(j)=t1QI+t2QI
           jjb=j+jb
           c(jjb)=t3QR-t4QI
           d(jjb)=t3QI+t4QR
           jjc=j+jc
           c(jjc)=t1QR-t2QR
           d(jjc)=t1QI-t2QI
           jjd=j+jd
           c(jjd)=t3QR+t4QI
           d(jjd)=t3QI-t4QR
           i=i+iskip
           j=j+iskip1
10      continue
          j=j+jump
         do 11 k=la,m-la,la
          k2=k+k
          kk=k2+1
          tr1QR=trigs(kk)
          tr1QI=trigs(kk+1)
          kk=kk+k2
          tr2QR=trigs(kk)
          tr2QI=trigs(kk+1)
          kk=kk+k2
          tr3QR=trigs(kk)
          tr3QI=trigs(kk+1)
          do 12 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           aibQR=a(iib)
           aibQI=b(iib)
           iic=i+ic
           aicQR=a(iic)
           aicQI=b(iic)
           iid=i+id
           aidQR=a(iid)
           aidQI=b(iid)
           t1QR=aiQR+aicQR
           t1QI=aiQI+aicQI
           t2QR=aibQR+aidQR
           t2QI=aibQI+aidQI
           t3QR=aiQR-aicQR
           t3QI=aiQI-aicQI
           t4QR=aibQR-aidQR
           t4QI=aibQI-aidQI
           c(j)=t1QR+t2QR
           d(j)=t1QI+t2QI
           x1QR=t3QR-t4QI
           x1QI=t3QI+t4QR
           x2QR=t1QR-t2QR
           x2QI=t1QI-t2QI
           x3QR=t3QR+t4QI
           x3QI=t3QI-t4QR
           jjb=j+jb
           c(jjb)=tr1QR*x1QR-tr1QI*x1QI
           d(jjb)=tr1QI*x1QR+tr1QR*x1QI
           jjc=j+jc
           c(jjc)=tr2QR*x2QR-tr2QI*x2QI
           d(jjc)=tr2QI*x2QR+tr2QR*x2QI
           jjd=j+jd
           c(jjd)=tr3QR*x3QR-tr3QI*x3QI
           d(jjd)=tr3QI*x3QR+tr3QR*x3QI
           i=i+iskip
           j=j+iskip1
12      continue
          j=j+jump
11      continue
        else
          do 13 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           aibQR=a(iib)
           aibQI=b(iib)
           iic=i+ic
           aicQR=a(iic)
           aicQI=b(iic)
           iid=i+id
           aidQR=a(iid)
           aidQI=b(iid)
           t1QR=aiQR+aicQR
           t1QI=aiQI+aicQI
           t2QR=aibQR+aidQR
           t2QI=aibQI+aidQI
           t3QR=aiQR-aicQR
           t3QI=aiQI-aicQI
           t4QR=aibQR-aidQR
           t4QI=aibQI-aidQI
           c(j)=t1QR+t2QR
           d(j)=t1QI+t2QI
           jjb=j+jb
           c(jjb)=t3QR+t4QI
           d(jjb)=t3QI-t4QR
           jjc=j+jc
           c(jjc)=t1QR-t2QR
           d(jjc)=t1QI-t2QI
           jjd=j+jd
           c(jjd)=t3QR-t4QI
           d(jjd)=t3QI+t4QR
           i=i+iskip
           j=j+iskip1
13      continue
          j=j+jump
         do 14 k=la,m-la,la
          k2=k+k
          kk=k2+1
          tr1QR=trigs(kk)
          tr1QI=-trigs(kk+1)
          kk=kk+k2
          tr2QR=trigs(kk)
          tr2QI=-trigs(kk+1)
          kk=kk+k2
          tr3QR=trigs(kk)
          tr3QI=-trigs(kk+1)
          do 15 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           aibQR=a(iib)
           aibQI=b(iib)
           iic=i+ic
           aicQR=a(iic)
           aicQI=b(iic)
           iid=i+id
           aidQR=a(iid)
           aidQI=b(iid)
           t1QR=aiQR+aicQR
           t1QI=aiQI+aicQI
           t2QR=aibQR+aidQR
           t2QI=aibQI+aidQI
           t3QR=aiQR-aicQR
           t3QI=aiQI-aicQI
           t4QR=aibQR-aidQR
           t4QI=aibQI-aidQI
           c(j)=t1QR+t2QR
           d(j)=t1QI+t2QI
           x1QR=t3QR+t4QI
           x1QI=t3QI-t4QR
           x2QR=t1QR-t2QR
           x2QI=t1QI-t2QI
           x3QR=t3QR-t4QI
           x3QI=t3QI+t4QR
           jjb=j+jb
           c(jjb)=tr1QR*x1QR-tr1QI*x1QI
           d(jjb)=tr1QI*x1QR+tr1QR*x1QI
           jjc=j+jc
           c(jjc)=tr2QR*x2QR-tr2QI*x2QI
           d(jjc)=tr2QI*x2QR+tr2QR*x2QI
           jjd=j+jd
           c(jjd)=tr3QR*x3QR-tr3QI*x3QI
           d(jjd)=tr3QI*x3QR+tr3QR*x3QI
           i=i+iskip
           j=j+iskip1
15      continue
          j=j+jump
14      continue
        endif
        return
500     continue
c ------------  factor 5 ---------------------------
        ic=ib+mskip
        id=ic+mskip
        ie=id+mskip
        jc=jb+laskip
        jd=jc+laskip
        je=jd+laskip
         if(isign.gt.0) then
          do 16 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           aibQR=a(iib)
           aibQI=b(iib)
           iic=i+ic
           aicQR=a(iic)
           aicQI=b(iic)
           iid=i+id
           aidQR=a(iid)
           aidQI=b(iid)
           iie=i+ie
           aieQR=a(iie)
           aieQI=b(iie)
           t1QR=aibQR+aieQR
           t1QI=aibQI+aieQI
           t2QR=aicQR+aidQR
           t2QI=aicQI+aidQI
           t3QR=sin72*(aibQR-aieQR)
           t3QI=sin72*(aibQI-aieQI)
           t4QR=sin72*(aicQR-aidQR)
           t4QI=sin72*(aicQI-aidQI)
           t5QR=t1QR+t2QR
           t5QI=t1QI+t2QI
           t6QR=sq54*(t1QR-t2QR)
           t6QI=sq54*(t1QI-t2QI)
           t7QR=aiQR-t5QR*0.25
           t7QI=aiQI-t5QI*0.25
           t8QR=t7QR+t6QR
           t8QI=t7QI+t6QI
           t9QR=t7QR-t6QR
           t9QI=t7QI-t6QI
           t10QR=t3QR+sin3672*t4QR
           t10QI=t3QI+sin3672*t4QI
           t11QR=sin3672*t3QR-t4QR
           t11QI=sin3672*t3QI-t4QI
           c(j)=aiQR+t5QR
           d(j)=aiQI+t5QI
           jjb=j+jb
           c(jjb)=t8QR-t10QI
           d(jjb)=t8QI+t10QR
           jjc=j+jc
           c(jjc)=t9QR-t11QI
           d(jjc)=t9QI+t11QR
           jjd=j+jd
           c(jjd)=t9QR+t11QI
           d(jjd)=t9QI-t11QR
           jje=j+je
           c(jje)=t8QR+t10QI
           d(jje)=t8QI-t10QR
           i=i+iskip
           j=j+iskip1
16        continue
          j=j+jump
         do 17 k=la,m-la,la
          k2=k+k
          kk=k2+1
          tr1QR=trigs(kk)
          tr1QI=trigs(kk+1)
          kk=kk+k2
          tr2QR=trigs(kk)
          tr2QI=trigs(kk+1)
          kk=kk+k2
          tr3QR=trigs(kk)
          tr3QI=trigs(kk+1)
          kk=kk+k2
          tr4QR=trigs(kk)
          tr4QI=trigs(kk+1)
          do 18 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           aibQR=a(iib)
           aibQI=b(iib)
           iic=i+ic
           aicQR=a(iic)
           aicQI=b(iic)
           iid=i+id
           aidQR=a(iid)
           aidQI=b(iid)
           iie=i+ie
           aieQR=a(iie)
           aieQI=b(iie)
           t1QR=aibQR+aieQR
           t1QI=aibQI+aieQI
           t2QR=aicQR+aidQR
           t2QI=aicQI+aidQI
           t3QR=sin72*(aibQR-aieQR)
           t3QI=sin72*(aibQI-aieQI)
           t4QR=sin72*(aicQR-aidQR)
           t4QI=sin72*(aicQI-aidQI)
           t5QR=t1QR+t2QR
           t5QI=t1QI+t2QI
           t6QR=sq54*(t1QR-t2QR)
           t6QI=sq54*(t1QI-t2QI)
           t7QR=aiQR-t5QR*0.25
           t7QI=aiQI-t5QI*0.25
           t8QR=t7QR+t6QR
           t8QI=t7QI+t6QI
           t9QR=t7QR-t6QR
           t9QI=t7QI-t6QI
           t10QR=t3QR+sin3672*t4QR
           t10QI=t3QI+sin3672*t4QI
           t11QR=sin3672*t3QR-t4QR
           t11QI=sin3672*t3QI-t4QI
           c(j)=aiQR+t5QR
           d(j)=aiQI+t5QI
           x1QR=t8QR-t10QI
           x1QI=t8QI+t10QR
           x2QR=t9QR-t11QI
           x2QI=t9QI+t11QR
           x3QR=t9QR+t11QI
           x3QI=t9QI-t11QR
           x4QR=t8QR+t10QI
           x4QI=t8QI-t10QR
           jjb=j+jb
           jjc=j+jc
           jjd=j+jd
           jje=j+je
           c(jjb)=tr1QR*x1QR-tr1QI*x1QI
           d(jjb)=tr1QI*x1QR+tr1QR*x1QI
           c(jjc)=tr2QR*x2QR-tr2QI*x2QI
           d(jjc)=tr2QI*x2QR+tr2QR*x2QI
           c(jjd)=tr3QR*x3QR-tr3QI*x3QI
           d(jjd)=tr3QI*x3QR+tr3QR*x3QI
           c(jje)=tr4QR*x4QR-tr4QI*x4QI
           d(jje)=tr4QI*x4QR+tr4QR*x4QI
           i=i+iskip
           j=j+iskip1
18      continue
          j=j+jump
17      continue
        else
          do 20 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           iic=i+ic
           iid=i+id
           iie=i+ie
           aibQR=a(iib)
           aibQI=b(iib)
           aicQR=a(iic)
           aicQI=b(iic)
           aidQR=a(iid)
           aidQI=b(iid)
           aieQR=a(iie)
           aieQI=b(iie)
           t1QR=aibQR+aieQR
           t1QI=aibQI+aieQI
           t2QR=aicQR+aidQR
           t2QI=aicQI+aidQI
           t3QR=sin72*(aibQR-aieQR)
           t3QI=sin72*(aibQI-aieQI)
           t4QR=sin72*(aicQR-aidQR)
           t4QI=sin72*(aicQI-aidQI)
           t5QR=t1QR+t2QR
           t5QI=t1QI+t2QI
           t6QR=sq54*(t1QR-t2QR)
           t6QI=sq54*(t1QI-t2QI)
           t7QR=aiQR-t5QR*0.25
           t7QI=aiQI-t5QI*0.25
           t8QR=t7QR+t6QR
           t8QI=t7QI+t6QI
           t9QR=t7QR-t6QR
           t9QI=t7QI-t6QI
           t10QR=t3QR+sin3672*t4QR
           t10QI=t3QI+sin3672*t4QI
           t11QR=sin3672*t3QR-t4QR
           t11QI=sin3672*t3QI-t4QI
           c(j)=aiQR+t5QR
           d(j)=aiQI+t5QI
           jjb=j+jb
           jjc=j+jc
           jjd=j+jd
           jje=j+je
           c(jjb)=t8QR+t10QI
           d(jjb)=t8QI-t10QR
           c(jjc)=t9QR+t11QI
           d(jjc)=t9QI-t11QR
           c(jjd)=t9QR-t11QI
           d(jjd)=t9QI+t11QR
           c(jje)=t8QR-t10QI
           d(jje)=t8QI+t10QR
           i=i+iskip
           j=j+iskip1
20      continue
          j=j+jump
         do 21 k=la,m-la,la
          k2=k+k
          kk=k2+1
          tr1QR=trigs(kk)
          tr1QI=-trigs(kk+1)
          kk=kk+k2
          tr2QR=trigs(kk)
          tr2QI=-trigs(kk+1)
          kk=kk+k2
          tr3QR=trigs(kk)
          tr3QI=-trigs(kk+1)
          kk=kk+k2
          tr4QR=trigs(kk)
          tr4QI=-trigs(kk+1)
          do 22 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           iic=i+ic
           iid=i+id
           iie=i+ie
           aibQR=a(iib)
           aibQI=b(iib)
           aicQR=a(iic)
           aicQI=b(iic)
           aidQR=a(iid)
           aidQI=b(iid)
           aieQR=a(iie)
           aieQI=b(iie)
           t1QR=aibQR+aieQR
           t1QI=aibQI+aieQI
           t2QR=aicQR+aidQR
           t2QI=aicQI+aidQI
           t3QR=sin72*(aibQR-aieQR)
           t3QI=sin72*(aibQI-aieQI)
           t4QR=sin72*(aicQR-aidQR)
           t4QI=sin72*(aicQI-aidQI)
           t5QR=t1QR+t2QR
           t5QI=t1QI+t2QI
           t6QR=sq54*(t1QR-t2QR)
           t6QI=sq54*(t1QI-t2QI)
           t7QR=aiQR-t5QR*0.25
           t7QI=aiQI-t5QI*0.25
           t8QR=t7QR+t6QR
           t8QI=t7QI+t6QI
           t9QR=t7QR-t6QR
           t9QI=t7QI-t6QI
           t10QR=t3QR+sin3672*t4QR
           t10QI=t3QI+sin3672*t4QI
           t11QR=sin3672*t3QR-t4QR
           t11QI=sin3672*t3QI-t4QI
           c(j)=aiQR+t5QR
           d(j)=aiQI+t5QI
           x1QR=t8QR+t10QI
           x1QI=t8QI-t10QR
           x2QR=t9QR+t11QI
           x2QI=t9QI-t11QR
           x3QR=t9QR-t11QI
           x3QI=t9QI+t11QR
           x4QR=t8QR-t10QI
           x4QI=t8QI+t10QR
           jjb=j+jb
           jjc=j+jc
           jjd=j+jd
           jje=j+je
           c(jjb)=tr1QR*x1QR-tr1QI*x1QI
           d(jjb)=tr1QI*x1QR+tr1QR*x1QI
           c(jjc)=tr2QR*x2QR-tr2QI*x2QI
           d(jjc)=tr2QI*x2QR+tr2QR*x2QI
           c(jjd)=tr3QR*x3QR-tr3QI*x3QI
           d(jjd)=tr3QI*x3QR+tr3QR*x3QI
           c(jje)=tr4QR*x4QR-tr4QI*x4QI
           d(jje)=tr4QI*x4QR+tr4QR*x4QI
           i=i+iskip
           j=j+iskip1
22      continue        
          j=j+jump
21      continue
        endif
        return
600     continue
c ------------  factor 6 ---------------------------
        ic=ib+mskip
        id=ic+mskip
        ie=id+mskip
        ig=ie+mskip
        jc=jb+laskip
        jd=jc+laskip
        je=jd+laskip
        jg=je+laskip
         if(isign.gt.0) then
          do 23 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           iic=i+ic
           iid=i+id
           iie=i+ie
           iig=i+ig
           aibQR=a(iib)
           aibQI=b(iib)
           aicQR=a(iic)
           aicQI=b(iic)
           aidQR=a(iid)
           aidQI=b(iid)
           aieQR=a(iie)
           aieQI=b(iie)
           aigQR=a(iig)
           aigQI=b(iig)
           t1QR=aicQR+aieQR
           t1QI=aicQI+aieQI
           t2QR=aiQR-0.5*t1QR
           t2QI=aiQI-0.5*t1QI
           t3QR=sin60*(aicQR-aieQR)
           t3QI=sin60*(aicQI-aieQI)
           y0QR=aiQR+t1QR
           y0QI=aiQI+t1QI
           y4QR=t2QR-t3QI
           y4QI=t2QI+t3QR
           y2QR=t2QR+t3QI
           y2QI=t2QI-t3QR
           t1QR=aigQR+aibQR
           t1QI=aigQI+aibQI
           t2QR=aidQR-0.5*t1QR
           t2QI=aidQI-0.5*t1QI
           t3QR=sin60*(aigQR-aibQR)
           t3QI=sin60*(aigQI-aibQI)
           y3QR=aidQR+t1QR
           y3QI=aidQI+t1QI
           y1QR=t2QR-t3QI
           y1QI=t2QI+t3QR
           y5QR=t2QR+t3QI
           y5QI=t2QI-t3QR
           x0QR=y0QR+y3QR
           x0QI=y0QI+y3QI
           x4QR=y4QR+y1QR
           x4QI=y4QI+y1QI
           x2QR=y2QR+y5QR
           x2QI=y2QI+y5QI
           x3QR=y0QR-y3QR
           x3QI=y0QI-y3QI
           x1QR=y4QR-y1QR
           x1QI=y4QI-y1QI
           x5QR=y2QR-y5QR
           x5QI=y2QI-y5QI
           c(j)=x0QR
           d(j)=x0QI
           jjb=j+jb
           jjc=j+jc
           jjd=j+jd
           jje=j+je
           jjg=j+jg
           c(jjb)=x1QR
           d(jjb)=x1QI
           c(jjc)=x2QR
           d(jjc)=x2QI
           c(jjd)=x3QR
           d(jjd)=x3QI
           c(jje)=x4QR
           d(jje)=x4QI
           c(jjg)=x5QR
           d(jjg)=x5QI
           i=i+iskip
           j=j+iskip1
23      continue
          j=j+jump
         do 24 k=la,m-la,la
          k2=k+k
          kk=k2+1
          tr1QR=trigs(kk)
          tr1QI=trigs(kk+1)
          kk=kk+k2
          tr2QR=trigs(kk)
          tr2QI=trigs(kk+1)
          kk=kk+k2
          tr3QR=trigs(kk)
          tr3QI=trigs(kk+1)
          kk=kk+k2
          tr4QR=trigs(kk)
          tr4QI=trigs(kk+1)
          kk=kk+k2
          tr5QR=trigs(kk)
          tr5QI=trigs(kk+1)
          do 25 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           iic=i+ic
           iid=i+id
           iie=i+ie
           iig=i+ig
           aibQR=a(iib)
           aibQI=b(iib)
           aicQR=a(iic)
           aicQI=b(iic)
           aidQR=a(iid)
           aidQI=b(iid)
           aieQR=a(iie)
           aieQI=b(iie)
           aigQR=a(iig)
           aigQI=b(iig)
           t1QR=aicQR+aieQR
           t1QI=aicQI+aieQI
           t2QR=aiQR-0.5*t1QR
           t2QI=aiQI-0.5*t1QI
           t3QR=sin60*(aicQR-aieQR)
           t3QI=sin60*(aicQI-aieQI)
           y0QR=aiQR+t1QR
           y0QI=aiQI+t1QI
           y4QR=t2QR-t3QI
           y4QI=t2QI+t3QR
           y2QR=t2QR+t3QI
           y2QI=t2QI-t3QR
           t1QR=aigQR+aibQR
           t1QI=aigQI+aibQI
           t2QR=aidQR-0.5*t1QR
           t2QI=aidQI-0.5*t1QI
           t3QR=sin60*(aigQR-aibQR)
           t3QI=sin60*(aigQI-aibQI)
           y3QR=aidQR+t1QR
           y3QI=aidQI+t1QI
           y1QR=t2QR-t3QI
           y1QI=t2QI+t3QR
           y5QR=t2QR+t3QI
           y5QI=t2QI-t3QR
           x0QR=y0QR+y3QR
           x0QI=y0QI+y3QI
           x4QR=y4QR+y1QR
           x4QI=y4QI+y1QI
           x2QR=y2QR+y5QR
           x2QI=y2QI+y5QI
           x3QR=y0QR-y3QR
           x3QI=y0QI-y3QI
           x1QR=y4QR-y1QR
           x1QI=y4QI-y1QI
           x5QR=y2QR-y5QR
           x5QI=y2QI-y5QI
           c(j)=x0QR
           d(j)=x0QI
           jjb=j+jb
           jjc=j+jc
           jjd=j+jd
           jje=j+je
           jjg=j+jg
           c(jjb)=tr1QR*x1QR-tr1QI*x1QI
           d(jjb)=tr1QI*x1QR+tr1QR*x1QI
           c(jjc)=tr2QR*x2QR-tr2QI*x2QI
           d(jjc)=tr2QI*x2QR+tr2QR*x2QI
           c(jjd)=tr3QR*x3QR-tr3QI*x3QI
           d(jjd)=tr3QI*x3QR+tr3QR*x3QI
           c(jje)=tr4QR*x4QR-tr4QI*x4QI
           d(jje)=tr4QI*x4QR+tr4QR*x4QI
           c(jjg)=tr5QR*x5QR-tr5QI*x5QI
           d(jjg)=tr5QI*x5QR+tr5QR*x5QI
           i=i+iskip
           j=j+iskip1
25      continue
          j=j+jump
24      continue
        else
          do 26 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           iic=i+ic
           iid=i+id
           iie=i+ie
           iig=i+ig
           aibQR=a(iib)
           aibQI=b(iib)
           aicQR=a(iic)
           aicQI=b(iic)
           aidQR=a(iid)
           aidQI=b(iid)
           aieQR=a(iie)
           aieQI=b(iie)
           aigQR=a(iig)
           aigQI=b(iig)
           t1QR=aicQR+aieQR
           t1QI=aicQI+aieQI
           t2QR=aiQR-0.5*t1QR
           t2QI=aiQI-0.5*t1QI
           t3QR=sin60*(aicQR-aieQR)
           t3QI=sin60*(aicQI-aieQI)
           y0QR=aiQR+t1QR
           y0QI=aiQI+t1QI
           y4QR=t2QR+t3QI
           y4QI=t2QI-t3QR
           y2QR=t2QR-t3QI
           y2QI=t2QI+t3QR
           t1QR=aigQR+aibQR
           t1QI=aigQI+aibQI
           t2QR=aidQR-0.5*t1QR
           t2QI=aidQI-0.5*t1QI
           t3QR=sin60*(aigQR-aibQR)
           t3QI=sin60*(aigQI-aibQI)
           y3QR=aidQR+t1QR
           y3QI=aidQI+t1QI
           y1QR=t2QR+t3QI
           y1QI=t2QI-t3QR
           y5QR=t2QR-t3QI
           y5QI=t2QI+t3QR
           x0QR=y0QR+y3QR
           x0QI=y0QI+y3QI
           x4QR=y4QR+y1QR
           x4QI=y4QI+y1QI
           x2QR=y2QR+y5QR
           x2QI=y2QI+y5QI
           x3QR=y0QR-y3QR
           x3QI=y0QI-y3QI
           x1QR=y4QR-y1QR
           x1QI=y4QI-y1QI
           x5QR=y2QR-y5QR
           x5QI=y2QI-y5QI
           c(j)=x0QR
           d(j)=x0QI
           jjb=j+jb
           jjc=j+jc
           jjd=j+jd
           jje=j+je
           jjg=j+jg
           c(jjb)=x1QR
           d(jjb)=x1QI
           c(jjc)=x2QR
           d(jjc)=x2QI
           c(jjd)=x3QR
           d(jjd)=x3QI
           c(jje)=x4QR
           d(jje)=x4QI
           c(jjg)=x5QR
           d(jjg)=x5QI
           i=i+iskip
           j=j+iskip1
26      continue
          j=j+jump
         do 27 k=la,m-la,la
          k2=k+k
          kk=k2+1
          tr1QR=trigs(kk)
          tr1QI=-trigs(kk+1)
          kk=kk+k2
          tr2QR=trigs(kk)
          tr2QI=-trigs(kk+1)
          kk=kk+k2
          tr3QR=trigs(kk)
          tr3QI=-trigs(kk+1)
          kk=kk+k2
          tr4QR=trigs(kk)
          tr4QI=-trigs(kk+1)
          kk=kk+k2
          tr5QR=trigs(kk)
          tr5QI=-trigs(kk+1)
          do 28 l=1,la
           aiQR=a(i)
           aiQI=b(i)
           iib=i+ib
           iic=i+ic
           iid=i+id
           iie=i+ie
           iig=i+ig
           aibQR=a(iib)
           aibQI=b(iib)
           aicQR=a(iic)
           aicQI=b(iic)
           aidQR=a(iid)
           aidQI=b(iid)
           aieQR=a(iie)
           aieQI=b(iie)
           aigQR=a(iig)
           aigQI=b(iig)
           t1QR=aicQR+aieQR
           t1QI=aicQI+aieQI
           t2QR=aiQR-0.5*t1QR
           t2QI=aiQI-0.5*t1QI
           t3QR=sin60*(aicQR-aieQR)
           t3QI=sin60*(aicQI-aieQI)
           y0QR=aiQR+t1QR
           y0QI=aiQI+t1QI
           y4QR=t2QR+t3QI
           y4QI=t2QI-t3QR
           y2QR=t2QR-t3QI
           y2QI=t2QI+t3QR
           t1QR=aigQR+aibQR
           t1QI=aigQI+aibQI
           t2QR=aidQR-0.5*t1QR
           t2QI=aidQI-0.5*t1QI
           t3QR=sin60*(aigQR-aibQR)
           t3QI=sin60*(aigQI-aibQI)
           y3QR=aidQR+t1QR
           y3QI=aidQI+t1QI
           y1QR=t2QR+t3QI
           y1QI=t2QI-t3QR
           y5QR=t2QR-t3QI
           y5QI=t2QI+t3QR
           x0QR=y0QR+y3QR
           x0QI=y0QI+y3QI
           x4QR=y4QR+y1QR
           x4QI=y4QI+y1QI
           x2QR=y2QR+y5QR
           x2QI=y2QI+y5QI
           x3QR=y0QR-y3QR
           x3QI=y0QI-y3QI
           x1QR=y4QR-y1QR
           x1QI=y4QI-y1QI
           x5QR=y2QR-y5QR
           x5QI=y2QI-y5QI
           c(j)=x0QR
           d(j)=x0QI
           jjb=j+jb
           jjc=j+jc
           jjd=j+jd
           jje=j+je
           jjg=j+jg
           c(jjb)=tr1QR*x1QR-tr1QI*x1QI
           d(jjb)=tr1QI*x1QR+tr1QR*x1QI
           c(jjc)=tr2QR*x2QR-tr2QI*x2QI
           d(jjc)=tr2QI*x2QR+tr2QR*x2QI
           c(jjd)=tr3QR*x3QR-tr3QI*x3QI
           d(jjd)=tr3QI*x3QR+tr3QR*x3QI
           c(jje)=tr4QR*x4QR-tr4QI*x4QI
           d(jje)=tr4QI*x4QR+tr4QR*x4QI
           c(jjg)=tr5QR*x5QR-tr5QI*x5QI
           d(jjg)=tr5QI*x5QR+tr5QR*x5QI
           i=i+iskip
           j=j+iskip1
28      continue
          j=j+jump
27      continue
        endif
        return
        end
        
   
      FUNCTION FACTRL(N)
      implicit real*8(a-h,p-y)
      dimension a(33)
      data ntop,a(1)/0,1./
      if (n.lt.0) then
         pause 'negative factorial'
      else if (n.le.ntop) then
         factrl=a(n+1)
      else if (N.le.32) then
         do 11 J=ntop+1,n
            A(j+1)=j*a(j)
11       continue
         ntop=n
         factrl=a(n+1)
      else
         factrl=dexp(GAMMLN(n+1.))
      endif
      return
      end
            
c     FUNCTION HERMITE(noa,q)
c     implicit real*8(a-h,o-z)
c     if (noa.eq.-1) then 
c         hermite=0.
c     else if (noa.eq.0) then 
c         hermite=1.
c     else 
c         hermite=2.*q*hermite(noa-1,q)-2.*(noa-1)*hermite(noa-2,q)
c     endif
c     return
c     end
  
        
        FUNCTION FM0(sf,y)
	implicit real*8(a-h,o-y)      
	COMMON/OSCILATOR/W1,D,BETA
	facto=dexp(GAMMLN(sf))
        fm0log=0.5*(dlog(beta)-GAMMLN(sf)-y+sf*dlog(y))
c	fM0=dsqrt(beta/facto*dexp(-y)*y**sf)
        fM0=dexp(fm0log)
	return
	end

cFUNCTION fMORSE(nmor,sf,y)
c	implicit real*8(a-h,o-y)
c	facto=dsqrt(nmor*(sf+nmor))
c	if (nmor.eq.-1) then
c	   fMORSE=0.
c 	else if (nmor.eq.0) then
c           fMORSE=fM0(sf,y)
c 	else 
c	    var1=2*nmor+sf-1-y
c	    var2=dsqrt((nmor-1)*(nmor+sf-1))
c            fMORSE=(var1*fMORSE(nmor-1,sf,y)-var2
c     .       *fMORSE(nmor-2,sf,y))/facto
c        end if
c	return
c	end

c 
C*********************************************************
C**   FUNCION QUE CALCULA EL LOGARITMO DE LA FUNCION    **
C**   GAMAFACTORIAL UN ARGUMENTO REAL XX .              **   
C**   NUMERICAL RECIPIES                                **
C*********************************************************
C
       FUNCTION GAMMLN(Xx)
       implicit  real*8(a-h,o-y)
       ! STP,HALF,ONE,FPF,X,TMP,SER
       dimension cof(6)
       DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *      -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
       DATA HALF,ONE,FPF/.5D0,1.D0,5.5D0/
       X=Xx-ONE
       !write(*,*)'***',x,xx
       TMP=X+FPF
       TMP=(X+HALF)*dLOG(TMP)-TMP
       SER=ONE
       DO  J=1,6
         X=X+ONE 
         SER=SER+COF(J)/X
         !write(*,*),x,xx,ser
       END DO
       GAMMLN=TMP+dLOG(STP*SER) 
       !write(*,*)tmp,stp
       !write(*,*)     
       RETURN
       END



