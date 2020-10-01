c---------------------------- vch3he.f ---------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
c      collection of subroutines for the 3D van der Waals potential
c      between CH3 and He derived from ab initio calculations and
c      fitted to an expansion into a HFD like form with coefficients
c      expanded over real spherical harmonics.
c                                                  m. lewerenz feb/11
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c----------------------------- vch3he ----------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine vch3he(r,costh,phi,v,nv,iout)
c
c      3D potential surface model for CH3-He. input is in polar
c      coordinates relative to the center of mass of CH3. theta
c      is measured from the positive z-axis pointing away from
c      the hydrogens, phi=0 is a plane containing a hydrogen.
c      input distance r in angstroem units, phi in radian, output
c      v is in cm-1
c
c      r     : vector of length nv with polar distances in AA units; input
c      costh : vector of length nv with cos(theta); input
c      phi   : vector of length nv with azimuthal angles in radian; input
c      v     : vector of length nv for potential values in cm-1; output
c      nv    : length of vectors r, costh, phi, and v; input
c      iout  : unit for print output, iout.le.0 suppresses output; input
c
c      during the first call a parameter file is read which is just
c      a short bit of the least squares fit output file. these values
c      are directly injected into the actual potential routine in
c      all subsequent calls. currently the name of this parameter file
c      is ch3he.p
c                                                  m. lewerenz feb/11
c
      include 'real.inc'
      parameter(mp=50,ipdata=10,model=1)
      dimension r(nv),costh(nv),phi(nv),v(nv),p(mp)
      data icall /0/
      save
c
      if(icall.eq.0) then
        icall=1
        call getpar(ipdata,'ch3he.p',p,mp,model,iout)
      end if
      call ch3he(r,costh,phi,v,nv,p,np,iout)
c
      return
      end
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c------------------------------ getpar ---------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine getpar(ipdata,pfile,p,mp,mod,iout)
c
c      reads a file with starting parameters for the fit or to
c      initialise the potential evaluation. the parameters are in the
c      format as they are written in the fitting routine output file.
c      just cut out the relevant lines and save them into a file.
c      a dummy call to ch3he is made in order to retrieve information
c      about the model (lmax,np).
c                                                 m. lewerenz feb/11
c
      include 'real.inc'
      dimension p(mp)
      character pfile*(*),intext*15
      common /fitpar/ model,lmax,np
c
      model=mod
      call ch3he(r,costh,phi,v,0,p,np,iout)
c
c     call openf(ipdata,pfile,'o','k','f')
      open(unit=ipdata,file=pfile,status='unknown')
      rewind (unit=ipdata)
      call getstr(pfile,lpf,0)
      write(iout,'(3a)')
     &        '  surface parameter input from file ',pfile(1:lpf),' : '
      npp=0
   10 continue
      call nextin(ipdata,0,items,0)
      if(items.le.0) goto 20
      read(ipdata,'(a15,e18.10)',end=20) intext,pnew
      write(iout,'(a15,e18.10)') intext,pnew
      npp=npp+1
      if(npp.le.mp) then
        p(npp)=pnew
      else
        call errprt(iout,'getpar','parameter buffer too small',-1)
      end if
      goto 10
   20 continue
      if(np.ne.npp) then
        call errprt(iout,'getpar','inconsistent parameter file',-1)
c     else
c       do i=1,np
c         write(iout,'(a,i2,e15.7)') '  parameter ',i,p(i) 
c       end do
      end if
c
      return
      end
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c------------------------------ ch3he ----------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine ch3he(r,costh,phi,v,nv,p,np,iout)
c
c      3D potential surface model for CH3-He.
c
c      r is the distance from He to the c.o.m of CH3 in angstroem,
c      costh the cosine of the polar angle relative to the symmetry
c      axis of CH3, phi is the azimuthal angle, phi=0 is a hydrogen.
c
c      V = A0*exp[-b*(r-re)]-C6/r**6-C8/R**8
c      
c      X(theta,phi) = sum_i x_lm T_lm(cos(theta),phi), X=re, b, C6, (C8) 
c      with symmetry restrictions on the allowed l and m values due to
c      the periodicity of the potential (m=0,3k; l-m even at alpha=90).
c                                                    m. lewerenz feb/11
c
      include 'real.inc'
      parameter (maxlm=50)
      dimension r(nv),costh(nv),phi(nv),v(nv),p(np),
     &          lm(2,maxlm),tlm(maxlm)
      common /fitpar/ model,lmax,mp
      data icall/0/
      save 
c
      rc=1.8

      open(110,file='vch3he.out')
      open(111,file='coef.out')
c
      if(icall.eq.0) then
        icall=1
        lmax=6
        call lmlist(0,lmax,lm,maxlm,nlm,0)
        if(model.eq.0) then
          write(iout,'(/a)')
     &         '  T_lm expansion up to l=6 for R_e, b, and C6 '
          mp=3*nlm+2
        else if(model.eq.1) then
          write(iout,'(/a)')
     &         '  T_lm expansion up to l=6 for R_e, b, C6, and C8 '
          mp=4*nlm+1
        else
          call errprt(iout,'ch3he','undefined potential model',-1)
        end if
        rfp=four*sqrt(atan(one))
      end if
c
      a0=p(1)
      do i=1,nv
        npbase=2
        call tlmtab(costh(i),phi(i),lm,tlm,nlm,iout)
        call vecdot(p(npbase),tlm,re,nlm)
        npbase=npbase+nlm
        call vecdot(p(npbase),tlm,b,nlm)
        npbase=npbase+nlm
        call vecdot(p(npbase),tlm,c6,nlm)
        npbase=npbase+nlm
        if(model.eq.0) then
          c8=p(npbase)
        else if(model.eq.1) then
          call vecdot(p(npbase),tlm,c8,nlm)
          c8=c8*rfp
        end if
        re=re*rfp
        b=b*rfp
        c6=c6*rfp


c
        if (r(i).lt.rc) then 
        
        rsqi=one/rc**2
        rsqi3=rsqi*rsqi*rsqi
        vdisp=-(c8*rsqi+c6)*rsqi3
        vrep=a0*exp(-b*(rc-re))
        vrc=vrep+vdisp

        vdispd=(6*c6+8*c8*rsqi)*rsqi3/rc
        vrcd=-b*vrep+vdispd

        b1=-vrcd/vrc
        a1=exp(b1*rc)*vrc

        v(i)=a1*exp(-b1*r(i))
        write(110,*) v(i)
        write(111,*) costh,phi


c
       write(iout,'(i6,a,4f10.7,2e14.7)')
     &    i,'  costh, phi, re, b, C6 : ',costh(i),phi(i),re,b,c6,c8
        write(iout,'(10f12.5)') (tlm(k),k=1,maxtlm)
c        write(iout,*) a0
c        write(iout,*) vrep,vdisp

       else
        
        rsqi=one/r(i)**2
        rsqi3=rsqi*rsqi*rsqi
        vdisp=-(c8*rsqi+c6)*rsqi3
        vrep=a0*exp(-b*(r(i)-re))
        v(i)=vrep+vdisp
     
        write(iout,'(i6,a,4f10.7,2e14.7)')
     &    i,'  costh, phi, re, b, C6 : ',costh(i),phi(i),re,b,c6,c8
        write(iout,'(10f12.5)') (tlm(k),k=1,maxtlm)      
       write(110,*) v(i)
       endif

      end do
      return
      end
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c------------------------------ tlmtab ---------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine tlmtab(costh,phi,lm,tlm,nlm,iout)
c
c      creates a table of the required tesseral harmonics
c
c                                                    m. lewerenz feb/11
      include 'real.inc'
      dimension lm(2,nlm),tlm(nlm)
c
      do i=1,nlm
        l=lm(1,i)
        m=lm(2,i)
        tlm(i)=tesslm(l,m,costh,phi,iout)
      end do
      return
      end
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c------------------------------ tesslm ---------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      function tesslm(l,m,costh,phi,iout)
c
c      tesseral harmonic T_lm(z,phi) , z=cos(theta)
c      real version of spherical harmonics according to
c
c      T_lm = sqrt(2) N_lm P_lm(z) cos(m*phi)    m>0
c      T_lm = sqrt(2) N_lm P_l-m(z) sin(m*phi)   m<0
c      T_lm = Y_l0 = N_l0 P_l0(z)                m=0
c      N_lm = sqrt((2l+1)/(4*pi) (l-|m|)!/(l+|m|)!)
c                                                    m. lewerenz feb/11
c
      include 'real.inc'
      rfpi=one/(four*sqrt(atan(one)))
      if(m.eq.0) then
        tlm=rfpi*dpleg(l,0,costh,iout)
      else
        ma=abs(m)
        rnlm=rfpi*sqrt(two*l+one)*
     &       exp(half*(dgamln((one+l-ma))-dgamln((one+l+ma))))
        plm=dpleg(l,ma,costh,iout)
        tlm=sqrt(two)*rnlm*plm
c       write(iout,'(a,2i3,4e15.8)') 'tesslm',l,m,costh,phi,rnlm,plm
        if(m.gt.0) then
          tlm=tlm*cos(ma*phi)
        else
          tlm=tlm*sin(ma*phi)
        end if
      end if
      tesslm=tlm
c     write(iout,'(a,i1,i2,a,2f8.5,a,e15.8)')
c    &      ' t_',l,m,'(',costh,phi,') = ',tesslm
      return
      end
c
c-----------------------------------------------------------------------
c
      function dpleg(ll,m,x,iout)
c
c      !!!! may be unreliable in a parallel calling situation !!!!
c      associated legendre polynomials using recursion similar to
c      arfken eq. (12.17a) using gradshteyn 8.731 2.; see also ams55
c      8.5.3 and introduction p. xiii; convention for negative m as
c      arfken eq. (12.81).         dpleg=P(l,m;x)
c
c      ll   : lower index l of legendre polynomial; input
c      m    : upper index m of legendre polynomial; input
c      x    : function argument; input
c      iout : unit number for messages, silent for iout.le.0; input
c      subroutines called: dlogam,errprt          m. lewerenz 13/sep/89
c
      implicit real*8 (a-h,o-z)
      parameter (maxn=100,zero=0.d0,one=1.d0,two=2.d0)
      dimension pleg(maxn)
      save xold,mold,imax
      data xold,mold,imax/zero,-1,-1/
      zgamln(xx)=dlogam(xx,iout)
c
c      check for bad arguments
c
      ma=iabs(m)
      if(ma.gt.ll.or.ll.lt.0) then
        call errprt(iout,'dpleg','illegal argument(s)',1)
        dpleg=zero
        return
      end if
c
      indx=ll-ma
      if(x.eq.xold.and.ma.eq.mold.and.imax.gt.0) goto 15
      mold=ma
      xold=x
      imax=-1
c
c      get starting values for recursion for each new value of m
c
      pmm=one
      if(ma.eq.0) goto 10
      somx2=sqrt((one-x)*(one+x))
      mfac=1
      do 5 i=1,ma
      pmm=-pmm*somx2*mfac
      mfac=mfac+2
    5 continue
   10 pleg(1)=pmm
      pleg(2)=pmm*x*(2*ma+1)
      imax=2
c
c      is this value already in the table ?
c
   15 if(indx.lt.imax) goto 30
      minl=imax
      maxl=indx
      if(indx.ge.maxn) maxl=maxn-1
      onemtm=1-2*ma
      twox=two*x
c
c      now recursion for increasing l-m to fill the table
c
      if(minl.eq.maxn) goto 40
      do 20 l=minl,maxl
      pleg(l+1)=twox*pleg(l)-pleg(l-1)-(x*pleg(l)-pleg(l-1))*onemtm/l
   20 continue
      imax=maxl+1
      if(indx.ge.maxn) goto 40
   30 dpleg=pleg(indx+1)
      goto 60
c
c      values outside range of table generated from last table entries
c
   40 plm1=pleg(maxn)
      plm2=pleg(maxn-1)
      do 50 lmm=maxn,indx
      dpleg=twox*plm1-plm2-(x*plm1-plm2)*onemtm/lmm
      plm2=plm1
      plm1=dpleg
   50 continue
   60 if(m.ge.0) return
      fac=exp(zgamln(one+indx)-zgamln(one+ll+ma))
      if(mod(ma,2).ne.0) fac=-fac
      dpleg=dpleg*fac
      return
      end
c
c
c----------------------------------------------------------------------
c
      FUNCTION dgamln(XX,iout)
c
C      LOGARITHMIC GAMMA FUNCTION FROM "NUMERICAL RECIPES" CHAP. 6.1
C      WITH UPGRADES FOR XX<1. IOUT IS A UNIT FOR ERROR MESSAGES
c      METHOD: lanczos' formula with ERROR < 2*10**-10 WITH RECURSION
c              AND REFLEXION FORMULA FOR XX<1.
c      iout is a fortran unit for error messages
c      SUBROUTINES NEEDED: errprt                  m. lewerenz 16.9.89
c
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,half=0.5d0,one=1.d0,fpf=5.5d0,thrsh=1.d-10,
     #           pi=3.141592653589793d0,stp=2.50662827465d0)
      dimension cof(6)
      save cof
      DATA COF/ 76.18009173D0,-86.50532033D0,24.01409822D0,
     *         -1.231739516D0,  .120858003D-2,  -.536382D-5/
      DATA SINLOG/ZERO/
c
      X=XX-ONE
      if(xx.lt.one) then
        x=one-xx
        sinx=sin(pi*x)/pi
        if(sinx.gt.thrsh) then
          sinlog=log(sinx)
          x=x-one
        else
          call errprt(iout,'dgamln','negative gamma function',1)
          dgamln=zero
          return
        end if
      end if
c
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
c$dir scalar
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      dgamln=TMP+LOG(STP*SER)
      if(xx.lt.one) dgamln=-dgamln-sinlog
      RETURN
      END
c
c----------------------------------------------------------------------
c
      FUNCTION dlogam(X,iout)
C
c     i changed the original name which was gamln, because there is
c     another routine with this name in this package. this one came
c     with the besj package.
c     see also: w.j. cody, k.e. hillstrom, chebyshev approximations for
c               the natural logarithm of the gamma function,
c               math. comp. 21, 198 (1967)
c                                                    m. lewerenz oct/89
c
c -----------------------------------------------------------------
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES                     *
C  *                   A PRIME CONTRACTOR TO THE                       *
C  *                UNITED STATES DEPARTMENT OF ENERGY                 *
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE   *
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE      *
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,    *
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES  *
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL      *
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR     *
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS      *
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE          *
C  * OWNED RIGHTS.                                                     *
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS     *
C  * PART IS SAND77-1441.                                              *
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C     WRITTEN BY D. E. AMOS, SEPTEMBER, 1977.
C
C     REFERENCES
C         SAND-77-1518
C
C         COMPUTER APPROXIMATIONS BY J.F.HART, ET.AL., SIAM SERIES IN
C         APPLIED MATHEMATICS, WILEY, 1968, P.135-136.
C
C         NBS HANDBOOK OF MATHEMATICAL FUNCTIONS, AMS 55, BY
C         M. ABRAMOWITZ AND I.A. STEGUN, DECEMBER. 1955, P.257.
C
C     ABSTRACT
C         GAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
C         X.GT.0. A RATIONAL CHEBYSHEV APPROXIMATION IS USED ON
C         8.LT.X.LT.1000., THE ASYMPTOTIC EXPANSION FOR X.GE.1000. AND
C         A RATIONAL CHEBYSHEV APPROXIMATION ON 2.LT.X.LT.3. FOR
C         0.LT.X.LT.8. AND X NON-INTEGRAL, FORWARD OR BACKWARD
C         RECURSION FILLS IN THE INTERVALS  0.LT.X.LT.2 AND
C         3.LT.X.LT.8. FOR X=1.,2.,...,100., GAMLN IS SET TO
C         NATURAL LOGS OF FACTORIALS.
C
C     DESCRIPTION OF ARGUMENTS
C
C         INPUT
C           X      - X.GT.0
C
C         OUTPUT
C           GAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT X
C
C     ERROR CONDITIONS
C         IMPROPER INPUT ARGUMENT - A FATAL ERROR
c
c      subroutines called :errprt
c
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,two=2.d0,half=0.5d0)
      DIMENSION GLN(100),P(5),Q(2),PCOE(9),QCOE(4)
C
      SAVE XLIM1,XLIM2,RTWPIL,P,Q,PCOE,QCOE,GLN
      DATA XLIM1,XLIM2,RTWPIL/    8.d0,  1000.d0 , 9.18938533204673d-01/
      DATA  P              / 7.66345188000000d-04,-5.94095610520000d-04,
     1 7.93643110484500d-04,-2.77777775657725d-03, 8.33333333333169d-02/
      DATA  Q              /-2.77777777777778d-03, 8.33333333333333d-02/
      DATA PCOE            / 2.97378664481017d-03, 9.23819455902760d-03,
     1 1.09311595671044d-01, 3.98067131020357d-01, 2.15994312846059d+00,
     2 6.33806799938727d+00, 2.07824725317921d+01, 3.60367725300248d+01,
     3 6.20038380071273d+01/
      DATA QCOE            / 1.00000000000000d+00,-8.90601665949746d+00,
     1 9.82252110471399d+00, 6.20038380071270d+01/
      DATA(GLN(I),I=1,60)  /         2*zero      , 6.93147180559945d-01,
     1 1.79175946922806d+00, 3.17805383034795d+00, 4.78749174278205d+00,
     2 6.57925121201010d+00, 8.52516136106541d+00, 1.06046029027453d+01,
     3 1.28018274800815d+01, 1.51044125730755d+01, 1.75023078458739d+01,
     4 1.99872144956619d+01, 2.25521638531234d+01, 2.51912211827387d+01,
     5 2.78992713838409d+01, 3.06718601060807d+01, 3.35050734501369d+01,
     6 3.63954452080331d+01, 3.93398841871995d+01, 4.23356164607535d+01,
     7 4.53801388984769d+01, 4.84711813518352d+01, 5.16066755677644d+01,
     8 5.47847293981123d+01, 5.80036052229805d+01, 6.12617017610020d+01,
     9 6.45575386270063d+01, 6.78897431371815d+01, 7.12570389671680d+01,
     A 7.46582363488302d+01, 7.80922235533153d+01, 8.15579594561150d+01,
     B 8.50544670175815d+01, 8.85808275421977d+01, 9.21361756036871d+01,
     C 9.57196945421432d+01, 9.93306124547874d+01, 1.02968198614514d+02,
     D 1.06631760260643d+02, 1.10320639714757d+02, 1.14034211781462d+02,
     E 1.17771881399745d+02, 1.21533081515439d+02, 1.25317271149357d+02,
     F 1.29123933639127d+02, 1.32952575035616d+02, 1.36802722637326d+02,
     G 1.40673923648234d+02, 1.44565743946345d+02, 1.48477766951773d+02,
     H 1.52409592584497d+02, 1.56360836303079d+02, 1.60331128216631d+02,
     I 1.64320112263195d+02, 1.68327445448428d+02, 1.72352797139163d+02,
     J 1.76395848406997d+02, 1.80456291417544d+02, 1.84533828861449d+02/
c
      DATA(GLN(I),I=61,100)/ 1.88628173423672d+02, 1.92739047287845d+02,
     1 1.96866181672890d+02, 2.01009316399282d+02, 2.05168199482641d+02,
     2 2.09342586752537d+02, 2.13532241494563d+02, 2.17736934113954d+02,
     3 2.21956441819130d+02, 2.26190548323728d+02, 2.30439043565777d+02,
     4 2.34701723442818d+02, 2.38978389561834d+02, 2.43268849002983d+02,
     5 2.47572914096187d+02, 2.51890402209723d+02, 2.56221135550010d+02,
     6 2.60564940971863d+02, 2.64921649798553d+02, 2.69291097651020d+02,
     7 2.73673124285694d+02, 2.78067573440366d+02, 2.82474292687630d+02,
     8 2.86893133295427d+02, 2.91323950094270d+02, 2.95766601350761d+02,
     9 3.00220948647014d+02, 3.04686856765669d+02, 3.09164193580147d+02,
     A 3.13652829949879d+02, 3.18152639620209d+02, 3.22663499126726d+02,
     B 3.27185287703775d+02, 3.31717887196928d+02, 3.36261181979198d+02,
     C 3.40815058870799d+02, 3.45379407062267d+02, 3.49954118040770d+02,
     D 3.54539085519441d+02, 3.59134205369575d+02/
C
      IF(X) 90,90,5
   5  NDX=X
      T=X-NDX
      IF(T.EQ.zero) GO TO 51
      DX=XLIM1-X
      IF(DX.LT.zero) GO TO 40
C
C     RATIONAL CHEBYSHEV APPROXIMATION ON 2.LT.X.LT.3 FOR GAMMA(X)
C
      NXM=NDX-2
      PX=PCOE(1)
      DO 10 I=2,9
   10 PX=T*PX+PCOE(I)
      QX=QCOE(1)
      DO 15 I=2,4
   15 QX=T*QX+QCOE(I)
      DGAM=PX/QX
      IF(NXM.GT.0) GO TO 22
      IF(NXM.EQ.0) GO TO 25
C
C     BACKWARD RECURSION FOR 0.LT.X.LT.2
C
      DGAM=DGAM/(one+T)
      IF(NXM.EQ.-1) GO TO 25
      DGAM=DGAM/T
      dlogam=LOG(DGAM)
      RETURN
C
C     FORWARD RECURSION FOR 3.LT.X.LT.8
C
   22 XX=two+T
      DO 24 I=1,NXM
      DGAM=DGAM*XX
   24 XX=XX+one
   25 dlogam=LOG(DGAM)
      RETURN
C
C     X.GT.XLIM1
C
   40 RX=one/X
      RXX=RX*RX
      IF((X-XLIM2).LT.zero) GO TO 41
      PX=Q(1)*RXX+Q(2)
      dlogam=PX*RX+(X-half)*LOG(X)-X+RTWPIL
      RETURN
C
C     X.LT.XLIM2
C
   41 PX=P(1)
      SUM=(X-half)*LOG(X)-X
      DO 42 I=2,5
      PX=PX*RXX+P(I)
   42 CONTINUE
      dlogam=PX*RX+SUM+RTWPIL
      RETURN
C
C     TABLE LOOK UP FOR INTEGER ARGUMENTS LESS THAN OR EQUAL 100.
C
   51 IF(NDX.GT.100) GO TO 40
      dlogam=GLN(NDX)
      RETURN
 90   CONTINUE
      call errprt(iout,'dlogam','zero or negative argument',1)
      dlogam=zero
      RETURN
      END
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c-----------------------------------------------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine lmlist(isym,lmax,lm,maxlm,nlm,iout)
c
c      generates a list of (l,m) indices for the basis,
c      m is 0,3k due to the ch3 symmetry and the T_lm functions
c      have to be even in theta if isym.ne.0 (planar ch3).
c      the additional condition V(phi)=V(-phi) eliminates the sin terms
c      corresponding to m < 0
c                                                    m. lewerenz feb/11
c
      dimension lm(2,maxlm)
c
      ifail=0
      nlm=0
c
      do ll=1,lmax+1
        l=ll-1
        do mm=1,ll
          m=mm-1
          if(mod(m,3).eq.0) then
            lmok=1
            if(isym.ne.0) then
              lmok=0
              if(mod(l-m,2).eq.0) lmok=1
            end if
            if(lmok.eq.1) then
              nlm=nlm+1
              if(nlm.le.maxlm) then
                lm(1,nlm)=l
                lm(2,nlm)=m
              else
                ifail=1
              end if
            end if
          end if
        end do
      end do
c
      if(iout.gt.0) then
        write(iout,'(/a)') '  (l,m) basis indices from lmlist:'
        do i=1,nlm
          write(iout,'(2x,2i3)') lm(1,i),lm(2,i)
        end do
      end if
c
      if(ifail.ne.0) call errprt(iout,'lmlist','basis too large',-1)
      return
      end
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c---------------------------- last line --------------------------------
