c_____f77 pp.f fft.f
      program  pp     
      implicit none

      character*4 title

c_____variables for fft__________________________________________________
      integer    nexp1,nexp2,nexp3
      parameter (nexp1=8)                  ! <-- modify only this: exponent such that nr=2^nexp1
      parameter (nexp2=0)                  ! 
      parameter (nexp3=0)                  ! 
c_____variables for fft: automatically set, no not modify________________
      integer    ni                        ! size of sine table ...
      parameter (ni=2**(nexp1-2))          ! which has to be ni=n/4 (of the largest n for multi-D)
      integer    ierr,inv(ni),mmm(3)       ! some arrays
      real*8     s(ni)                     ! sine tables
c_____end variables for fft:_____________________________________________

      integer    nr,nj,nl,no
      integer    nall
      integer    nro
      integer    nth
      integer    nt

      parameter (nr=2**nexp1)

      parameter (nj=2)
      parameter (nl=2)
      parameter (no=2)
      parameter (nth=80)

      parameter (nall=nr*(nj+1)*(nl+1)*(2*no+1))


c-----reduced mass  and rot constant
      real*8     mass,B 
      parameter (mass=4*1836.d0)
      parameter (B=0.25021562/219474.d0)
    
c_____potential matrices
      real*8     pot_He,pot_Co,pot_La
      complex*16 EV_He(0:nl,0:nl,-no:no,nr)    
      complex*16 EV_Co(-no:no,-no:no,0:nl,0:nj)

      real*8     VV_La(0:nj,0:nj,-no:no)  
      complex*16 EV_La(0:nj,0:nj,-no:no)  

      real*8     coef1,coef2

      real*8     hh(4,4)
      complex*16 tout(4,4)

      real*8     POP_C ,POP_O
      real*8     POP_L ,POP_J      
      real*8     POP_T
      real*8     V_rot
      real*8     threej      
      real*8     NORM
      real*8     COS_2,COS_3,COS_X
      real*8     P_th

      logical    relax
      integer    i,j,k,l,o,jp,lp,lpp,op,opp
      integer    iter,count

      real*8     timestep,gauss,pi
      real*8     norm1,norm2

c_____wavefunction
      complex*16 zz(nr,0:nj,0:nl,-no:no)
      complex*16 zzin(nr,0:nj,0:nl,-no:no)
      complex*16 kk(nr)
      complex*16 z1(0:nl)
      complex*16 zsolap

c-----observables
      real*8     jjmat(0:nj,0:nj),cc(0:nj,0:nj)
      real*8     cosav,jav,cos2

      real*8     o1

c_____variablables for laser pulse
      real*8     e0,time,et,Tpul
c     parameter (e0=0.09000,Tpul=100.d0,f=0.057d0)
c     parameter (e0=0.0004400d0,Tpul=600.d0/0.024189d0)  !  this worked 
c     parameter (e0=0.0008800d0,Tpul=600.d0/0.024189d0)
      parameter (Tpul=450.d0/0.024189d0)
c     parameter (e0=0.00d0,Tpul=800.d0,f=0.057d0)

c-----variables for Ylmb
      real*8     rvec(3)
      real*8     Y(nl*nl)
      complex*16 Ylo(0:nth,0:nl,-no:no)
c_____variablables for radial grid   
      real*8     x0,length,dx,x(nr),kx(nr)
      parameter (x0=4.1d0,length=16.d0)
c_____variables for absorber
      integer    icut
      real*8     abx(nr)

c_____variablables for initial state and gaussian  
      integer    jinit,minit
      real*8     xeq,xwidth
      parameter (xeq=7.0d0,xwidth=2.d0)    ! initial state

      data       mmm /nexp1,nexp2,nexp3/   ! array of exponents, up to 3D

      real*8    phi,th,rr
      real*8    timeps
      integer    interv

c     title='test'

      read(5,'(a)') title
      read(5,*)    jinit,minit  
      read(5,*)    e0, interv

c_____write parameters
      open (33,file=title//'.param')
      write(33,*) x0,'  x0'
      write(33,*) length,'  length'
      write(33,*) nexp1,'  nexpth'
      write(33,*) nj,'  nj'
      write(33,*) e0,'  e0'
      write(33,*) Tpul,'  Tpul'
      write(33,*) jinit,minit,'  jinit and minit'
      close(33)



c     goto 999

      dx=length/dble(nr)
      call fft(zz,mmm,inv,s,0,ierr)

c_____create absorber
      icut=nr-nr/8
      abx=1.d0
      do i=1,nr
      if (i.gt.icut) abx(i)=dexp(-1.5d0*(dble(i-icut)/dble(nr-icut))**2)
      enddo
c
c-----calculate spherical harmonics

      pi=dacos(-1.d0)
      phi=0.d0
      nt=0
      

      do j=0,nth

      th= 2.d0*pi*dble(j)/dble(nth)

      rvec(1)=sin(th)*cos(phi)
      rvec(2)=sin(th)*sin(phi)
      rvec(3)=cos(th)

      call Ylmb(rvec,nl,Y)

      count=0

      do l=0,nl
       do o=-l,l
       count=count+1
       Ylo(j,l,o)=Y(count)
      enddo
      enddo

      enddo




c____some testing________________________________________
c     o1=-99.d0
c     do i=1,nr
c       x(i)=x0+dble(i-1)*dx
c       write(31,*) x(i)
c    .    , pot_He(0,0,1,x(i),mass,B),pot_He(1,1,1,x(i),mass,B)
c    .    , pot_He(2,2,1,x(i),mass,B),pot_He(3,3,1,x(i),mass,B)
c    .    , pot_He(4,4,1,x(i),mass,B),pot_He(5,5,1,x(i),mass,B)
c     enddo
c_____end testing________________________________________

     
       

  
        
c_____create laser interaction matrix element

      VV_La=0.d0

      do o=-no,no
      do j=abs(o),nj
      do jp=abs(o),nj
        VV_La(j,jp,o)=2.d0/3.d0*pot_La(j,jp,o,minit)
      enddo
      enddo
      enddo

c_____end create laser interaction matrix element
c     goto 999


      zz=(0.d0,0.d0)
      
      open(41,file="wfin.dat")
c_____create initial state________________________________
      do i=1,nr
       read(41,*) rr,(zz(i,jinit,l,0), l=1,nl)
      enddo

      
      call NORMALIZE(nall,norm1,zz)
      write(*,*)  POP_C(nr,nj,nl,no,zz,jinit,0,0)
c_____end create initial state____________________________





c==============================================================
c_____relaxation part
      open (32,file=title//'.relax')

      timestep=5.0d0
      relax=.true.

      call CRKIN(relax,0.5d0*timestep,mass,length,nr,kx,kk)
      call CR_He(relax,      timestep,mass,B,x0,dx,nr,nl,no,EV_He)
      call CR_Co(relax,0.5d0*timestep,B,nj,nl,no,EV_Co)

      do iter=1,15000
      write(*,*) 'iter relax',iter

        call NORMALIZE(nall,norm1,zz)

        call PROP__Co(nr,nj,nl,no,zz,EV_Co) 

        call PROP_KIN(nr,nj,nl,no,ni,zz,kk,inv,s,mmm)  
        call PROP__He(nr,nj,nl,no,zz,EV_He) 
        call PROP_KIN(nr,nj,nl,no,ni,zz,kk,inv,s,mmm)  

        call PROP__Co(nr,nj,nl,no,zz,EV_Co) 

        call NORMALIZE(nall,norm2,zz)

        write(32,80) iter,-0.5d0/timestep*dlog(norm2)
     .              ,(POP_L(nr,nj,nl,no,zz,l),l=0,nl)

        norm1=norm2

      enddo  ! end of relaxation loop
      call NORMALIZE(nall,norm2,zz)

      close(32)
c_____end__relaxation part
c==============================================================

c      goto 999




c==============================================================
c_____propoagation part with laser pulse
      open (84,file=title//'.J')
      open (85,file=title//'.O')
      open (86,file=title//'.L')
      open (88,file=title//'.cos2')
c      open (40,file=title//'.Pth')
 
c      open (100,file=title//'.Vrot')
c      open (101,file=title//'.autocorr')
c      open (102,file=title//'.kin')

      timestep=1.0/0.024189d0
      relax=.false.

      call CRKIN(relax,0.5d0*timestep,mass,length,nr,kx,kk)
      call CR_He(relax,      timestep,mass,B,x0,dx,nr,nl,no,EV_He)
      call CR_Co(relax,0.5d0*timestep,B,nj,nl,no,EV_Co)

      time = 0.d0

      do iter=1,450
      write(*,*) 'iter laser',iter

        time=dble(iter)*timestep

        et = e0*dsin(pi*time/Tpul)**4    ! field envelope

        call CR_La(relax,0.25d0*timestep,et,nj,no,VV_La,EV_La)

        call PROP__La(nr,nj,nl,no,zz,EV_La) 
        call PROP__Co(nr,nj,nl,no,zz,EV_Co) 
        call PROP__La(nr,nj,nl,no,zz,EV_La) 

        call PROP_KIN(nr,nj,nl,no,ni,zz,kk,inv,s,mmm)  
        call PROP__He(nr,nj,nl,no,zz,EV_He) 
        call PROP_KIN(nr,nj,nl,no,ni,zz,kk,inv,s,mmm)  

        call PROP__La(nr,nj,nl,no,zz,EV_La) 
        call PROP__Co(nr,nj,nl,no,zz,EV_Co) 
        call PROP__La(nr,nj,nl,no,zz,EV_La) 

        call APPLYABS(nr,nj,nl,no,zz,abx) 
        
        timeps=time/40000.d0
       
       if (mod(iter,interv).eq.0) then

       do k=0,nth

       th= 2.d0*pi*dble(k)/dble(nth)
       write(40,80) timeps,th,P_th(k,nth,nr,nj,nl,no,Ylo,zz)
      
       enddo  
        write(40,80)
       

c_____calculate observalbles here:
        write(84,80) timeps,et,(POP_J(nr,nj,nl,no,zz,j),j=0,nj)
c        write(85,80) timeps,et,(POP_O(nr,nj,nl,no,zz,o),o=-no,no)
        write(86,80) timeps,et,(POP_L(nr,nj,nl,no,zz,l),l=0,nl)
c        write(88,80)  timeps,COS_3(nr,nj,nl,no,zz,VV_La),NORM(nall,zz)
c     .                   ,COS_X(nr,nj,nl,no,zz,VV_La)
        
c        write(100,80) time,V_rot(nr,nj,nl,no,zz,mass,B,x0,dx)
c        call solap(zsolap,nr,nj,nl,no,zzin,zz)
c        write(101,80) time,cdabs(zsolap)**2
        
c_____end observalbles here:
       endif
      enddo  ! end of time loop


c==============================================================
c_____propoagation part without laser pulse

      do iter=451,120000
c      write(*,*) 'iter prop',iter

        time=dble(iter)*timestep

        call PROP__Co(nr,nj,nl,no,zz,EV_Co)
        call PROP_KIN(nr,nj,nl,no,ni,zz,kk,inv,s,mmm)
        call PROP__He(nr,nj,nl,no,zz,EV_He)
        call PROP_KIN(nr,nj,nl,no,ni,zz,kk,inv,s,mmm)
        call PROP__Co(nr,nj,nl,no,zz,EV_Co)

        timeps=time/40000.d0
         
         if (mod(iter,interv).eq.0) then
c         write(*,*) timeps
c_____calculate observalbles here:
        write(84,80) timeps,et,(POP_J(nr,nj,nl,no,zz,j),j=0,nj)
c        write(85,80) timeps,et,(POP_O(nr,nj,nl,no,zz,o),o=-no,no)
        write(86,80) time,et,(POP_L(nr,nj,nl,no,zz,l),l=0,nl)
c        write(88,80)  time,COS_3(nr,nj,nl,no,zz,VV_La),NORM(nall,zz)
c     .                   ,COS_X(nr,nj,nl,no,zz,VV_La)


c        write(100,80) time,V_rot(nr,nj,nl,no,zz,mass,B,x0,dx)
c        call solap(zsolap,nr,nj,nl,no,zzin,zz)
c        write(101,80) time,cdabs(zsolap)**2

        do k=0,nth

        th= 2.d0*pi*dble(k)/dble(nth)
         write(40,80) timeps,th,P_th(k,nth,nr,nj,nl,no,Ylo,zz)

       enddo 
        
         write(40,80)
c_____end observalbles here:
      endif
      enddo   ! end of tiome loop

c_____end__propagation part
      close(84)
      close(85)
      close(86)
      close(88)
      close(40)
c==============================================================


999   continue



c==============================================================
c_____final analysis 

c     do i=1,nr
c     write(89,80) x(i),(cdabs(zz(i,2,k,-1))**2,k=0,nl)
c     enddo


c_____end final analysis 
c==============================================================







c_____closing time.....     
      close(20)
      close(21)
      close(86)


80    format(80G16.8)
      end
c==============================================================
c==============================================================





c==============================================================
      function POP_C(nr,nj,nl,no,zz,j,l,o)
      implicit none

      integer    nr,nj,nl,no,i,j,l,o
      complex*16 zz(nr,0:nj,0:nl,-no:no)

      real*8     POP_C,o1

      o1=0.d0

      do i=1,nr
      o1=o1+cdabs(zz(i,j,l,o))**2
      enddo

      POP_C=o1

      return
      end
c==============================================================

      function POP_L(nr,nj,nl,no,zz,l)
      implicit none

      integer    nr,nj,nl,no,i,j,l,o
      complex*16 zz(nr,0:nj,0:nl,-no:no)

      real*8     POP_L,o1

      o1=0.d0

      do j=0,nj
      do o=-no,no
      do i=1,nr
      o1=o1+cdabs(zz(i,j,l,o))**2
      enddo
      enddo
      enddo

      POP_L=o1

      return
      end

c==============================================================

      function POP_O(nr,nj,nl,no,zz,o)
      implicit none

      integer    nr,nj,nl,no,i,j,l,o
      complex*16 zz(nr,0:nj,0:nl,-no:no)

      real*8     POP_O,o1

      o1=0.d0

      do j=0,nj
      do l=0,nl
      do i=1,nr
      o1=o1+cdabs(zz(i,j,l,o))**2
      enddo
      enddo
      enddo

      POP_O=o1

      return
      end


c==============================================================


      function POP_T(nr,nj,nl,no,zz)
      implicit none

      integer    nr,nj,nl,no,i,j,l,o
      complex*16 zz(nr,0:nj,0:nl,-no:no)

      real*8     POP_T,o1

      o1=0.d0

      do j=0,nj
      do l=0,nl
      do o=-no,no
      do i=1,nr
      o1=o1+cdabs(zz(i,j,l,o))**2
      enddo
      enddo
      enddo
      enddo

      POP_T=o1

      return
      end


c==============================================================


      function POP_J(nr,nj,nl,no,zz,j)
      implicit none

      integer    nr,nj,nl,no,i,j,l,o
      complex*16 zz(nr,0:nj,0:nl,-no:no)

      real*8     POP_J,o1

      o1=0.d0

      do l=0,nl
      do o=-no,no
      do i=1,nr
      o1=o1+cdabs(zz(i,j,l,o))**2
      enddo
      enddo
      enddo

      POP_J=o1

      return
      end

c==============================================================

      function COS_2(nr,nj,nl,no,zz,VV_La)
      implicit none

      integer     nr,nj,nl,no,nall
      integer     i,j,jp,l,o
      complex*16  zz(nr,0:nj,0:nl,-no:no)
      complex*16  zh(nr,0:nj,0:nl,-no:no)
      real*8      VV_La(0:nj,0:nj,-no:no)
      complex*16  c1,SCAL
      real*8      COS_2

      nall=nr*(nj+1)*(nl+1)*(2*no+1)
      

c_____cos^2

      c1=dcmplx(0.d0,0.d0)
      zh=(0.d0,0.d0)

      do l=0,0     ! diagonal
c     do l=0,nl    ! diagonal
c     do o=-no,no  ! diagonal
      do o=0,0     ! diagonal
      do i=1,nr    ! diagonal
      do j=abs(o),nj
      zh(i,j,l,o)=(0.d0,0.d0)
      do jp=abs(o),nj
      zh(i,j,l,o)=zh(i,j,l,o)+VV_La(j,jp,o)*zz(i,jp,l,o)
      enddo
      enddo
      enddo
      enddo
      enddo
 
      c1=SCAL(nall,zh,zz) 

      COS_2=(2.d0*dreal(c1)+1.d0)/3.d0

      return
      end


c==============================================================


      function COS_3(nr,nj,nl,no,zz,VV_La)
      implicit none

      integer     nr,nj,nl,no
      integer     i,j,jp,l,o
      complex*16  zz(nr,0:nj,0:nl,-no:no)
      real*8      VV_La(0:nj,0:nj,-no:no)
      complex*16  c1
      real*8      COS_3

c_____cos^2

      c1=dcmplx(0.d0,0.d0)

      do l=0,nl    ! diagonal
      do o=-no,no  ! diagonal
      do j=abs(o),nj
      do jp=abs(o),nj
      do i=1,nr    ! diagonal
      c1=c1+dconjg(zz(i,j,l,o))*VV_La(j,jp,o)*zz(i,jp,l,o)
      enddo
      enddo
      enddo
      enddo
      enddo

      COS_3=dreal(c1)+1.d0/3.d0

      return
      end

      function COS_X(nr,nj,nl,no,zz,VV_La)
      implicit none

      integer     nr,nj,nl,no
      integer     i,j,jp,l,o
      complex*16  zz(nr,0:nj,0:nl,-no:no)
      real*8      VV_La(0:nj,0:nj,-no:no)
      complex*16  c1
      real*8      COS_X

c_____cos^2

      c1=dcmplx(0.d0,0.d0)

      do l=0,nl    ! diagonal
      do o=-no,no  ! diagonal
      do j=abs(o),nj
      do i=1,nr    ! diagonal
      c1=c1+dconjg(zz(i,j,l,o))*VV_La(j,j,o)*zz(i,j,l,o)
      enddo
      enddo
      enddo
      enddo

      COS_X=dreal(c1)+1.d0/3.d0

      return
      end


c---------------------------------------------
      function V_rot(nr,nj,nl,no,zz,mass,B,x0,dx)
      implicit none

      integer     nr,nj,nl,no
      integer     i,j,jp,l,o
      complex*16  zz(nr,0:nj,0:nl,-no:no)
      complex*16  c1
      real*8      V_rot
      real*8      x(nr)
      real*8      x0,dx
      real*8      B,mass



       c1=dcmplx(0.d0,0.d0)

       do l=0,nl    ! diagonal
       do o=-no,no  ! diagonal
       do j=abs(o),nj
       do i=1,nr    ! diagonal
       x(i)=x0+dble(i-1)*dx
       c1=c1+dconjg(zz(i,j,l,o))*(dble(l)*dble(l+1)/2.d0/mass/x(i)/x(i)
     .  +B*dble(l)*dble(l+1))*zz(i,j,l,o)
       enddo
       enddo
       enddo
       enddo
  
       V_rot=dreal(c1)


      return

      end


c______________________________________________
      subroutine APPLYABS(nr,nj,nl,no,zz,abx)
      implicit none

      integer     nr,nj,nl,no
      integer     i,j,l,lp,o
      complex*16  zz(nr,0:nj,0:nl,-no:no)
      real*8      abx(nr)                      

c_____He
      do i=1,nr    ! diagonal
      do j=0,nj    ! diagonal
      do l=0,nl    ! diagonal
      do o=-no,no  ! diagonal
         zz(i,j,l,o)=abx(i)*zz(i,j,l,o)
      enddo
      enddo
      enddo
      enddo

      return
      end





c==============================================================

c______________________________________________
      subroutine PROP__He(nr,nj,nl,no,zz,EV_He)
      implicit none

      integer     nr,nj,nl,no
      integer     i,j,l,lp,o
      complex*16  zz(nr,0:nj,0:nl,-no:no)
      complex*16  EV_He(0:nl,0:nl,-no:no,nr)
      complex*16  c1(0:nl)

c_____He
      do i=1,nr    ! diagonal
      do j=0,nj    ! diagonal
      do o=-no,no    ! diagonal

      do l=0,nl
         c1(l)=dcmplx(0.d0,0.d0)
         do lp=0,nl
         c1(l)=c1(l)+EV_He(l,lp,o,i)*zz(i,j,lp,o)
         enddo
      enddo

      do l=0,nl
         zz(i,j,l,o)=c1(l)
      enddo

      enddo
      enddo
      enddo

      return
      end


c==============================================================

c______________________________________________
      subroutine PROP__Co(nr,nj,nl,no,zz,EV_Co)
      implicit none

      integer     nr,nj,nl,no
      integer     i,j,l,o,op
      complex*16  zz(nr,0:nj,0:nl,-no:no)
      complex*16  EV_Co(-no:no,-no:no,0:nl,0:nj)
      complex*16  c1(-no:no)

c_____Co
      do i=1,nr    ! diagonal
      do j=0,nj    ! diagonal
      do l=0,nl    ! diagonal

      do o=-no,no   
         c1(o)=dcmplx(0.d0,0.d0)
         do op=-no,no 
         c1(o)=c1(o)+EV_Co(op,o,l,j)*zz(i,j,l,op)
         enddo
      enddo

      do o=-no,no   
         zz(i,j,l,o)=c1(o)
      enddo

      enddo      ! loop over l
      enddo      ! loop over j
      enddo      ! loop over i

      return
      end


c==============================================================

c______________________________________________
      subroutine PROP__La(nr,nj,nl,no,zz,EV_La)
      implicit none

      integer     nr,nj,nl,no
      integer     i,j,jp,l,o
      complex*16  zz(nr,0:nj,0:nl,-no:no)
      complex*16  EV_La(0:nj,0:nj,-no:no)
      complex*16  c1(0:nj)

c_____He
      do i=1,nr    ! diagonal
      do l=0,nl    ! diagonal
      do o=-no,no  ! diagonal

      do j=0,nj
         c1(j)=dcmplx(0.d0,0.d0)
         do jp=0,nj
         c1(j)=c1(j)+EV_La(j,jp,o)*zz(i,jp,l,o)
         enddo
      enddo

      do j=0,nj
         zz(i,j,l,o)=c1(j)
      enddo

      enddo
      enddo
      enddo

      return
      end


c==============================================================

c_______________________________________________________________________
      subroutine PROP_KIN(nr,nj,nl,no,ni,zz,kk,inv,s,mmm)
      implicit none

c_____variables for fft: automatically set, no not modify________________
      integer    ni                        ! size of sine table ...
      integer    ierr,inv(ni),mmm(3)       ! some arrays
      real*8     s(ni)                     ! sine tables
c_____end variables for fft:_____________________________________________

      integer     nr,nj,nl,no,i,j,l,o
      complex*16  kk(nr)
      complex*16  zz(nr,0:nj,0:nl,-no:no)

      do j=0,nj
      do l=0,nl
      do o=-no,no
        call fft(zz(1,j,l,o),mmm,inv,s,-2,ierr)
        do i=1,nr
          zz(i,j,l,o)=kk(i)*zz(i,j,l,o)
        enddo
        call fft(zz(1,j,l,o),mmm,inv,s, 2,ierr)
      enddo  ! o
      enddo  ! l
      enddo  ! j

      return
      end

      
c==============================================================


c_______________________________________________________________________
      function CHECK_HERM(n,aa)
      implicit none

      integer n,i,j
      complex*16 c1,aa(n,n)
      complex*16 CHECK_HERM

      c1=(0.d0,0d0)
      do i=1,n
      do j=1,n
         c1=c1+dconjg(aa(j,i))*aa(i,j)
      enddo
      enddo
      CHECK_HERM=c1

      return
      end



c==============================================================

      


c_______________________________________________________________________
      subroutine NORMALIZE2(n,norm,z)
      implicit none

      integer    i,n
      real*8     norm
      complex*16 z(n)

      norm=0.d0
      do i=1,n
          norm=norm+cdabs(z(i))**2
      enddo

      if(norm .NE. 0) then
      do i=1,n
            z(i)=z(i)/dsqrt(norm)
      enddo
      endif

      return 
      end

c==============================================================
c_______________________________________________________________________
      function SCAL(n,z1,z2)
      implicit none
      integer    i,n
      complex*16 cc,SCAL   
      complex*16 z1(n),z2(n)
      cc=(0.d0,0.d0)
      do i=1,n
        cc=cc+dconjg(z1(i))*z2(i)
      enddo
      SCAL=cc
      return
      end

 
c_______________________________________________________________________
      function NORM(n,z)
      implicit none
      integer    i,n
      real*8     nn,NORM
      complex*16 z(n)
      nn=0.d0
      do i=1,n
        nn=nn+cdabs(z(i))**2
      enddo
      NORM=nn
      return
      end


c_______________________________________________________________________
      subroutine NORMALIZE(n,norm,z)
      implicit none

      integer    i,n
      real*8     norm
      complex*16 z(n)

      norm=0.d0
      do i=1,n
        norm=norm+cdabs(z(i))**2
      enddo

      if(norm .NE. 0) then
        do i=1,n
          z(i)=z(i)/dsqrt(norm)
        enddo
      endif

      return 
      end



c_______________________________________________________________________
      subroutine CR_He(relax,timestep,mass,B,x0,dx,nr,nl,no,EV_He)
      implicit none

      logical    relax
      integer    nr,nl,no
      integer    i,l,lp,o
      real*8     x0,dx,timestep,mass,B
      real*8     pot_He
      real*8     x(nr)
      real*8     hh(nl+1,nl+1)
      complex*16 EV_He(0:nl,0:nl,-no:no,nr)  
      complex*16 tout(nl+1,nl+1)


c______create EV_He operator

      EV_He=(0.d0,0d0)

      do i=1,nr
      do o=-no,no 

      hh=0.d0

      do l=abs(o),nl 
      do lp=abs(o),nl 
        x(i)=x0+dble(i-1)*dx
        hh(l+1,lp+1)=pot_He(l,lp,o,x(i),mass,B)
      enddo
      enddo

      do l=0,nl
      do lp=0,nl
        write(120,*) i,l,lp,hh(l+1,lp+1)
      enddo
      enddo

      call MATEXP(nl+1,relax,timestep,hh,tout)

      do l =abs(o),nl 
      do lp=abs(o),nl 
        EV_He(l,lp,o,i)=tout(l+1,lp+1)
      enddo
      enddo

      enddo   ! loop over o
      enddo   ! loop over i

      return
      end


c______________________________________________________________________
      subroutine CR_Co(relax,timestep,B,nj,nl,no,EV_Co)
      implicit none

      logical    relax
      integer    nj,nl,no
      integer    j,l,o,op
      real*8     timestep
      real*8     B,pot_Co
      complex*16 EV_Co(-no:no,-no:no,0:nl,0:nj)

      real*8     hh(2*no+1,2*no+1)
      complex*16 tt(2*no+1,2*no+1)

c______create EV_Co operator

      EV_Co=(0.d0,0d0)

      do j=0,nj    ! diagonal 
      do l=0,nl    ! diagonal 

      hh=0.d0

      do o= -no,no
      do op=-no,no
         if (abs( o).gt.j) goto 333
         if (abs(op).gt.j) goto 333
         if (abs( o).gt.l) goto 333
         if (abs(op).gt.l) goto 333
         hh(o+no+1,op+no+1)= pot_Co(o,op,l,j,B) 
333      continue
      enddo
      enddo

      call MATEXP(2*no+1,relax,timestep,hh,tt)

      do o= -no,no
      do op=-no,no
         EV_Co(o,op,l,j)=tt(o+no+1,op+no+1) 
      enddo
      enddo

      enddo   ! loop over l
      enddo   ! loop over j

      return
      end




c__________________________________________________________
      subroutine CR_La(relax,timestep,et,nj,no,VV_La,EV_La)
      implicit none

      logical    relax
      integer    nj,no
      integer    j,jp,o   
      real*8     timestep
      real*8     pot_la,et
      real*8     VV_La(0:nj,0:nj,-no:no)
      complex*16 EV_La(0:nj,0:nj,-no:no) 

      real*8     hh(nj+1,nj+1)
      complex*16 tt(nj+1,nj+1)


c_____create EV_La operator

      EV_La=(0.d0,0d0)

      do o=-no,no

      hh=0.d0

      do j=abs(o),nj
      do jp=abs(o),nj
        hh(j+1,jp+1)=-et*VV_La(j,jp,o)
      enddo
      enddo

      call MATEXP(nj+1,relax,timestep,hh,tt)

      do j =abs(o),nj
      do jp=abs(o),nj
        EV_La(j,jp,o)=tt(j+1,jp+1)
      enddo
      enddo

      enddo   ! loop over o

      return
      end



c_______________________________________________________________________
      subroutine CRKIN(relax,timestep,mass,length,n,kx,kk)
      implicit none
      integer    n,j
      logical    relax
      real*8     timestep,kx(n)
      complex*16 kk(n)
      real*8     mass,length
      real*8     pi,twopi


      pi=dacos(-1.d0)
      twopi=2.d0*pi

c    * creating kk file: e^{- timestep p^2/2m}

      do j=1,n
        if (j.le.n/2+1) kx(j)= dble(j-1)  *(twopi/length)
        if (j.gt.n/2+1) kx(j)= dble(j-1-n)*(twopi/length)
      enddo

      do j=1,n
         kk(j)=cdexp(-dcmplx(0.d0,timestep*0.5d0*kx(j)**2/mass))
      if (relax)
     .   kk(j)=cdexp(-dcmplx(timestep*0.5d0*kx(j)**2/mass,0.d0))
      enddo

      return
      end
c_______________________________________________________________________



c_______________________________________________________________________
      function gauss(xc,sig,x)
      implicit none
      real*8 gauss,xc,sig,x
      gauss=dexp(-(x-xc)**2/(2.d0*sig**2))
      return
      end
c_______________________________________________________________________



c_______________________________________________________________________
      subroutine MATEXP(n,relax,timestep,hh,tout)
c     hh is REAL
      implicit none

      logical    relax
      integer    n,j,k,info,lwork

      real*8     ew(n)
      real*8     work(3*n)             ! dim=3*n-2
      real*8     hh(n,n)
      real*8     timestep
      complex*16 te(n,n),tout(n,n),ex(n)

      call DSYEV('V','U',n,hh,n,ew,work,3*n,info)

c______________________________________________________________________
      do 10 j=1,n
                    ex(j)=cdexp(dcmplx(0.d0,-timestep*ew(j)))   ! real time prop
         if (relax) ex(j)=cdexp(dcmplx(-timestep*ew(j),0.d0))   ! imag time prog

         do 12 k=1,n
            te(j,k)=hh(k,j)*ex(j)
12          continue
10    continue
      call mm_r(n,hh,te,tout)

      end
c______________________________________________________________________




c______________________________________________________________________
      subroutine mm_r(n,a,b,c)   ! a is REAL, b / c COMPLEX
      implicit none
      integer    n,i,j,k
      real*8     a(n,n)
      complex*16 b(n,n),c(n,n)
      do k=1,n
      do i=1,n
      c(i,k)=dcmplx(0.d0,0.d0)
      do j=1,n
      c(i,k)=c(i,k)+a(i,j)*b(j,k)
      enddo
      enddo
      enddo
      return
      end
c     Computes the vector product C = A.B
c                                 =   = =
c     where A is a (nlA,ncA) matrix
c           B is a (ncA,ncB) matrix
c           C is a (nlA,ncB) matrix
c           stcX is the stride between row  elements of X-matrix
c           stlX is the stride between line elements of X-matrix
c
c     The product is computed as (in order to vectorize) :
c     for k=1,ncB
c     for j=1,ncA
c     for i=1,nlB
c     C(i,k) = C(i,k) + A(i,j).B(j,k)
c     =================================================================


   
      subroutine solap(zsolap,nr,nj,nl,no,zw1,zw2)
      
      implicit none
      integer nr,nj,nl,no,i,j,l,o     
      complex*16 zsum,zsolap
      complex*16 zw1(nr,0:nj,0:nl,-no:no)
      complex*16 zw2(nr,0:nj,0:nl,-no:no)


      zsum=(0.,0.)

      do j=0,nj
      do l=0,nl
      do o=-no,no
      do i=1,nr
        zsum=zsum+dconjg(zw1(i,j,l,o))*zw2(i,j,l,o)
      end do
      enddo
      enddo
      enddo

      zsolap=zsum

      return
      end


c ----------------------------------------------
c ----------------------------------------------

       function P_th(indk,nth,nr,nj,nl,no,Ylo,zz)

       implicit none
       integer icount,jcount,indk
       integer nr,nj,nl,no,nth,i,j,k,l,lp,o,th
       complex*16 zz(nr,0:nj,0:nl,-no:no)
       complex*16 Ylo(0:nth,0:nl,-no:no)
       real*8 Pth,P_th

       Pth=0.d0
        do i=1,nr
          do j=0,nj
          do o=-no,no
           do l=abs(o),nl
           icount=o+l
           do lp=abs(o),nl
            jcount=o+lp
c           write(*,*) icount, jcount
            Pth=Pth+dconjg(zz(i,j,lp,o))*zz(i,j,l,o)
     .      *Ylo(indk,l,o)*Ylo(indk,lp,o)
           enddo 
         enddo
       enddo
       enddo
       enddo   
       P_th=Pth
       
       return
       end
