c
c-----------------------------------------------------------------------
c
      subroutine hepair(ipot,irunit,ivunit,iout)
c
c      initializes the common blocks holding pair potential parameters
c      for helium
c
c      iout   : message unit number, silent for iout.le.0
c      subroutines called: none                     m. lewerenz jan/93
c
      include 'real.inc'
      character hfdref*40,hfdcrf*40,ttref*40,ttmref*40,ljref*40,
     #          ttkref*40,ttmrf2*40,mref*40,text*2
      data text/'he'/,ljref/' '/
c
c ---------------------- helium parameters ---------------------------
c   hfd   :    aziz in inert gases, m.l. klein ed. springer 1984
c
c   lennard-jones parameters
      data eps,re/7.6093d0,2.963d0/
c   hfd-b parameters
      data c6b,c8b,c10b,astarb,rrmb,db,alsb,betsb,gammab,epsb
     #             /1.461d0,14.11d0,183.5d0,184431.01d0,2.963d0,
     #              1.4826d0,10.43329537d0,-2.27965105d0,zero,7.6093d0/,
     #     hfdref/'mol. phys. 61, 1487 (1987)'/
c   hfd parameters
      data c6c,c8c,c10c,astarc,rrmc,dc,alsc,betsc,gammac,epsc
     #             /1.461d0,14.11d0,183.5d0,544850.4d0,2.9673d0,
     #              1.241314d0,13.353384d0,zero,zero,7.5064d0/,
     #     hfdcrf/'inert gases, m.l. klein ed., 1984'/
c  tang-toennies data, aziz et al., z. phys. d, 21, 251 (1991)
      data ahett/25.42830068d0/,betahe/2.41734082d0/,
     #     cc6,cc8,cc10,cc12,cc14,cc16/1.461,14.11,183.5,3213.5,
     #     75779.4,2406339./,
     #     ttref/'z. phys. d, 21, 251 (1991)'/
c  modified tang-toennies data
      data ahettm/6.42140197d0/,b1mhe,b2mhe,b3mhe/1.83899568d0,
     #     -0.08401d0,0.003761d0/,xcut,rettm/4.5d0,2.9673d0/,
     #     ttmref/'z. phys. d, 21, 251 (1991)'/
c  modified tang-toennies data 1997 (retardation left out)
c       !!!!!partially junk data!!!!! currently not used
c     data ahettm/6.42140197d0/,b1mhe,b2mhe,b3mhe/1.83899568d0,
c    #     -0.08401d0,0.003761d0/,xcut,rettm/4.5d0,2.9673d0/,
      data ccm6,ccm8,ccm10,ccm12,ccm14,ccm16/1.4609778,14.117855,
     #     183.69125,3265.,76440.,2275000./,
     #     ttmrf2/'j. chem. phys. 107, 914 (1997)'/
c   modified tang-toennies data, cybulski 1999, all in a.u.
      data cm6,cm8,cm10,cm12,cm14,cm16,attm,beta1,beta2,bdamp
     &             /1.46098d0,14.1179d0,183.691d0,
     &              3265.27d0,7.64399d4,2.27472d6,
     &              6.62002d0,1.88553d0,0.0639819d0,1.85822d0/,
     #     ttkref/'j. chem. phys. 111, 10520 (1999)(av5z+)'/
c  morse potential
      data dem,rem,am/7.61d0,2.963d0,2.126d0/,
     #     mref/'j. chem. phys. 94, 2817 (1991)'/
c
C-NH      if(iabs(ipot).eq.0) then
C-NH        call ljset(eps,re,text,ljref,irunit,ivunit,iout)
C-NH      else if(iabs(ipot).eq.1) then
        call hfdset(epsb,rrmb,c6b,c8b,c10b,zero,zero,zero,
     #              astarb,alsb,betsb,
     #              gammab,db,text,'hfd-b',hfdref,irunit,ivunit,iout)
C-NH      else if(iabs(ipot).eq.2) then
C-NH        call hfdset(epsc,rrmc,c6c,c8c,c10c,zero,zero,zero,
C-NH     #              astarc,alsc,betsc,
C-NH     #              gammac,dc,text,'hfd-c',hfdcrf,irunit,ivunit,iout)
C-NH      else if(iabs(ipot).eq.3) then
C-NH        call ttset(cc6,cc8,cc10,ahett,betahe,text,ttref,
C-NH     #             irunit,ivunit,iout)
C-NH      else if(iabs(ipot).eq.4) then
C-NH        call ttmset(cc6,cc8,cc10,ahettm,b1mhe,b2mhe,b3mhe,xcut,rettm,
C-NH     #              text,ttmref,irunit,ivunit,iout)
C-NH      else if(iabs(ipot).eq.5) then
C-NH        call morset(dem,rem,am,text,mref,irunit,ivunit,iout)
C-NH      else if(iabs(ipot).eq.6) then
C-NH        if(iout.gt.0) write(iout,'(/a)')
C-NH     #   '  tang-toennies-yiu potential for helium; prl 74, 1546 (1995)'
C-NH      else if(iabs(ipot).eq.8) then
C-NH        if(iout.gt.0) write(iout,'(/a)')
C-NH     #   '  lm2m2 potential for helium; jcp 94, 8047 (1991)'
C-NH      else if(iabs(ipot).eq.16) then
C-NH        call ttkset(cm6,cm8,cm10,cm12,cm14,cm16,attm,beta1,beta2,bdamp,
C-NH     &              text,ttmrf2,irunit,ivunit,iout)
C-NH      else if(iabs(ipot).eq.19) then
C-NH        if(iout.gt.0) write(iout,'(/a)')
C-NH     #   '  jeziorska potential for helium; jcp 127, 124303 (2007)'
c
c      some manipulated potentials
c      
C-NH      else if(iabs(ipot).eq.100) then
C-NH        xeps=0.5d0
C-NH        xre=3.8d0
C-NH        call ljset(xeps,xre,text,ljref,irunit,ivunit,iout)
C-NH      else if(iabs(ipot).eq.200) then
C-NH        xeps=0.5d0
C-NH        xre=5.0d0
C-NH        call ljset(xeps,xre,text,ljref,irunit,ivunit,iout)
C-NH      else if(iabs(ipot).eq.300) then
C-NH        xeps=0.8d0
C-NH        xre=3.9d0
C-NH        call ljset(xeps,xre,text,ljref,irunit,ivunit,iout)
C-NH      else if(iabs(ipot).eq.105) then
C-NH        xde=0.5d0
C-NH        xre=3.8d0
C-NH        xa=1.d0
C-NH        call morset(xde,xre,xa,text,mref,irunit,ivunit,iout)
C-NH      else
C-NH        call errprt(iout,'hepair','undefined potential',1)
C-NH      end if
      return
      end
c
c----------------- common block initialization routines ----------------
c-----------------------------------------------------------------------
c
      subroutine hfdset(eps,rrm,cc6,cc8,cc10,cc12,cc14,cc16,
     #                  astar,alstar,bstar,
     #                  gam,dd,name,type,ref,irunit,ivunit,iout)
c
c      converts input data to the desired units and initializes 
c      common block /hfddat/ for hfd potential routines 
c      subroutines called : r2aa,e2cm,strlen         m. lewerenz feb/93
c
      include 'real.inc'
      character name*(*),type*(*),ref*(*)
      common /hfddat/ c(6),a,rm,d,alpha,beta,gamma
C     -NH
      common /units/bohr2Angs,ehtocm
c
c -NH
c      call r2aa(0,a0tor,iout)
c      call e2cm(0,ehtocm,iout)
c      call runit(irunit,rscale,iout)
c     call vunit(ivunit,vscale,iout)
      a0tor = 1.d0
      rscale = 1.d0
      vscale = 1.d0
c
      c(1)=cc6*vscale*rscale**6
      c(2)=cc8*vscale*rscale**8
      c(3)=cc10*vscale*rscale**10
      c(4)=cc12*vscale*rscale**12
      c(5)=cc14*vscale*rscale**14
      c(6)=cc16*vscale*rscale**16
      a=astar*eps/ehtocm*vscale
      rm=rrm/a0tor*rscale
      d=dd
      alpha=alstar/rm
      beta=bstar/rm**2
      gamma=gam
      if(iout.gt.0) then
        call strlen(name,lname,iout)
        call strlen(type,ltype,iout)
        call strlen(ref,lref,iout)
        write(iout,'(/2x,5a)') name(1:lname),
     #       ' pair potential initialized for ',type(1:ltype),', ',
     #       ref(1:lref)
      end if
      return
      end
c
c
c-----------------------------------------------------------------------
c----------------------- backup of common blocks -----------------------
c-----------------------------------------------------------------------
c
      subroutine ptsave(ipot,io,save,nsave,iout)
c
c      saves contents of pair potential common blocks in a vector or
c      copies vector contents into a common block for pair potential
c      routines
c
c      ipot  : potential identifier. 0 -> lennard jones (6/12),
c              1 -> hfd potential, 2 -> hfd-c
c              3 -> tang-toennies, 4 -> modified tt, 5 -> morse
c      io    : read/write control: io > 0 : common block -> save vector,
c              io < 0 : save vector -> common block, else -> no action;
c              input
c      save  : vector of length nsave for data; input/output
c      nsave : length of vector save, must be > 20; input
c      iout  : unit number for messages, silent for iout.le.0
c      subroutines called : errprt                   m. lewerenz feb/93
c
      include 'real.inc'
      parameter (maxc=16,maxpar=200,maxspl=200,maxchb=128)
      dimension save(nsave)
      common /mrdata/ remors,amors,demors
      common /ljdata/ epslj,rlj
      common /hfddat/ c(6),a,rm,d,alpha,beta,gamma
      common /ttdata/ ctt(maxc),att,beta1,beta2,beta3,rcut
      common /coul00/ rep
      common /prdata/ params(maxpar)
      common /prspl/ pspl(maxspl)
      common /chbdat/ chrmin,ipower,nc,chbc(maxchb)
c
c......lennard-jones potential data, 0 => 6/12, 9 => 4/8
c
      if(ipot.eq.0.or.ipot.eq.9.or.ipot.eq.100.
     &             or.ipot.eq.200.or.ipot.eq.300) then
        minsav=2
        if(nsave.lt.minsav.and.iout.gt.0) then
          write(iout,'(/a,i3,a/)')
     #     '  *** error, vector is shorter than ',minsav,
     #     ' in ptsave, return without action *** '
        else
          if(io.gt.0) then
            save(1)=epslj
            save(2)=rlj
          else if(io.lt.0) then
            epslj=save(1)
            rlj=save(2)
          end if
        end if
c
c......hfd-b and hfd-c potential data
c
      else if(ipot.eq.1.or.ipot.eq.2) then
        minsav=12
        if(nsave.lt.minsav.and.iout.gt.0) then
          write(iout,'(/a,i3,a/)')
     #     '  *** error, vector is shorter than ',minsav,
     #     ' in ptsave, return without action *** '
        else
          if(io.gt.0) then
            save(1)=c(1)
            save(2)=c(2)
            save(3)=c(3)
            save(4)=c(4)
            save(5)=c(5)
            save(6)=c(6)
            save(7)=a
            save(8)=rm
            save(9)=d
            save(10)=alpha
            save(11)=beta
            save(12)=gamma
          else if(io.lt.0) then
            c(1)=save(1)
            c(2)=save(2)
            c(3)=save(3)
            c(4)=save(4)
            c(5)=save(5)
            c(6)=save(6)
            a=save(7)
            rm=save(8)
            d=save(9)
            alpha=save(10)
            beta=save(11)
            gamma=save(12)
          end if
        end if
c
c......tang-toennies and modified tang-toennies potential data
c
      else if(ipot.eq.3.or.ipot.eq.4.or.ipot.eq.11.or.ipot.eq.12
     &                                            .or.ipot.eq.15) then
        minsav=maxc+5
        if(nsave.lt.minsav.and.iout.gt.0) then
          write(iout,'(/a,i3,a/)')
     #     '  *** error, vector is shorter than ',minsav,
     #     ' in ptsave, return without action *** '
        else
          if(io.gt.0) then
            do 10 i=1,maxc
            save(i)=ctt(i)
  10        continue
            save(maxc+1)=att
            save(maxc+2)=beta1
            save(maxc+3)=beta2
            save(maxc+4)=beta3
            save(maxc+5)=rcut
          else if(io.lt.0) then
            do 20 i=1,maxc
            ctt(i)=save(i)
  20        continue
            att=save(maxc+1)
            beta1=save(maxc+2)
            beta2=save(maxc+3)
            beta3=save(maxc+4)
            rcut=save(maxc+5)
          end if
        end if
c
c......morse potential data
c
      else if(ipot.eq.5.or.ipot.eq.105) then
        minsav=3
        if(nsave.lt.minsav.and.iout.gt.0) then
          write(iout,'(/a,i3,a/)')
     #     '  *** error, vector is shorter than ',minsav,
     #     ' in ptsave, return without action *** '
        else
          if(io.gt.0) then
            save(1)=remors
            save(2)=amors
            save(3)=demors
          else if(io.lt.0) then
            remors=save(1)
            amors=save(2)
            demors=save(3)
          end if
        end if
c
c......tang-toennies-yiu potential
c
      else if(ipot.eq.6) then
        minsav=0
c
c......constant potential
c
      else if(ipot.eq.7) then
        minsav=1
        if(nsave.lt.minsav.and.iout.gt.0) then
          write(iout,'(/a,i3,a/)')
     #     '  *** error, vector is shorter than ',minsav,
     #     ' in ptsave, return without action *** '
        else
          if(io.gt.0) then
            call vcopy(params,save,minsav)
          else if(io.lt.0) then
            call vcopy(save,params,minsav)
          end if
        end if
c
c......special lm2m2 potential for he_2
c
      else if(ipot.eq.8) then
        minsav=0
c
c      ipot=9 is defined above with ipot=0 !!!
c
c......coulomb repulsion
c
      else if(ipot.eq.10) then
        minsav=1
        if(nsave.lt.minsav.and.iout.gt.0) then
          write(iout,'(/a,i3,a/)')
     #     '  *** error, vector is shorter than ',minsav,
     #     ' in ptsave, return without action *** '
        else
          if(io.gt.0) then
            save(1)=rep
          else if(io.lt.0) then
            rep=save(1)
          end if
        end if
c
c......general exponential/damped dispersion potential(13) and
c      joined morse/pade potential(14)
c      both use common block /prdata/
c
      else if(ipot.eq.13.or.ipot.eq.14.) then
        minsav=maxpar
        if(nsave.lt.minsav) call errprt(iout,'ptsave',
     &                       'vector buffer too short',0)
c         write(iout,'(/a,i3,a/)')
c    #     '  *** error, vector is shorter than ',minsav,
c    #     ' in ptsave, return without action *** '
c       else
          if(io.gt.0) then
            call vcopy(params,save,minsav)
          else if(io.lt.0) then
            call vcopy(save,params,minsav)
          end if
c       end if
c
c......generalised lennard-jones potential(18)
      else if(ipot.eq.18) then
        minsav=min(30,maxpar)
        if(nsave.lt.minsav) call errprt(iout,'ptsave',
     &                       'vector buffer too short',0)
        if(io.gt.0) then
          call vcopy(params,save,minsav)
        else if(io.lt.0) then
          call vcopy(save,params,minsav)
        end if
c
c......special vhe2je potential for he_2
c
      else if(ipot.eq.19) then
        minsav=0
c
c......chebyshev approximation in r**-k
c
      else if(ipot.eq.20) then
        if(io.gt.0) then
          ndata=nc+3
          if(ndata.le.nsave) then
            save(1)=chrmin
            save(2)=ipower
            save(3)=nc
            call vcopy(chbc,save(4),nc)
          else
            call errprt(iout,'ptsave',
     &           'vector buffer too small for chebyshev data',1)
          end if
        else if(io.lt.0) then
          chrmin=save(1)
          ipower=save(2)
          nc=save(3)
          if(nc.le.maxchb) then
            call vcopy(save(4),chbc,nc)
          else
            call errprt(iout,'ptsave',
     &           'common block too short for chebyshev data',1)
          end if
        end if
c
c......generalised exponential + inverse power potential, vxion
c
      else if(ipot.eq.21.or.ipot.eq.22) then
        minsav=maxc+2
        if(nsave.lt.minsav.and.iout.gt.0) then
          write(iout,'(/a,i3,a/)')
     #     '  *** error, vector is shorter than ',minsav,
     #     ' in ptsave, return without action *** '
        else
          if(io.gt.0) then
            do i=1,maxc
              save(i)=ctt(i)
            end do
            save(maxc+1)=att
            save(maxc+2)=beta1
          else if(io.lt.0) then
            do i=1,maxc
              ctt(i)=save(i)
            end do
            att=save(maxc+1)
            beta1=save(maxc+2)
          end if
        end if
c
c......general spline potential; uses common block /prspl/
      else if(ipot.eq.99.or.ipot.eq.98) then
        minsav=maxspl
        if(io.gt.0) then
          nsp=pspl(1)
          ndata=3*nsp+1
          if(ndata.le.nsave) then
            call vcopy(pspl,save,ndata)
          else
            call errprt(iout,'ptsave',
     &                       'vector buffer too short for spline',1)
          end if
        else if(io.lt.0) then
          nsp=save(1)
          ndata=3*nsp+1
          if(ndata.le.maxspl) then
            call vcopy(save,pspl,ndata)
          else
            call errprt(iout,'ptsave',
     &                       'common buffer too short for spline',1)
          end if
        end if
c
c......error exit
c
      else
        write(iout,*) ' ipot = ',ipot
        call errprt(iout,'ptsave','undefined potential type',-1)
      end if
      return
      end
c
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c---------------------- extracted from DMC code ------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine vgpair(rx,ry,rz,mpsi,npsi,v,natoms,ngroup,
     #                  potpar,mptpar,ipot,iout)
c
c      pair potential evaluation for groupwise organized particles
c      unrolled version of 05/dec/08
c
c      we need (3*natom+3)*8 bytes per walker.
c      optimal block length is cache-size/(3*natom+3)*8
c      2 MB/natom=500: mblock=128
c
c      rx,ry,rz : arrays mpsi*natom with cartesian coordinates for
c               npsi replicas and natom particles; input
c      mpsi   : column dimension of rx,ry,rz as declared in calling
c               routine; input
c      npsi   : actual number of replicas; input
c      v      : vector of length npsi returning potential for npsi
c               replicas; output
c      natoms : vector of length ngroup with numbers of identical atoms
c               in each group; input
c      ngroup : number of distinct particle groups; input
c      potpar : array mptpar*(ngroup*(ngroup+1)/2) with pair potential
c               parameters for all possible pair types; input
c      mptpar : column dimension of potpar; input
c      ipot   : vector of length ngroup*(ngroup+1)/2 with potential
c               type identifiers; input
c      iout   : unit number for messages, silent for iout.le.0; input
c      subroutines called : errprt, ivsum, vecset, ptsave, pairpt
c                                     m. lewerenz nov/94, OpenMP aug/07
c
      include 'real.inc'
cparallel processing:
c     parameter (irsq=1,mblock=32)
      parameter (irsq=1,mblock=128)
cvector:
c     parameter (irsq=1,mblock=511)
      dimension rx(mpsi,*),ry(mpsi,*),rz(mpsi,*),v(npsi),
     &          natoms(ngroup),potpar(mptpar,*),ipot(*)
      dimension rij(mblock*4),work(mblock*4)
c
      if(npsi.le.0.or.mpsi.le.0.or.npsi.gt.mpsi.or.ngroup.le.0) then
        call errprt(iout,'vgpair','illegal array dimension(s)',1)
      else
        call ivsum(natoms,ngroup,natom)
        call vecset(v,npsi,zero)
        nloop=(npsi+mblock-1)/mblock
c
c      loop over groups of identical particles
c
        index=0
        ih=0
        do 650 ii=1,ngroup
        il=ih+1
        ih=ih+natoms(ii)
        jh=il-1
        do 600 jj=ii,ngroup
        jl=jh+1
        jh=jh+natoms(jj)
        index=index+1
c
c      initialize potential data and loop over all atoms in the
c      i and j groups in strips of max. length mblock
c
          call ptsave(ipot(index),-1,potpar(1,index),mptpar,iout)
c$doacross local(kloop,nstart,nblock,i,j,n,dx,dy,dz,rsq,rij,work)
cc$&        ,MP_SCHEDTYPE=RUNTIME
!$OMP DO PRIVATE(kloop,nstart,nblock,i,j,n,dx,dy,dz,rsq,rij,work)
          do 550 kloop=1,nloop
          nstart=(kloop-1)*mblock
          nblock=min(npsi-nstart,mblock)
c
          do i=il,min(ih,natom)
            jmin=max(i+1,jl)
            jmax=min(jh,natom)
            jmod=jmax-jmin+1 
c--  turn off unrolling by commenting out the next line
            jmod=mod(jmax-jmin+1,4) 
c
c      first clean up loop over atoms
c
            do j=jmin,jmin+jmod-1
              do n=1,nblock
                dx=rx(n+nstart,i)-rx(n+nstart,j)
                rsq=dx*dx
                dy=ry(n+nstart,i)-ry(n+nstart,j)
                rsq=rsq+dy*dy
                dz=rz(n+nstart,i)-rz(n+nstart,j)
                rij(n)=rsq+dz*dz
              end do
cdir$ inline always pairpt
c -NH           call pairpt(ipot(index),rij,work,work,nblock,irsq,0,iout) 
              call hfdpot(rij,work,work,nblock,irsq,0,iout)
              do n=1,nblock
                v(nstart+n)=v(nstart+n)+work(n)
              end do
            end do

            do j=jmin+jmod,jmax,4
              do n=1,nblock
                dx0=rx(n+nstart,i)-rx(n+nstart,j)
                dx1=rx(n+nstart,i)-rx(n+nstart,j+1)
                dx2=rx(n+nstart,i)-rx(n+nstart,j+2)
                dx3=rx(n+nstart,i)-rx(n+nstart,j+3)
                rsq0=dx0*dx0
                rsq1=dx1*dx1
                rsq2=dx2*dx2
                rsq3=dx3*dx3
                dy0=ry(n+nstart,i)-ry(n+nstart,j)
                dy1=ry(n+nstart,i)-ry(n+nstart,j+1)
                dy2=ry(n+nstart,i)-ry(n+nstart,j+2)
                dy3=ry(n+nstart,i)-ry(n+nstart,j+3)
                rsq0=rsq0+dy0*dy0
                rsq1=rsq1+dy1*dy1
                rsq2=rsq2+dy2*dy2
                rsq3=rsq3+dy3*dy3
                dz0=rz(n+nstart,i)-rz(n+nstart,j)
                dz1=rz(n+nstart,i)-rz(n+nstart,j+1)
                dz2=rz(n+nstart,i)-rz(n+nstart,j+2)
                dz3=rz(n+nstart,i)-rz(n+nstart,j+3)
                rij(n         )=rsq0+dz0*dz0
                rij(n+  nblock)=rsq1+dz1*dz1
                rij(n+2*nblock)=rsq2+dz2*dz2
                rij(n+3*nblock)=rsq3+dz3*dz3
              end do
cdir$ inline always pairpt
C-NH              call pairpt(ipot(index),rij,work,work,4*nblock,irsq,0,
C-NH     &             iout)
              call hfdpot(rij,work,work,4*nblock,irsq,0,iout)
              
              do n=1,nblock
                v(nstart+n)=v(nstart+n)+work(n)+work(n+nblock)+
     &                         work(n+2*nblock)+work(n+3*nblock)
              end do
            end do
          end do
  550     continue
  600   continue
  650   continue
      end if
      return
      end
c-----------------------------------------------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine hfdpot(rsq,v,dvr,nr,irsq,ideriv,iout)
c
c      hfd (hartree fock dispersion) potential routine compatible with
c      the sum over pair routines for mc, md, and dmc calculations.
c      potential parameters are passed through common block /hfddat/
c
c            v(r)=a*exp(-alpha*r+beta*r**2)
c                 -switch(r)*(c6/r**6+c8/r**8+c10/r**10)
c      switch(r)=1, r > rm*d, switch(r)=exp(-(rm*d/r-1)**2) elsewhere
c      beta=0 is the original hfd-potential, otherwise it is hfd-b type
c
c      rsq   : vector of length nr with r or r**2 values; input
c      v     : pair potential at each rsq entry; output
c      dvr   : first derivative of pair potential wrt rsq; output
c              only computed for ideriv=1, may coincide with v if
c              ideriv=0
c      nr    : number of points= vector length of rsq,v,dvr; input
c      irsq  : 0 -> input in rsq is distance, else rsq = r**2; input
c      ideriv: 0 -> only v-evaluation
c              1 -> v + dvr evaluation; input
c      iout  : unit number for messages, silent for iout.le.0; input
c
c      ca. 43 ns per interaction on xeon5130/2ghz/ifort
c      ca. 103 ns on power6 4.7 ghz, xlf90 -O4
c      ca. 52 ns on power6 4.7 ghz, xlf90 -O4 -lmass
c      subroutines called: none                 m. lewerenz dec/92
c
      include 'real.inc'
      dimension rsq(nr),v(nr),dvr(nr)
      common /hfddat/ c(6),a,rm,d,alpha,beta,gamma
c
      rmd=rm*d
      if(ideriv.eq.0) then
        if(irsq.eq.0) then
          do i=1,nr
            if(rsq(i).ge.rmd) then
              switch=one
            else
              switch=exp(-(rmd/rsq(i)-one)**2)
            end if
            ri2=one/rsq(i)**2
            poly=switch*(((((c(6)*ri2+c(5))*ri2+c(4))*ri2+c(3))
     &                           *ri2+c(2))*ri2+c(1))*(ri2*ri2*ri2)
            v(i)=a*exp((beta*rsq(i)-alpha)*rsq(i))-poly
          end do
        else
          do i=1,nr
            rr=sqrt(rsq(i))
            if(rr.ge.rmd) then
              switch=one
            else
              switch=exp(-(rmd/rr-one)**2)
            end if
            ri2=one/rsq(i)
            poly=switch*(((((c(6)*ri2+c(5))*ri2+c(4))*ri2+c(3))
     &                           *ri2+c(2))*ri2+c(1))*(ri2*ri2*ri2)
            v(i)=a*exp((beta*rr-alpha)*rr)-poly
          end do  
        end if
      else
        call errprt(iout,'hfdpot','ideriv.ne.0 not implemented',1)
      end if
      return
      end
c
