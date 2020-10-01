c
c---------------------------- vch3ihen.f -------------------------------

c
      subroutine vch3ihen(x,y,z,natom,vso,iout)
c
c      intermolecular potential for CH3I-He composed of a sum over
c      CH3-He, I-He and He-He terms. CH3 has to have C_3v symmetry
c      and is currently restricted to be completely rigid. 
c
c      we assume that atom 1 is I, 2-4 are H, 5 is C, 6-natom are He;
c      geometry input in a0, energy output in Eh.
c
c      currently we use the HFD-B He-He potential of Aziz 1987.
c      other preprogrammed potentials can be invoked by setting
c      ihepot to other values.
c      the I-He_n part is fitted to ab initio CCSD(T) calculations
c      and uses the experimental spin-orbit splitting constant.
c      the CH3-He part is a 3D fit to CCSD(T) data for rigid CH3.
c
c      we use subroutines from the dmc package (vgpair, vch3hen, viheso)
c      for each of these contributions to the overall interaction.
c
c      x     : vector of length natom with x-coordinates of all atoms;
c              input
c      y     : vector of length natom with y-coordinates of all atoms;
c              input
c      z     : vector of length natom with z-coordinates of all atoms;
c              input
c      natom : total number of atoms; input
c      vso   : vector of length 6 with potential energy for all
c              spin-orbit states; output
c      iout  : fortran unit for error messages; input
c                                                  m. lewerenz feb/11
c
      include 'real.inc'
      parameter (clight=299792458.d0,planck=6.6260755d-34,
     &     emass=9.1093897d-31, bohr=0.529177249d-10)
      parameter (mhepar=20,ihepot=1,ngroup=2,mblock=128,itest=1) 
      parameter (ii=1,ih1=2,ih2=3,ih3=4,ic=5,ndop=5)
      dimension x(natom),y(natom),z(natom),vso(6)
      dimension heparm(mhepar),natoms(ngroup)
      allocatable xihe(:),yihe(:),zihe(:)
      common /hsodat/ ev1(mblock),ev2(mblock),ev3(mblock),
     &                ev4(mblock),ev5(mblock),ev6(mblock)
C     -NH
      common /units/bohr2Angs,ehtocm
      data icall/0/
      save
c
c      perform a few initialisation actions during the first call
c
      if(icall.eq.0) then
C-NH
         pi = dacos(-1.d0)
         ehtocm=0.01d0*planck/((4*pi**2)*clight*(emass*bohr**2))
         bohr2Angs = bohr*1.d10

        call hepair(ihepot,0,0,iout)
        call ptsave(ihepot,1,heparm,mhepar,iout)
        natoms(1)=1
        ev1(1)=zero
        ev2(1)=zero
        ev3(1)=zero
        ev4(1)=zero
        ev5(1)=zero
        ev6(1)=zero
        icall=1
      end if
c
c      total interaction between helium atoms (sum over pairs)
c
      nhe=natom-ndop
      vtotal=zero

      if (nhe.gt.1) then
      call vgpair(x(ndop+1),y(ndop+1),z(ndop+1),1,1,vhehe,nhe,1,
     &            heparm,mhepar,ihepot,iout)
      vtotal=vtotal+vhehe
      if(itest.ne.0) then
        write(iout,'(/a,e15.8)') '  vgpair output; vhehe: ',vhehe
      end if

      end if
c
c      total interaction between CH3 and helium atoms (sum over pairs)
c
      vmehe=zero
      call vch3hen(x(2),y(2),z(2),1,1,nhe+4,vmehe,iout)
      if(itest.ne.0) then
        write(iout,'(/a,e15.8)') '  vch3hen output; vmehe: ',vmehe
      end if
      vtotal=vtotal+vmehe
c
c      finally add the spin-orbit eigenvalues of the I-He_n part:
c      copy I and He coordinates into local vectors and call viheso.
c      lowest eigenvalue is returned through calling list.
c      complete eigenvalue set is in common block /hsodat/
c
      natoms(2)=nhe
      allocate (xihe(nhe+1),yihe(nhe+1),zihe(nhe+1))
      xihe(1)=x(ii)
      yihe(1)=y(ii)
      zihe(1)=z(ii)
      do i=1,nhe
        xihe(1+i)=x(ndop+i)
        yihe(1+i)=y(ndop+i)
        zihe(1+i)=z(ndop+i)
      end do
!      call viheso(xihe,yihe,zihe,1,1,vihen,natoms,ngroup,iout)
      deallocate (xihe,yihe,zihe)
      vso(1)=vtotal
      vso(2)=vtotal
      vso(3)=vtotal
      vso(4)=vtotal
      vso(5)=vtotal
      vso(6)=vtotal
      if(itest.ne.0) then
        write(iout,'(/a,e15.8)') '  viheso output; vihen: ',vihen
        do i=1,6
          write(iout,'(a,e15.8)') '  after viheso; vso : ',vso(i)
        end do
      end if
c
      return
      end
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c----------------------------- vch3hen ---------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine vch3hen(rx,ry,rz,mpsi,npsi,natom,v,iout)
c
c      CH3-He_n intercation potential evaluation compatible with dmc.
c      we assume that atom 1 is C, 2-4 are H, 5-natom are He.
c      geometry input in a0, energy output in Eh.
c
c      rx    : matrix with x-coordinates of all atoms in all walkers;
c              input
c      ry    : matrix with y-coordinates of all atoms in all walkers;
c              input
c      rz    : matrix with z-coordinates of all atoms in all walkers;
c              input
c      mpsi  : maximum number of walkers; input
c      npsi  : number of walkers; input
c      natom : total number of atoms; input
c      v     : vector of length npsi for potential energy; output
c      iout  : fortran unit for error messages; input
c                                                  m. lewerenz feb/11
c
      include 'real.inc'
      parameter (ih1=1,ih2=2,ih3=3,ic=4,ndop=4,itest=0)
      parameter (third=one/three,rflat=1.d-7)
      dimension rx(mpsi,natom),ry(mpsi,natom),rz(mpsi,natom),v(*)
      data icall/0/
      save

C     -NH
      common /units/bohr2Angs,ehtocm
c
c      perform a few initialisation actions during the first call
c
      if(icall.eq.0) then
c-NH     call e2eh(1,escale,iout)
         escale = 1.d0/ehtocm
c-NH     call r2aa(0,rscale,iout)
         rscale = bohr2Angs
c     cmass=atmass(12,6,iout)
        cmass = 12.d0
c     hmass=atmass(1,1,iout)
        hmass = 1.007825035d0
        tmassi=one/(cmass+three*hmass)
        icall=1
      end if
c
      do j=1,npsi
        xc=rx(j,ic)
        xh1=rx(j,ih1)
        xh2=rx(j,ih2)
        xh3=rx(j,ih3)
        yc=ry(j,ic)
        yh1=ry(j,ih1)
        yh2=ry(j,ih2)
        yh3=ry(j,ih3)
        zc=rz(j,ic)
        zh1=rz(j,ih1)
        zh2=rz(j,ih2)
        zh3=rz(j,ih3)
c
c      summation over individual CH3-He contributions.
c      first we need the CH3 center of mass coordinates.
c
        xcom=(cmass*xc+hmass*(xh1+xh2+xh3))*tmassi
        ycom=(cmass*yc+hmass*(yh1+yh2+yh3))*tmassi
        zcom=(cmass*zc+hmass*(zh1+zh2+zh3))*tmassi
c
c      we also need a unit vector which falls into the symmetry axis
c      of CH3 to define the polar coordinates of He relative to CH3.
c      we first try the vector pointing from the H center of mass
c      to C or, if CH3 is planar, we take the vector product of 
c      two CH bond vectors.
c      both versions assume that at least C3v symmetry is present,
c      which underlies the construction of the CH3-He potential surface.
c      we also compute the vector product between the symmetry axis
c      vector and the C->H1 vector for the phi definition below. 
c
        xhsum=(xh1+xh2+xh3)*third
        yhsum=(yh1+yh2+yh3)*third
        zhsum=(zh1+zh2+zh3)*third
        dx=xc-xhsum
        dy=yc-yhsum
        dz=zc-zhsum
        rc=sqrt(dx**2+dy**2+dz**2)
        if(itest.ne.0) write(iout,'(a,f8.5)') '  rc : ',rc
        dx1=xh1-xc
        dy1=yh1-yc
        dz1=zh1-zc
        rchi=one/sqrt(dx1**2+dy1**2+dz1**2)
        if(rc.lt.rflat) then
          if(itest.ne.0) write(iout,'(a)') '  ch3 is planar '
          dx2=xh2-xc
          dy2=yh2-yc
          dz2=zh2-zc
          dx=dy1*dz2-dz1*dy2
          dy=dz1*dx2-dx1*dz2
          dz=dx1*dy2-dy1*dx2
          rc=sqrt(dx**2+dy**2+dz**2)
        end if
        rci=one/rc
c
c      unit reference vectors
c
        xref=dx*rci
        yref=dy*rci
        zref=dz*rci
        xphi=(yref*dz1-zref*dy1)
        yphi=(zref*dx1-xref*dz1)
        zphi=(xref*dy1-yref*dx1)
        rphi=one/sqrt(xphi**2+yphi**2+zphi**2)
        xphi=xphi*rphi
        yphi=yphi*rphi
        zphi=zphi*rphi
        if(itest.ne.0) then
          write(iout,'(a,3f8.5)') '  C_3v axis  : ',xref,yref,zref
          write(iout,'(a,3f8.5)') '  phi vector : ',xphi,yphi,zphi
        end if
c
c      for each helium atom we compute the polar distance,
c      the dot product with the symmetry axis (costh) and
c      the angle between the planes spanned by the symmetry
c      axis and C->H1 and by the symmetry axis and the 
c      com->He vector. there is no need to figure out if
c      phi is above or below pi due to the symmetry of the
c      potential surface V(phi)=V(2*pi-phi).
c      if He is on the C_3v axis (theta=0) the He reference vector 
c      has a vanishing norm but phi is irrelevant and set to zero
c
        vmehe=zero
        do i=ndop+1,natom
          dxhe=rx(j,i)-xcom
          dyhe=ry(j,i)-ycom
          dzhe=rz(j,i)-zcom
          rhe=sqrt(dxhe**2+dyhe**2+dzhe**2)
          costh=(xref*dxhe+yref*dyhe+zref*dzhe)/rhe
          xhe=(yref*dzhe-zref*dyhe)
          yhe=(zref*dxhe-xref*dzhe)
          zhe=(xref*dyhe-yref*dxhe)
          rnorm=sqrt(xhe**2+yhe**2+zhe**2)
          if(rnorm.ne.zero) then
            cosphi=(xhe*xphi+yhe*yphi+zhe*zphi)/rnorm
            phi=acos(cosphi)
          else
            phi=zero
          end if
          rjac=rhe*rscale
          call vch3he(rjac,costh,phi,vmhe,1,iout)
          vmhe=vmhe*escale
          if(itest.ne.0) then
            write(iout,'(a,3f8.5)') ' H1 : ',xh1,yh1,zh1
            write(iout,'(a,3f8.5)') ' H2 : ',xh2,yh2,zh2
            write(iout,'(a,3f8.5)') ' H3 : ',xh3,yh3,zh3
            write(iout,'(a,3f8.5)') '  C : ',xc,yc,zc
            write(iout,'(a,3f8.5)') ' He : ',rx(j,i),ry(j,i),rz(j,i)
            write(iout,'(a,3f8.5)') ' (He,z)  vector : ',xhe,yhe,zhe
            write(iout,'(a,3f8.5)')
     &                '  rjac, costh, phi : ',rjac,costh,phi
            write(iout,'(a,i3,a,e15.8)')
     &                '  vch3he output; vmhe for He',i-ndop,' : ',vmhe
          end if
          vmehe=vmehe+vmhe
        end do
        v(j)=vmehe
      end do
c
      return
      end
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c----------------------------- iheso.f ---------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
c    potential for I in he_n, particle 1 is I, all others He.
c    cleaned up version of iheso_ji.f; m. lewerenz 03/feb/11
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c------------------------------ viheso ---------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine viheso(rx,ry,rz,mpsi,npsi,v,natoms,ngroup,iout)
c
c      interaction potential for I-He_n including spin-orbit
c      parameters for ^2\Sigma and ^2\Pi states fitted to ext.T.T series.
c      DIM style matrix for mixing of sigma and pi states and additional
c      spin-orbit coupling matrix.
c      v0=(1/3)*(2*v_pi+v_sigma), v2=(5/3)*(v_sigma-v_pi)
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
c      iout   : unit number for messages, silent for iout.le.0; input
c      subroutines called : errprt, ivsum, vecset, vcopy, vecsum
c                           vxion, vind, pbpots, hso
c
c                   m. lewerenz may/09 (original Ar+), J. Jiang jul/10
c                                     cleaned up by m. lewerenz feb/11
c
      include 'real.inc'
      parameter (itype=1,irsq=0,maxc=16,maxatm=400)
c----------------
c   explicit strip mining for cache optimisation or vector processors
cparallel processing:
      parameter (mblock=128)
cvector:
c     parameter (mblock=511)
c----------------
c      I spin orbit splitting:
      parameter (delta=-0.03464250d0)
c----------------
      dimension rx(mpsi,*),ry(mpsi,*),rz(mpsi,*),v(npsi),natoms(ngroup)
      dimension dxn(2*mblock),dyn(2*mblock),dzn(2*mblock),rij(2*mblock),
     &          work(2*mblock),vsigma(2*mblock),vpi(2*mblock)
      dimension dipx(mblock,maxatm),dipy(mblock,maxatm),
     &          dipz(mblock,maxatm)
      dimension csigma(maxc),cpi(maxc)
      common /ttdata/ cn(maxc),ax,bx,cx,dum2,dum3
      common /hsodat/ hr11(mblock),hr22(mblock),hr12(mblock),
     &                hr13(mblock),hi12(mblock),hi13(mblock)
      data icall /0/
      save
c
      if(npsi.le.0.or.mpsi.le.0.or.npsi.gt.mpsi.or.ngroup.le.0) then
        call errprt(iout,'viheso','illegal array dimension(s)',1)
      else
        if(icall.eq.0) then
          if(natoms(1).ne.1) 
     &      call errprt(iout,'viheso','too many atoms in group 1',-1)
          call ivsum(natoms,ngroup,natom)
          if(natom.gt.maxatm)
     &            call errprt(iout,'viheso','too many atoms',-1)
          call ihepot(itype,asigma,bsigma,cexps,csigma,
     &                      api,bpi,cexppi,cpi,maxc,iout)
          root2=sqrt(two)
          aso=delta/three
          write(iout,'(a,f8.6)') '  Atomic SO-splitting (E_h): ',delta
          icall=1
        end if
c
        call ivsum(natoms,ngroup,natom)
        call vecset(v,npsi,zero)
        nloop=(npsi+mblock-1)/mblock
c----top of block loop
        do kloop=1,nloop
          nstart=(kloop-1)*mblock
          nblock=min(npsi-nstart,mblock)
          call vecset(hr11,nblock,zero)
          call vecset(hr12,nblock,zero)
          call vecset(hr13,nblock,zero)
          call vecset(hr22,nblock,zero)
          call vecset(hi12,nblock,zero)
          call vecset(hi13,nblock,zero)
c
c      first compute i-he distances and induced he-dipoles
c
          jleft=natom-1
c--- comment out the next line to turn off unrolling
          jleft=mod(jleft,2)
          do j=2,1+jleft
            do n=1,nblock
              dx=rx(n+nstart,j)-rx(n+nstart,1)
              rsq=dx*dx
              dy=ry(n+nstart,j)-ry(n+nstart,1)
              rsq=rsq+dy*dy
              dz=rz(n+nstart,j)-rz(n+nstart,1)
              rsq=rsq+dz*dz
              rij(n)=sqrt(rsq)
              dxn(n)=dx
              dyn(n)=dy
              dzn(n)=dz
            end do
c
c     get the sigma and pi energies
c
            call vcopy(csigma,cn,maxc)
            ax=asigma
            bx=bsigma
            cx=cexps
            call ttxpot(rij,vsigma,work,nblock,irsq,0,iout)
            call vcopy(cpi,cn,maxc)
            ax=api
            bx=bpi
            cx=cexppi
            call ttxpot(rij,vpi,work,nblock,irsq,0,iout)
c
c      update the complex coupling matrix in p_-1, p_0, p_1 basis
c
            do n=1,nblock
              v0=(one/three)*(two*vpi(n)+vsigma(n))
              v2r=(half/(three*rij(n)**2))*(vsigma(n)-vpi(n))
              v2rz=v2r*three*root2*dzn(n)
              v2rzz=v2r*(rij(n)**2-three*dzn(n)**2)
              hr11(n)=hr11(n)+v0+v2rzz
              hr12(n)=hr12(n)+v2rz*dxn(n)
              hi12(n)=hi12(n)+v2rz*dyn(n)
              hr13(n)=hr13(n)+v2r*three*(dyn(n)**2-dxn(n)**2)
              hi13(n)=hi13(n)-v2r*three*two*dxn(n)*dyn(n)
              hr22(n)=hr22(n)+v0-two*v2rzz
            end do
          end do
c
c      unrolled version for speed
c
          do j=2+jleft,natom,2
            do n=1,nblock
              dx0=rx(n+nstart,j  )-rx(n+nstart,1)
              dx1=rx(n+nstart,j+1)-rx(n+nstart,1)
              rsq0=dx0*dx0
              rsq1=dx1*dx1
              dy0=ry(n+nstart,j  )-ry(n+nstart,1)
              dy1=ry(n+nstart,j+1)-ry(n+nstart,1)
              rsq0=rsq0+dy0*dy0
              rsq1=rsq1+dy1*dy1
              dz0=rz(n+nstart,j  )-rz(n+nstart,1)
              dz1=rz(n+nstart,j+1)-rz(n+nstart,1)
              rsq0=rsq0+dz0*dz0
              rsq1=rsq1+dz1*dz1
              rij0=sqrt(rsq0)
              rij1=sqrt(rsq1)
              rij(n)=rij0
              dxn(n)=dx0
              dyn(n)=dy0
              dzn(n)=dz0
              rij(n+nblock)=rij1
              dxn(n+nblock)=dx1
              dyn(n+nblock)=dy1
              dzn(n+nblock)=dz1
            end do
c
c     get the sigma and pi energies
c
            call vcopy(csigma,cn,maxc)
            ax=asigma
            bx=bsigma
            cx=cexps
            call ttxpot(rij,vsigma,work,2*nblock,irsq,0,iout)
            call vcopy(cpi,cn,maxc)
            ax=api
            bx=bpi
            cx=cexppi
            call ttxpot(rij,vpi,work,2*nblock,irsq,0,iout) 
c
c      update the complex coupling matrix in p_-1, p_0, p_1 basis
c
            do n=1,nblock
              nn=n+nblock
              v0=(one/three)*(two*vpi(n)+vsigma(n))
              v2r=(half/(three*rij(n)**2))*(vsigma(n)-vpi(n))
              v2rz=v2r*three*root2*dzn(n)
              v2rzz=v2r*(rij(n)**2-three*dzn(n)**2)
              hr11(n)=hr11(n)+v0+v2rzz
              hr12(n)=hr12(n)+v2rz*dxn(n)
              hi12(n)=hi12(n)+v2rz*dyn(n)
              hr13(n)=hr13(n)+v2r*three*(dyn(n)**2-dxn(n)**2)
              hi13(n)=hi13(n)-v2r*three*two*dxn(n)*dyn(n)
              hr22(n)=hr22(n)+v0-two*v2rzz
              v0=(one/three)*(two*vpi(nn)+vsigma(nn))
              v2r=(half/(three*rij(nn)**2))*(vsigma(nn)-vpi(nn))
              v2rz=v2r*three*root2*dzn(nn)
              v2rzz=v2r*(rij(nn)**2-three*dzn(nn)**2)
              hr11(n)=hr11(n)+v0+v2rzz
              hr12(n)=hr12(n)+v2rz*dxn(nn)
              hi12(n)=hi12(n)+v2rz*dyn(nn)
              hr13(n)=hr13(n)+v2r*three*(dyn(nn)**2-dxn(nn)**2)
              hi13(n)=hi13(n)-v2r*three*two*dxn(nn)*dyn(nn)
              hr22(n)=hr22(n)+v0-two*v2rzz
            end do
          end do
c
c      build the spin-orbit matrices and diagonalize,
c      then optionally add the induced dipole interactions
c
          call hso(hr11,hr22,hr12,hr13,hi12,hi13,v(nstart+1),nblock,
     &             aso,iout)
c----end of block loop
        end do
      end if
      return
      end
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c-------------------------------- hso ----------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine hso(hr11,hr22,hr12,hr13,hi12,hi13,v,nv,aso,iout)
c
c      builds complex 6x6 matrices and finds their lowest eigenvalue.
c      uses lapack subroutines zhetrd and zsteqr.
c      complete set of eigenvalues is returned in the six input vectors
c      for special purposes.
c                                                    m. lewerenz feb/09
c
      include 'real.inc'
      dimension hr11(nv),hr12(nv),hr13(nv),hr22(nv),hi12(nv),hi13(nv),
     &          v(nv)
      dimension hsor(6,6),hsoi(6,6),vecr(6,6),veci(6,6),
     &          vr(6),vi(6),aux(6),int(6)
      complex*16 zhso(6,6),ztau(6),zwork(6*20)
c
c      expand to full matrix, add spin-orbit coupling and diagonalize
c
      root2=sqrt(two)
      do n=1,nv
        do k=1,6
          do l=k,6
            hsor(k,l)=zero
            hsoi(k,l)=zero
          end do
        end do
        hsor(1,1)=hr11(n)-aso
        hsor(4,4)=hr11(n)+aso
        hsor(1,2)=hr12(n)
        hsor(4,5)=hr12(n)
        hsor(1,3)=hr13(n)
        hsor(4,6)=hr13(n)
        hsor(2,2)=hr22(n)
        hsor(5,5)=hr22(n)
        hsor(2,3)=-hr12(n)
        hsor(5,6)=-hr12(n)
        hsor(3,3)=hr11(n)+aso
        hsor(6,6)=hr11(n)-aso
        hsoi(1,2)=hi12(n)
        hsoi(4,5)=hi12(n)
        hsoi(1,3)=hi13(n)
        hsoi(4,6)=hi13(n)
        hsoi(2,3)=-hi12(n)
        hsoi(5,6)=-hi12(n)
        hsor(1,5)=aso*root2
        hsor(2,6)=aso*root2
        do k=1,6
          do l=k+1,6
            hsor(l,k)=hsor(k,l)
            hsoi(l,k)=-hsoi(k,l)
          end do
        end do
        do k=1,6
          do l=1,6
            zhso(k,l)=cmplx(hsor(k,l),hsoi(k,l))
          end do
        end do
c
c  lapack/mkl etc.,
        call zhetrd('u',6,zhso,6,vr,vi,ztau,zwork,6*20,info)
        call zsteqr('n',6,vr,vi,zhso,6,zwork,info)
c       write(iout,'(a/6e14.6)')
c    &             ' eigenvalues: ',(vr(k),k=1,6)
c
c    eiscg1 seems to give correct eigenvalues, but unordered
c       call eiscg1(6,6,hsor,hsoi,vr,vi,vecr,veci,ierr,aux,int)
c       write(iout,'(a/a,6e14.6/a,6e14.6)')
c    &    ' eigenvalues: ',' re: ',(vr(k),k=1,6),' im: ',(vi(k),k=1,6)
c
c      establish increasing order to be on the safe side
c
c       call indexx(6,vr,int,iout)
c       write(iout,*) (int(k),k=1,6)
c       call incseq(vr,int,aux,6,iout)
        v(n)=vr(1)-aso
        hr11(n)=vr(1)-aso
        hr22(n)=vr(2)-aso
        hr12(n)=vr(3)-aso
        hr13(n)=vr(4)-aso
        hi12(n)=vr(5)-aso
        hi13(n)=vr(6)-aso
      end do
      return
      end
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c------------------------------ ihepot ---------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine ihepot(itype,asigma,bsigma,cexps,csigma,
     &                        api,bpi,cexppi,cpi,nc,iout)
c
c      initializes parameters for extended tang-toennies model
c      for ihe interactions fitted to ab initio calculations by
c      jiang ji mar/10.
c
c      itype=1: extrapolated CCSD(T) data
c                               m. lewerenz oct/09, j. jiang jul/10
c
      include 'real.inc'
      dimension csigma(nc),cpi(nc)
c
      call vecset(csigma,nc,zero)
      call vecset(cpi,nc,zero)
      if(itype.eq.1) then
        asigma   = 0.2911528d+02
        bsigma   = 0.1569085d+01
        cexps    = 0.1963841d-01
        csigma(4)= 0.0000000d+00
        csigma(6)= 0.1730855d+02
        csigma(8)= 0.8515043d+03
        api      = 0.2787605d+02
        bpi      = 0.1407436d+01
        cexppi   = 0.2519996d-01
        cpi(4)   = 0.0000000d+00
        cpi(6)   = 0.1730855d+02
        cpi(8)   = 0.1143582d+04
        write(iout,'(/2a)') '  I-He Sigma/Pi + SO,',
     &    ' extrapolated Q56 CCSD(T), ext. tang-toennies fit, mar/09 '
        write(iout,'(a/1x,6f10.5)')
     &    '  Sigma state: A, b, c, c4, c6, c8: ',
     &     asigma,bsigma,cexps,csigma(4),csigma(6),csigma(8)
        write(iout,'(a/1x,6f10.5)')
     &    '  Pi state: A, b, c, c4, c6, c8: ',
     &     api,bpi,cexppi,cpi(4),cpi(6),cpi(8)
      else
        call errprt(iout,'ihepot','undefined potential type',-1)
      end if
      return
      end
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c-----------------------------------------------------------------------
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine ttxpot(rsq,v,dvr,nr,irsq,ideriv,iout)
c
c      extended tang-toennies potential compatible with the sum over
c      pair routines for mc, md, and dmc calculations.
c
c      v(r)=a*exp(-beta(r)*r)-sum_k f_k c_k/r**k,  k=kmin,kmax
c      beta(r)=beta1+beta2*r ; f_k= tang-toennies damping function
c                                   with b(r)=beta1+2*beta2*r
c
c      potential parameters are passed through common block /ttdata/
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
c      subroutines called: ttdamp, vscale, errprt    m. lewerenz nov/97
c
      include 'real.inc'
      parameter (mr=511,mc=16)
      dimension rsq(nr),v(nr),dvr(nr)
      dimension rlocal(mr),rinv(mr),vbeta(mr),damp(mr,mc),rexp(mr)
      common /ttdata/ cn(mc),a,beta1,beta2,beta3,rcut
c
c      first determine kmin and kmax by checking cn-coefficients
c
      kmax=mc
      do 50 i=mc,1,-1
      if(cn(i).ne.zero) goto 55
      kmax=kmax-1
   50 continue
   55 if(kmax.eq.0) then
        call errprt(iout,'ttxpot','dispersion coefficients are zero',0)
      end if
      kmin=1
      do 60 i=1,mc
      if(cn(i).ne.zero) goto 65
      kmin=kmin+1
   60 continue
   65 continue
c
c      block loop over local buffer lengths
c
      if(ideriv.eq.0) then
        nstart=0
        nleft=nr
  120   nblock=min(mr,nleft)
c
c      set local beta, r, and r**-1 vectors
c
        if(irsq.eq.0) then
          do 150 i=1,nblock
          rlocal(i)=rsq(nstart+i)
          rinv(i)=one/rlocal(i)
          vbeta(i)=beta1+(two*beta2)*rlocal(i)
  150     continue
        else
          do 160 i=1,nblock
          rlocal(i)=sqrt(rsq(nstart+i))
          rinv(i)=one/rlocal(i)
          vbeta(i)=beta1+(two*beta2)*rlocal(i)
  160     continue
        end if
c
c      calculate damping functions, then polynomial in r**-n
c      skipping zero coefficients. finally add the exponential term
c
        call ttdamp(rlocal,vbeta,rexp,damp,damp,mr,nblock,
     #              kmin,kmax,0,ideriv,iout)
        call vscale(damp(1,kmax),v(nstart+1),nblock,cn(kmax))
        do 400 n=kmax-1,kmin,-1
        if(cn(n).eq.zero) then
          do 300 i=1,nblock
          v(nstart+i)=v(nstart+i)*rinv(i)
  300     continue
        else
          do 350 i=1,nblock
          v(nstart+i)=v(nstart+i)*rinv(i)+cn(n)*damp(i,n)
  350     continue
        end if
  400   continue
        do 500 i=1,nblock
        beta=beta1+beta2*rlocal(i)
        v(nstart+i)=a*exp(-beta*rlocal(i))-v(nstart+i)*rinv(i)**kmin
  500   continue
c
        nstart=nstart+nblock
        nleft=nleft-nblock
        if(nleft.gt.0) goto 120
      else
        call errprt(iout,'ttxpot','ideriv.ne.0 not implemented',1)
      end if
      return
      end
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine ttdamp(r,beta,work,damp,ddamp,mr,nr,kmin,kmax,
     #                  irsq,ideriv,iout)
c
c      tang-toennies damping function and optional first derivative.
c      each r may have a different beta value associated with it.
c      !!derivatives are correct only if beta is a constant!!
c
c        f(n;r)=1-exp(-beta(r)*r)*sum(k=0,n)(beta(r)*r)**k/k!
c
c      r      : vector of distances, length nr; input
c      beta   : vector of range scaling parameters, length nr; input
c      work   : work space of length nr, holding exp(-beta*r) on return
c      damp   : matrix mr*kmax for damping function values; output
c               f(n;r(i)) is stored in damp(i,n)
c      ddamp  : matrix mr*kmax for first derivative values; output
c               df(n;r)/dr is stored in ddamp(i,n)
c               only used if ideriv.ne.0
c      mr     : max. column size of damp and ddamp; input
c      nr     : actual column size of damp and ddamp and length of
c               vector r; input
c      kmin   : smallest n for which f(n;r) is desired; input
c      kmax   : largest n for which f(n;r) is desired; input
c      irsq   : currently not active
c      ideriv : 0 -> only function evaluation. in this case damp and
c               ddamp may coincide in memory, else -> also ddamp
c               is calculated; input
c      iout   : unit number for messages, silent for iout.le.0; input
c      subroutines called: errprt                 m. lewerenz nov/92
c
      include 'real.inc'
      dimension r(nr),beta(nr),work(nr),damp(mr,kmax),ddamp(mr,kmax)
c
      if(nr.le.0.or.mr.le.0.or.nr.gt.mr.
     #           or.kmin.gt.kmax.or.kmax.lt.0) then
        call errprt(iout,'ttdamp','illegal parameter(s)',1)
      else
        if(irsq.eq.0) then
c
c      special case kmax=1 is dealt with first
c
          if(kmax.eq.1) then
            if(ideriv.eq.0) then
              do i=1,nr
                br=beta(i)*r(i)
                work(i)=exp(-br)
                damp(i,1)=one-work(i)*(br+one)
              end do
            else
              do i=1,nr
                br=beta(i)*r(i)
                work(i)=exp(-br)
                damp(i,1)=one-(work(i)*br)-work(i)
                ddamp(i,1)=beta(i)*(work(i)*br)
              end do
            end if
          else
c
c      generate all sum(k)(beta*r)**k/k! using unrolled loops
c      damp(i,kmax) serves as work space
c
            if(mod(kmax,2).eq.0) then
              do i=1,nr
                work(i)=beta(i)*r(i)
                damp(i,kmax)=work(i)
                damp(i,1)=work(i)+one
              end do
              kstart=2
            else
              do i=1,nr
                work(i)=beta(i)*r(i)
                damp(i,kmax)=work(i)*work(i)*half
                damp(i,1)=work(i)+one
                damp(i,2)=damp(i,1)+damp(i,kmax)
              end do
              kstart=3
            end if
c
            do k=kstart,kmax-1,2
              rki=one/k
              rki2=one/(k+1)
              do i=1,nr
                temp=damp(i,kmax)*work(i)*rki
                damp(i,k)=damp(i,k-1)+temp
                damp(i,kmax)=temp*work(i)*rki2
                damp(i,k+1)=damp(i,k)+damp(i,kmax)
              end do
            end do
c
c      clean up loop and work=exp(beta*r)
c
            rki=one/kmax
            do i=1,nr
              damp(i,kmax)=damp(i,k-1)+damp(i,kmax)*work(i)*rki
              work(i)=exp(-work(i))
            end do
c
c      optional first derivative beta*exp(-beta*r)*(damp(n)-damp(n-1))
c      this part assumes constant beta!!!!!
c
            if(ideriv.ne.0) then
              do i=1,nr
                ddamp(i,kmax)=beta(i)*work(i)
                ddamp(i,1)=ddamp(i,kmax)*r(i)*beta(i)
              end do
              do k=max(kmin,2),kmax
                do i=1,nr
                  ddamp(i,k)=ddamp(i,kmax)*(damp(i,k)-damp(i,k-1))
                end do
              end do
            end if
c
c      multiply with exponential and subtract from 1, starting at kmin
c
            kstart=max(kmin,1)
            kover=mod(kmax-kstart+1,2)
            if(kover.eq.1) then
              do i=1,nr
                damp(i,kstart)=one-damp(i,kstart)*work(i)
              end do
            end if
            do k=kstart+kover,kmax,2
              do i=1,nr
                damp(i,k)=one-damp(i,k)*work(i)
                damp(i,k+1)=one-damp(i,k+1)*work(i)
              end do
            end do
          end if
        else
          call errprt(iout,'ttdamp','irsq.ne.0 not implemented',1)
        end if
      end if
      return
      end
c
c$---------------------------------------------------------------------
c
      subroutine vecset(v,nv,const)
c
c      set all elements of vector v of length nv to const
c      subroutines called : none                 m. lewerenz 17.9.89
c
      implicit real*8 (a-h,o-z)
      dimension v(nv)
c
      if(nv.gt.0) then
        do 10 i=1,nv
        v(i)=const
   10   continue
      end if
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine vcopy(v1,v2,nv)
c
c      copies vector v1 into v2.
c      subroutines called : none                           m. lewerenz 
c
      implicit real*8 (a-h,o-z)
      dimension v1(nv),v2(nv)
c
      if(nv.gt.0) then
        do 10 i=1,nv
        v2(i)=v1(i)
   10   continue
      end if
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine vecdot(v1,v2,vdot,n)
c
c      returns dot product of v1 and v2 in vdot. v1 and v2 may coincide
c      subroutines called: none                   m. lewerenz 23.9.89
c      unrolling allows overlapping of add and multiply on rs6000
c
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,mroll=4)
      dimension v1(n),v2(n)
      common /doctrl/ nroll
c
      s1=zero
      if(nroll.le.1) then
        do 10 i=1,n
        s1=s1+v1(i)*v2(i)
   10   continue
      else
        nover=mod(n,mroll)
        do 100 i=1,nover
        s1=s1+v1(i)*v2(i)
  100   continue
        s2=zero
        s3=zero
        s4=zero
        do 200 i=nover+1,n,mroll
        s1=s1+v1(i)*v2(i)
        s2=s2+v1(i+1)*v2(i+1)
        s3=s3+v1(i+2)*v2(i+2)
        s4=s4+v1(i+3)*v2(i+3)
  200   continue
        s1=s1+s2+s3+s4
      end if
      vdot=s1
      return
      end
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
      subroutine vscale(v1,v2,nv,scale)
c
c      multiplies all elements of vector v1 with scalar scale and
c      stores result in v2. v1 and v2 may coincide.
c      subroutines called : vecset, veccop        m. lewerenz   sep/89
c
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.d0,mroll=4)
      dimension v1(nv),v2(nv)
      common /doctrl/ nroll
c
      if(nv.gt.0) then
        if(scale.eq.zero) then
          call vecset(v2,nv,zero)
        else if(scale.eq.one) then
          call veccop(v1,v2,nv)
        else
          if(nroll.le.1) then
            do 10 i=1,nv
            v2(i)=v1(i)*scale
   10       continue
          else
            nover=mod(nv,mroll)
            do 20 i=1,nover
            v2(i)=v1(i)*scale
   20       continue
            do 30 i=nover+1,nv,mroll
            v2(i)=v1(i)*scale
            v2(i+1)=v1(i+1)*scale
            v2(i+2)=v1(i+2)*scale
            v2(i+3)=v1(i+3)*scale
   30       continue
          end if
        end if
      end if
      return
      end
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
      subroutine veccop(v1,v2,nv)
c
c      copies vector v1 into v2.
c      subroutines called : none               m. lewerenz 17/sep/89
c
      implicit real*8 (a-h,o-z)
      dimension v1(nv),v2(nv)
c
      if(nv.gt.0) then
        do 10 i=1,nv
        v2(i)=v1(i)
   10   continue
      end if
      return
      end
c
c---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
c
c---------------------------- last line --------------------------------
