c         SPO-FGM method for propagation in single PES
     
      parameter (npmax1=1024,npmax2=1024)
      parameter (npulsemax=5)
      IMPLICIT REAL*8(A-H,O-y)
      real*8 Tcin(0:npmax1-1,2)
      real*8 LC
      complex*16 pot
      complex*16 zsolap,zproj
      complex*16 im
      complex*16 zwf(0:npmax1-1,0:npmax2-1)
      complex*16 zwfin(0:npmax1-1,0:npmax2-1)
      complex*16 zwfout(0:npmax1-1,0:npmax2-1)
      complex*16 zwftarg(0:npmax1-1,0:npmax2-1)
      complex*16 zwfn(0:npmax1-1,0:npmax2-1)
      complex*16 expk(0:npmax1-1,2)
      complex*16 expV2(0:npmax1-1,0:npmax2-1)
      complex*16 expdip(0:npmax1-1,0:npmax2-1)
      common/FFTprep/trigx(0:20000,2)
      common/FFTfac/ifax(30,2),nfax(2)
      common/param/np(2),npp(2)
      common/potential/pot(0:npmax1-1,0:npmax2-1)
      common/rlim/dr(2)
      common/zgrid/z(0:npmax1-1)
      common/xgrid/x(0:npmax2-1)
      common/mass/xmu(2)
      common/dipole/dipole(0:npmax1-1,0:npmax2-1)
      common/tfinal/nfin
      common/dt/dt
      common/npulses/npulses
      dimension pulseshape(npulsemax),freq(npulsemax),pulse(npulsemax)
      dimension eta(npulsemax)
     
      integer nt
 
      nt=0
       
      open(1,file='propagation')
      open(70,file='autocorr.dat')
      open(115,file='proj.dat')
                                
      CALL input       
      call initstate(zwfin) 
              
!      open(50,file='wftarg.dat')
  
!      do i=0,npp(1)
!       do j=0,npp(2)
!        read(50,*) bb,cc,wf1,wfim1
!        zwftarg(i,j)=dcmplx(wf1,wfim1)
!       enddo
!      enddo
        
      read(1,*)nfich
      read(1,*)ires
      read(1,*)icpaq
 
C---  time propagation of the wave packet -------
 
      time=0.

      suma = 0.
      do ind=0,npp(1)
       do jnd=0,npp(2)
          zwf(ind,jnd)=zwfin(ind,jnd)
          suma = suma + cdabs(zwf(ind,jnd))**2
       end do
      end do 

c     Prepare arrays for fourier transform

      do i = 1,2
         call KEXP(expk(:,i),np(i),dr(i),xmu(i),dt)
         call KCIN(Tcin(:,i),np(i),dr(i),xmu(i))
         call prepr1d(np(i),trigx(:,i),ifax(:,i),nfax(i))
      enddo
    
      im = dcmplx(0.d0,1.d0)

c     suma = 0.       
      do i=0,npp(1)
         do j=0,npp(2)
          expV2(i,j) =  cdexp(-im*pot(i,j)*dt/2.)
         enddo
      enddo

c       diss(1)=0.
c       diss(2)=0.
c       diss(3)=0.
 
      call symmetry(zwf,symm)   
      call ENERGIA(zwf,pot,Tcin,time,symm)
      call Defwp(zwf,time)
      call pintafuncion(zwf,time,nt)
      nt=nt+1
      call solap(zsolap,zwfin,zwf)
c      call solap(zsolaptarg,zwftarg,zwf)
      call projection (zwf,zwftarg,time,LC,zproj)
      write(70,*)time/40.,cdabs(zsolap)**2,real(zsolap),dimag(zsolap)
      write (115,*) time/40.,LC,(cdabs(zproj))**2
      

      DO it=1,nfin
        
77      FORMAT(F10.5,4G15.7)
      
        call projection(zwf,zwftarg,time,LC,zproj)
          write (115,*) time/40.,LC,(cdabs(zproj))**2
        call solap(zsolap,zwfin,zwf)
        call solap(zsolaptarg,zwftarg,zwf)
        write(70,*)time/40.,cdabs(zsolap)**2,real(zsolap),dimag(zsolap)
        call Lasershape(time,LC,pulseshape)

        call Frequency(time,freq,eta)
        do k = 1, npulses
           pulse(k) = pulseshape(k) * dcos ( eta(k) )
        enddo
c        if wanted: print the pulses and/or other information
           if(mod(it-1,ires).eq.0.) then
               write(25,77)time/40.,(pulseshape(k),k=1,npulses)
               write(26,77)time/40.,(pulse(k),k=1,npulses)
           endif
          

        call Coupling(dipole,pulse,dt,expdip)
        call spo(zwf,expdip,expV2,expk,dt,zwfout)
   
        do ind = 0,npp(1)
        do jnd = 0,npp(2)
           zwf(ind,jnd) = zwfout(ind,jnd)
        enddo
        enddo
 
         call symmetry(zwf,symm)

c        call renormwf(zwf,zwfn)

c        do ind = 0,npp(1)
c         do jnd = 0,npp(2)
c           zwf(ind,jnd) = zwfn(ind,jnd)
c         enddo
c        enddo

         

         time=time+dt

          if(mod(it,ires).eq.0.) then
            call ENERGIA(zwf,pot,Tcin,time,symm)
c            call RESULTADOS(zwf,pot,Tcin,time,nfich,symm)
            call Defwp(zwf,time)
          
          end if
     
         if(mod(it,icpaq).eq.0.) then
c      ...  wave packets
            call pintafuncion(zwf,time,nt)
c            call pintafunction2(zfwf)
            nt=nt+1
          end if
    
  
       END DO

       open(60,file = 'wfend.dat')
       do i = 0,npp(1)
       do j = 0,npp(2)     
          write(60,*)z(i),x(j),real(zwf(i,j)),dimag(zwf(i,j))
       enddo
        write (60,*)
       enddo
       
      end



c-------------------------------------------------
      subroutine Coupling(dipole,pulse,dt,expdip)
      parameter (npmax1=1024,npmax2=1024,npulsemax=5)
      IMPLICIT REAL*8(A-H,O-y)
      real*8 dipole(0:npmax1-1,0:npmax2-1)
      complex*16 im
      complex*16 expdip(0:npmax1-1,0:npmax2-1)
      common/param/np(2),npp(2)
      common/zgrid/z(0:npmax1-1)
      common/xgrid/x(0:npmax2-1)
      common/npulses/npulses
      dimension pulse(npulsemax)

      im=dcmplx(0.d0,1.d0)
         
      do i=0, npp(1)
        do j=0 , npp(2)

         pulsetot = 0.d0
         do k = 1,npulses
            pulsetot = pulsetot + pulse(k)
         enddo
    
         expdip(i,j) = exp(-im*dipole(i,j)*pulsetot*dt/2.)
       
        enddo
      enddo   

      return
      end 
c ------------------------------------------------
      subroutine symmetry(zwf,symm)
      IMPLICIT REAL*8(A-H,O-y)
      parameter (npmax1=1024,npmax2=1024)
      common/rlim/dr(2)
      common/zgrid/z(0:npmax1-1)
      common/xgrid/x(0:npmax2-1)
      common/param/np(2),npp(2)
      complex*16 zwf(0:npmax1-1,0:npmax2-1)
      
      symm=0.
      do i=0,npp(1)
         do j=0,npp(2)
            symm = symm + dconjg(zwf(i,j))*zwf(j,i)
         enddo
      enddo
      symm=symm*dr(1)*dr(2)

      return
      end 
c-------------------------------------------------
      subroutine spo(zwin,expdip,expV2,expk,dt,zwout)
      IMPLICIT REAL*8(A-H,O-y)
      parameter (npmax1=1024,npmax2=1024)
      common/rlim/dr(2)
      common/zgrid/z(0:npmax1-1)
      common/xgrid/x(0:npmax2-1)
      common/param/np(2),npp(2)     
      common/FFTprep/trigx(0:20000,2)
      common/FFTfac/ifax(30,2),nfax(2)    
      complex*16 zwin(0:npmax1-1,0:npmax2-1)
     .          ,zwout(0:npmax1-1,0:npmax2-1)
     .          ,expk(0:npmax1-1,2)
      complex*16 expV2(0:npmax1-1,0:npmax2-1)
      complex*16 expdip(0:npmax1-1,0:npmax2-1)
c     complex*16 im
      complex*16 aux(0:20000)

c     im = dcmplx(0.d0,1.d0)

c   #  Create exp(-i V dt/2)
      do i = 0, npp(1)
        do j = 0, npp(2)
c            zwout(i,j) = expV2(i,j) * zwin(i,j)
           zwout(i,j) =  expdip(i,j) * expV2(i,j) * zwin(i,j)
        enddo
      enddo
  
c     second derivative with respect to first index 
c     the procedure is:
c     FFT followed by product to exp(-ik_11^2/2m) followed by FFT^-1

        do j=0,npp(2)     
           call tfft(zwout(:,j),aux,np(1),trigx(:,1)
     .               ,ifax(:,1),nfax(1),1,1)
           do i=0,npp(1)
              zwout(i,j)=zwout(i,j)*expk(i,1)/float(np(1))
           enddo

           call tfft(zwout(:,j),aux,np(1),trigx(:,1)
     .               ,ifax(:,1),nfax(1),1,-1)
        enddo
      
   
c     second derivative with respect to second index 
c     the procedure is:
c     FFT followed by product to exp(-ik_22^2/2m) followed by FFT^-1

        do i=0,npp(1)
           call tfft(zwout(i,:),aux,np(2),trigx(:,2)
     .              ,ifax(:,2),nfax(2),1,1)
           do j=0,npp(2)
              zwout(i,j)=zwout(i,j)*expk(j,2)/float(np(2))
           enddo

           call tfft(zwout(i,:),aux,np(2),trigx(:,2)
     .              ,ifax(:,2),nfax(2),1,-1)
        enddo


      do i = 0, npp(1)
      do j = 0, npp(2)
c         zwout(i,j) = expV2(i,j) * zwout(i,j)
         zwout(i,j) = expdip(i,j) * expV2(i,j) * zwout(i,j)
      enddo
      enddo

      return
      end

c------------------------------------------------------------
      SUBROUTINE Lasershape(time,LC,pulseshape)
      parameter(npulsemax=5)
      IMPLICIT REAL*8(A-H,O-y)
      implicit complex*16 (z)
      real*8 LC
      common/npulses/npulses
      common/radpar/rp(5,npulsemax)
      common/nparam/nparam(3)            !!! 3 types of pulse shapes?
      common/ltype/ltype(npulsemax)

      dimension pulseshape(npulses)

c      solaptarg=(cdabs(zstarg))**2   
 
!      open(116,file='input.laser')
!      read(116,*) tt,pulseLC
   
      npar = 1
      do k = 1, npulses

      if (ltype(k).eq.1) pulseshape(k) = pulseG(time,rp(1,k),nparam(1))
      if (ltype(k).eq.2) pulseshape(k) = pulseS2(time,rp(1,k),nparam(2))
      npar = npar + nparam(ltype(k))
      if (ltype(k).eq.3) pulseshape(k) = pulseS2(time,rp(1,k),nparam(3))
     . *LC
      if (ltype(k).eq.4) pulseshape(k) = pulseLC
      
      npar = npar + nparam(ltype(k))
      
!      open (110,file='pruebaLC.dat')

!      write(110,*) time/40., pulseLC

      enddo

      return
      end

c------------------------------------------------------------
      SUBROUTINE Frequency(time,freq,eta)
      parameter(npulsemax=5)
      IMPLICIT REAL*8(A-H,O-y)
      common/npulses/npulses
      common/lfreqvar/lfreqvar(npulsemax)
      common/freqpar/freqvar(npulsemax,4)
      common/nfpar/nfpar(4)

      dimension freq(npulses),eta(npulses)

c     freq(t) are the pulse frequencies
c     eta(t) are the pulse optical phase
c     by definition freq(t) = d eta(t) / dt
c     we first write what freq we want and obtain eta(t) = int freq(t')dt'

c     we choose freqvar(k,1) as t0 for all frequency modulations
      pi = dacos(-1.d0)

      do k = 1, npulses

      dt = time - freqvar(k,1)

c     -------> a transformed limited field:
         if (lfreqvar(k).eq.1) then
             freq(k) = freqvar(k,2)
             eta(k)  = freqvar(k,2) * dt
          endif
c     -------> quadratic chirp:
         if (lfreqvar(k).eq.2) then
             freq(k) = freqvar(k,2) +freqvar(k,3)*dt +freqvar(k,4)*dt**2
             eta(k)  = freqvar(k,2)*dt + freqvar(k,3)*dt**2/2.
     .               + freqvar(k,4)*dt**3/3.
         endif
c     -------> tanh chirp:  w(t) = w0 + DW/2 * tanh( 5*(t-t0)/DT ) 
c                          The pulse phase is eta(t) = int_0^t w(t') dt'
c                          = w0 (t-t0) + DW*DT/10 ln( cosh( 5*(t-t0)/DT ) )
         if (lfreqvar(k).eq.3) then
             freq(k) = freqvar(k,2) + freqvar(k,3)/2.
     .               * dtanh(5.*dt/freqvar(k,4))
             eta(k)  = freqvar(k,2)*dt + freqvar(k,3)*freqvar(k,4)/10.
     .               * dlog(cosh(5*dt/freqvar(k,4)))
         endif
c     --------> linear chirp with sign flip
         if (lfreqvar(k).eq.4) then
             if (time.le.freqvar(k,1)) then
                                              
                   freq(k) = freqvar(k,2) - freqvar(k,3) * dt
                   eta(k)  = -freqvar(k,2)*dt + freqvar(k,3)*dt**2/2.
              else
                   freq(k) = freqvar(k,2) + freqvar(k,3) * dt
                   eta(k)  = freqvar(k,2)*dt + freqvar(k,3)*dt**2/2.
             endif
         endif
c     --------> cos chirp : w(t) = dw * cos(2pi*dt/T) + wmin
c                           where dt = t-t0; dw = wmax - wmin
c                          (k=1, t0; k=2, wmin; k=3, dw; k=4, T)
         if (lfreqvar(k).eq.5) then
             freq(k) = freqvar(k,3)*dcos(2.*pi*dt/freqvar(k,4))
     .               + freqvar(k,2)
             eta(k)  = freqvar(k,3)*freqvar(k,4)/2./pi * dsin(
     .                 2.*pi*dt/freqvar(k,4)) + freqvar(k,2)*dt
         endif

      enddo

      return
      end

c--------------------------------------------------------
      Subroutine Defwp(zwp,time)
      IMPLICIT REAL*8(A-H,O-y)
      parameter(npmax1=1024,npmax2=1024)
      complex*16 zwp(0:npmax1-1,0:npmax2-1)
      common/param/np(2),npp(2)
      common/zgrid/z(0:npmax1-1)
      common/xgrid/x(0:npmax2-1)
      common/rlim/dr(2)
     
         avx1=0.       
         sdx1=0.
         xn=0.
         avx2=0.
         sdx2=0.
         avxR1=0.
         avxR2=0.
         do ind = 0,npp(1)
            do jnd =0,npp(2)
               density = cdabs(zwp(ind,jnd))**2
             
               xn = xn + density
             
               avx1 = avx1 + z(ind)*density
               avx2 = avx2 + x(jnd)*density

               avxR1 = avxR1 + density*(z(ind)-x(jnd))**2
               avxR2 = avxR2 + density*(z(ind)+x(jnd))**2

               sdx1 = sdx1 + density*z(ind)**2
               sdx2 = sdx2 + density*x(jnd)**2
            enddo
         enddo
      
c         xn=xn*dr(1)*dr(2)

         if(xn.gt.1.d-6) then
            avx1=avx1/xn
            sdx1=dsqrt(sdx1/xn-avx1**2)
c            avxR1=avxR1/xn
c            avxR2=avxR2/xn
            avxR1=1./2.*dsqrt(avxR1/xn)
            avxR2=1./2.*dsqrt(avxR2/xn)
            avx2=avx2/xn
            sdx2=dsqrt(sdx2/xn-avx2**2)
         else
            avx1=0.
            sdx1=0.
            avx2=0.
            sdx2=0.
            avxR1=0.
            avxR2=0.
         endif
     
      write(32,13)time/40.,avx1,sdx1,avx2,sdx2
      write(33,13)time/40.,avxR1,avxR2
13    FORMAT(G12.5,4F14.5)

      return
      end

c-----------------------------------------------
      function pulseG(time,rp,npar)
      IMPLICIT REAL*8 (a-h,o-z)
      dimension rp(npar)

c         arg1 = (time-rp(2))**2/2./rp(3)**2
c         arg1 = 4.d0*dlog(2.d0)*((time-rp(2))/rp(3))**2
         arg1 = ((time-rp(2))/rp(3))**2
         pulseG = rp(1) * dexp(-arg1)

      end

c------------------------------------------------------
      function pulseS2(time,rp,npar)
      IMPLICIT REAL*8 (a-h,o-z)
      dimension rp(npar)

      pi = dacos(-1.d0)

      en1arg = (time-rp(2))/2./(rp(3)-rp(2))
      ap1arg = (time-rp(4))/2./(rp(5)-rp(4))

c  #  start:
      if(time.lt.rp(2)) pulseS2=0.
c  #  switching on function:
         if(time.ge.rp(2).and.time.le.rp(3))
     .    pulseS2 = rp(1) * dsin(pi*en1arg)**2
c  #  constant field
         if(time.gt.rp(3).and.time.lt.rp(4)) pulseS2 = rp(1)
c  #  switching off function:
         if(time.ge.rp(4).and.time.le.rp(5))
     .    pulseS2 = rp(1) * dcos(pi*ap1arg)**2
c  #  off 
         if(time.gt.rp(5)) pulseS2 = 0.

      end


                                                     
c-----------------------------------------------------

      subroutine ENERGIA(zwf,V,Tcin,time,symm)
      IMPLICIT REAL*8(A-H,O-y)
      parameter (npmax1=1024,npmax2=1024)
      complex*16 zwf,zkin11,zkin22,zconj 
      complex*16 V
      common/zgrid/z(0:npmax1-1)
      common/xgrid/x(0:npmax2-1)
      common/param/np(2),npp(2)
      common/rlim/dr(2)
      common/dt/dt
      common/mass/xmu(2)
      common/FFTprep/trigx(0:20000,2)
      common/FFTfac/ifax(30,2),nfax(2)
      dimension zwf(0:npmax1-1,0:npmax2-1)
      dimension zkin11(0:npmax1-1,0:npmax2-1)
      dimension zkin22(0:npmax1-1,0:npmax2-1)
      dimension zconj(0:npmax1-1,0:npmax2-1)
      dimension V(0:npmax1-1,0:npmax2-1), sumK(2)
      dimension Tcin(0:npmax1-1,2)
      complex*16 aux(0:20000)
 
      do i=0,npp(1)
        do j=0,npp(2)
         zconj(i,j)=dconjg(zwf(i,j))
         zkin11(i,j)= zwf(i,j)
         zkin22(i,j)= zwf(i,j)
       enddo
      end do
    
      do j=0,npp(2)
           call tfft(zkin11(:,j),aux,np(1),trigx(:,1)
     .              ,ifax(:,1),nfax(1),1,1)
           do i=0,npp(1)
              zkin11(i,j)=zkin11(i,j)*Tcin(i,1)/dble(np(1))
           enddo
           call tfft(zkin11(:,j),aux,np(1),trigx(:,1)
     .              ,ifax(:,1),nfax(1),1,-1)
        enddo
    
c     second derivative with respect to second index 
c     the procedure is:
c     FFT followed by product to exp(-ik_22^2/2m) followed by FFT^-1

        do i=0,npp(1)
           call tfft(zkin22(i,:),aux,np(2),trigx(:,2)
     .              ,ifax(:,2),nfax(2),1,1)
           do j=0,npp(2)
              zkin22(i,j)=zkin22(i,j)*Tcin(j,2)/dble(np(2))
           enddo
           call tfft(zkin22(i,:),aux,np(2),trigx(:,2)
     .              ,ifax(:,2),nfax(2),1,-1)
        enddo
       
         sumK(1)=0.
         sumK(2)=0.
         sumP=0.
         sumN=0.
         sumNl=0.
         sumNr=0.

         do i=0,(npp(1)+1)/2-1
          do j=1,npp(2)+1
          sumNl=sumNl + cdabs(zwf(i,j))**2
          enddo
         enddo

         do i=(npp(1)+1)/2,npp(1)
          do j=0,npp(2)
          sumNr=sumNr + cdabs(zwf(i,j))**2
          enddo
         enddo
      
       do ind=0,npp(1)
        do jnd=0,npp(2)
         sumP=sumP + cdabs(zwf(ind,jnd))**2 *V(ind,jnd)
         sumK(1)=sumK(1)+dble(zconj(ind,jnd)*zkin11(ind,jnd))
         sumK(2)=sumK(2)+dble(zconj(ind,jnd)*zkin22(ind,jnd))
         sumN=sumN + cdabs(zwf(ind,jnd))**2 

        enddo
      end do
      
       dtau = dr(1)*dr(2)
       xnorm=sumN*dtau
       xnorml=sumNl*dtau
       xnormr=sumNr*dtau
       energy=(sumP+sumK(1)+sumK(2))*dtau
       Ekin=(sumK(1)+sumK(2))*dtau

      if (xnorm.eq.0.) then
          energy=0.
          symm=0.
          Ekin=0.
      else
          energy = energy / xnorm
          Ekin=Ekin / xnorm
          symm=symm / xnorm
      endif

      write(22,66) time/40.,xnorm,xnorml,xnormr
      write(23,66) time/40.,energy, symm
66    format(f14.5,5F14.6)


      write(4,62)time/40,sumP*dtau,sumK(1)*dtau,sumK(2)*dtau,sumN*dtau
62    format(4G15.5)

      return
      end
     

c-------------------------------------------- 
      SUBROUTINE INITSTATE(zwfin)

      IMPLICIT REAL*8(A-H,O-y)
      PARAMETER(npmax1=1024,npmax2=1024)
      COMPLEX*16 re
      complex*16 zwfin(0:npmax1-1,0:npmax2-1)
      common/param/np(2),npp(2)
c     common/wpin/zwfin(0:npmax1-1,0:npmax2-1)
   
      open(52,file='wfinit.dat')

c     The subroutine is not general: we are asuming that the
c     wave function is initially on a single potential.
c     It is not difficult to generalized the subroutine although
c     it will probably be rarely needed

       re=dcmplx(1.d0,0.d0)

       do ind=0,npp(1)
         do jnd=0,npp(2)
            read(52,*)r,y,wf,wfim
            zwfin(ind,jnd)=dcmplx(wf,wfim)
         enddo
        read(52,*)
       end do
       close(52)

      RETURN
      END

c--------------------------------------------       
      SUBROUTINE INPUT
      
      IMPLICIT REAL*8(A-H,O-y)
      parameter (npmax1=1024,npmax2=1024,npulsemax=5)
      real*8 zmin,zmax,xmin,xmax
      real*8 zi,zd,xd,xi
      complex*16 pot
      common/param/np(2),npp(2)
      common/rlim/dr(2) 
      common/zgrid/z(0:npmax1-1)
      common/xgrid/x(0:npmax2-1)
      common/mass/xmu(2)
      common/tfinal/nfin 
      common/dt/dt
      common/potential/pot(0:npmax1-1,0:npmax2-1)
      common/radpar/rp(5,npulsemax)
      common/npulses/npulses
      common/nparam/nparam(3)            !!! 3 types of pulse shapes?
      common/ltype/ltype(npulsemax)
      common/nfpar/nfpar(4)              !!! 3 types of freq variations?
      common/lfreqvar/lfreqvar(npulsemax)
      common/freqpar/freqvar(npulsemax,4)
      common/dipole/dipole(0:npmax1-1,0:npmax2-1)
      common/barrier/fint(0:npmax1-1,3)

      dimension xrp(10000)

C--------->grid data
      open(10,file='grid')
      read(10,*) np(1),np(2)
      npp(1)=np(1)-1
      npp(2)=np(2)-1
      read(10,*)zmin,zmax
      read(10,*)xmin,xmax
      dr(1) = (zmax - zmin) / dble(npp(1))
      dr(2) = (xmax - xmin) / dble(npp(2))
      close(10)


c--------> MASS
c     xmu=21112.48
      open(11,file='mass')
      read(11,*)xmu(1),xmu(2)
      close(11)

c--------> ELECTRONIC POTENTIALS      
      open(14,file='pes.dat')
      open(16,file='pesim.dat')
      open(78,file='zpot.dat')
      
      do i=0,npp(1)
         do j=0, npp(2)
           read(14,*)z(i),x(j),rpot
           read(16,*)z(i),x(j),potim
           pot(i,j)=dcmplx(rpot,potim)
         write(78,*) z(i),x(j),dreal(pot(i,j)),dimag(pot(i,j))

         enddo
           read(14,*)
           read(16,*)
           write(78,*)
      enddo
      close(14) 

C--------->TIME STEP AND NUMBER OF STEPS 
      write(*,*) 'introduzca dt (u.a.)'
      read(1,*) dt      
      write(*,*)'dt= ',dt
      write(*,*) 'numero de ciclos?'
      read(1,*) nfin
      write(*,*) 'numero de ciclos = ',nfin
     
      

c  Field parameters
      nparam(1) = 3     !! # parameters of pulsetype 1 -> Gaussian
      nparam(2) = 5     !! # parameters of pulsetype 2 -> sinesquare
                        !! several others can be implemented
      nparam(3) = 5

      open(55,file='param.dat')
      read(55,*)npulses
      nnp = 10000
      do 51 i=1,nnp
           read(55,*,end=61)xrp(i)
51    continue
61    continue
      close(55)

      do i = 1,5           !!! max # of parameters of pulse shape
      do j = 1,npulsemax
         rp(i,j) = 0.
      enddo
      enddo

      ipar = 1
      do i = 1, npulses
         ltype(i) = idint(xrp(ipar))       !!!! real to integer
         do j = 1,nparam(ltype(i))
            rp(j,i) = xrp(ipar+j)
         enddo
         ipar = ipar + nparam(ltype(i))+1
      enddo

c  chirp parameters
      nfpar(1) = 2     !! # parameters of freqvar 1 -> constant freq.
      nfpar(2) = 4     !! # parameters of freqvar 2 -> quadratic chirp
      nfpar(3) = 4     !! # parameters of freqvar 3 -> tanh chirp
      nfpar(4) = 3     !! # parameters of freqvar 4 -> linear chirp with
c                                                      sign flip

      open(55,file='freq.dat')
      nnp = 10000
      do 71 i=1,nnp
           read(55,*,end=81)xrp(i)
71    continue
81    continue
      close(55)
      
      do i = 1,npulsemax
      do j = 1,4             !!! max. # of parameters of freq changes
         freqvar(i,j) = 0.
      enddo
      enddo

      ipar = 1
      do i = 1, npulses
         lfreqvar(i) = idint(xrp(ipar))
         do j = 1,nfpar(lfreqvar(i))
            freqvar(i,j) = xrp(ipar+j)
         enddo
         ipar = ipar + nfpar(lfreqvar(i))+1
      enddo

c     write(*,*)(xrp(i),i=1,5)
      write(*,*)lfreqvar(1),(freqvar(1,j),j=1,nfpar(lfreqvar(1)))

c--------> DIPOLE COUPLINGs: Dipole functions and FC factors

      open(15,file='dipole.dat')
      do ind = 0,npp(1)
        do jnd =0,npp(2)
         read(15,*)z(ind),x(jnd),dipole(ind,jnd)
        enddo
        read(15,*)
      enddo
      close(15)


      return
      end

c--------------------------------------------------
      subroutine solap(zsolap,zw1,zw2)
      parameter(npmax1=1024,npmax2=1024)
      implicit real*8(a-h,o-y)
      complex*16 zsum,zsolap
      complex*16 zw1(0:npmax1-1,0:npmax2-1),zw2(0:npmax1-1,0:npmax2-1)
      common/param/np(2),npp(2)
      common/rlim/dr(2)
      
      zsum = dcmplx(0.,0.)

      do i=0,npp(1)
         do j=0,npp(2)
            zsum = zsum + dconjg(zw1(i,j))*zw2(i,j)
         enddo
      end do
      zsolap = zsum*dr(1)*dr(2)

      return
      end       


c-----------------------------------------------------
        subroutine KEXP(ker,n,delta,mass,dt)
 
        IMPLICIT REAL*8(A-H,O-y)
        complex*16 ker(0:n-1)
        real*8 mass

        pi=dacos(-1.d0)
        anorm=2.*pi/n/delta
        ker(0)=(1.,0.)
        do 10 j=1,n/2
                temp=(j*anorm)**2
                ker(j)=cdexp(dcmplx(0.d0,-1.*
     .                   temp*dt/2./mass))
                ker(n-j)=+ker(j)
10      end do
        continue

        return
        end

c----------------------------------------------------
        subroutine KCIN(ker,n,delta,mass)

        IMPLICIT REAL*8(A-H,O-y)
        real*8 ker(0:n-1)
        real*8 mass

        pi=dacos(-1.d0)
        anorm=2*pi/n/delta
        ker(0)= 0.

        do  j=1,n/2
                
                temp=(j*anorm)**2        
                ker(n-j)= temp/2./mass
                ker(j)=ker(n-j)

10      end do
              
        continue
 
        return
        end

c----------------------------------------------------
        subroutine pintafuncion(zwf,time,nt)
        parameter(npmax1=1024,npmax2=1024)
        implicit real*8(a-h,o-y)
        complex*16 zwf(0:npmax1-1,0:npmax2-1)
        common/param/np(2),npp(2)
        common/rlim/dr(2)
        common/dt/dt
        common/avener/etot
        common/zgrid/z(0:npmax1-1)
        common/xgrid/x(0:npmax2-1)
      
        integer nt
        character*6 fich

        write(fich,"('wf',I4.4)") nt
        open(99,file=fich)
        
        do i=0,npp(1)
         do j=0,npp(2)
           zconj=dconjg(zwf(i,j))
           prob=zconj*zwf(i,j)
c          write(nfich,*)i,prob+dt*float(itime)/80000. !!!,dimag(zwf(i-1))
           write(99,38)time/40.,z(i),x(j),prob    !!!+etot
c          i=400
c          write(nfich,38)time,zwf(i-1)
        enddo
        write(99,*)
       enddo
         
        close(99)
 
38      FORMAT(F12.5,3F14.7)
 
        return
                                                         
        end
c------------------------------------------------
      subroutine projection(zwf,zwftarg,time,LC,zproj)
      parameter(npmax1=1024,npmax2=1024)
      implicit real*8(a-h,o-y)
      implicit complex*16(z)
      real*8 LC
      common/param/np(2),npp(2)
      common/rlim/dr(2)
      common/dipole/dipole(0:npmax1-1,0:npmax2-1)
      dimension zwf(0:npmax1-1,0:npmax2-1)
      dimension zwftarg(0:npmax1-1,0:npmax2-1)


      zint1=(0.,0.)
      zint2=(0.,0.)

      do i=0,npp(1)
       do j=0,npp(2)
       zint1=zint1+(dconjg(zwftarg(i,j))*zwf(i,j))
       zint2=zint2+(dconjg(zwf(i,j))*zwftarg(i,j))*dipole(i,j)
       enddo
      enddo


      zproj=zint1*dr(1)*dr(2)
      zint2=zint2*dr(1)*dr(2)
      LC=(dimag(zproj*zint2))

      return
      end

c---------------------------------------------------------------
      subroutine renormwf(zwf,zwfn)
      parameter(npmax1=1024,npmax2=1024)
      implicit real*8(a-h,o-y)
      implicit complex*16(z)
      common/param/np(2),npp(2)
      common/rlim/dr(2)
      dimension zwf(0:npmax1-1,0:npmax2-1)
      dimension zwfn(0:npmax1-1,0:npmax2-1)


      suma=0.0
      do i=0,npp(1)
       do j=0,npp(2)
        suma=suma+cdabs(zwf(i,j))**2
       enddo
      enddo
 
      suma=suma*dr(1)*dr(2)
  
      do i=0,npp(1)
       do j=0,npp(2)
        zwfn(i,j)=zwf(i,j)/dsqrt(suma)
       enddo
      enddo

      end
