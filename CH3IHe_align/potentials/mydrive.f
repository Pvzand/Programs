       program mydrive
       implicit none

       character*2 c_def(100)
       real*8     bohr2Angs,pi
       parameter (bohr2Angs=0.529177249d0)
       parameter (pi=3.14156235d0)

       integer    i,iout
       integer    natom,nat
       parameter (natom=6)
       parameter (nat=natom+1)

       integer    type(nat)
       real*8     x(natom),y(natom),z(natom)
       real*8     rx(nat),ry(nat),rz(nat)
       real*8     xcom,ycom,zcom
       real*8     m(natom)
       real*8     vso(6)

       real*8     r,th,ph
       real*8     rCI,rCH,alpha
       parameter (rCH=1.084d0/bohr2Angs)
       parameter (rCI=2.136d0/bohr2Angs)
       parameter (alpha=107.47d0/360.d0*2.d0*pi)

       real*8     mtotal
       real*8     mp,mH,mC,mI,mHe
       parameter (mp=1836.d0)
       parameter (mH=mp,mC=12*mp,mI=127*mp,mHe=4*mp)

       c_def(1)= ' H'
       c_def(2)= 'HE'
       c_def(6)= ' C'
       c_def(53)=' I'
       c_def(100)='XX'

       r=3.d0
       th=0.5d0
       ph=0.d0


c______position I :
       type(1)=53
       m(1)= mI
       x(1)=0.d0
       y(1)=0.d0
       z(1)=rCI

c______position H :
       type(2)=1
       type(3)=1
       type(4)=1
       m(2)= mH
       m(3)= mH
       m(4)= mH

       x(2)= rCH*sin(alpha)  
       y(2)= 0  
       z(2)= rCH*cos(alpha)

       x(3)=-rCH*sin(alpha)/2  
       y(3)= (dsqrt(3.d0)/2.d0)*rCH*sin(alpha) 
       z(3)= rCH*cos(alpha)

       x(4)=-rCH*sin(alpha)/2  
       y(4)=-(dsqrt(3.d0)/2.d0)*rCH*sin(alpha)
       z(4)= rCH*cos(alpha)

c______position C :
       type(5)=6
       m(5)= mC
       x(5)= 0.d0
       y(5)= 0.d0
       z(5)= 0.d0

c______detemine com
       xcom=0.d0
       ycom=0.d0
       zcom=0.d0
       mtotal=0.d0
       do i=1,5
          mtotal=mtotal+m(i)
          xcom=xcom+m(i)*x(i)
          ycom=ycom+m(i)*y(i)
          zcom=zcom+m(i)*z(i)
       enddo
       xcom=xcom/mtotal
       ycom=ycom/mtotal
       zcom=zcom/mtotal

       write(*,*) xcom,ycom,zcom*bohr2Angs

       
       do i=1,5
       rx(i)=x(i)-xcom
       ry(i)=y(i)-ycom
       rz(i)=z(i)-zcom
       write(60,*) i,rx(i)*bohr2Angs,ry(i)*bohr2Angs,rz(i)*bohr2Angs
c      rx(i)=x(i)
c      ry(i)=y(i)
c      rz(i)=z(i)
       enddo


c______position He
       type(6)=2
       m(2)=mHe
       rx(6)=r*sin(th)*cos(ph)
       ry(6)=r*sin(th)*sin(ph)
       rz(6)=r*cos(th)

c______add com as pseud-atom
       type(7)=100
       rx(7)=0.d0
       ry(7)=0.d0
       rz(7)=0.d0

c______here geometry done

       call write_xbs(56,nat,type,c_def,rx,ry,rz)

       iout=0
       call vch3ihen(rx,ry,rz,natom,vso,iout)


       do i=1,6
       write(*,*) i,vso(i)
       enddo


       return
       end



       subroutine write_xbs(unit,nat,type,c_def,x,y,z)
       implicit none

       integer unit,nat,i
       integer type(nat)
       real*8  x(nat),y(nat),z(nat)
       character*2 c_def(53)

c______write xbs file
       do i=1,nat
          write(unit,99) 'atom  ',c_def(type(i)),x(i),y(i),z(i)
       enddo

       write(unit,*) ' spec    I  0.5   Yellow'
       write(unit,*) ' spec    C  0.5   Brown'
       write(unit,*) ' spec    H  0.2   Gray'
       write(unit,*) ' spec   HE  0.5   Red'
       write(unit,*) ' spec   XX  0.2   Black'
       write(unit,*) ' bonds  I C 1.0 8.0 0.1 0.0'
       write(unit,*) ' bonds  H C 1.0 3.0 0.1 0.0'
       write(unit,*) ' bonds  C O 1.0 3.0 0.1 0.5'
       write(unit,*) ' bonds  XX HE 1.0 8.0 0.1 0.0'

99     format(A6,A2,3X,F12.8,F12.8,F12.8)

       return
       end


