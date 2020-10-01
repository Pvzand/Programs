      program COULOMB_POTENTIAL
      implicit none
      real*8 z,zmin,zmax,dz,r0,Q,V,z0,Ve,Vn
      real*8 Rmin, Rmax,dR, R
      integer i,j,npe,npn

      open(10,file='grid')
      open(15,file='input.coulomb')       
      open(20,file='pes.dat')
!      open(21,file='PES.dat')
      

      read(10,*)npe,npn
      read(10,*)zmin,zmax
      read(10,*)Rmin,Rmax
      
      dz = (zmax - zmin) / dble(npe-1)
      dR = (Rmax - Rmin) / dble(npn-1)
      
c    Soft Coulomb Potencial H2+
c    R internuclear distance
c    z electronic coordinate

      read(15,*)Q
      read(15,*)r0

      z = zmin

      do i = 1,npe
       R=Rmin
       do j = 1,npn
        
      Ve = -1/dsqrt(1.d0+(z-R/2)**2)-1/dsqrt(1.d0+(z+R/2)**2)
      Vn = 1/dsqrt(0.03d0+R**2)
        
        V=Ve+Vn

        write(20,*) z,R,V
        R=R+dR
                
      end do

      write(20,*)
      z = z+dz

      enddo
      
      end
