c     program test-ch3ihenpot

c     first compile with makefile -f Makefiletest
ccc   then move test-ch3ihenpot to data/
ccc   then cd data/
ccc   choose which set of parameters:
ccc        ch3he_90.p (CH3 planar)
ccc     or ch3he_1095.p  (CH3 pyramidal)
ccc   and copy it to ch3he.p
ccc   (coefficients for the expansion of the potential parameters in Tesseral harmonics)
ccc   then ./test-ch3ihenpot

      
      implicit real*8 (a-h,o-z)

      character*80    atitle

      parameter (nhemax=10,natmax=nhemax+5)

      character*3     aatom(natmax)
      dimension xang(natmax), yang(natmax), zang(natmax)
      dimension xau(natmax), yau(natmax), zau(natmax)
      dimension vso(6)

      parameter(bohr2Angs=0.529177249d0)

chm   open(20,file='./ch3ihe_90.xyz',status='old')
      open(20,file='./geom.xyz',status='old')
      read(20,*) natom
      read(20,'(a)') atitle
      print *, atitle
      natom = 0
      print *
      print *,' reading input coords (in Angs.) from file 20 '
      do iat = 1, natmax
         read(20,*,end=99000) aatom(iat),xang(iat),yang(iat),zang(iat)
         print *,aatom(iat),xang(iat),yang(iat),zang(iat)
         xau(iat) = xang(iat)/bohr2Angs
         yau(iat) = yang(iat)/bohr2Angs
         zau(iat) = zang(iat)/bohr2Angs
         natom = natom + 1
      enddo
99000 continue
      print *,natom, ' input coordinates read '

      iout = 6
      call vch3ihen(xau,yau,zau,natom,vso,iout)
      

      
      stop
      end
      
