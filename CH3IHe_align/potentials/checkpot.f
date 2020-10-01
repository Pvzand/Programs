      program checkpot
      implicit none


      integer i,n
      real*8 j,r,v00,v10,v20,v30,v40
      real*8 radpot


      open(50,file='fort.50')
      open(60,file='pc.dat')

      do i=1,60
    
      read(50,*) j,r,v00,v10,v20,v30,v40

      write(60,*) j,r
     .  ,v00,radpot(0,r)
     .  ,v10,radpot(1,r)
     .  ,v20,radpot(2,r)
     .  ,v30,radpot(3,r)
     .  ,v40,radpot(4,r)

      enddo

      close(50)
      close(60)

      return
      end

