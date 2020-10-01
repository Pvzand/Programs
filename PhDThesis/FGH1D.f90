program FGH1D
implicit none
integer n
complex*16,allocatable :: wf(:)
real*8, allocatable :: PES(:),H(:,:),r(:),autoval(:)
real*8 dr,dk,kk
integer i,j,k
integer pot

integer nfunc
character*50 fich
real*8 mass,masa

write(6,*) "Number of points and potential"
read(5,*) n,pot
write(6,*) n,pot
allocate( wf(n),PES(n),H(n,n),r(n),autoval(n) )
write(6,*) "Number of eigenfunctions"
read(5,*) nfunc
write(6,*) "I will write the first ",nfunc

open(1,file="PES.dat")
do i=1,n
 read(1,*) r(i),PES(i)
enddo
close(1)
dr=r(2)-r(1)
dk=2.d0*dacos(-1.d0)/n/dr
write(6,*) "PES.dat read"

do i=1,n
 do j=1,n
  H(i,j)=0.d0
 enddo
 H(i,i)=PES(i)
enddo
write(6,*) "Potential done"

do j=1,n
!Preparo para FFT
!write(2,*) "Funcion ",j
 do k=1,n
  wf(k)=dcmplx(0.d0,0.d0)
  if (k.eq.j) wf(k)=dcmplx(1.d0,0.d0)
 enddo
 call tfft(wf,n,n,1)
 do k=1,n
  if (k.le.n/2) then
   kk=k-1
  else
   kk=k-n-1
  endif
  wf(k)=(kk*dk)**2*wf(k)
 enddo
 call tfft(wf,n,n,-1)
 do k=1,n
  mass=masa(pot,r(k))
  wf(k)=wf(k)/2./mass
! write(2,*) wf(k)
 enddo
 do i=1,n
   H(i,j)=H(i,j)+dreal(wf(i))
 enddo
enddo

write(6,*) "Kinetic part ready"

call diag(H,n,autoval)

!reorder
do i=1,n-1
 do j=i+1,n
  if (autoval(i).ge.autoval(j)) then
   kk=autoval(i)
   autoval(i)=autoval(j)
   autoval(j)=kk
   do k=1,n
    kk=H(k,i)
    H(k,i)=H(k,j)
    H(k,j)=kk
   enddo
  endif
 enddo
enddo

write(fich,"('val',I3.3,'.dat')") pot
open(1,file=fich)
do i=1,n
 write(1,900) autoval(i)
enddo
close(1)
 
open(10,file='wf.dat')
  do j=1,n
   write(10,900)r(j),(H(j,i),i=1,nfunc)
  enddo
  close(10)



900 format(1000(1X,E20.10e3))
end
