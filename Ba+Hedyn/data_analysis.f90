PROGRAM analyse
  IMPLICIT NONE
  logical :: lterminate
  !real     (kind=8), dimension(6)  ::  caract, caract_old
  integer :: traj,traj2,tmax,nbtraj,i,itt,l
  integer :: nHe,param, nbhop, ioldhop,ihop,num_traj
  double precision ::  bohr,ang2bohr,ligne_old,rdist,ligne,pico,lecture_energ, distance,ligne2,t,maxn
  real(kind=8), dimension(:), allocatable  ::  caract,caract_old,r,r_old,energie
  real(kind=8), dimension(:), allocatable  :: distance_instantannee,distance_instantannee_old,tab
  double precision, dimension(501) :: temps
 integer :: nbdissoc, compteur, k,kk
  integer, dimension(501) :: nbdissoc_tot  

  open(33,file="exploitation.in")
  read(33,*) nHe
  read(33,*) rdist
  read(33,*) nbtraj
  param=3*nHe+3
  compteur=0
  k=2
  
  allocate (caract(param),caract_old(param))
  allocate (r(nHe),r_old(nHe),distance_instantannee(nHe))
  allocate (distance_instantannee_old(nHe),tab(nHe))
  allocate(energie(nHe))

  bohr = 0.529177249d0
  pico=2.418884327d-5
  
  
  tmax=tmax/pico
  rdist=rdist/bohr
  
  
  open(20,file='./3ihop=1/coord3',STATUS='UNKNOWN')
  open(21,file='./3ihop=1/energies3',STATUS='UNKNOWN')
  
 



  
      write(101,'(1(1x,g25.16))') '#nombre de trajectoire :'
      write(101,'(1(1x,g25.16))') nbtraj
      write(101,'(1(1x,g25.16))') '#Distance temps initial'
 
      write(100,'(1(1x,g25.16))') '#nombre de trajectoire :'
      write(100,'(1(1x,g25.16))') nbtraj
      write(100,'(1(1x,g25.16))')  '#Distance temps final'
 
      write(102,'(1(1x,g25.16))') '#nombre de trajectoire :'
      write(102,'(1(1x,g25.16))') nbtraj
      write(102,'(1(1x,g25.16))')'#energie totale par trajectoire'
 
      write(103,'(1(1x,g25.16))')'#nombre de trajectoire :'
      write(103,'(1(1x,g25.16))') nbtraj
      write(103,*)'#Distance instantannee Ba+He (pour He non &
                                  dissocies, les He dissocies sont fixes à 0)'
      write(104,'(1(1x,g25.16))')'#nombre de trajectoire :'
      write(104,'(1(1x,g25.16))') nbtraj 
      write(104,*) '#A t final : Cas avec tout les He dissociees'
      write(104,*) '#Numéro trajectoire  Nombre de saut  Surface_finale'

      write(105,'(1(1x,g25.16))')'#nombre de trajectoire :'
      write(105,'(1(1x,g25.16))') nbtraj 
      write(105,*) '#A t final : Cas où certains He sont liee '
      write(105,*) '#Numéro de trajectoire','Nombre de He dissocies','distance des He liee'
        
      write(106,*) '#nombre de trajectoires'
      write(106,*) nbtraj
      write(106,*) 'Nombre de He dissocies en fonction du temps'
     
      write(107,*) '#nombre de trajectoires'
      write(107,*) nbtraj
      write(107,*) '#Positions finales'





      ligne_old=0.d0
      traj2=1
      read(20,*,end=99)ligne,caract(1:3+3*nHe)
      read(21,*,end=99)ligne,lecture_energ
      do i=1,nHe
        r(i)=sqrt((caract(1)-caract(3*i+1))**2+(caract(2)-caract(3*i+2))**2+(caract(3)-caract(3*i+3))**2)
      end do
      write(102,*) traj2, lecture_energ
      write(101,*) traj2, r  !ajout de la première valeur du temps initial
      write(103,*) ligne, r 
      write(106,*) ligne,0
      nbdissoc_tot(1)=0
      kk=2
10    continue
     
      read(20,*,end=99)ligne,caract(1:3+3*nHe) 
      read(21,*)ligne2,lecture_energ 
      do i=1,nHe
         distance=sqrt((caract(1)-caract(3*i+1))**2+(caract(2)-caract(3*i+2))**2+&
                   (caract(3)-caract(3*i+3))**2)
         if (distance<rdist) then
            distance_instantannee(i)=distance
         else
            distance_instantannee(i)=0
         end if
      end do

     
      nbdissoc=0      
      do i=1,nHe
          if (distance_instantannee(i)<1.d-15) then
              nbdissoc=nbdissoc+1
          end if 
      end do
      
     write(106,*) ligne, nbdissoc

      l=0
      if ( ligne-ligne_old<-1.d-15) then
         kk=1
        ! print *,ligne,caract(1:3+3*nHe)
         write(103,*)
            

         do i=1,nHe 
            r(i)=sqrt((caract(1)-caract(3*i+1))**2+(caract(2)-caract(3*i+2))**2+(caract(3)-caract(3*i+3))**2)
            r_old(i)=sqrt((caract_old(1)-caract_old(3*i+1))**2+(caract_old(2)-caract_old(3*i+2))**2+&
                     (caract_old(3)-caract_old(3*i+3))**2)
         end do
       write(100,*) traj2, r_old(1:nHe)    ! temps final
       write(107,*) caract_old(1:param)
       write(107,*) '#Traj suivante' 
     
      
         if ( sum (distance_instantannee_old(1:nHe))<1.d-10) then
                open(22,file='./3ihop=1/saut_fin3')

             do i=1,compteur     
                 read(22,*) 
             end do
                
             do i=1,10  
                  read(22,*) itt,lterminate,t,ihop,nbhop
                  if (lterminate .eqv. .false.) then
                        
                         if ( itt==traj2-compteur) then
                               write(104,*) traj2, nbhop,  ihop


                         end if
                  end if
            end do               
         
         
          else
            do i=1,nHe
                if (distance_instantannee_old(i)>1.d-12) then 
                   tab(i-l)=distance_instantannee_old(i)  
                else 
                   l=l+1
                end if
            end do
            write(105,*) traj2,l,tab(1:nHe-l)            
          end if
      

      
       if (modulo(traj2,10)==0) then
         compteur=compteur+10
       end if  
        traj2=traj2+1
       write(102,*) traj2,lecture_energ
       write(101,*) traj2,r(1:nHe)        ! temps initial
       
      end if
      
      if (dabs(ligne-ligne_old )>1000) then
         if (k<501) then
        temps(k)=ligne
      k=k+1
      end if

          nbdissoc_tot(kk)=nbdissoc_tot(kk)+nbdissoc
          print *, kk,ligne,nbdissoc_tot(kk),nbdissoc
          kk=kk+1
      end if
     

      write(103,*) ligne, distance_instantannee(1:nHe)
      ligne_old=ligne
      caract_old=caract
     
   
      distance_instantannee_old=distance_instantannee
      close(22)      

      goto 10
      
99  continue
    
!    write(107,*) caract_old(1:param)
        
    do i=1,nHe   
       r_old(i)=sqrt((caract(1)-caract(3*i+1))**2+(caract(2)-caract(3*i+2))**2+(caract(3)-caract(3*i+3))**2)
    end do
    write(100,*) traj2, r_old(1:nHe)    ! ajout de la dernière valeur du temps final 
    
    l=0
!    if ( sum (distance_instantannee_old(1:nHe))<1d-10) then
!                open(22,file='./results/saut_fin')
!                print *,distance_instantannee_old(2)
!                do i=1,nbtraj
!                 read(22,*) itt,lterminate,t,ihop,nbhop
!                  if (lterminate .eqv. .false.) then
!
!                         if ( itt==traj2) then
!                               write(104,*) itt, nbhop,  ihop
!

 !                        end if
 !                 end if
 !           end do
if ( sum (distance_instantannee_old(1:nHe))<1.d-10) then
                open(22,file='./3ihop=1/saut_fin3')

             do i=1,compteur
                 read(22,*)
             end do

             do i=1,10
                  read(22,*) itt,lterminate,t,ihop,nbhop
                  if (lterminate .eqv. .false.) then

                         if ( itt==traj2-compteur) then
                               write(104,*) traj2, nbhop,  ihop


                         end if
                  end if
            end do


      else
            do i=1,nHe
                if (distance_instantannee_old(i)>1.d-8) then
                   tab(i-l)=distance_instantannee_old(i)
                else
                   l=l+1
                end if
            end do
            write(105,*) traj2,l,tab(1:nHe-l)
      end if
      
      maxn=0
      do i=1,501
          !print *,'transititon'
          !print *, nbdissoc_tot(i)
          if ( nbdissoc_tot(i)<maxn) then 
             write(108,*) temps(i),1.d0*maxn/100
         else
             write(108,*) temps(i),1.d0*nbdissoc_tot(i)/100
             maxn=nbdissoc_tot(i)
          end if
      
      end do
close(20)
close(21)
close(22)
close(100)
close(101)
close(102)
close(103)
close(104)
close(105)
end PROGRAM
