
       subroutine write_xbs(unit,nat,type,c_def,cc)
       implicit none

       integer unit,nat,i
       integer type(nat)
       real*8  cc(3,nat)
       character*2 c_def(26)

c______write xbs file
       do i=1,nat
          write(unit,99) 'atom  ',c_def(type(i)),cc(1,i),cc(2,i),cc(3,i)
       enddo

       write(unit,*) ' spec   FE  0.5   Yellow'
       write(unit,*) ' spec    C  0.5   Brown'
       write(unit,*) ' spec    O  0.5   IndianRed'
       write(unit,*) ' spec    N  0.5   Blue'
       write(unit,*) ' spec    H  0.2   Gray'
       write(unit,*) ' spec   XX  0.5   Black'
       write(unit,*) ' bonds FE C 1.0 8.0 0.1 0.0'
       write(unit,*) ' bonds FE N 1.0 4.0 0.1 0.0'
       write(unit,*) ' bonds  C C 1.0 3.0 0.1 0.0'
       write(unit,*) ' bonds  C N 1.0 3.0 0.1 0.0'
       write(unit,*) ' bonds  H C 1.0 3.0 0.1 0.0'
       write(unit,*) ' bonds  H N 1.0 3.0 0.1 0.0'
       write(unit,*) ' bonds  C O 1.0 3.0 0.1 0.5'
       write(unit,*) ' bonds  XX XX 1.0 3.0 0.1 0.0'
       write(unit,*) ' bonds  XX C 1.0 3.0 0.1 0.0'
       write(unit,*) ' bonds  XX N 1.0 3.0 0.1 0.0'
       write(unit,*) ' bonds  H XX 1.0 3.0 0.1 0.0'
       write(unit,*) ' bonds  XX N 1.0 3.0 0.1 0.0'
       write(unit,*) ' bonds  XX O 1.0 3.0 0.1 0.5'

99     format(A6,A2,3X,F12.8,F12.8,F12.8)

       return
       end


