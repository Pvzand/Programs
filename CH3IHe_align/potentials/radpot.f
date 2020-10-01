      function radpot(l,r)
      implicit none
c
c     the He-CH3I interaction potential iin grd state, only m=0
c


      integer l
      real*8  radpot,r

      real*8  d,b,c


      select case (l) 

c  fort.50 2:3
      case (0)
      d               = 0.000185445   !  +/- 2.003e-06    (1.08%)
      b               = 0.782913      !  +/- 0.005996     (0.7659%)
      c               = 9.5243        !  +/- 0.01072      (0.1125%)

c  fort.50 2:4
      case (1)
      d               = -4.51657e-05  !  +/- 5.099e-07    (1.129%)
      b               = 0.837569      !  +/- 0.004355     (0.52%)
      c               = 10.3772       !  +/- 0.007319     (0.07053%)

c  fort.50 2:5
      case (2)
      d               = 2.53688e-05   !  +/- 2.933e-07    (1.156%)
      b               = 0.849398      !  +/- 0.002787     (0.3281%)
      c               = 10.6988       !  +/- 0.004548     (0.04251%)

c  fort.50 2:6
      case (3)
      d               = -6.51825e-06  !  +/- 4.938e-08    (0.7576%)
      b               = 0.779254      !  +/- 0.0061       (0.7828%)
      c               = 11.454        !  +/- 0.007429     (0.06486%)

c  fort.50 2:7
      case (4)
      d               = 1.5032e-06    !  +/- 2.002e-08    (1.332%)
      b               = 0.787553      !  +/- 0.0054       (0.6856%)
      c               = 12.236        !  +/- 0.01018      (0.08318%)

      case default
      stop 'in radpot: l out of range'
      end select

c     Morse function
      radpot=d*(1.d0-dexp(-b*(r-c)))**2-d

      return 
      end




