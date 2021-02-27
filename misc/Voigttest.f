      program voigttest
      implicit none
      real*8 voigt
      real*8 alpha(401), v(1001)
      real*8 H(401, 1001)
      real*8 x
      integer clock_start, clock_stop, clock_rate
      real e_time
      integer i, j

      do 10 i = 1,401
        alpha(i) = (10_8**((i-1_8)*0.01_8 - 3_8))/2_8
10    continue

      do 20 i = 1, 1001
        v(i) = (i-1_8)*0.01_8
20    continue
            
      call system_clock(count_rate=clock_rate) !Find the time rate
      call system_clock(count=clock_start)     !Start Timer
      do 30 i = 1, 401
        do 40 j = 1, 1001
          H(i,j) = voigt(alpha(i), v(j))
40      continue
30    continue
      call system_clock(count=clock_stop)      ! Stop Timer
      e_time = real(clock_stop-clock_start)/real(clock_rate)
      write (*,*) e_time
      write (*,*)

      write (*,*) 'verification values'
      do 50 i = 1, 400, 50
        write (*,*) H(i,i)
50    continue

      stop
      end
