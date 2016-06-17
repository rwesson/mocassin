program makenuryd
  implicit none

  integer :: ios, i
  real :: lambda,c=2.9979250e10,ryd=3.28984e15, nu(105)


  open(file='dustysed.dat', unit=10, iostat=ios, status = 'old')
  open(file='nuRyd.dat', unit=11, iostat=ios, status = 'unknown')
  
  do i = 1, 4
     read(10, *)
  end do

  do i = 1, 102
     read(10, *) lambda
     nu(i) = c/(lambda*1.e-4*ryd)
  end do
  
  do i = 102, 1, -1
     write(11, *) nu(i)
  end do          

  close(10)
  close(11)

end program makenuryd
     
     
