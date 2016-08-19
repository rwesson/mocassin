program makenk
  implicit none

  character(len=50) :: fn, fnout, rubbish
  real,pointer :: lam(:), n(:), k(:)
  real :: skip
  integer :: ios,nwav,i,nh

  print*, 'filename in?'
  read*, fn
  print*, 'filename out?'
  read*, fnout
  print*, 'how many lines in the header?'
  read*, nh

  open(unit=10,file=fn,status='old')
  do i = 1,nh
     read(10,*) rubbish
  end do
  
  do i = 1,1000000
     read(unit=10,fmt=*,iostat=ios)
     if (ios/=0) then
        nwav = i-1
        print*, 'lines of data = ', nwav
        exit
     end if
  end do

  rewind(10)

  allocate(lam(nwav))
  allocate(n(nwav))
  allocate(k(nwav))
  
  do i = 1,nh
     read(10,*) rubbish
  end do
  
  do i = 1,nwav
     read(unit=10,fmt=*,iostat=ios) lam(i), skip, skip, n(i),k(i)
     if (ios/=0) then
        print*, 'problem reading'
        stop
     end if
  end do
  
  close(10)


  open (unit=11,file=fnout, status='unknown')
  do i = nwav, 1, -1
     write(11,fmt='(3e12.5)') lam(i), n(i)+1, k(i)
  end do

  close(11)

end program makenk
