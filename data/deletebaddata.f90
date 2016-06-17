program dbd
  implicit none

  integer :: i, lc, j, nlev

  character(len =30) :: fn, skip

  open(file = 'fileNames.dat', unit = 10)

  do i = 1, 219

     read(10,*) fn


     open(file = fn, unit = 11)

     read(11,*) lc

     do j = 1, lc
        read(11, *) skip
     end do
     read(11,*) nlev

     if (nLev == 1 .or. nLev == 0) print*, fn, nlev
     close(11)

  end do
  close(10)

end program dbd
