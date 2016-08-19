! this program produces sets of random locations 
! distributed within a sphere of a given radius
! and one twice the size of the original one
! normalised to 1.
program randomStars
  implicit none

  integer :: nsets, nstars, i, j

  real :: r, x, y, z, ran1, ran2, ran3, rs

  print*, 'please enter number of sets, number of sources and radius of the sphere normalised to 1.'
  read*, nsets, nstars, r
 

  open(file='randomstars1.out', unit=10)
  open(file='randomstars2.out', unit=11)

  call random_seed()

  write(10,*) 'radius of the sphere ', r
  write(11,*) 'radius of the sphere ', 2.*r

  do j =1, nsets
     write(10,*) 'set ', j
     write(11,*) 'set ', j
     
     do i = 1, nstars

        do 

           call random_number(ran1)
           call random_number(ran2)
           call random_number(ran3)
        
           rs = sqrt(ran1*ran1+ran2*ran2+ran3*ran3)
           
           if (rs<r) exit
           
        end do
        
        write(10,*) ran1, ran2, ran3
        write(11,*) 2.*ran1, 2.*ran2, 2.*ran3
        
     end do
  end do

  close(10)
  close(11)


end program randomStars

     
