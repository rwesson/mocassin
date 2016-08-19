program irasConvolve
  implicit none

  integer :: nbins, i,j,k,ii,jj,kk,ios

  real :: lambda(700), Fnu(700), dnu(4000), nlam(4),band(3,4000),&
       & nu(700),D,Pi = 3.1416, TotF(4)=0., normalization, nuRyd(4000)

  character(len=10) :: filein(4)


  nlam = (/851,776,601,3751/)

  nbins = 700

  ! distance in cm
  D = 4.6e21

  ! the mocassin file
  open(unit=10, file='SED.out', iostat=ios)
  
  do i = 1, 4
     read(10,*)
  end do
  do i = 1, nbins
     read(10,*) nu(i), lambda(i), Fnu(i)
  end do

  close(10)

  Fnu = Fnu*Pi

  ! calculate Fnu in Jy (now you have nu Fnu/ dnu (ryd) )
  do i = 1, nbins

     Fnu(i) = Fnu(i)*1.e23*8.*1.e36/(nu(i)*3.28984e15*4.*Pi*D*D)

  end do
  i=nbins

  filein(1) = '12.dat'
  filein(2) = '25.dat'
  filein(3) = '60.dat'
  filein(4) = '100.dat'

  ! the transmission curves
  open(unit=11, file=filein(1), iostat=ios)  
  open(unit=12, file=filein(2), iostat=ios)  
  open(unit=13, file=filein(3), iostat=ios)  
  open(unit=14, file=filein(4), iostat=ios)  


  print*, '        band      F[Jy]     normalization'
  do i = 1, 4
     
     ii = 10+i


     ! A to ryd
     do j = 1, nlam(i)
        read(ii, *) band(1, j), band(2, j)

        nuRyd(j) = 2.9978e10/(band(1,j)*1.e-8*3.28984e15)

     end do

     ! reverse order
     do j = 1, nlam(i)
        band(1,j) = nuRyd(nlam(i)+1-j)

     end do

     ! find dnu in ryd 
     do j = 1,nlam(i)-1
        dnu(j) = band(1,j+1)-band(1,j)

     end do
     dnu(nlam(i)) = dnu(nlam(i)-1)

     normalization = 0.
     do j = 1, nlam(i)
        
        ! locate lambda_filter on nu array        
        kk=0
        do k = 1, nbins
           if (band(1,j) > nu(k) ) then
              if (band(1,j) < nu(k+1) ) then
                 kk = k
                 if (kk<1) then
                    print*, 'problem 1', band(1,j), nu(1), nu(700)
                    stop
                 end if
                 exit
              end if
           end if
        end do

        if (Fnu(kk+1)>=Fnu(kk)) then
           band(3,j) = Fnu(kk) + ((Fnu(kk+1)-Fnu(kk))*(band(1,j)-nu(kk))/(nu(kk+1)-nu(kk)))
        else
           band(3,j) = Fnu(kk+1) - ( (Fnu(kk+1)-Fnu(kk))*(nu(kk+1)-band(1,j))/(nu(kk+1)-nu(kk)))
        end if

        TotF(i) = TotF(i) + band(3,j)*band(2,j)*dnu(j)

        normalization = normalization+dnu(j)*band(2,j)
     end do

     TotF(i) = TotF(i)/normalization

     print*, i, filein(i), TotF(i), normalization
    
  end do

end program irasConvolve
