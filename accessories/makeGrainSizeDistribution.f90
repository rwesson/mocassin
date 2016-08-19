! 
! this program greates a grain size distribution file for use 
! by the MOCASSIN code
! it creates nsizes equally spaced bins in log10 space from 
! amin and amax
! 
! B. Ercolano 23/09/04
program makeGrainSizeDistribution
  implicit none
  
  character(len=50) :: outfile
  integer :: nsizes, i
  real    :: slope, amin, amax, astep
  
  print*, 'how many bins?'
  read*, nsizes
  print*, 'what is the exponential giving the slope of the distribution? (p)'
  read*, slope
  print*, 'what are amin and amax in microns?'
  read*, amin, amax
  print*, 'what should the output file be called?'
  read*, outfile

  astep = (log10(amax)-log10(amin))/real(nsizes-1)

  open(unit=10,file=outfile,status='unknown')
  
  write(10,*) nsizes, ' sizes ' 
  do i=1,nsizes
     write(10,*) i, 10.**(log10(amin)+(i-1)*astep), (10.**(log10(amin)+(i-1)*astep))**(-slope)
  end do
  
  write(10,*)
  write(10,*)
  write(10,*) 'slope (p): ', slope

  
end program makeGrainSizeDistribution
  
  

  
  
