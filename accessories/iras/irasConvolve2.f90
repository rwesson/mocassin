program irascolvolve2
  implicit none

  real, pointer :: nuArray(:), SED(:)
  real, parameter :: c = 2.9979250e10            ! speed of light [cm/s]
  real, parameter :: fr1Ryd = 3.28984e15         ! frequency a 1 Ryd [Hz]
  real, parameter :: pi = 3.141592654

  real :: firas(4), skip, distance
  
  integer :: nbins, i

  character(len=10) :: charskip

  logical :: lgSymmetricXYZ

  nbins = 700
  lgSymmetricXYZ = .true.

  print*, 'please enter distance in cm'
!  read*, distance
distance=4.26e21
  print*, 'merci'

  allocate(nuArray(nbins))
  allocate(SED(nbins))

  open(file='SED.out',unit=10)

  do i = 1, 5
     read(10,*) charskip
  end do
  do i =1, nbins
     read(10,*) nuArray(i), skip, SED(i)
  end do

  close(10)
  
  call  irasColours(SED,firas)

  print*, 'firas: ', firas

  contains
    
    ! calculates iras colours - adapted from dusty (ivezic, 1997)
    ! Filters data from Neugebauer et al, 1984, ApJ, 278, L1.
    ! Procedure described in Bedijn, 1987, A&A, 186, 136.
    ! sedin must be [erg/sec/str^-1]
    subroutine irasColours(sedin, colours)
      implicit none

      real, dimension(nbins), intent(in)   :: sedin
      real, dimension(4)    , intent(out)  :: colours

      real, dimension(nbins) :: fnu,faux,lambdaum
      real, dimension(4,7) :: phi, al
      real, dimension(4) :: cl, prz

      real :: wmid, fac

      integer :: i,j

      ! Data for 4 IRAS filters
      ! wavelengths
      al(1,:)=(/7.55, 8.0, 10.3, 11.5, 13.4, 14.7, 15.5/)
      al(2,:)=(/16.6, 22.5, 25.6, 26.8, 27.5, 29.3, 31.1/)
      al(3,:)=(/30.5, 40.1, 40.2, 65.2, 74.2, 83.8, 83.9/)
      al(4,:)=(/72.7, 95.4, 111.2, 116.6, 137.4, 137.5, 137.6/)
      ! transmittivities
      phi(1,:)=(/0.000, 0.618, 0.940, 0.750, 1.022, 0.906, 0.000/)
      phi(2,:)=(/0.235, 0.939, 0.939, 0.745, 0.847, 0.847, 0.000/)
      phi(3,:)=(/0.000, 0.102, 0.260, 1.026, 0.842, 0.001, 0.000/)
      phi(4,:)=(/0.000, 0.910, 1.000, 0.330, 0.002, 0.001, 0.000/)
      
      ! first get rid of the str^-1 and transform to Jy
      fnu=sedin*Pi
      if (lgSymmetricXYZ) fnu = fnu*8.
      do i =1 , nbins
         fnu(i) = fnu(i)*1.e23*1.e36/(nuArray(i)*fr1Ryd*4.*Pi*distance*distance)
         print*, i, c*1.e4/(nuArray(i)*fr1ryd),fnu(i)
      end do
      faux=fnu

      ! express nuArray and SED in lambda[um]
      do i = 1, nbins
         lambdaum(nbins-i+1) = c*1.e4/(nuArray(i)*fr1ryd)
         fnu(nbins-1+1) = faux(i)
      end do

      ! find IRAS colors
      do j = 2, nbins
         ! mid wav
         wmid = 0.5*(lambdaum(j-1)+lambdaum(j))
         ! interpolate filters
         call phiLam(wmid,prz,al,phi)

         ! add contribution to integral
         do i = 1, 4
            colours(i) = colours(i) + prz(i)*0.5*(fnu(j-1)+fnu(j))*&
                 & c/((-lambdaum(j-1)+lambdaum(j))*1.e4*fr1ryd)
            cl(i) = cl(i)+prz(i)*(-lambdaum(j-1)+lambdaum(j))
         end do
      end do
      do i = 1, 4
!         colours(i) = colours(i) / cl(i)
      end do
      
    end subroutine irasColours

    subroutine phiLam(lin,tout,lam,tran)
      implicit none
      
      real, intent(in)  :: lin
      real, intent(out) :: tout(4)
      
      real, dimension(4,7) :: lam,tran
      real :: a, b
      
      integer :: iloc, im
      
      do iloc = 1, 4
         tout(iloc) = 0.0
         im = 0
         do
            im = im + 1
            if ( (lin-lam(iloc,im))*(lin-lam(i,im+1)) <= 0.) then
               a = (tran(iloc,im+1)-tran(iloc,im))&
                    &/log10(lam(iloc,im+1)/lam(iloc,im))
               b = tran(iloc,im) - a*log10(lam(iloc,im))
               tout(iloc) = a*log10(lin) + b
               if ( tout(iloc) > 1.0 ) tout(iloc) = 1.0
            end if
            if (im == 6) exit
         end do
      end do
    end subroutine phiLam   

  end program irascolvolve2
