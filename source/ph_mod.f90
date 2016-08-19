! Copyright (C) 2005 Barbara Ercolano 
!
! Version 2.02
module xSec_mod

    use common_mod
    use constants_mod
    use elements_mod
    use interpolation_mod
    implicit none

        
    real, private, allocatable, dimension(:) :: xSecArrayTemp ! temporary xSec array


    public

    ! common variables
    integer, save                      :: bremsXSecP  ! pointer to the brems x section in xSecArray
    integer, save, dimension(nHlevel)  :: HlevXSecP   ! pointer to the nth H level in xSecArray
    integer, save, dimension(9)        :: HeISingXSecP! pointer to the HeI singlet in xSecArray
    integer, save, dimension(9)        :: HeIIXSecP   ! pointer to the HeII bound free in xSecArray
    integer, save, dimension(7)        :: level
    integer, save, dimension(30)       :: nInn
    integer, save, dimension(30)       :: nTot
    integer, save                      :: xSecTop     ! pointer to the top entry of xSecArray


    contains

    subroutine initGammaCont()
      implicit none

      integer :: err ! allocation, i/o error
      integer :: ntkold, itk, i

      print*, 'in initGammaCont'

      close(21)
      open(file='data/gammaHI.dat', unit=21, status='old', iostat = err)
      if (err /= 0) then
         print*, "! initGammaCont: can't open file: data/gammaHI"
         stop
      end if
      close(22)
      open(file='data/gammaHeI.dat', unit=22, status='old', iostat = err)
      if (err /= 0) then
         print*, "! initGammaCont: can't open file: data/gammaHeI"
         stop
      end if
      close(23)
      open(file='data/gammaHeII.dat', unit=23, status='old', iostat = err)
      if (err /= 0) then
         print*, "! initGammaCont: can't open file: data/gammaHeII"
         stop
      end if

      read(21, *) nTkGamma, nlimGammaHI
      ntkold=nTkGamma

      read(22, *) nTkGamma, nlimGammaHeI
      if (ntkold /= ntkGamma) then
         print*, '! initGammaCont: the number of temperature points must &
              & be the same  for all gammas'
         stop
      end if

      read(23, *) nTkGamma, nlimGammaHeII
      if (ntkold /= ntkGamma) then
         print*, '! initGammaCont: the number of temperature points must &
              & be the same  for all gammas'
         stop
      end if

      allocate(tkGamma(nTkGamma), stat=err)
      if (err /= 0) then
         print*, "! initGammaCont: Cannot allocate grid array -0"
         stop
      end if 

      allocate(nuGammaHI(nlimGammaHI), stat=err)
      if (err /= 0) then
         print*, "! initGammaCont: Cannot allocate grid array -01"
         stop
      end if 
      allocate(nuGammaHeI(nlimGammaHeI), stat=err)
      if (err /= 0) then
         print*, "! initGammaCont: Cannot allocate grid array -02"
         stop
      end if 
      allocate(nuGammaHeII(nlimGammaHeII), stat=err)
      if (err /= 0) then
         print*, "! initGammaCont: Cannot allocate grid array -03"
         stop
      end if 
      allocate(HINuEdgeP(nlimGammaHI), stat=err)
      if (err /= 0) then
         print*, "! initGammaCont: Cannot allocate grid array -03"
         stop
      end if 
      allocate(HeINuEdgeP(nlimGammaHeI), stat=err)
      if (err /= 0) then
         print*, "! initGammaCont: Cannot allocate grid array -03"
         stop
      end if 
      allocate(HeIINuEdgeP(nlimGammaHeII), stat=err)
      if (err /= 0) then
         print*, "! initGammaCont: Cannot allocate grid array -03"
         stop
      end if 

      read (21,*) (tkGamma(i), i=1,ntkgamma)
      read (22,*) (tkGamma(i), i=1,ntkgamma)
      read (23,*) (tkGamma(i), i=1,ntkgamma)



      ! allocate logGammaHI, logGammaHeI, logGammaHeII
      allocate(logGammaHI(nTkGamma, nlimGammaHI), stat=err)
      if (err /= 0) then
         print*, "! initGammaCont: Cannot allocate grid array -1"
         stop
      end if
      allocate(logGammaHeI(nTkGamma, nlimGammaHeI), stat=err)
      if (err /= 0) then
         print*, "! initGammaCont: Cannot allocate grid array -2"
         stop
      end if
      allocate(logGammaHeII(nTkGamma, nlimGammaHeII), stat=err)
      if (err /= 0) then
         print*, "! initGammaCont: Cannot allocate grid array -3"
         stop
      end if

      do i = 1, nlimGammaHI
         read(21,*) nuGammaHI(i), (logGammaHI(itk, i), itk=1, nTkGamma)         
         call locate(nuArray(1:nbins), nuGammaHI(i), HINuEdgeP(i))
         if (HINuEdgeP(i)<nbins) then
            if (nuGammaHI(i) > (nuArray(HINuEdgeP(i))+nuArray(HINuEdgeP(i)+1))/2.) &
                 & HINuEdgeP(i) = HINuEdgeP(i)+1
         end if
      end do

      do itk = 1, ntkgamma
         do i = 1, nlimGammaHI
            if (logGammaHI(itk,i)>0.) then         
               logGammaHI(itk,i) = log10(logGammaHI(itk,i))
            else
               logGammaHI(itk,i) = 0.
            end if
         end do
      end do

      do i = 1, nlimGammaHeI
         read(22,*) nuGammaHeI(i), (logGammaHeI(itk,i), itk=1, nTkGamma)         
         call locate(nuArray(1:nbins), nuGammaHeI(i), HeINuEdgeP(i))
         if (HeINuEdgeP(i)<nbins) then
            if (nuGammaHeI(i) > (nuArray(HeINuEdgeP(i))+nuArray(HeINuEdgeP(i)+1))/2.) &
                 & HeINuEdgeP(i) = HeINuEdgeP(i)+1
         end if
      end do

      do itk = 1, ntkgamma
         do i = 1, nlimGammaHeI
            if (logGammaHeI(itk,i)>0.) then         
               logGammaHeI(itk,i) = log10(logGammaHeI(itk,i))
            else
               logGammaHeI(itk,i) = 0.
            end if
         end do
      end do

      do i = 1, nlimGammaHeII
         read(23,*) nuGammaHeII(i), (logGammaHeII(itk,i), itk=1, nTkGamma)         
         call locate(nuArray(1:nbins), nuGammaHeII(i), HeIINuEdgeP(i))
         if (HeIINuEdgeP(i)<nbins) then
            if (nuGammaHeII(i) > (nuArray(HeIINuEdgeP(i))+nuArray(HeIINuEdgeP(i)+1))/2.) &
                 & HeIINuEdgeP(i) = HeIINuEdgeP(i)+1
         end if
      end do

      do itk = 1, ntkgamma
         do i = 1, nlimGammaHeII
            if (logGammaHeII(itk,i)>0.) then         
               logGammaHeII(itk,i) = log10(logGammaHeII(itk,i))
            else
               logGammaHeII(itk,i) = 0.
            end if
         end do
      end do

      close(21)
      close(22)
      close(23)

      HINuEdgeP(1) = 1      
      HeINuEdgeP(1) = 1
      HeIINuEdgeP(1) = 1

      print*, 'out initGammaCont'

    end subroutine initGammaCont

    subroutine initXSecArray()
        implicit none

        ! local variables
        integer :: err                        ! allocation error status
        integer :: i                          ! counter
        integer :: n                          ! principle QN
       

        real :: thomson = 6.65e-25            ! thomson x section
        real :: thres                         ! energy threshold [eV]
        real :: xSec                          ! x Section
        real :: x                             ! log10(W/W0) (for use in phFitHIon)
        real :: z                             ! Z

        print*, 'in initXSecArray'  
        
        ! allocate memory for xSecArrayTemp
        allocate (xSecArrayTemp(1:1000000), stat = err)
        if (err /= 0) then
            print*, "! initXSecArray: Cannot allocate grid array"
            stop
        end if


        ! zero out arrays
        xSecArrayTemp = 0.

        ! initialize top of stack pointer, xSecTop
        xSecTop = 0

        if (lgGas) then

           ! Hydrogen Lyman continuum

           HlevXSecP(1) = xSecTop + 1                    ! set pointer to the Lyman continuum        

           ! set pointers on the energy array nuArray
           call setPointers()

           do i = HlevNuP(1), nbins
              thres = max(nuArray(i)*RydToeV, ph1(1, 1, 1, 1)) 
              call phFitEl(1, 1, 1, thres, xSec)
              xSecArrayTemp(i-HlevNuP(1)+HlevXSecP(1)) = xSec*1e-18   
           end do
           xSecTop = xSecTop + nbins - HlevNuP(1) + 1

           ! Balmer, Paschen etc. continua (up to n= 10)
           do n = 2, nHlevel        
              HlevXSecP(n) = xSecTop + 1               ! set the pointers to the various continua
              do i = HlevNuP(n), HlevNuP(1)
                 x = log10(nuArray(i)) - log10(nuArray(HlevNuP(n)))           
                 xSecArrayTemp(i-HlevNuP(n) + HlevXSecP(n)) = phFitHIon(x, n, 1.)
              end do
              xSecTop = xSecTop + HlevNuP(1) - HlevNuP(n) +1
              
           end do

           ! free-free opacity (bremstrahlung)
           bremsXSecP = xSecTop+1
           do i = 1, nbins
              ! factor of 1e-20 missing to avoid underflow
              ! free-free opacity needs g(ff)*(1-exp(hn/kT))/sqrt(T)*1e-20
              xSecArrayTemp(i-1+bremsXSecP) = 1.03680e-18/(nuArray(i)**3)
           end do
           xSecTop = xSecTop + nbins -1 +1
        
           ! HeI singlet neutral Helium ground
           HeISingXSecP(1) = xSecTop + 1
           do i = HeIlevNuP(1), nbins
              call phFitEl(2, 2, 1, nuArray(i)*RydToeV, xSec)
              xSecArrayTemp(i-HeIlevNuP(1)+HeISingXSecP(1)) = xSec*1e-18
           end do
           xSecTop = xSecTop + nbins - HeIlevNuP(1) + 1      
           
           ! HeI neutral He 21S from StewartJPhysB 11, L431
           call powLawXSec(HeIlevNuP(2), HeIlevNuP(1), 0.4*8.7e-18, 1.5, HeISingXSecP(2))

           ! HeI singlet bound free excited states
           z = 1.
           do n = 3, nHeIlevel
              xSec = 7.906e-18*float(n) / (z*z)
              call powLawXSec(HeIlevNuP(n), HeIlevNuP(1), xSec, 3.0, HeISingXSecP(n))
           end do
           
           ! HeII ionized Helium bound free ground state
           HeIIXSecP(1)  = xSecTop + 1
           do i = HeIIlevNuP(1), nbins
              thres = max(nuArray(i)*RydToeV, ph1(1,2,1,1))
              call phFitEl(2, 1, 1, thres, xSec)
              xSecArrayTemp(i-HeIIlevNuP(1)+HeIIXSecP(1)) = xSec*1e-18
           end do
           xSecTop = xSecTop + nbins - HeIIlevNuP(1) +1
           
           ! HeII ionized Helium bound free excited state
           z = 2.
           do n = 2, nHeIIlevel
              xSec = 7.906e-18*float(n) / (z*z)
              HeIIXSecP(n) = xSecTop + 1
              do i = HeIIlevNuP(n), HeIIlevNuP(1)
                 x = log10(nuArray(i)) - log10(nuArray(HeIIlevNuP(n)))
                 xSecArrayTemp(i-HeIIlevNuP(n)+HeIIXSecP(n)) = phFitHIon(x, n, 2.)
              end do
              xSecTop = xSecTop + HeIIlevNuP(1) - HeIIlevNuP(n) + 1
           end do

           do i = 3, nElements
              if (lgElementOn(i)) call makeOpacity(i)
           end do
           
           ! set up dust constants and calculate dust extinction cross-sections
           if (lgDust) call makeDustXsec()

        else if (lgDust) then

           call makeDustXsec()

        else

           print*, '! initXSecArray: no gas or dust! the grid is empty.'
           stop

        end if
           

        ! allocate just enough space to XSecArray
        allocate (xSecArray(1:XSecTop), stat = err)
        if (err /= 0) then
           print*, "! initXSecArray: Cannot allocate grid array"
           stop 
        end if

        ! assign the values to XSecArray
        do i=1, XSecTop
           xSecArray(i) = xSecArrayTemp(i)
        end do

        ! free the space occupated by XSecArrayTemp
        if (allocated(XSecArrayTemp)) deallocate(XSecArrayTemp)

        print*, 'out initXSecArray' 

      end subroutine initXSecArray
        

       ! this subroutines reads the ph1 and ph2 parameters needed by the subroutine phFitEl
       ! from the data/ph1.dat and the data/ph2.dat files. also it the level, nInn and 
       ! nTot arrays are initialised here 
       subroutine phInit()
       
         ! local variables 
         integer :: ios             ! I/O error status
         integer :: i, j, k , l     ! counters

         ! initialise arrays
         level = (/0, 0, 1, 0, 1, 2, 0/)
         nInn  = (/0,0,1,1,1,1,1,1,1,1,3,3, &
              &  3,3,3,3,3,3,5,5,5,5,5,5,5,5,5,5,5,5/)
         nTot  = (/1,1,2,2,3,3,3,3,3,3,4,4, &
              &  5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7/)

         ! open data/ph1.dat file and check for errors in opening 
         close(10)
         open(unit = 10, file = "data/ph1.dat", status = "old", position = "rewind", iostat = ios)
         if (ios /= 0) then
            print*, "! phInit: can't open file: data/ph1.dat"
            stop
         end if

         ! open data/ph2.dat file and check for errors in opening  
         close(20)
         open(unit = 20, file = "data/ph2.dat", status = "old", position = "rewind", iostat = ios)
         if (ios /= 0) then
            print*, "! phInit: can't open file: data/ph2.dat"
            stop
         end if

         ! read data/ph1.dat file and check for errors in reading
         do l = 1, 7
            do k = 1, 30
               do j = 1, 30
                  read(unit = 10, fmt = *) (ph1(i, j, k, l), i=1, 6)
               end do
            end do
         end do

         ! read data/ph2.dat file and check for errors in reading    
         do k = 1, 30
            do j = 1, 30 
               read(unit = 20, fmt = *)  (ph2(i, j, k), i=1, 7)
            end do
         end do

         close(10)
         close(20)

       end subroutine phInit

       ! this subroutine calculates partial photoionization cross sections for
       ! all ionization stages of all atoms from H to Zn (Z=30) using the following
       ! fit parameters:
       ! outer shell of the Opacity Project (OP) elements: Verner, Ferland, Korista
       !    Yakovlev, 1996, ApJ, 465, 487
       ! inner shells of all elements, and outer shells of the non OP elements:
       !    Verner and Yakovlev, 1995, A&AS, 109, 125
       ! shell numbers:
       ! 1 - 1s; 2 - 2s; 3 - 2p; 4 - 3s; 5 - 3p; 6 - 3d; 7 - 4s  
       subroutine phFitEl(nz, ne, shell, photEn, xSec)
         implicit none
         
         integer, intent(in) :: nz                ! atomic number from 1 to 30
         integer, intent(in) :: ne                ! # of electrons from 1 to nz
         integer, intent(in) :: shell             ! shell number
     
         real, intent(in) :: photEn               ! photon energy [eV]
         real, intent(out) :: xSec                ! cross section [Mb]

         ! local variables
         integer :: nOut, nInt                    ! # of outer or internal shell
 
         real :: eInn, p1, y, q, a, b, x, z

         ! initialize the cross section to 0.0
         xSec = 0.0
    
         ! xSec = 0.0 if atomic number < 1 or > 30
         if( (nz < 1) .or. (nz > 30) ) return
         ! xSec = 0.0 if # of electrons  < 1 or > than atomic number
         if( (ne < 1) .or. (ne > nz) ) return         

         ! set outer shell number
         nOut = nTot(ne)
         if ( (nz == ne) .and. (nz > 18) ) nOut = 7 
         if ( (nz == (ne+1)) .and. ( (nz == 20) .or. (nz == 21) .or. &
&             (nz == 22) .or. (nz == 25) .or. (nz == 26) ) ) nOut = 7
         if ( shell > nOut) return
         if ( (shell == 6) .and. ( (nz == 20) .or. (nz == 19) ) .and. &
&             (ne >= 19) ) return
         if ( photEn < ph1(1, nz, ne, shell) ) return

         ! set internal shell number
         nInt = nInn(ne)
    
         ! set eInn
         if ( (nz == 15) .or. (nz == 17) .or. (nz == 19) .or. &
&             ( (nz > 20) .and. (nz /= 26) ) ) then 
            eInn = 0.0
         else
            if (ne < 3) then 
               eInn = 1.e30
            else
               eInn = ph1(1, nz, ne, nInt)
            end if
         end if

         ! calculate the cross section
         if( (shell <= nInt) .or. (photEn >= eInn) ) then
            p1 = -ph1(5, nz, ne, shell)
            y = photEn / ph1(2, nz, ne, shell)
            q = -0.5*p1 - level(shell) - 5.5
            a = ph1(3, nz, ne, shell) * ( (y-1.0)**2 + ph1(6, nz, nz, shell)**2 )
            b = sqrt(y/ph1(4, nz, ne, shell)) + 1.
            xSec = a * y**q * b**p1
         else
            if( (shell < nOut) .and. (shell > nInt) .and. (photEn < eInn) ) return
            p1 = -ph2(4, nz, ne)
            q = -0.5*p1 - 5.5
            x = photEn / ph2(1, nz, ne) - ph2(6, nz, ne)
            z = sqrt(x*x + ph2(7, nz, ne)*ph2(7, nz, ne))
            a = ph2(2, nz, ne) * ((x-1.)*(x-1.) + ph2(5, nz, ne)*ph2(5, nz, ne)) 
            b = sqrt(z/ph2(3, nz, ne)) + 1.
            xSec = a * z**q * b**p1
         end if
    
       end subroutine phFitEl


       ! this function calculates photoionisation cross sections for hydrogenic ions
       ! routine was adapted from cloudy and originally obtained from Anil Pradhan, Nov90
       ! W = photon energy ; W0 = threshold photon energy
       function phFitHIon(x, n, z)
         implicit none
         
         real :: phFitHIon              ! cross section [1/cm^2]
         
         integer :: n                   ! principal QN
         
         real :: x                      ! log10(W/W0)
         real :: z                      ! nuclear charge

         ! local variables

         real, dimension(10) ::  a, &   ! fit parameters
&             b, c, d, e, f                               

         ! assign the the fit parameters

         a = (/-17.2004,-16.8584,-16.6670,-16.5339,-16.4319,-16.3491, &
&              -16.2795,-16.2194,-16.1666,-16.1195/)
         b = (/-2.6671, -2.8068, -2.8549, -2.8812, -2.8983, -2.9103, &
&              -2.9193, -2.9264, -2.9321, -2.9369/)
         c = (/-0.3228, -0.1323, -0.0931, -0.0735, -0.0612, -0.0527, &
&              -0.0465, -0.0418, -0.0380, -0.0350/)
         d = (/ 0.0608,  0.0224,  0.0224,  0.0204,  0.0181,  0.0163, &
&               0.0147,  0.0135,  0.0125,  0.0117/)
         e = (/-16.9991,-16.7709,-16.6188,-16.5012,-16.4069,-16.3289, &
&              -16.2624,-16.2046,-16.1536,-16.1078/)
         f = (/-3.1304, -3.0041, -2.9738, -2.9671, -2.9662, -2.9669, &
&              -2.9681, -2.9694, -2.9707, -2.9719/)

         ! check n is within the range [1-10]
         if ( (n > 10) .or. (n<1) ) then
            print*, "! phFitHIon: quantum number n out of range [1,10]"
            stop
         end if
  
         ! calculate cross sections
         if (x <= 1) then
            phFitHIon = a(n) + x*(b(n) + x*(c(n) + x*d(n)))
         else
            phFitHIon = e(n) + x*f(n)
         end if

         phFitHIon = 10.**phFitHIon
         phFitHIon = phFitHIon / (z*z)

       end function phFitHIon


       ! this subroutine calculates the photoionisation cross-section for
       ! excited  hydrogenic states
       subroutine HFit(h0, pow, n)
         implicit none
         
         integer, intent(in) :: n
         
         real, intent(out) :: h0        ! crossection
         real, intent(out) :: pow       ! power

         ! better than 1% fit for n>=30
         h0 = (1.3962688e-10 + 2.7479352e-9 * sqrt(float(n)) ) * &
&             (1.3962688e-10 + 2.7479352e-9 * sqrt(float(n)) ) 

         ! better than 1% fit for n>=30
         if(n<1000) then
            pow = 2.99377 + 1.11784e-5*n - 0.3382905/n - &
&            6.72051e-9*n*n
         else
            pow = 3.
         end if
       end subroutine hFit

       subroutine H2plusXSec(low, high, opP)
         implicit none

         integer, intent(in) :: low, high
         integer, intent(out) :: opP                       ! opacity pointer
         
         ! local variables
         integer, parameter :: nData = 10
         integer :: i, ixSec, j                            ! counters
         integer :: nPairs                                 ! # of data pair
         
         real ::  slope 				          ! slope for linear interpolation
         real, dimension(nData) :: H2plusAr
         real, dimension(nData/2) :: enerData,  &          ! energy and cross 
&                                    xSecData              ! section arrays

         opP = xSecTop + 1
         nPairs = nData/2

         ! the array H2plusAr has ordered pairs of elements.
         ! the first is the energy in eV and the second is the cross
         ! section in Mb
         H2plusAr = (/6.75,0.24 , 8.68,2.5, 10.54,7.1, 12.46,6.0, &
&                      14.28,2.7/)
 
         do i = 1, nPairs
            enerData(i) = H2plusAr( (i-1)*2 + 1) / 13.6
            xSecData(i) = H2plusAr(i*2) * 1e-18
         end do

         if ( enerData(1) > nuArray(low) ) then
            print*, "! H2plusOp: the entered opacity energy bandwidth &
&                       is not largee enough [low]"
            stop
        end if
        
        slope = ( xSecData(2) - xSecData(1) ) / ( enerData(2) - enerData(1) )
    
        ixSec = 1

        ! fill in the opacities. use linear interpolation
        do i = low, high
           if ( (nuArray(i) > enerData(ixSec)) .and. &
&                 (nuArray(i) <= enerData(ixSec+1)) ) then
              xSecArrayTemp(i - low + opP) = xSecData(ixSec) + &
&                                           slope * (nuArray(i) - enerData(ixSec))
           else
              ixSec = ixSec + 1
              if (ixSec + 1 > nPairs) then
                 print*, " H2plusOp: the entered opacity energy bandwidth &
&                              is not large enough [high]"
                 stop
              end if
              slope = ( xSecData(ixSec+1) - xSecData(ixSec) ) / &
&                         ( enerData(ixSec+1) - enerData(ixSec) )
              if ( (nuArray(i) > enerData(ixSec)) .and. &
&                     (nuArray(i) <= enerData(ixSec+1)) ) then
                 xSecArrayTemp(i - low + opP) = xSecData(ixSec) + &
&                                               slope * (nuArray(i) - enerData(ixSec))
              else
                 print*, "! H2plusOp: internal logical error.&
&                               the desired energy is not within the next energy bound"
                 stop
              end if
           end if
        end do

        ! set the pointer to the top of the xSec array
        xSecTop = xSecTop + high - low + 1
         
        ! check that xSecTop doesn't exceed xSecMax
        if ( xSecTop > xSecMax) then
            print*, "! H2plusOp: xSecMax exceeded"
            stop
        end if

      end subroutine H2plusXSec
 
      ! this subroutine generates an array of cross 
      ! sections using a simple power law fit
      subroutine powLawXSec(low, high, cross, s, xSecP)
        implicit none
    
        integer, intent(in)  :: low, high ! pointers to boundaries of freq region
        integer, intent(out) :: xSecP     ! pointer to xSec array in xSecArray 

        real, intent(in)     :: cross, s

        ! local variables
        integer :: i                      ! counter
        
        real :: thres                     ! threshold

        xSecP = xSecTop + 1               ! set x Sec pointer
        thres =  nuArray(low)             ! set threshold

        ! generate array of x section
        do i = low, high
            xSecArrayTemp(i-low+xSecP) = &
&                                       cross * (nuArray(i) / thres)**(-s)
        end do

        xSecTop = xSecTop + high - low + 1
        
        if ( xSecTop > xSecMax )  then
            print*, "! powLawXSec: xSecMax exceeded"
            stop
        end if
      end subroutine powLawXSec

      ! this subroutine generates ionic subshell opacities
      ! by calling the procedure phFitEl
      subroutine makeOpacity(nElem)
        implicit none

        integer, intent(in) :: nElem

        ! local variables
        integer :: binsNeed              ! number of continuum bins needed to store opacity
        integer :: i                     ! counter
        integer :: nElec                 ! number of bound electrons
        integer :: nIon                  ! ionic number
        integer :: shell                 ! shell(suit) number
    
        real    :: energy                ! energy [eV]
        real    :: xSec                  ! xSec [1/cm^2]
        

        do nIon = 1, min(nElem, nstages)
           ! find number of bound electrons
           nElec = nElem - nIon +1

           do shell = 1, nShells(nElem, nIon)
              
              elementP(nElem,nIon , shell, 3) = xSecTop+1
              ! find how many cont bins are needed to store opacity
              binsNeed = elementP(nElem,nIon , shell, 2) - &
&                           elementP(nElem,nIon , shell, 1) + 1

              ! check that there are enough cells left in xSecArrayTemp to store opacities
              if ( (xSecTop + binsNeed) > xSecMax) then
                 print*, "! makeOpacity: xSecMax exceeded"
                 stop
              end if
                 
              if ( elementP(nElem,nIon , shell, 1)>elementP(nElem,nIon , shell, 2) .and. &
&                    ( .not.(elementP(nElem,nIon , shell, 1)==2) .or. &
&                     .not.(elementP(nElem,nIon , shell, 2)==1) ) ) then

                 print*, "! makeOpacity: upper energy limit for this shell is smaller &
&                              than threshold [elem,ion,shell, highlim,thres]", nElem,nIon, &
&                              shell,nuArray(elementP(nElem,nIon,shell,2)),&
&                              nuArray(elementP(nElem,nIon,shell,1))
                 stop
              end if

              do i = elementP(nElem,nIon , shell, 1), elementP(nElem,nIon , shell, 2)
                 energy = max( nuArray(i)*RydToeV, ph1(1, nElem, nElec, shell) )
                 call phFitEl(nElem, nElec, shell, energy, xSec)
                 xSecArrayTemp(i-elementP(nElem,nIon , shell, 1)+&
&                                  elementP(nElem,nIon , shell, 3)) = xSec*1e-18
              end do

              if (binsneed>=1) then       
                 xSecTop = xSecTop + elementP(nElem,nIon , shell, 2) - &
&                               elementP(nElem,nIon , shell, 1) +1
              end if

           end do
        end do

      end subroutine makeOpacity

      ! optical constants data files needed
      ! note that the number of wavelength points must be the same for all files
      ! makes dust xsections [cm^2] ( pi a^2 Q )
      subroutine makeDustXsec()

        character(len = 15) :: textString, dustFileType
        character(len=50)   :: extinctionFile

        integer :: err ! allocation error status
        integer :: i, j, iwav, n, iSize ! counter
        integer :: ios ! I/O error status
        integer :: iP ! array pointer
        integer :: nSpec, ai ! counter
        integer :: nn, iskip ! counter
        integer :: Nwav, NwavOld ! # of wavelength points used
        integer :: nRadii ! # of radii in Qfile

        real :: intValue ! interpolated value
        real :: normWeight ! normalization constant for grain size weighting
        real :: value ! general value variable
        real, pointer :: da(:)
        real, pointer :: tmp1(:), tmp2(:), tmp3(:), agrain(:), tmp11(:,:), tmp22(:,:),&
             & tmp33(:,:), tmpWav(:), QaTemp(:), QsTemp(:), gTemp(:), temp1nbins(:,:), &
             & temp2nbins(:,:), temp3nbins(:,:)
        real, pointer :: wav(:) ! wavelength in [um]
        real, pointer :: gCos(:,:,:) ! phase function parameter
        real, pointer :: Cabs(:,:,:) ! abs cross-section [um^2] for each grain species and size
        real, pointer :: Csca(:,:,:) ! sca cross-section [um^2] for each grain species and size
        real, pointer :: CTabs(:) ! total abs cross-section [um^2] for grain mixture
        real, pointer :: CTsca(:) ! total sca cross-section [um^2] for grain mixture
        real, pointer :: Ere(:), Eim(:) 
        real, pointer :: temp(:)

        print*, 'in makeDustXSec'

        open (unit=10, file=dustFile(1), iostat = ios, status = 'old', position = 'rewind')
        if(ios/=0) then
           print*, '! makeDustXSec: cannot open file for reading ', dustFile(1)
           stop
        end if

        ! check how many species are used
        read (10, *) nSpecies        

        ! allocate abundances array for dust
        allocate (absOpacSpecies(1:nSpecies, 1:nbins), stat=err)
        if (err/=0) then
           print*, "! makeDustXsec: error allocation memory for absOpacSpecies array"
           stop
        end if
        absOpacSpecies=0.

        ! allocate abundances array for dust
        allocate (grainAbun(1:nSpecies), stat=err)
        if (err/=0) then
           print*, "! makeDustXsec: error allocation memory for grainAbun array"
           stop
        end if
        grainAbun=0.

        allocate(rho(1:nSpecies), stat = err)
        if (err /= 0) then
           print*, "! makeDustXsec: can't allocate rho memory"
           stop
        end if
        rho=0.
        allocate(grainVn(1:nSpecies), stat = err)
        if (err /= 0) then
           print*, "! makeDustXsec: can't allocate grainVn memory"
           stop
        end if
        grainVn=0.
        allocate(MsurfAtom(1:nSpecies), stat = err)
        if (err /= 0) then
           print*, "! makeDustXsec: can't allocate surfAtom memory"
           stop
        end if
        MsurfAtom=0

        allocate (grainLabel(1:nSpecies), stat=err)
        if (err/=0) then
           print*, "! makeDustXsec: error allocation memory for grainLabel array"
           stop
        end if
        grainLabel=""
        allocate (TdustSublime(1:nSpecies), stat=err)
        if (err/=0) then
           print*, "! makeDustXsec: error allocation memory for TdustSublime array"
           stop
        end if
        TdustSublime=0.        

        open (unit=11, file=dustFile(2), iostat = ios, status = 'old', position = 'rewind')
        if(ios/=0) then
           print*, '! makeDustXSec: cannot open file for reading ', dustFile(2)
           stop
        end if
        ! chekc how many sizes are used
        read (11, *) nSizes
        allocate (grainRadius(1:nSizes), stat=err)
        if (err/=0) then
           print*, "! makeDustXsec: error allocation memory for grainRadius array"
           stop
        end if
        grainRadius=0.
        allocate (da(1:nSizes), stat=err)
        if (err/=0) then
           print*, "! makeDustXsec: error allocation memory for da array"
           stop
        end if
        da=0.

        allocate (grainWeight(1:nSizes), stat=err)
        if (err/=0) then
           print*, "! makeDustXsec: error allocation memory for grainWeight array"
           stop
        end if
        grainWeight=0.
        normWeight =0.

        allocate (dustScaXSecP(0:nSpecies,1:nSizes), stat=err)
        if (err/=0) then
           print*, "! makeDustXsec: error allocation memory for dustScaXSecP array"
           stop
        end if
        dustScaXSecP = -1
        allocate (dustAbsXSecP(0:nSpecies,1:nSizes), stat=err)
        if (err/=0) then
           print*, "! makeDustXsec: error allocation memory for dustAbsXSecP array"
           stop
        end if
        dustAbsXSecP = -1

        do ai = 1, nSizes
           read(11,*) iskip, grainRadius(ai), grainWeight(ai)       
        end do

        if (nSizes>1) then
           da(1) = grainRadius(2)-grainRadius(1)                    
           do ai = 2, nSizes-1
              da(ai) = (grainRadius(ai+1)-grainRadius(ai-1))/2.              
           end do
           da(nSizes) = grainRadius(nSizes)-grainRadius(nSizes-1)
        end if
        normWeight=  0.
        do ai = 1, nSizes
           normWeight = normWeight+grainWeight(ai)*da(ai) 
        end do
        if (nSizes>1) then
           do ai = 1, nSizes
              grainWeight(ai) = (grainWeight(ai)*da(ai))/normWeight
              if (.not.grainWeight(ai)>=0.) then
                 print*, '! makeDustXSec : Invalid grain weight ', grainWeight(ai), ai
                 stop
              end if
           end do
        else if (nSizes==1) then
           grainWeight(1) = 1.
        else
           print*, '! makeDustXSec : Invalid nSizes', nSizes
           stop
        end if

        if (taskid == 0)  then
           print*, '! makeDustXSec : Size Distribution ' 
           print*, ' index, a [um], da [um], weight '
           do ai = 1, nSizes
              print*, ai, grainRadius(ai), da(ai), grainWeight(ai)
           end do
        end if

        do nSpec = 1, nSpecies
           
           read(10, *) extinctionFile, grainAbun(nSpec)
           
           open (unit=20, file=extinctionFile, iostat = ios, status = 'old', position = 'rewind')
           if(ios/=0) then
              print*, '! makeDustXSec: cannot open file for reading ', extinctionFile
              stop
           end if

           read (20, *) dustFileType

           select case (dustFileType)
           case ('nk')

              read (unit=20, fmt=*, iostat=ios) textString

              ! count number of frequency points
              Nwav = 0
              do
                 read (unit=20, fmt=*, iostat=ios) value
                 if (ios /= 0 ) exit
                 Nwav = Nwav+1
              end do
              rewind 20
              
              if (nSpec>1 .and. nWav/=nWavOld) then
                 print*, '! makeDustXSec: [warning] extinction files &
                      & do not have the same number of frequency points'
              end if

              allocate (wav(1:nWav), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for wav array"
                 stop
              end if
              allocate (tmp1(1:nWav), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for tmp1 array"
                 stop
              end if
              allocate (tmp2(1:nWav), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for tmp2 array"
                 stop
              end if
              allocate (tmp3(1:nWav), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for tmp3 array"
                 stop
              end if
           
              ! initialise arrays
              tmp1 = 0.
              tmp2 = 0.
              tmp3 = 0.
              wav = 0.
              
              if (nSpec == 1) then
                 ! allocate the pointers' memory
                 allocate (Csca(1:nSpecies,0:nSizes,1:nbins), stat=err)
                 if (err/=0) then
                    print*, "! makeDustXsec: error allocation memory for Csca array"
                    stop
                 end if
                 allocate (Cabs(1:nSpecies,0:nSizes,1:nbins), stat=err)
                 if (err/=0) then
                    print*, "! makeDustXsec: error allocation memory for Cabs array"
                    stop
                 end if
                 allocate (CTsca(1:nbins), stat=err)
                 if (err/=0) then
                    print*, "! makeDustXsec: error allocation memory for CTsca array"
                    stop
                 end if
                 allocate (CTabs(1:nbins), stat=err)
                 if (err/=0) then
                    print*, "! makeDustXsec: error allocation memory for CTabs array"
                    stop
                 end if
                 Cabs = 0.
                 Csca = 0.
                 CTabs = 0.
                 CTsca = 0.
              end if

              allocate (Ere(1:nWav), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for Ere array"
                 stop
              end if
              allocate (Eim(1:nWav), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for Eim array"
                 stop
              end if
              Ere = 0.
              Eim = 0.
              
              
              read(20, *) dustFileType
              read(20, *) grainLabel(nSpec), TdustSublime(nSpec),&
                & rho(nSpec), grainVn(nSpec), MsurfAtom(nSpec)
           
              do i = 1, Nwav
                 read (20,*) wav(i), Ere(i), Eim(i) 
              end do
              
              close(20)
           
              nWavOld = nWav
              
              ! reverse order of wav, Ere and Eim and map onto mocassin's grid
              do i = nWav, 1, -1
                 ! convert wav into energy [Ryd]
                 wav(i) = c/(wav(i)*fr1Ryd*1.e-4)
                 tmp1(nWav-i+1) = Ere(i)
                 tmp2(nWav-i+1) = Eim(i)
                 tmp3(nWav-i+1) = wav(i)
              end do
              wav=tmp3
                            
              if (associated(Ere)) deallocate(Ere)
              if (associated(Eim)) deallocate(Eim)
              
              allocate (Ere(1:nbins), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for Ere array - 2"
                 stop
              end if
              allocate (Eim(1:nbins), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for Eim array - 2"
                 stop
              end if
              Ere = 0.
              Eim = 0.
              
              ! map onto mocassin grid
              call linearMap(tmp1, wav, nwav, Ere, nuArray, nbins)
              call linearMap(tmp2, wav, nwav, Eim, nuArray, nbins)

              ! calculate efficiencies
              call getQs(Ere,Eim,Cabs(nSpec,1:nSizes,1:nbins),Csca(nSpec,1:nSizes,1:nbins))
              
              if (associated(wav)) deallocate(wav)
              if (associated(Ere)) deallocate(Ere)
              if (associated(Eim)) deallocate(Eim)
              if (associated(tmp1)) deallocate(tmp1)
              if (associated(tmp2)) deallocate(tmp2)
              if (associated(tmp3)) deallocate(tmp3)

              
           case ('Q')

              read (unit=20, fmt=*, iostat=ios) textString
              do i = 1, 2 
                 read(20,*) textString
              end do
              read(20,*) nRadii
              read(20,*) nwav
              read(20,*) textString
              
              allocate (wav(1:nWav), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for wav array"
                 stop
              end if
              allocate (agrain(1:nRadii), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for agrain array"
                 stop
              end if
              allocate (QaTemp(1:nWav), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for QaTemp array"
                 stop
              end if
              allocate (QsTemp(1:nWav), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for QaTemp array"
                 stop
              end if
              allocate (gTemp(1:nWav), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for QaTemp array"
                 stop
              end if
              allocate (temp1nbins(1:nRadii,1:nbins), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for temp1nbins array"
                 stop
              end if
              allocate (temp2nbins(1:nRadii,1:nbins), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for temp1nbins array"
                 stop
              end if
              allocate (temp3nbins(1:nRadii,1:nbins), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for temp1nbins array"
                 stop
              end if
              allocate (tmp11(1:nRadii,1:nWav), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for tmp11 array"
                 stop
              end if
              allocate (tmp22(1:nRadii,1:nWav), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for tmp22 array"
                 stop
              end if
              allocate (tmp33(1:nRadii,1:nWav), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for tmp33 array"
                 stop
              end if
              allocate (tmpWav(1:nWav), stat=err)
              if (err/=0) then
                 print*, "! makeDustXsec: error allocation memory for tmpWav array"
                 stop
              end if
                            
              ! initialise arrays
              QaTemp = 0.
              QsTemp = 0.
              gTemp = 0.
              temp1nbins = 0.
              temp2nbins = 0.
              temp3nbins = 0.
              tmp11 = 0.
              tmp22 = 0.
              tmp33 = 0.
              tmpWav = 0.
              wav = 0.
              agrain = 0.
              
              if (nSpec == 1) then
                 ! allocate the pointers' memory
                 allocate (gCos(1:nSpecies,0:nSizes,1:nbins), stat=err)
                 if (err/=0) then
                    print*, "! makeDustXsec: error allocation memory for Cabs array"
                    stop
                 end if
                 allocate (Csca(1:nSpecies,0:nSizes,1:nbins), stat=err)
                 if (err/=0) then
                    print*, "! makeDustXsec: error allocation memory for Csca array"
                    stop
                 end if
                 allocate (Cabs(1:nSpecies,0:nSizes,1:nbins), stat=err)
                 if (err/=0) then
                    print*, "! makeDustXsec: error allocation memory for Cabs array"
                    stop
                 end if
                 allocate (CTsca(1:nbins), stat=err)
                 if (err/=0) then
                    print*, "! makeDustXsec: error allocation memory for CTsca array"
                    stop
                 end if
                 allocate (CTabs(1:nbins), stat=err)
                 if (err/=0) then
                    print*, "! makeDustXsec: error allocation memory for CTabs array"
                    stop
                 end if
                 gCos = 0.
                 Cabs = 0.
                 Csca = 0.
                 CTabs = 0.
                 CTsca = 0.
              end if
           
              do i = 1, nRadii
                 do j = 1, nWav
                    read(20,*) agrain(i), wav(j), tmp11(i,j), tmp22(i,j), tmp33(i,j)
                 end do

                 ! reverse order of wav, tmp11, tmp22 and tmp33 and map onto mocassin's grid
                 do iwav = nWav, 1, -1
                    ! convert wav into energy [Ryd]
                    tmpWav(iwav) = c/(wav(iwav)*fr1Ryd*1.e-4)
                    QaTemp(iwav) = tmp11(i,iwav)
                    QsTemp(iwav) = tmp22(i,iwav)
                    gTemp(iwav) = tmp33(i,iwav)
                 end do
                 wav=tmpWav

                 ! map  data onto mocassin's nu grid 

                 call linearMap(QaTemp, wav, nwav, temp1nbins(i,:), nuArray, nbins)
                 call linearMap(QsTemp, wav, nwav, temp2nbins(i,:), nuArray, nbins)
                 call linearMap(gTemp, wav, nwav, temp3nbins(i,:), nuArray, nbins)                 

              end do

              close(20)

              if (associated(QaTemp)) deallocate(QaTemp) 
              if (associated(QsTemp)) deallocate(QsTemp) 
              if (associated(gTemp)) deallocate(gTemp) 
              if (associated(tmp11)) deallocate(tmp11)
              if (associated(tmp22)) deallocate(tmp22)
              if (associated(tmp33)) deallocate(tmp33)           
              if (associated(wav)) deallocate(wav)
              if (associated(tmpWav)) deallocate(tmpWav)           

              ! interpolate over the grain size distribution
              do i = 1, nSizes

                 call locate(agrain, grainRadius(i), iSize)

                 if (grainRadius(i)<agrain(1)) then
                    print*, '! makeDustXsec: [warning] Size distribution &
                         & extends to radii smaller than & 
                         & the smallest radius listed in the external Qfile provided'
                    print*, 'Value of the smallest radius was assigned'
                    Cabs(nSpec,i,:) =  temp1nbins(1,:)
                    Csca(nSpec,i,:) =  temp2nbins(1,:)
                    gCos(nSpec,i,:) =  temp3nbins(1,:)                 
                 else if (grainRadius(i)>agrain(nRadii)) then
                    print*, '! makeDustXsec: [warning] Size distribution & 
                         & extends to radii larger than & 
                         & the largest radius listed in the external Qfile provided'
                    print*, 'Value of the largest radius was assigned'

                    Cabs(nSpec,i,:) =  temp1nbins(nRadii,:)
                    Csca(nSpec,i,:) =  temp2nbins(nRadii,:)
                    gCos(nSpec,i,:) =  temp3nbins(nRadii,:)

                 else
                    do j =1, nbins

                       Cabs(nSpec,i,j) = temp1nbins(iSize,j)+&
                            & (grainRadius(i)-agrain(iSize))*&
                            & (temp1nbins(iSize+1,j)-temp1nbins(iSize,j))/&
                            & (agrain(iSize+1)-agrain(iSize))
                       
                       Csca(nSpec,i,j) = temp2nbins(iSize,j)+&
                            & (grainRadius(i)-agrain(iSize))*&
                            & (temp2nbins(iSize+1,j)-temp2nbins(iSize,j))/&
                            & (agrain(iSize+1)-agrain(iSize))
                       
                       gCos(nSpec,i,j) = temp3nbins(iSize,j)+&
                            & (grainRadius(i)-agrain(iSize))*&
                            & (temp3nbins(iSize+1,j)-temp3nbins(iSize,j))/&
                            & (agrain(iSize+1)-agrain(iSize))
                    end do
                 end if

              end do              

              if (associated(agrain)) deallocate(agrain)
              if (associated(temp1nbins)) deallocate(temp1nbins)
              if (associated(temp2nbins)) deallocate(temp2nbins)
              if (associated(temp3nbins)) deallocate(temp3nbins)           


           case default
              print*, "! makeDustXsec: invalid dustFileType", dustFileType,extinctionFile
              stop
           end select

        end do

        close(10)                

        print*, "! makeDustXsec: Grain Abundances: " 
        print*, "(index, label, abundance, sublimation T)"        
        do nSpec = 1, nSpecies
           print*, nSpec, grainLabel(nSpec),grainAbun(nSpec),TdustSublime(nSpec)
        end do
        print*, ""
        
        ! calculate cross sections [um^2]        
        ! combine individual species/sizes cross-sections into total cross-section for
        ! the grain mixture
        do i = 1, nbins
           do nSpec = 1, nSpecies
              do ai = 1, nSizes

                 Csca(nSpec,ai,i) = Csca(nSpec,ai,i)*Pi*grainRadius(ai)*grainRadius(ai)*1.e-8
                 Cabs(nSpec,ai,i) = Cabs(nSpec,ai,i)*Pi*grainRadius(ai)*grainRadius(ai)*1.e-8

                 CTsca(i) = CTsca(i) + grainAbun(nSpec)*Csca(nSpec,ai,i)*grainWeight(ai)
                 CTabs(i) = CTabs(i) + grainAbun(nSpec)*Cabs(nSpec,ai,i)*grainWeight(ai)

              end do
           end do
        end do


        do i = 1, nbins
           xSecArrayTemp(xSecTop+i) = CTsca(i)
           xSecArrayTemp(xSecTop+i+nbins) = CTabs(i)

           nn=2
           do nSpec = 1, nSpecies
              do ai = 1, nSizes

                 xSecArrayTemp(xSecTop+i+nn*nbins) = Csca(nSpec, ai, i)
                 nn=nn+1

                 xSecArrayTemp(xSecTop+i+nn*nbins) = Cabs(nSpec, ai, i)
                 nn=nn+1

              end do
           end do

        end do

        ! set up dust cross-section pointers
        ! individual species
        dustScaXsecP(0,:) = xSecTop+1
        dustAbsXsecP(0,:) = xSecTop+1+nbins

        nn = 2
        do n = 1, nSpecies
           do ai = 1, nSizes
              dustScaXsecP(n,ai) = xSecTop+1+nn*nbins

              nn = nn+1
              dustAbsXsecP(n,ai) = xSecTop+1+nn*nbins

              nn=nn+1
           end do
        end do
        ! update value of xSecTop
        xSecTop = xSecTop + 2*nbins*(nSpecies+1)*nSizes

        if (associated(Csca)) deallocate(Csca)
        if (associated(Cabs)) deallocate(Cabs)
        if (associated(CTsca)) deallocate(CTsca)
        if (associated(CTabs)) deallocate(CTabs)
        if (associated(gCos)) deallocate(gCos)

        do n = 1, nSpecies
           do i = 1, nbins
              do ai = 1, nSizes
                 absOpacSpecies(n, i) = absOpacSpecies(n, i) + xSecArrayTemp(dustAbsXsecP(n,ai)+i-1)*&
                      & grainWeight(ai)


              end do
           end do
        end do


      end subroutine makeDustXsec

      subroutine getQs(Ere_in,Eim_in,Qabs,Qsca)
        implicit none

        real, intent(in) :: Ere_in(*), Eim_in(*)
        real, intent(out) :: Qabs(nsizes,nbins), Qsca(nSizes, nbins)      
      
        complex :: refIndex
        real    :: sizeParam, Qback(nSizes, nbins)
        
        integer :: i, ai 
        
        do i = 1, nbins

           ! complex refraction index
           refIndex = cmplx(Ere_in(i),Eim_in(i))
           
           do ai = 1, nSizes

              ! size parameter
              sizeParam=2.0*3.14159265*grainRadius(ai)/&
                   &( 2.9979250e14/(nuArray(i)*fr1Ryd) )

              ! if size parameter > 100 use 100 (geometrical optics)
              if (sizeParam > 100.) sizeParam=100.

              ! now calculate the efficiencies
              call BHmie(sizeParam,refIndex,Qabs(ai,i),Qsca(ai,i),Qback(ai,i))
              Qabs(ai,i) = Qabs(ai,i) - Qsca(ai,i)
            
              if (.not.lgDustScattering) Qsca(ai,i)=0.

           end do
           
        end do

      end subroutine getQs

      ! ***********************************************************************
      ! This subroutine was originally obtained from prof. P. Menguc, Dept. of
      ! Mechanical
      ! Eng., University of Kentucky, by Z.Ivezic., Aug. 1996
      ! Ported to F90 and integrated into MOCASSIN by B. Ercolano, Sep 2004
      ! -----------------------------------------------------------------------
      !     __________________________________________________________________
      !
      !     SUBROUTINE BHMIE CALCULATES AMPLITUDE SCATTERING MATRIX ELEMENTS
      !     & EFFICIENCIES FOR EXTINCTION, TOTAL SCATTERING AND BACSCATTERING,
      !     FOR A GIrVEN SIZE PARAMETER AND RELATIVE REFRACTIVE INDEX
      !     __________________________________________________________________
      
      subroutine BHmie (x,refrel,qext,qsca,qback)
        implicit none
        
        complex, intent(in) :: refrel
        real, intent(in)    :: x
        
        real, intent(out)   :: qext, qsca, qback
        
        real, dimension(100) ::  amu, theta,pii,tau,pi0,pi1
        
        complex              :: d(3000),y,xi,xi0,xi1,an,bn,s1(200),s2(200)
        real                 ::  psi0,psi1,psi,dn,dx,xstop,ymod,dang,chi0,&
             &chi1,apsi0,apsi1,apsi,fn,chi,p,t
        
        integer              :: nang,j,nstop,nmx,i,n,nn,rn,jj
        
        nang = 2
        dx=x
        y=x*refrel
        
        ! series terminated after nstop terms
        ! ___________________________________________________________________
        xstop=x+4.*x**.3333 +2.0
        nstop=int(xstop)
        ymod=cabs(y)
        nmx=int(amax1(xstop,ymod)) + 15
        dang=1.570796327/float(nang-1)
        do j = 1,nang
           theta(j)= (float(j)-1.)*dang
           amu(j)=cos(theta(j))
        end do
        ! __________________________________________________________________
        ! logarithmic derivative d(j) calculated by downward recurrence
        ! beginning with initial value 0.0 + i*0.0 at j = nmx
        !__________________________________________________________________
        d(nmx)=cmplx(0.0,0.0)
        nn=nmx-1
        do n=1,nn
           rn=nmx-n+1
           d(nmx-n)=(rn/y)-(1./(d(nmx-n+1)+rn/y))
        end do
        do j=1,nang
           pi0(j)=0.0
           pi1(j)=1.0
        end do
      nn=2*nang-1
      
      do j=1,nn
         s1(j)=cmplx(0.0,0.0)
         s2(j)=cmplx(0.0,0.0)
      end do
      ! __________________________________________________________________
      ! riccati bessel functions with real argument x calculated by upward
      ! recurrence
      ! __________________________________________________________________
      !
      psi0=cos(dx)
      psi1=sin(dx)
      chi0=-sin(x)
      chi1=cos(x)
      apsi0=psi0
      apsi1=psi1
      xi0=cmplx(apsi0,-chi0)
      xi1=cmplx(apsi1,-chi1)
      qsca=0.0
      n=1
      do 
         dn=n
         rn=n
         fn=(2.*rn+1.)/(rn*(rn+1.))
         psi=(2.*dn-1.)*psi1/dx-psi0
         apsi=psi
         chi=(2.*rn-1.)*chi1/x -  chi0
         xi = cmplx(apsi,-chi)
         an=(d(n)/refrel+rn/x)*apsi - apsi1
         an=an/((d(n)/refrel+rn/x)*xi - xi1)
         bn=(refrel *d(n)+rn/x)*apsi - apsi1
         bn=bn/((refrel*d(n)+rn/x)*xi - xi1)
         qsca=qsca+(2.*rn+1.)*(cabs(an)*cabs(an)+cabs(bn)*cabs(bn))
         do j=1,nang
            jj=2*nang-j
            pii(j)=pi1(j)
            tau(j)=rn*amu(j)*pii(j) - (rn+1.)*pi0(j)
            p=(-1.)**(n-1)
            s1(j)=s1(j)+fn*(an*pii(j)+bn*tau(j))
            t=(-1.)**n
            s2(j)=s2(j) + fn*(an*tau(j)+bn*pii(j))
            if (j == jj) exit
            s1(jj)=s1(jj) + fn*(an*pii(j)*p + bn*tau(j)*t)
            s2(jj)=s2(jj) + fn*(an*tau(j)*t + bn*pii(j)*p)
         end do
         psi0=psi1
         psi1=psi
         apsi1=psi1
         chi0=chi1
         chi1=chi
         xi1=cmplx(apsi1,-chi1)
         n=n+1
         rn=n
         do j=1,nang
            pi1(j)=((2.*rn-1.)/(rn-1.))*amu(j)*pii(j)
            pi1(j)=pi1(j) - rn*pi0(j)/(rn-1.)
            pi0(j) = pii(j)
         end do
         if (n-1-nstop>=0) exit
      end do

      qsca=(2./(x*x))*qsca
      qext=(4./(x*x))*real(s1(1))
      qback=(4./(x*x))*cabs(s1(2*nang -1))*cabs(s1(2*nang -1))

    end subroutine BHmie

end module xSec_mod
       



