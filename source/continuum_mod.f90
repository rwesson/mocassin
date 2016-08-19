! Copyright (C) 2005 Barbara Ercolano 
!
! Version 2.02
module continuum_mod
    use constants_mod
    use common_mod
    use interpolation_mod

    ! common variables
    real, parameter      :: hcRyd_k = 157893.94    ! constant: h*cRyd/k (Ryd at inf used) [K]
    real, parameter      :: hcRyd = 2.1799153e-11  ! constant: h*c*Ryd (Ryd at inf used) [erg]
    real, parameter      :: cRyd = 3.2898423e15    ! constant: c*Ryd (Ryd at inf used) [Hz]


    real, pointer, save :: inSpectrumErg(:)        ! input specrum energy distribution [erg/(cm^2*s*Hz*sr)] 
    real, pointer, save :: inSpectrumPhot(:)       ! input specrum energy distribution [phot/(cm^2*s*Hz*sr)]
    real, pointer, save :: inSpectrumProbDen(:,:)  ! probability density for input spectrum (nstars,nbins)

    real, save          :: normConstantErg = 0.    ! normalization constant (area beyond input spectrum)
    real, save          :: normConstantPhot= 0.    ! normalization constant (area beyond input spectrum)    
    real                :: correctionPhot
    real,save :: RStar
               
    integer, parameter  :: maxLim = 1000000        ! max number of rows in the input spectrum file
    integer             :: nuP                     ! frequency pointer
    integer             :: lymanP                  ! frequency pointer to Lyman limit

    contains


    subroutine setContinuum()
        implicit none

        ! local variables

        character(len=30) :: filein          ! input file name

        integer :: enP                       ! pointer in enArray
        integer :: err                       ! allocation error status
        integer :: i, j, k, iStar, iloop     ! counters
        integer :: ios                       ! I/O error status
        integer :: numLam
        integer, parameter :: sb99nuLim=1221

        integer,dimension(nbins) :: lamCount
        real    :: SStar, skip, time         ! stellar surface [e36 cm^2]
        real,dimension(maxLim)    :: tmp1, tmp2

        real, dimension(maxLim) :: enArray  ! freq array as read from input spectrum file [Hz]
        real, dimension(maxLim)  :: Hflux    ! flux array as read from input spectrum file [erg/cm^2/s/Hz/sr]

  
        print*, 'in setContinuum', contShape

        ios = 0

        lamCount=0

        ! initialize arrays 
        allocate(inSpectrumErg(nbins), stat = err)
        if (err /= 0) then
            print*, "! setContinuum: can't allocate grid memory"
            stop
        end if
        allocate(inSpectrumPhot(nbins), stat = err)                
        if (err /= 0) then    
            print*, "setContinuum: can't allocate grid memory"    
            stop    
        end if      
        allocate(inSpectrumProbDen(nStars,nbins), stat = err)
        if (err /= 0) then
            print*, "setContinuum: can't allocate grid memory"
            stop
        end if

        ! find the Lyman limit
        call locate(nuArray, 1., lymanP)

        do iStar=1, nStars

           ! initialise arrays
           enArray           = 0.
           Hflux             = 0.
           inSpectrumErg     = 0.
           inSpectrumPhot    = 0.
           inSpectrumProbDen(iStar,:) = 0.          
           
           filein = contShape(iStar)           
                 
           ! open file for reading
           close(12)
           open(unit = 12, file=filein, status = 'old', position = 'rewind', iostat=err)

           if (err /= 0) then     

              ! then maybe a set continuum has been chosen
              print*, 'contShape ', contShape(iStar)
              do i = 1, nbins
                 inSpectrumErg(i) = getFlux(nuArray(i), Tstellar(iStar), contShape(iStar))
                 inSpectrumPhot(i) = inSpectrumErg(i)/ (nuArray(i)*hcRyd) 

              end do


           else

              select case(spID(iStar)) 
              case ('sb99')

                 ! Starburst99 input
                 ! NOTE: the file must be in the format of 5 columns
                 !       the first containing the time step
                 !       the second the lambda points (A) 
                 !       and the fourth containing the stellar spectrum points 
                 !       (erg/s/A/sr). The file must be in ascending
                 !       wavenegth order. 
                 !       First 7 lines are comments
                 ! check if the end of the file has been reached 
                 do j = 1, 7                   
                    read(unit=12, fmt=*) 
                 end do

                 do j = 1, maxLim                 

                    read(unit=12, fmt=*, iostat=ios) time
                    if (ios < 0) exit ! end of file reached

                    if (time == tStep(iStar)) then
                       backspace(12)
                       do k = 1, sb99nuLim                          
                          read(12, *) time, enArray(k), skip, Hflux(k)

                          ! change Log L(lambda) [erg/s/A] into L(nu) [erg/s/Hz]                         
                          tmp1(k) = (10.**(Hflux(k))*(1.e-8*enArray(k)**2.))/c
                    
                          ! change to Ryd 
                          tmp2(k) = (c/(enArray(k)*1.e-8))/cRyd

                       end do
                       exit
                    end if
                    
                 end do

                 close(12)

                 do i = 1,sb99nuLim
                    
                    Hflux(i) = tmp1(sb99nuLim-i+1)
                    enArray(i) = tmp2(sb99nuLim-i+1)
                    
                 end do
             
                 do j = 1, sb99nuLim
                    
                    if (enArray(j) >= nuArray(1)-widflx(1)/2. .and.&
                         & enArray(j) < nuArray(1)+widflx(1)/2.) then
                       inSpectrumErg(1) = inSpectrumErg(1)+Hflux(j)
                       lamCount(1) = lamCount(1)+1
                    end if
                    do i = 2, nbins-1
                       if(enArray(j)>=(nuArray(i-1)+nuArray(i))/2. .and.&
                            & enArray(j)<(nuArray(i+1)+nuArray(i))/2.) then
                          inSpectrumErg(i) = inSpectrumErg(i)+Hflux(j)
                          lamCount(i) = lamCount(i)+1
                       end if
                    end do
                    if(enArray(j)>=nuArray(nbins)-widflx(nbins)/2. .and.&
                         & enArray(j) < nuArray(nbins)+widflx(nbins)/2.) then
                       inSpectrumErg(nbins) = inSpectrumErg(nbins)+Hflux(j)
                       lamCount(nbins) = lamCount(nbins)+1
                    end if
                    
                 end do
                 

                 do i = 1, nbins
                    if (lamCount(i)>0) then

                       inSpectrumErg(i) = inSpectrumErg(i)/real(lamCount(i))
                    end if

                    inSpectrumPhot(i) = inSpectrumErg(i)/ (nuArray(i)*hcRyd)
                 end do
                 
                 do i = 2, nbins-1
                    if (inSpectrumErg(i) <= 0. .and. inSpectrumErg(i-1)>0.) then
                       do j = i+1, nbins
                          if (inSpectrumErg(j) > 0.) then
                             inSpectrumErg(i) = inSpectrumErg(i-1)-&
                                  & ( (nuArray(i)-nuArray(i-1))*&
                                  & (inSpectrumErg(i-1)-inSpectrumErg(j))/&
                                  & (nuArray(j)-nuArray(i-1)))
                             exit
                          end if
                       end do
                    end if
                 end do

              case ('rauch')
                 do j = 1, maxLim

                    ! check if the end of file has been reached
                    ! NOTE: the file must be in the format of two columns
                    !       the first containing the lambda points (A) 
                    !       and the second containing the input spectrum points 
                    !       (erg/cm^2/s/A/sr). The file must be in ascending
                    !       wavenegth order. See e.g. Thomas Rauch's tables 
                    !       http://astro.uni-tuebingen.de/~rauch/
                    
                    ! check if the end of the file has been reached 
                    read(unit=12, fmt=*, iostat=ios) 
                    if (ios < 0) exit ! end of file reached
                    
                    backspace(12)
                    read(12, *) enArray(j), Hflux(j)
                    
                    ! change Hflux(lambda) [erg/cm^"/s/A/sr] into Hflux(nu) [erg/cm^2/s/Hz]
                    
                    tmp1(j) = (Hflux(j)*(1.e-8*enArray(j)**2.))/c/4.
                    
                    ! change to Ryd 
                    tmp2(j) = (c/(enArray(j)*1.e-8))/cRyd
                    
                 end do
             
                 close(12)

                 do i = 1,j-1
                    
                    Hflux(i) = tmp1(j-1-i+1)
                    enArray(i) = tmp2(j-1-i+1)
                    
                 end do
             
                 numLam = j-1
                 
                 do j = 1, numLam
                    
                    if (enArray(j) >= nuArray(1)-widflx(1)/2. .and.&
                         & enArray(j) < nuArray(1)+widflx(1)/2.) then
                       inSpectrumErg(1) = inSpectrumErg(1)+Hflux(j)
                       lamCount(1) = lamCount(1)+1
                    end if
                    do i = 2, nbins-1
                       if(enArray(j)>=(nuArray(i-1)+nuArray(i))/2. .and.&
                            & enArray(j)<(nuArray(i+1)+nuArray(i))/2.) then
                          inSpectrumErg(i) = inSpectrumErg(i)+Hflux(j)
                          lamCount(i) = lamCount(i)+1
                       end if
                    end do
                    if(enArray(j)>=nuArray(nbins)-widflx(nbins)/2. .and.&
                         & enArray(j) < nuArray(nbins)+widflx(nbins)/2.) then
                       inSpectrumErg(nbins) = inSpectrumErg(nbins)+Hflux(j)
                       lamCount(nbins) = lamCount(nbins)+1
                    end if
                    
                 end do
                 
                 do i = 1, nbins
                    if (lamCount(i)>0) then
                       inSpectrumErg(i) = inSpectrumErg(i)/real(lamCount(i))
                    end if
                    inSpectrumPhot(i) = inSpectrumErg(i)/ (nuArray(i)*hcRyd)
                 end do
                 
                 ! get physical flux
                 inSpectrumErg = fourPi*inSpectrumErg
                 inSpectrumPhot = fourPi*inSpectrumPhot

              case default

                 print*, '! setContinuum: unsupported spectrum ID', i, spID(i)
                 stop

              end select
                 
           end if

          ! set the input spectrum probability density distribution
          call setProbDen(iStar)


       end do
        
       if (nStars==1) then
          if (LPhot > 0. ) then  ! calculate Lstar [e36 erg/s]
             
             RStar = sqrt(Lphot / (fourPi*normConstantPhot*cRyd))
              
              LStar(1) = fourPi*RStar*RStar*sigma*TStellar(1)**4.

              deltaE(1)=Lstar(1)/nPhotons(1)
           
              print*, "Q(H) = ", LPhot, " [e36 phot/s]"
              print*, "LStar= ", LStar(1), " [e36 erg/s]"
              print*, "RStar = ", RStar, " [e18 cm]"
              print*, "deltaE = ", deltaE, " [e36 erg/s]"

           else if (Lstar(1) > 0.) then ! calculate LPhot [e36 phot/s]
              
              deltaE(1)=Lstar(1)/nPhotons(1)

              RStar = sqrt(Lstar(1) / (fourPi*sigma*TStellar(1)**4.))
              
              LPhot = fourPi*RStar*RStar*normConstantPhot*cRyd
              
              print*, "RStar = ", RStar, " [e18 cm]"
              print*, "Q(H) = ", LPhot, " [e36 phot/s]"
              print*, "LStar = ", LStar(1), " [e36erg/s]"
              print*, "deltaE = ", deltaE, " [e36 erg/s]"              

           end if
           
        end if

        if (associated(inSpectrumErg)) deallocate(inSpectrumErg)
        if (associated(inSpectrumPhot)) deallocate(inSpectrumPhot)

        print*, 'out setContinuum'

      end subroutine setContinuum


    ! this function returns energy distribution [erg/(cm^2*s*Hz*sr)] for the 
    ! continuum specified by contShape at the specified cell  
    function getFlux(energy, temperature, cShape)

        real, intent(in) :: energy     ! energy [Ryd]
        real, intent(in) :: temperature! temperature [K]

        real :: constant               ! units constant
        real :: denominator            ! denominator
        real :: getFlux                ! erg/(cm^2*s*Hz*sr)

        character(len=50), intent(in) :: cShape


        ! calculate the input spectrum
        ! this is the Planck function divided by h (to avoid underflow)
        select case (cShape)
        case ("blackbody")
            constant = 0.5250229/hPlanck ! 2*fr1Ryd^3 / c^2 [erg/cm^2]

            ! Planck's law: B(T) = (2 h nu^3/c^2) / (exp(h nu / (k T)) - 1)
            ! transformaion of frequency into Ryd: nu = Ryd c nu [Ryd]

            if (hcRyd_k*energy/temperature > 86.) then 

               ! Wien distribution
               getFlux = constant*energy*energy*energy*&
                    & exp(-hcRyd_k*energy/temperature)
               return               
            end if
            
            denominator = (exp(hcRyd_k*energy/temperature)-1.)

            if (denominator <= 0.) then

               ! Rayleigh-Jeans distribution
               ! 2*fr1Ryd2*k/c^2 = 3.32154e-6
               getFlux = 3.32154e-6*energy*energy*temperature/hPlanck
               return
            end if

            getFlux = constant*energy*energy*energy/&
                    & denominator

            return

        case default
            print*, "! getFlux: invalid continuum shape", cShape
            stop
        end select 

    end function getFlux

    subroutine setProbDen(iS)
        implicit none

        ! local variables
        integer, intent(in) :: iS         ! central star index

        integer :: err                    ! allocation error status
        integer :: enP                    ! nuArray index 
        integer :: i                      ! counter
        integer :: nuP                    ! frequency pointer

        real :: aFit, bFit, cFit          ! fit parameters
        real :: term                      ! general calculations term
        real :: delNu                     ! frequency step in nuArray
        real :: RStar                     ! stellar radius [e18 cm]
        real :: SStar                     ! stellar surface [e36 cm^2]

        real, pointer :: inSpSumErg(:)    ! partial input spectrum sum [erg/s]
        real, pointer :: inSpSumPhot(:)   ! partial input spectrum sum [phot/s]


        allocate(inSpSumPhot(nbins), stat = err)
        if (err /= 0) then
            print*, "! setProbDens: can't allocate grid memory"
            stop
        end if
        inSpSumPhot = 0.

        allocate(inSpSumErg(nbins), stat = err)
        if (err /= 0) then
            print*, "! setProbDens: can't allocate grid memory"
            stop
        end if
        inSpSumErg = 0.

        do i = 1, nbins
           inSpSumErg(i)  =  inSpectrumErg(i)
           inSpSumPhot(i) =  inSpectrumPhot(i)
        end do

        ! calculate normalization constants 
        do i = 1, nbins
            normConstantErg   = normConstantErg  + inSpSumErg(i)*widFlx(i)   
        end do

        if (lgDust) then
           do i = 1, nbins
              normConstantPhot  = normConstantPhot + inSpSumPhot(i)*widFlx(i)
           end do
        else
           do i = lymanP, nbins
              normConstantPhot  = normConstantPhot + inSpSumPhot(i)*widFlx(i)
           end do
        end if

        ! normalise spectrum and calculate PDF
        inSpectrumProbDen(iS,1) = inSpSumErg(1)*widFlx(1)/normConstantErg        
        do i = 2, nbins
           inSpectrumProbDen(iS,i) = inSpectrumProbDen(iS,i-1) + &
&                  inSpSumErg(i)*widFlx(i)/normConstantErg 

        end do

        inSpectrumProbDen(iS,nbins) = 1.

        if (contShape(iS)=='blackbody') then
           normConstantErg = Pi*normConstantErg*hPlanck
           normConstantPhot = Pi*normConstantPhot*hPlanck
        end if

        if (associated(inSpSumErg)) deallocate(inSpSumErg)
        if (associated(inSpSumPhot)) deallocate(inSpSumPhot)

      end subroutine setProbDen

end module continuum_mod

