! Copyright (C) 2005 Barbara Ercolano 
! 
! Version 2.02
module photon_mod
    
    use common_mod
    use constants_mod
    use continuum_mod
    use grid_mod
    use interpolation_mod
    use pathIntegration_mod
    use vector_mod

    ! common variables 

    real :: Qphot = 0.

    integer     , parameter :: safeLim = 10000          ! safety limit for the loops

    type(vector), parameter :: origin=vector(0.,0.,0.)  ! origin of the cartesian grid axes

    contains

    subroutine energyPacketDriver(iStar, n, grid, plot)
        implicit none
        
        integer, intent(in)            :: n           ! number of energy packets 
        integer, intent(in)            :: iStar       ! central star index 
       
        type(plot_type), intent(inout), optional &
             & :: plot                                ! only used in the mocassinPlot version
        type(grid_type), dimension(:), intent(inout) :: grid        ! the grid(s)
        
        type(vector)                   :: posVector   ! initial position vector for dust emi

        ! local variables
        integer                        :: ian         ! angle counter
        integer                        :: ifreq       ! freq counter
        integer                        :: iview       ! viewing angle counter
        integer                        :: freqP       ! pointer to frequency
        integer                        :: i,j,k,iG,ii ! counters        
        integer                        :: iCell       ! cell counter
        integer                        :: igrid,ix,iy,iz ! location indeces
        integer                        :: ierr        ! allocation error status
        integer                        :: iPhot       ! counter
        integer                        :: plotNum     ! counter
        integer                        :: seedSize    ! pseudo random number generator seed 
        integer, dimension(nGrids)     :: inX,inY,inZ ! initial position indeces
        integer, pointer               :: seed(:)     ! seed array
        integer                        :: msec        ! millisecs of the sec
        integer                        :: dt(8)       ! date and time values

        integer                        :: totalEscaped  
        real                           :: absInt      ! total number of absorption events
        real                           :: scaInt      ! total number of scattering events
        real                           :: JDifTot     ! tot JDif
        real                           :: JsteTot     ! tot Jste
        real                           :: radius      ! radius
        


        character(len=7)               :: chTypeD     ! character type for driver

        type(vector)                   :: absPosition ! position of packet absorption

        print*, "in energyPacketDriver"

        ! initilize interactions stats counters
        if (iStar == 1) then
           absInt = 0.
           scaInt = 0.
        end if

        call date_and_time(values=dt)
        msec=dt(8)

        call random_seed(seedSize) 

        allocate(seed(1:seedSize), stat= ierr)
        if (ierr /= 0) then
            print*, "energyPacketDriver: can't allocate array memory: seed"
            stop
        end if

        seed = 0

        call random_seed(get = seed)
 
        seed = seed + msec + taskid

        call random_seed(put = seed)
        
        if (associated(seed)) deallocate(seed)

        Qphot = 0.
        
        do iPhot = 1, n
           chTypeD = "stellar"
           call energyPacketRun(chTypeD, starPosition(iStar))
        end do
        print*, 'Star: ', iStar
        print*, 'Qphot = ', Qphot

        if (lgDust.and.convPercent>=resLinesTransfer .and.&
             & .not.lgResLinesFirst&
             & .and. (.not.nIterateMC==1)) then

           print*, "! energyPacketDriver: starting resonance line packets transfer"


           iCell = 0
           do igrid = 1, nGrids

              do ix = 1, grid(igrid)%nx
                 do iy = 1, grid(igrid)%ny
                    do iz = 1, grid(igrid)%nz
                       iCell = iCell+1

                       if (mod(iCell-(taskid+1),numtasks)==0) then
                          if (grid(igrid)%active(ix,iy,iz)>0) then

                             do iPhot = 1, grid(igrid)%resLinePackets(grid(igrid)%active(ix,iy,iz))

                                chTypeD = "diffuse"
                                posVector%x = grid(igrid)%xAxis(ix)
                                posVector%y = grid(igrid)%yAxis(iy)
                                posVector%z = grid(igrid)%zAxis(iz)
                                
                                inX=-1
                                inY=-1
                                inZ=-1
                                inX(igrid)=ix
                                inY(igrid)=iy
                                inZ(igrid)=iz
                                
                                if (igrid>1) then
                                   ! check location on mother grid
                                   call locate(grid(grid(igrid)%motherP)%xAxis, &
                                        & posVector%x,inX(grid(igrid)%motherP))
                                   if (posVector%x > (grid(grid(igrid)%motherP)%xAxis(inX(grid(igrid)%motherP))+&
                                        & grid(grid(igrid)%motherP)%xAxis(inX(grid(igrid)%motherP)+1) )/2.) &
                                        & inX(grid(igrid)%motherP) = inX(grid(igrid)%motherP)+1                                   
                                   call locate(grid(grid(igrid)%motherP)%yAxis, &
                                        & posVector%y,inY(grid(igrid)%motherP))
                                   if (posVector%y > (grid(grid(igrid)%motherP)%yAxis(inY(grid(igrid)%motherP))+&
                                        & grid(grid(igrid)%motherP)%yAxis(inY(grid(igrid)%motherP)+1) )/2.) &
                                        & inY(grid(igrid)%motherP) = inY(grid(igrid)%motherP)+1                                   
                                   call locate(grid(grid(igrid)%motherP)%zAxis, &
                                        & posVector%z,inZ(grid(igrid)%motherP))
                                   if (posVector%z > (grid(grid(igrid)%motherP)%zAxis(inZ(grid(igrid)%motherP))+&
                                        & grid(grid(igrid)%motherP)%zAxis(inZ(grid(igrid)%motherP)+1) )/2.) &
                                        & inZ(grid(igrid)%motherP) = inZ(grid(igrid)%motherP)+1                                   
                                end if

                                call energyPacketRun(chTypeD,posVector,inX,inY,inZ,igrid)

                             end do
                          end if
                       end if

                    end do
                 end do
              end do

           end do

           print*, "! energyPacketDriver: ending resonance line packets transfer"

        end if


        print*, 'Qphot = ', Qphot

        ! evaluate Jste and Jdif
        ! NOTE : deltaE is in units of [E36 erg/s] however we also need a factor of
        ! 1.e-45 from the calculations of the volume of the cell hence these
        ! two factors cancel each other out giving units of [E-9erg/s] so we need to
        ! multiply by 1.E-9
        ! NOTE : Jste and Jdif calculated this way are in units of
        ! [erg sec^-1 cm^-2] -> no Hz^-1 as they are summed over separate bins (see
        ! Lucy A&A (1999)                                                                                                                                   

        print*, 'Lstar', Lstar(iStar)

        if (iStar==nStars) then
           totalEscaped=0        
           do iG = 1, nGrids
              do i = 0, grid(iG)%nCells
                 grid(iG)%Jste(i,:) = grid(iG)%Jste(i,:) * 1.e-9
                 
                 if (lgDebug) grid(iG)%Jdif(i,:) = grid(iG)%Jdif(i,:) * 1.e-9
                 do ifreq = 1, nbins
                    totalEscaped = totalEscaped+&
                         & grid(iG)%escapedPackets(i,ifreq, 0)
                    do ian = 0, nAngleBins
                       grid(iG)%escapedPackets(i,ifreq,ian) = grid(iG)%escapedPackets(i,ifreq,ian)                           
                    end do
                 end do
                 
                 if (lgSymmetricXYZ) then
                    grid(iG)%Jste(i,:) = grid(iG)%Jste(i,:)/8.
                    grid(iG)%escapedPackets(i,:,:) = grid(iG)%escapedPackets(i,:,:)/8.
                    if (lgDebug) grid(iG)%Jdif(i,:) = grid(iG)%Jdif(i,:)/8.
                 end if
                 
              end do
           end do
           
           

           if (lgDust) then
              print*, "! energyPacketDriver: [Interactions] : total -- abs -- sca: "
              print*, "! energyPacketDriver: [Interactions] ", absInt+scaInt, " -- ", &
                   &  absInt*100./(absInt+scaInt),"% -- ", &
                   &  scaInt*100./(scaInt+absInt),"%"
           end if

           print*, " total Escaped Packets :",  totalEscaped

        end if

        print*, "out energyPacketDriver"

        contains

        recursive subroutine energyPacketRun(chType, position, xP, yP, zP, gP)
            implicit none

            character(len=7), intent(in)     :: chType           ! stellar or diffuse?

            integer, optional, dimension(:), intent(in)    :: xP, yP, &
                 & zP                                            ! cartesian axes indeces 

            integer, optional, intent(in)    :: gP               ! grid index

            type(vector),intent(in), optional:: position         ! the position of the photon
        
            ! local variables

            type(photon_packet)              :: enPacket         ! the energu packet

            integer                          :: err              ! allocation error status
            integer                          :: i, j             ! counters  
            integer                          :: idirP, idirT     ! direction cosines

            real :: number 
            real, save :: ionPhot = 0.       

            ! create a new photon packet
            select case (chType)
            ! if the energy packet is stellar
            case ("stellar")
                ! check for errors in the sources position
                if (present(position) ) then
                    if( position /= starPosition(iStar) ) then
                        print*, "! energyPacketRun: stellar energy packet must&
                             & start at the stellar position"
                        stop
                    end if
                end if

                ! create the packet        
                enPacket = newPhotonPacket(chType)

            ! if the photon is diffuse
            case ("diffuse")

                ! check that the position has been specified
                if (.not.present(position)) then
                    print*, "! energyPacketRun: position of the new diffuse &
                         & energy packet has not been specified"
                    stop
                end if

                ! check also that axes indeces have been carried through (save time)
                if (.not.(present(xP).and.present(yP).and.present(zP))) then
                    print*, "! energyPacketRun: cartesian axes indeces of the new diffuse &
                         & energy packet has not been specified"
                    stop
                end if
                ! check also that grid index have been carried through 
                if (.not.(present(gP))) then
                    print*, "! energyPacketRun: cartesian axes indeces of the new diffuse &
                         & energy packet has not been specified"
                    stop
                end if

                ! create the packet

                enPacket = newPhotonPacket(chType, position, xP, yP, zP, gP)

            case ("dustEmi")
                 ! check that the position has been specified
                if (.not.present(position)) then
                    print*, "! energyPacketRun: position of the new dust emitted &
                         & energy packet has not been specified"
                    stop
                end if

                ! check also that axes indeces have been carried through (save time)
                if (.not.(present(xP).and.present(yP).and.present(zP))) then
                    print*, "! energyPacketRun: cartesian axes indeces of the new dust emitted &
                         &energy packet has not been specified"
                    stop
                end if
                ! check also that grud index have been carried through
                if (.not.(present(gP))) then
                    print*, "! energyPacketRun: cartesian axes indeces of the new dust emitted &
                         &energy packet has not been specified"
                    stop
                end if

                ! create the packet
                enPacket = newPhotonPacket(chType, position, xP, yP, zP, gP)

            end select 

            if (.not.lgDust .and. enPacket%nu < ionEdge(1) .and. .not.enPacket%lgLine) then

               ! the packet escapes without further interaction
                  
               idirT = int(acos(enPacket%direction%z)/dTheta)+1
               if (idirT>totangleBinsTheta) then
                  idirT=totangleBinsTheta
!                  print*, '! energyPacketRun: idir>totanglebins - error corrected', &
!                       & idir, totanglebins, enPacket%direction, dtheta
               end if
               if (idirT<1 .or. idirT>totAngleBinsTheta) then
                  print*, '! energyPacketRun: error in theta direction cosine assignment',&
                       &  idirT, enPacket, dTheta, totAngleBinsTheta
                  stop
               end if
              

               if (enPacket%direction%x<1.e-35) then
                  idirP = 0
               else
                  idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)             
               end if
               if (idirP<0) idirP=totAngleBinsPhi+idirP
               idirP=idirP+1

               if (idirP>totangleBinsPhi) then
                  idirP=totangleBinsPhi
!                  print*, '! energyPacketRun: idir>totanglebins - error corrected', &
!                       & idir, totanglebins, enPacket%direction, dtheta
               end if
             
               if (idirP<1 .or. idirP>totAngleBinsPhi) then
                  print*, '! energyPacketRun: error in phi direction cosine assignment',&
                       &  idirP, enPacket, dPhi, totAngleBinsPhi
                  stop
               end if

               if (nAngleBins>0) then
                  if (viewPointPtheta(idirT) == viewPointPphi(idirP).or. &
                       & (viewPointTheta(viewPointPphi(idirP))==viewPointTheta(viewPointPtheta(idirT))) .or. & 
                       & (viewPointPhi(viewPointPtheta(idirT))==viewPointPhi(viewPointPphi(idirP))) ) then
                     grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,viewPointPtheta(idirT)) = &
                          &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                          & enPacket%nuP,viewPointPtheta(idirT)) + deltaE(iStar)
                     if (viewPointPtheta(idirT)/=0) grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                          & enPacket%nuP,0) = &
                          & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                          & enPacket%nuP,0) +  deltaE(iStar)
                  else
                     grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                          & enPacket%nuP,0) = &
                          & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                          & enPacket%nuP,0) +  deltaE(iStar)
                  end if
                              
               else

                  grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                       & enPacket%nuP,0) = &
                       & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                       & enPacket%nuP,0) +  deltaE(iStar)
                  
               end if

               return
            end if

            ! if the new packet is capable of re-ionizing 
            if (.not.enPacket%lgLine) then

                ! compute the next segment of trajectory
                call pathSegment(enPacket)
                return

            else ! if the packet is a line packet 
                ! add to respective line packet bin

               if (lgDebug) &
                    & grid(gP)%linePackets(grid(gP)%active(enPacket%xP(gP), &
                    & enPacket%yP(gP), enPacket%zP(gP)), enPacket%nuP) = &
                    & grid(gP)%linePackets(grid(gP)%active(enPacket%xP(gP), &
                    & enPacket%yP(gP), enPacket%zP(gP)), enPacket%nuP) + deltaE(iStar)


            end if

        end subroutine energyPacketRun

        ! this function initializes a photon packet
        function initPhotonPacket(nuP,  position, lgLine, lgStellar, xP, yP, zP, gP)
            implicit none

            type(photon_packet)      :: initPhotonPacket  ! the photon packet
   
            real                     :: random            ! random number

            integer, intent(in)      :: nuP               ! the frequency of the photon packet
            integer, intent(in),dimension(:) :: xP, yP, &
                 & zP                                     ! indeces of position on the x, y and z axes            
            integer, intent(in)      :: gP                ! grid index


            logical, intent(in)      :: lgLine, lgStellar ! line, stellar packet?

            type(vector), intent(in) :: position          ! the position at which the photon
                                                          ! packet is created    
            ! local variables

            integer                  :: i                 ! counter

            initPhotonPacket%position = position

            initPhotonPacket%iG  = gP
            
            initPhotonPacket%nuP      = nuP       
           
            initPhotonPacket%lgStellar = lgStellar

            ! check if photon packen is line or continuum photon
            if ( lgLine ) then
                ! line photon
                initPhotonPacket%nu       = 0.
                initPhotonPacket%lgLine   = .true.
            else
                ! continuum photon
                initPhotonPacket%nu       = nuArray(nuP)
                initPhotonPacket%lgLine   = .false.          
            end if

            initPhotonPacket%xP  = -1 
            initPhotonPacket%yP  = -1
            initPhotonPacket%zP  = -1

            ! check is the photon is stellar or diffuse
            if (lgStellar) then
                
                initPhotonPacket%xP(gP)  = xP(1)
                initPhotonPacket%yP(gP)  = yP(1)
                initPhotonPacket%zP(gP)  = zP(1)

            else

                do i = 1, nGrids                                   
                   initPhotonPacket%xP(i)  = xP(i)
                   initPhotonPacket%yP(i)  = yP(i)
                   initPhotonPacket%zP(i)  = zP(i)
                end do

            end if

            ! cater for plane parallel ionization case
            if (initPhotonPacket%lgStellar .and. lgPlaneIonization) then
               
               ! get position

               ! x-direction
               call random_number(random)
               random = 1. - random
               initPhotonPacket%position%x = &
                    & -(grid(gP)%xAxis(2)-grid(gP)%xAxis(1))/2. + random*( &
                    & (grid(gP)%xAxis(2)-grid(gP)%xAxis(1))/2.+&
                    & (grid(gP)%xAxis(grid(gP)%nx)-grid(gP)%xAxis(grid(gP)%nx-1))/2.+&
                    & grid(gP)%xAxis(grid(gP)%nx))
               if (initPhotonPacket%position%x<grid(gP)%xAxis(1)) &
                    & initPhotonPacket%position%x=grid(gP)%xAxis(1)
               if (initPhotonPacket%position%x>grid(gP)%xAxis(grid(gP)%nx)) & 
                    initPhotonPacket%position%x=grid(gP)%xAxis(grid(gP)%nx)

               call locate(grid(gP)%xAxis, initPhotonPacket%position%x, initPhotonPacket%xP(gP))
               if (initPhotonPacket%xP(gP) < grid(gP)%nx) then
                  if (initPhotonPacket%position%x >= (grid(gP)%xAxis(initPhotonPacket%xP(gP))+&
                       & grid(gP)%xAxis(initPhotonPacket%xP(gP)+1))/2.) &
                       & initPhotonPacket%xP(gP) = initPhotonPacket%xP(gP)+1
               end if

               ! y-direction
               initPhotonPacket%position%y = 0.
               initPhotonPacket%yP(gP) = 1

               ! z-direction
               call random_number(random)
               random = 1. - random
               initPhotonPacket%position%z = -(grid(gP)%zAxis(2)-grid(gP)%zAxis(1))/2. + random*( &
                    & (grid(gP)%zAxis(2)-grid(gP)%zAxis(1))/2.+&
                    & (grid(gP)%zAxis(grid(gP)%nz)-grid(gP)%zAxis(grid(gP)%nz-1))/2.+&
                    & grid(gP)%zAxis(grid(gP)%nz))
               if (initPhotonPacket%position%z<grid(gP)%zAxis(1)) & 
                    & initPhotonPacket%position%z=grid(gP)%zAxis(1)
               if (initPhotonPacket%position%z>grid(gP)%zAxis(grid(gP)%nz)) & 
                    & initPhotonPacket%position%z=grid(gP)%zAxis(grid(gP)%nz)

              call locate(grid(gP)%zAxis, initPhotonPacket%position%z, initPhotonPacket%zP(gP))
               if (initPhotonPacket%zP(gP) < grid(gP)%nz) then               
                  if (initPhotonPacket%position%z >= (grid(gP)%xAxis(initPhotonPacket%zP(gP))+&
                       & grid(gP)%zAxis(initPhotonPacket%zP(gP)+1))/2.) initPhotonPacket%zP(gP) =& 
                       & initPhotonPacket%zP(gP)+1
               end if

               if (initPhotonPacket%xP(gP)<1) initPhotonPacket%xP(gP)=1             
               if (initPhotonPacket%zP(gP)<1) initPhotonPacket%zP(gP)=1

               ! direction is parallel to y-axis direction
               initPhotonPacket%direction%x = 0.
               initPhotonPacket%direction%y = 1.
               initPhotonPacket%direction%z = 0.

                planeIonDistribution(initPhotonPacket%xP(gP),initPhotonPacket%zP(gP)) = &
                     & planeIonDistribution(initPhotonPacket%xP(gP),initPhotonPacket%zP(gP)) + 1

             else

               ! get a random direction
               initPhotonPacket%direction = randomUnitVector()

            end if

            if ((lgSymmetricXYZ) .and. initPhotonPacket%lgStellar .and. nStars==1) then
                if (initPhotonPacket%direction%x<0.) &
                     & initPhotonPacket%direction%x = -initPhotonPacket%direction%x
                if (initPhotonPacket%direction%y<0.) &
                     & initPhotonPacket%direction%y = -initPhotonPacket%direction%y
                if (initPhotonPacket%direction%z<0.) &
                     & initPhotonPacket%direction%z = -initPhotonPacket%direction%z
            end if

            initPhotonPacket%origin(1) = gP
            initPhotonPacket%origin(2) = grid(gP)%active(initPhotonPacket%xP(gP),&
                 & initPhotonPacket%yP(gP), initPhotonPacket%zP(gP))

        end function initPhotonPacket

    
        ! this subroutine determines the frequency of a newly created photon packet
        ! according to the given probability density
        subroutine getNu(probDen, nuP)

            real, dimension(:), intent(in) :: probDen    ! probability density function
       
            integer, intent(out)           :: nuP         ! frequency index of the new

            ! local variables
            real                           :: random     ! random number

            ! get a random number
            call random_number(random)

            random = 1.-random

            ! see what frequency random corresponds to 
            call locate(probDen, random, nuP)
             if (nuP <= 0) nuP = 1

 !           if (probDen(nuP) /= random) then
 !              nuP = nuP+1               
 !           end if

             if (nuP<nbins) then
                if (random>=(probDen(nuP+1)+probDen(nuP))/2.) nuP=nuP+1
             end if

        end subroutine getNu

        
        ! this function creates a new photon packet
        function newPhotonPacket(chType, position, xP, yP, zP, gP)
    
            type(photon_packet)                :: newPhotonPacket! the photon packet to be created

            character(len=7), intent(in)       :: chType         ! stellar or diffuse?

            type(vector), intent(in), optional :: position       ! the position of the photon
                                                                 ! packet
            ! local variables
            integer                            :: nuP            ! the frequency index of the photon packet
            integer, dimension(1)              :: orX,orY,orZ    ! dummy
            integer, optional, dimension(:),intent(in) :: xP, yP, & 
                 & zP                                            ! cartesian axes indeces    
            integer, optional, intent(in)      :: gP             ! grid index
            logical                            :: lgLine_loc=.false.! line photon?

            real                               :: random         ! random number
                     
            select case (chType)

            ! if the photon is stellar
            case ("stellar")

                ! check for errors in the sources position
                if (present(position) ) then
                    if( position /= starPosition(iStar) ) then
                        print*, "! newPhotonPacket: stellar photon packet must&
                             & start at the stellar position"
                        stop
                    end if
                end if 
                if (present(gP) ) then
                    if( gP /= 1) then
                        print*, "! newPhotonPacket: stellar photon packet must&
                             & start from mother grid"
                        stop
                    end if
                end if 
                
                ! determine the frequency of the newly created photon packet
                call getNu(inSpectrumProbDen(iStar,:), nuP)
                
                if (nuP>nbins) then
                   print*, "! newPhotonPacket: insanity occured in stellar photon &
                        &nuP assignment (nuP,xP,yP,zP,activeP)", nuP, xP(gP),yP(gP),zP(gP), &
                        & grid(gP)%active(xP(gP),yP(gP),zP(gP))
                   print*, "inSpectrumProbDen: ",iStar,inSpectrumProbDen(iStar,:)
                   stop
                end if
                
                if (nuP < 1) then
                    print*, "! newPhotonPacket: insanity occured in stellar photon &
&                               nuP assignment"
                    stop
                end if

                ! initialize the new photon packet
                orX = starIndeces(iStar,1)
                orY = starIndeces(iStar,2)
                orZ = starIndeces(iStar,3)
                newPhotonPacket = initPhotonPacket(nuP, starPosition(iStar), .false., .true., orX,orY,orZ, 1)

                if (newPhotonPacket%nu>1.) then
                   Qphot = Qphot + deltaE(iStar)/(2.1799153e-11*newPhotonPacket%nu)
                end if

            ! if the photon is diffuse
            case ("diffuse")
               
                ! check that gas is present in the grid
                if (.not.lgGas) then
                   print*, "! newPhotonPacket: diffuse packet cannot be created in a no gas grid"
                   stop
                end if

                ! check that the position has been specified
                if (.not.present(position)) then
                    print*, "! newPhotonPacket: position of the new diffuse &
                         & photon packet has not been specified"
                    stop
                end if
            
                ! check that the position indeces have been specified
                if (.not.(present(xP).and.present(yP).and.present(zP))) then
                    print*, "! newPhotonPacket: position indeces of the new diffuse &
                         & photon packet have not been specified"
                    stop
                end if
                ! check that the grid indeces have been specified
                if (.not.(present(gP))) then
                    print*, "! newPhotonPacket: grid index of the new diffuse &
                         & photon packet has not been specified"
                    stop
                end if
 
                ! check that the position is not inside the inner region
                if (grid(gP)%active(xP(gP),yP(gP),zP(gP))<= 0) then                   
                    print*, "! newPhotonPacket: position of the new diffuse &
                         & photon packet is inside the inner region", xP(gP),yP(gP),zP(gP),gP
                    stop
                end if

                ! decide whether continuum or line photon
                call random_number(random)

                random = 1.-random

                if (random <= grid(gP)%totalLines(grid(gP)%active(xP(gP),yP(gP),zP(gP)))) then 
                   ! line photon
                   ! line photons escape so don't care which one it is unless debugging
                   if (lgDebug) then
                      call getNu( grid(gP)%linePDF(grid(gP)%active(xP(gP),yP(gP),zP(gP)),:), nuP )

                      if (nuP < 1) then
                         print*, "! newPhotonPacket: insanity occured in line photon &
                              & nuP assignment"
                         stop
                      end if
                           
                   else
                      
                      nuP = 0

                   end if

                   ! initialize the new photon packet
                   newPhotonPacket = initPhotonPacket(nuP, position, .true., .false., xP, yP, zP, gP)
                else 
                    ! continuum photon

                    call getNu(grid(gP)%recPDF(grid(gP)%active(xP(gP),yP(gP),zP(gP)),:), nuP)

                    if (nuP>nbins) then
                       print*, "! newPhotonPacket: insanity occured in diffuse photon &
                       & nuP assignment (nuP,xP,yP,zP,activeP)", nuP, xP(gP),yP(gP),zP(gP),&
                       &  grid(gP)%active(xP(gP),yP(gP),zP(gP))
                       print*, "recPDF: ", grid(gP)%recPDF(grid(gP)%active(xP(gP),yP(gP),zP(gP)),:)
                       stop
                    end if

                    if (nuP < 1) then
                        print*, "! newPhotonPacket: insanity occured in diffuse photon &
                             & nuP assignment", nuP, xP(gP),yP(gP),zP(gP), & 
                             & grid(gP)%active(xP(gP),yP(gP),zP(gP))
                       print*, "recPDF: ", grid(gP)%recPDF(grid(gP)%active(xP(gP),yP(gP),zP(gP)),:)
                        stop
                    end if

                    ! initialize the new photon packet
                    newPhotonPacket = initPhotonPacket(nuP, position, .false., .false., xP, yP, zP, gP)
                end if

            case ("dustEmi")
               ! check dust is present
               if (.not.lgDust) then
                  print*, "! newPhotonPacket: dust emitted packet cannot be created in a &
                       &no dust grid."
                  stop
               end if

               if (lgGas) then
                  print*, "! newPhotonPacket: dustEmi-type packet should be created in a &
                       & grid containing gas."
                  stop
               end if


               ! check that the position has been specified
               if (.not.present(position)) then
                  print*, "! newPhotonPacket: position of the new dust emitted &
                       &photon packet has not been specified"
                  stop
               end if

               ! check that the position indeces have been specified
               if (.not.(present(xP).and.present(yP).and.present(zP))) then
                  print*, "! newPhotonPacket: position indeces of the new dust emitted &
                       &photon packet have not been specified"
                  stop
               end if
               ! check that the position indeces have been specified
               if (.not.(present(gP))) then
                  print*, "! newPhotonPacket: grid index of the new dust emitted &
                       &photon packet has not been specified"
                  stop
               end if

               ! check that the position is not inside the inner region
               if (grid(gP)%active(xP(gP),yP(gP),zP(gP))<= 0) then
                  print*, "! newPhotonPacket: position of the new dust emitted &
                       &photon packet is inside the inner region", xP(gP),yP(gP),zP(gP)
                  stop
               end if

               call getNu(grid(gP)%dustPDF(grid(gP)%active(xP(gP),yP(gP),zP(gP)),:), nuP)
               
               if (nuP>nbins) then
                   print*, "! newPhotonPacket: insanity occured in dust emitted photon &
                       &nuP assignment (iphot, nuP,xP(gP),yP(gP),zP(gP),activeP)", iphot, &
                       & nuP, xP(gP),yP(gP),zP(gP), &
                       & grid(gP)%active(xP(gP),yP(gP),zP(gP))
                   print*, "dustPDF: ", grid(gP)%dustPDF(grid(gP)%active(xP(gP),yP(gP),zP(gP)),:)
                   print*, "grain T: ", grid(gP)%Tdust(:,0,grid(gP)%active(xP(gP),yP(gP),zP(gP)))
                   stop
                end if


               if (nuP < 1) then
                  print*, "! newPhotonPacket: insanity occured in dust emitted photon &
                       &nuP assignment", nuP,xP(gP),yP(gP),zP(gP), grid(gP)%active(xP(gP),yP(gP),zP(gP))
                  print*, "dustPDF: ", grid(gP)%dustPDF(grid(gP)%active(xP(gP),yP(gP),zP(gP)),:)
                  print*, "grain T: ", grid(gP)%Tdust(:, 0, grid(gP)%active(xP(gP),yP(gP),zP(gP)))
                  stop
               end if

               ! initialize the new photon packet
               newPhotonPacket = initPhotonPacket(nuP, position, .false., .false., xP, yP, zP, gP)


            ! if the photon packet type is wrong or missing
            case default
        
                print*, "! newPhotonPacket: wrong photon packet type - 'stellar', 'diffuse' and &
                     & dust emitted types allowed-"
                stop
            end select

        end function newPhotonPacket
            
        subroutine pathSegment(enPacket)
          implicit none

          type(photon_packet), intent(inout) :: enPacket ! the energy packet

          ! local variables
          type(vector)                    :: vHat     ! direction vector
          type(vector)                    :: rVec     ! position vector
          
          real                            :: absTau   ! optical depth
          real                            :: dlLoc    ! local displacement
          real                            :: dSx, dSy, dSz 
                                                      ! distances from x,y and z wall
          real                            :: dS       ! distance from nearest wall 
          real                            :: dV       ! lume of this cell
          real                            :: passProb ! prob of passing the next segment
          real                            :: probSca  ! prob that the packet scatters
          real                            :: radius   ! radius
          real                            :: random   ! random number
          real                            :: tauCell  ! local tau

          integer                         :: idirT,idirP ! direction cosine counters
          integer                         :: i, nS    ! counter
          integer                         :: xP,yP,zP ! cartesian axes indeces
          integer                         :: gP       ! grid index
          integer                         :: safeLimit =500000
                                                      ! safe limit for the loop

          character(len=7)                :: packetType ! stellar, diffuse, dustEmitted?

          logical                         :: lgScattered ! is the packet scattering with dust?
          logical                         :: lgReturn
          

          ! check that the input position is not outside the grid
          if ( (enPacket%iG <= 0).or.(enPacket%iG > nGrids) ) then   
             print*, "! pathSegment: starting position not in any defined gridhses",&
                  & enPacket
             stop
          else if ( (enPacket%xP(enPacket%iG) <= 0).or.&
               &(enPacket%xP(enPacket%iG) > grid(enPacket%iG)%nx) ) then
             print*, "! pathSegment: starting position in x is outside the grid",&
                  & enPacket
             stop
          else if ( (enPacket%yP(enPacket%iG) <= 0).or. & 
               & (enPacket%yP(enPacket%iG) > grid(enPacket%iG)%ny) ) then
             print*, "! pathSegment: starting position in y is outside the grid",&
                  & enPacket
             stop
          else if ( (enPacket%zP(enPacket%iG) <= 0).or.& 
               & (enPacket%zP(enPacket%iG) > grid(enPacket%iG)%nz) ) then   
             print*, "! pathSegment: starting position in z is outside the grid",&
                  & enPacket
             stop          
          end if

          ! define vHat and rVec
          rVec = enPacket%position
          vHat = enPacket%direction

          ! initialize xP, yP,zP
          xP = enPacket%xP(enPacket%iG)
          yP = enPacket%yP(enPacket%iG)
          zP = enPacket%zP(enPacket%iG) 
          gP = enPacket%iG

          ! initialise distance from walls
          dSx = 0.
          dSy = 0.
          dSz = 0.
          
          if (lg1D) then
             radius = 1.e10*sqrt((rVec%x/1.e10)*(rVec%x/1.e10) + &
                  &                               (rVec%y/1.e10)*(rVec%y/1.e10) + &
                  &                               (rVec%z/1.e10)*(rVec%z/1.e10))
             call locate(grid(1)%xAxis, radius, xP)
             if (nGrids > 1 .or. gP >1) then
                print*, " ! energyPacketRun: multiple grids are not allowed in a 1D simulation"
                stop
             end if
          end if

          ! initialize optical depth
          absTau = 0.

          ! get a random number
          call random_number(random)
                          
          ! calculate the probability 
          passProb = -log(1.-random)

          ! speed up photons that may be trapped
          if (lgPlaneIonization) then
             safeLimit=5000
          else
             safeLimit=500000
          end if

          do i = 1, safeLimit
             
             if (grid(gP)%active(xP,yP,zP)<0) then                

                ! packet is entering a subgrid
                enPacket%xP(gP) = xP
                enPacket%yP(gP) = yP
                enPacket%zP(gP) = zP

                gP = abs(grid(gP)%active(xP,yP,zP))

                ! where is the packet in the sub-grid?

                call locate(grid(gP)%xAxis, rVec%x, xP)
                if (xP==0) xP = xP+1
                if (xP< grid(gP)%nx) then
                   if (rVec%x >  (grid(gP)%xAxis(xP+1)+grid(gP)%xAxis(xP))/2.) &
                        & xP = xP + 1
                end if

                call locate(grid(gP)%yAxis, rVec%y, yP)
                if (yP==0) yP=yP+1
                if (yP< grid(gP)%ny) then
                   if (rVec%y >  (grid(gP)%yAxis(yP+1)+grid(gP)%yAxis(yP))/2.) &
                        & yP = yP + 1
                end if

                call locate(grid(gP)%zAxis, rVec%z, zP)
                if (zP==0) zP=zP+1
                if (zP< grid(gP)%nz) then
                   if (rVec%z >  (grid(gP)%zAxis(zP+1)+grid(gP)%zAxis(zP))/2.) &
                        & zP = zP + 1
                end if

             end if

             enPacket%iG = gP
             
             ! find distances from all walls

             if (vHat%x>0.) then
                if (xP<grid(gP)%nx) then
                   dSx = ( (grid(gP)%xAxis(xP+1)+grid(gP)%xAxis(xP))/2.-rVec%x)/vHat%x
                   if (abs(dSx)<1.e-10) then
                      rVec%x=(grid(gP)%xAxis(xP+1)+grid(gP)%xAxis(xP))/2.
                      xP = xP+1
                   end if
                else
                   dSx = ( grid(gP)%xAxis(grid(gP)%nx)-rVec%x)/vHat%x
                   if (abs(dSx)<1.e-10) then
                      rVec%x=grid(gP)%xAxis(grid(gP)%nx)
                      if (.not.lgPlaneIonization .and. gP==1) return
                   end if
                end if
             else if (vHat%x<0.) then
                if (xP>1) then
                   dSx = ( (grid(gP)%xAxis(xP)+grid(gP)%xAxis(xP-1))/2.-rVec%x)/vHat%x
                   if (abs(dSx)<1.e-10) then             
                      rVec%x=(grid(gP)%xAxis(xP)+grid(gP)%xAxis(xP-1))/2.
                      xP = xP-1
                   end if
                else
                   dSx = (grid(gP)%xAxis(1)-rVec%x)/vHat%x
                   if (abs(dSx)<1.e-10) then             
                      rVec%x=grid(gP)%xAxis(1)
                   end if
                end if
             else if (vHat%x==0.) then
                dSx = grid(gP)%xAxis(grid(gP)%nx)
             end if
             
             if (.not.lg1D) then 
                if (vHat%y>0.) then
                   if (yP<grid(gP)%ny) then
                      dSy = ( (grid(gP)%yAxis(yP+1)+grid(gP)%yAxis(yP))/2.-rVec%y)/vHat%y
                      if (abs(dSy)<1.e-10) then
                         rVec%y=(grid(gP)%yAxis(yP+1)+grid(gP)%yAxis(yP))/2.
                         yP = yP+1
                      end if
                   else
                      dSy = (  grid(gP)%yAxis(grid(gP)%ny)-rVec%y)/vHat%y
                      if (abs(dSy)<1.e-10) then
                         rVec%y=grid(gP)%yAxis(grid(gP)%ny)
                         if(gP==1) return
                      end if
                   end if
                else if (vHat%y<0.) then
                    if (yP>1) then
                       dSy = ( (grid(gP)%yAxis(yP)+grid(gP)%yAxis(yP-1))/2.-rVec%y)/vHat%y
                       if (abs(dSy)<1.e-10) then             
                          rVec%y=(grid(gP)%yAxis(yP)+grid(gP)%yAxis(yP-1))/2.
                          yP = yP-1
                       end if
                    else 
                       dSy = ( grid(gP)%yAxis(1)-rVec%y)/vHat%y
                       if (abs(dSy)<1.e-10) then             
                          rVec%y=grid(gP)%yAxis(1)
                       end if
                    end if
                 else if (vHat%y==0.) then
                    dSy = grid(gP)%yAxis(grid(gP)%ny)
                 end if

                 if (vHat%z>0.) then
                    if (zP<grid(gP)%nz) then
                       dSz = ( (grid(gP)%zAxis(zP+1)+grid(gP)%zAxis(zP))/2.-rVec%z)/vHat%z
                       if (abs(dSz)<1.e-10) then
                          rVec%z=(grid(gP)%zAxis(zP+1)+grid(gP)%zAxis(zP))/2.
                          zP = zP+1
                       end if
                    else
                       dSz = ( grid(gP)%zAxis(grid(gP)%nz)-rVec%z)/vHat%z
                       if (abs(dSz)<1.e-10) then
                          rVec%z=grid(gP)%zAxis(grid(gP)%nz)
                          if (.not.lgPlaneIonization .and. gP==1) return
                       end if
                    end if
                 else if (vHat%z<0.) then
                    if (zP>1) then             
                       dSz = ( (grid(gP)%zAxis(zP)+grid(gP)%zAxis(zP-1))/2.-rVec%z)/vHat%z
                       if (abs(dSz)<1.e-10) then             
                          rVec%z=(grid(gP)%zAxis(zP)+grid(gP)%zAxis(zP-1))/2.
                          zP = zP-1
                       end if
                    else
                       dSz = ( grid(gP)%zAxis(1)-rVec%z)/vHat%z
                       if (abs(dSz)<1.e-10) then             
                          rVec%z=grid(gP)%zAxis(1)
                       end if
                   end if
                else if (vHat%z==0.) then
                   dSz = grid(gP)%zAxis(grid(gP)%nz)
                end if

             end if

             ! cater for cells on cell wall
             if ( abs(dSx)<1.e-10 ) dSx = grid(gP)%xAxis(grid(gP)%nx)
             if ( abs(dSy)<1.e-10 ) dSy = grid(gP)%yAxis(grid(gP)%ny)
             if ( abs(dSz)<1.e-10 ) dSz = grid(gP)%zAxis(grid(gP)%nz)

             ! find the nearest wall
             dSx = abs(dSx)
             dSy = abs(dSy)
             dSz = abs(dSz)

             if (dSx<=0.) then
                print*, '! pathSegment: [warning] dSx <= 0.'
                dS = amin1(dSy, dSz)
             else if (dSy<=0.) then
                print*, '! pathSegment: [warning] dSy <= 0.'
                dS = amin1(dSx, dSz)
             else if (dSz<=0.) then
                print*, '! pathSegment: [warning] dSz <= 0.'
                dS = amin1(dSx, dSy)
             else
                dS = amin1(dSx,dSy,dSz)
             end if

             ! this should now never ever happen
             if (dS <= 0.) then
                print*, 'pathSegment: dS <= 0'
                stop
             end if

             ! calculate the optical depth to the next cell wall 
             tauCell = dS*grid(gP)%opacity(grid(gP)%active(xP,yP,zP), enPacket%nuP)

             ! find the volume of this cell
             dV = getVolume(grid(gP), xP,yP,zP)

             ! check if the packets interacts within this cell
             if ((absTau+tauCell > passProb) .and. (grid(gP)%active(xP,yP,zP)>0)) then

                ! packet interacts
                
                ! calculate where within this cell the packet is absorbed
                dlLoc = (passProb-absTau)/grid(gP)%opacity(grid(gP)%active(xP,yP,zP), enPacket%nuP)

                ! update packet's position
                rVec = rVec + dlLoc*vHat

                if (lgSymmetricXYZ .and. gP==1) then
                   if ( rVec%x <= grid(gP)%xAxis(1) ) then
                      if (vHat%x<0.) vHat%x = -vHat%x
                      rVec%x = grid(gP)%xAxis(1)
                   end if
                   if ( rVec%y <= grid(gP)%yAxis(1) ) then
                      if (vHat%y<0.) vHat%y = -vHat%y
                      rVec%y = grid(gP)%yAxis(1)
                   end if
                   if ( rVec%z <= grid(gP)%zAxis(1) ) then
                      if (vHat%z<0.) vHat%z = -vHat%z
                      rVec%z = grid(gP)%zAxis(1)
                   end if
                end if

                ! add contribution of the packet to the radiation field
                if (enPacket%lgStellar) then
                   grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                        grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dlLoc*deltaE(iStar) / dV
                else ! if the energy packet is diffuse
                   if (lgDebug) then
                      grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                           & grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dlLoc*deltaE(iStar) / dV
                   else
                      grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                           & grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dlLoc*deltaE(iStar) / dV
                   end if
                end if                  


                ! check if the position within the cell is still within the outer radius
                if ( sqrt( (rvec%x/1.e10)**2. + (rvec%y/1.e10)**2. + (rvec%z/1.e10)**2.)*1.e10 >= R_out &
                     & .and. R_out > 0.) then
                   
                   
                   ! the packet escapes without further interaction
                   
                   idirT = int(acos(enPacket%direction%z)/dTheta)+1
                   if (idirT>totangleBinsTheta) then
                      idirT=totangleBinsTheta
                   end if
                   if (idirT<1 .or. idirT>totAngleBinsTheta) then
                      print*, '! energyPacketRun: error in theta direction cosine assignment',&
                           &  idirT, enPacket, dTheta, totAngleBinsTheta
                      stop
                   end if

                   if (enPacket%direction%x<1.e-35) then
                      idirP = 0
                   else
                      idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)             
                   end if
                   if (idirP<0) idirP=totAngleBinsPhi+idirP
                   idirP=idirP+1
                   if (idirP>totangleBinsPhi) then
                      idirP=totangleBinsPhi
                   end if
                   
                   if (idirP<1 .or. idirP>totAngleBinsPhi) then
                      print*, '! energyPacketRun: error in Phi direction cosine assignment',&
                           &  idirP, enPacket, dPhi, totAngleBinsPhi
                      stop
                   end if
               
               
                   if (nAngleBins>0) then
                      if (viewPointPtheta(idirT) == viewPointPphi(idirP).or. &
                           & (viewPointTheta(viewPointPphi(idirP))==viewPointTheta(viewPointPtheta(idirT))) .or. & 
                           & (viewPointPhi(viewPointPtheta(idirT))==viewPointPhi(viewPointPphi(idirP))) ) then
                         grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,viewPointPtheta(idirT)) = &
                              &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              & enPacket%nuP,viewPointPtheta(idirT)) +  deltaE(iStar)
                         if (viewPointPtheta(idirT)/=0) grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                              & enPacket%nuP,0) = &
                              & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              & enPacket%nuP,0) +  deltaE(iStar)
                      else
                         grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                          & enPacket%nuP,0) = &
                          & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                          & enPacket%nuP,0) +  deltaE(iStar)
                      end if
                   else
                  

                      grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                           & enPacket%nuP,0) = &
                           & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                           & enPacket%nuP,0) +  deltaE(iStar)
                      
                   end if
                   
                   return
                end if
                
 
                ! check if the packet is absorbed or scattered 
                if (lgDust) then

                   probSca = grid(gP)%scaOpac(grid(gP)%active(xP,yP,zP),enPacket%nuP)/&
                        & (grid(gP)%opacity(grid(gP)%active(xP,yP,zP),enPacket%nuP))
                   
                   call random_number(random)
                   
                   random = 1.-random

                   if (random > probSca) then
                      lgScattered = .false.         
                   else if (random <= probSca) then
                      lgScattered = .true.         
                   else
                      print*, '! pathSegment: insanity occured and scattering/absorption &
                           & decision stage.'
                      stop
                   end if

                   if (.not. lgScattered) then

                      absInt = absInt + 1.
                            
                      if (.not.lgGas) then

                         ! packet is absobed by the dust
                         packetType = "dustEmi"
                         exit                         

                      else

                         ! packet is absobed by the dust+gas
                         packetType = "diffuse"
                         exit    

                      end if


                   else

                      scaInt = scaInt + 1.                           
                      
                      do nS = 1, nSpecies
                         if (grainabun(nS)>0. .and. grid(gP)%Tdust(nS, 0, & 
                              & grid(gP)%active(xP,yP,zP))<TdustSublime(nS)) exit
                      end do
                      if (nS>7) then
                         print*, "! pathSegment: packet scatters with dust at position where all &
                              &grains have sublimed."
                         print*, xP,yP,zP, grid(gP)%active(xP,yP,zP), tauCell, absTau, passProb
                         stop
                      end if                      

                      ! packet is scattered by the grain
                         
                      ! calculate new direction
                      ! for now assume scattering is isotropic, when phase
                      ! function is introduced the following must be changed                         
                         
                      enPacket%xP(gP) = xP
                      enPacket%yP(gP) = yP
                      enPacket%zP(gP) = zP            
                      
                      enPacket = initPhotonPacket(enPacket%nuP, rVec, .false., .false., enPacket%xP(1:nGrids), &
                           & enPacket%yP(1:nGrids), enPacket%zP(1:nGrids), gP)
                      
                      
                      vHat%x = enPacket%direction%x
                      vHat%y = enPacket%direction%y
                      vHat%z = enPacket%direction%z
                      
                      ! initialize optical depth
                      absTau = 0.
                      
                      ! get a random number
                      call random_number(random)
                      
                      ! calculate the probability 
                      passProb = -log(1.-random)
                      
                   end if
                      
                else
                   absInt = absInt + 1.
                   
                   if (.not.lgGas) then
                      print*, "! pathSegment: Insanity occured - no gas present when no dust interaction"
                      stop
                   end if
                   
                   ! packet interacts with gas
                   packetType = "diffuse"
                   exit
                   
                end if
                   
                
             else
                   
                ! the packet is not absorbed within this cell
                ! add contribution of the packet to the radiation field
                
                if (enPacket%lgStellar) then
                   grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                        grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dS*deltaE(iStar) / dV
                else ! if the energy packet is diffuse
                   if (lgDebug) then
                      grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                           & grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dS*deltaE(iStar) / dV
                   else
                      grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                           & grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dS*deltaE(iStar) / dV
                   end if
                end if
                
                ! update absTau
                absTau = absTau+tauCell
                
                ! update packet's position
                rVec = rVec+dS*vHat
                
                ! keep track of where you are on mother grid
                if (gP>1) then
                   if (vHat%x>0.) then
                      if ( enPacket%xP(grid(gP)%motherP) < grid(grid(gP)%motherP)%nx ) then
                         if ( rVec%x > (grid(grid(gP)%motherP)%xAxis(enPacket%xP(grid(gP)%motherP))+& 
                          & grid(grid(gP)%motherP)%xAxis(enPacket%xP(grid(gP)%motherP)+1))/2. ) then
                            enPacket%xP(grid(gP)%motherP) = enPacket%xP(grid(gP)%motherP)+1
                         end if
                      else
                         if ( rVec%x > grid(grid(gP)%motherP)%xAxis(enPacket%xP(grid(gP)%motherP))) then
                            print*, '! pathSegment: insanity occured at mother grid transfer (x axis +)', & 
                                 & rVec%x, gP, grid(gP)%motherP
                            stop
                         end if
                      end if
                   else
                      if ( enPacket%xP(grid(gP)%motherP) > 1 ) then
                         if ( rVec%x <= (grid(grid(gP)%motherP)%xAxis(enPacket%xP(grid(gP)%motherP)-1)+& 
                              & grid(grid(gP)%motherP)%xAxis(enPacket%xP(grid(gP)%motherP)))/2. ) then
                            enPacket%xP(grid(gP)%motherP) = enPacket%xP(grid(gP)%motherP)-1
                         end if
                      else
                         if (rVec%x < grid(grid(gP)%motherP)%xAxis(1)) then
                            print*, '! pathSegment: insanity occured at mother grid transfer (x axis-)',&  
                                 & rVec%x, gP, grid(gP)%motherP
                            stop
                         end if
                      end if
                   end if
                   if (vHat%y>0.) then
                      if (  enPacket%yP(grid(gP)%motherP) < grid(grid(gP)%motherP)%ny ) then
                         if ( rVec%y > (grid(grid(gP)%motherP)%yAxis( enPacket%yP(grid(gP)%motherP))+& 
                              & grid(grid(gP)%motherP)%yAxis( enPacket%yP(grid(gP)%motherP)+1))/2. ) then
                            enPacket%yP(grid(gP)%motherP) =  enPacket%yP(grid(gP)%motherP)+1
                         end if
                      else
                         if ( rVec%y > grid(grid(gP)%motherP)%yAxis( enPacket%yP(grid(gP)%motherP))) then
                            print*, '! pathSegment: insanity occured at mother grid transfer (y axis +)',&
                                 & rVec%y, gP, grid(gP)%motherP
                            stop
                         end if
                      end if
                   else
                      if (  enPacket%yP(grid(gP)%motherP) > 1 ) then
                         if ( rVec%y <= (grid(grid(gP)%motherP)%yAxis( enPacket%yP(grid(gP)%motherP)-1)+& 
                              & grid(grid(gP)%motherP)%yAxis( enPacket%yP(grid(gP)%motherP)))/2. ) then
                            enPacket%yP(grid(gP)%motherP) =  enPacket%yP(grid(gP)%motherP)-1
                         end if
                      else
                         if (rVec%y < grid(grid(gP)%motherP)%yAxis(1)) then
                            print*, '! pathSegment: insanity occured at mother grid transfer (y axis -)', & 
                                 & rVec%y, gP, grid(gP)%motherP
                            stop
                         end if
                      end if
                   end if
                   if (vHat%z>0.) then
                      if (  enPacket%zP(grid(gP)%motherP) < grid(grid(gP)%motherP)%nz ) then
                         if ( rVec%z > (grid(grid(gP)%motherP)%zAxis( enPacket%zP(grid(gP)%motherP))+&
                              & grid(grid(gP)%motherP)%zAxis( enPacket%zP(grid(gP)%motherP)+1))/2. ) then
                            enPacket%zP(grid(gP)%motherP) =  enPacket%zP(grid(gP)%motherP)+1
                         end if
                      else
                         if ( rVec%z > grid(grid(gP)%motherP)%zAxis( enPacket%zP(grid(gP)%motherP))) then
                            print*, '! pathSegment: insanity occured at mother grid transfer (z axis +)', &
                                 & rVec%z, gP, grid(gP)%motherP
                            stop
                         end if
                      end if
                   else
                  if (  enPacket%zP(grid(gP)%motherP) > 1 ) then
                     if ( rVec%z <= (grid(grid(gP)%motherP)%zAxis( enPacket%zP(grid(gP)%motherP)-1)+&
                          & grid(grid(gP)%motherP)%zAxis( enPacket%zP(grid(gP)%motherP)))/2. ) then
                        enPacket%zP(grid(gP)%motherP) =  enPacket%zP(grid(gP)%motherP)-1
                     end if
                  else
                     if (rVec%z < grid(grid(gP)%motherP)%zAxis(1)) then
                        print*, '! pathSegment: insanity occured at mother grid transfer (z axis -)', &
                             & rVec%z, gP, grid(gP)%motherP
                        stop
                     end if
                  end if
               end if
            end if

            if (.not.lg1D) then
               if ( (dS == dSx) .and. (vHat%x > 0.)  ) then
                  xP = xP+1
               else if ( (dS == dSx) .and. (vHat%x < 0.) ) then
                  xP = xP-1
               else if ( (dS == dSy) .and. (vHat%y > 0.) ) then
                  yP = yP+1
               else if ( (dS == dSy) .and. (vHat%y < 0.) ) then 
                  yP = yP-1
               else if ( (dS == dSz) .and. (vHat%z > 0.) ) then
                  zP = zP+1
               else if ( (dS == dSz) .and. (vHat%z < 0.) ) then
                  zP = zP-1
               else
                  print*, '! pathSegment: insanity occured in dS assignement', dS,dSx,dSy,dSz
               end if
            else
               radius = 1.e10*sqrt((rVec%x/1.e10)*(rVec%x/1.e10) + &
                    & (rVec%y/1.e10)*(rVec%y/1.e10) + &
                    & (rVec%z/1.e10)*(rVec%z/1.e10))
               call locate(grid(gP)%xAxis, radius , xP)
            end if

            if(lgPlaneIonization) then
               lgReturn=.false.
               
               if ( rVec%y <= grid(gP)%yAxis(1)-grid(gP)%geoCorrY .or. yP<1) then
                  
                  ! the energy packet escapes this grid
                  if (gP==1) then	
                     yP=1
                     lgReturn=.true.
                  else if (gP>1) then
                     xP = enPacket%xP(grid(gP)%motherP)
                     yP = enPacket%yP(grid(gP)%motherP)
                     zP = enPacket%zP(grid(gP)%motherP)
                     gP = grid(gP)%motherP
                  else
                     print*, '! pathSegment: insanity occured - invalid gP', gP
                     stop
                  end if
                  
               end if
               
               if (rVec%y > grid(gP)%yAxis(grid(gP)%ny)+grid(gP)%geoCorrY .or. yP>grid(gP)%ny) then

                  if (gP==1) then	
                     ! the energy packet escapes
                     yP = grid(gP)%ny
                     lgReturn=.true.
                  else if (gP>1) then
                     xP = enPacket%xP(grid(gP)%motherP)
                     yP =  enPacket%yP(grid(gP)%motherP)
                     zP =  enPacket%zP(grid(gP)%motherP)
                     gP = grid(gP)%motherP                      
                  else
                     print*, '! pathSegment: insanity occured - invalid gP', gP
                     stop
                  end if
                  
               end if
               
               if ( (rVec%x <= grid(gP)%xAxis(1) .or. xP<1) .and. gP==1) then
                  xP=1
                  rVec%x = grid(gP)%xAxis(1)
                  vHat%x = -vHat%x
                  
               end if
                   
               if ( (rVec%x <= grid(gP)%xAxis(1)-grid(gP)%geoCorrX .or. xP<1) &
                    & .and. gP>1) then
                  
                  xP = enPacket%xP(grid(gP)%motherP)
                  yP =  enPacket%yP(grid(gP)%motherP)
                  zP =  enPacket%zP(grid(gP)%motherP)
                  gP = grid(gP)%motherP
                   
               end if
               
               if ( (rVec%x >=  grid(gP)%xAxis(grid(gP)%nx) &
                    & .or. xP>grid(gP)%nx) .and. gP==1 )then
                  xP = grid(gP)%nx
                  rVec%x = grid(gP)%xAxis(grid(gP)%nx)
                  vHat%x = -vHat%x
                  
               end if
               
               if ( (rVec%x >=  grid(gP)%xAxis(grid(gP)%nx)+grid(gP)%geoCorrX&
                    & .or. xP>grid(gP)%nx) .and.  gP>1) then
                  
                  xP = enPacket%xP(grid(gP)%motherP)
                  yP =  enPacket%yP(grid(gP)%motherP)
                  zP =  enPacket%zP(grid(gP)%motherP)
                  gP = grid(gP)%motherP
                      
               end if
               
               if ( (rVec%z <= grid(gP)%zAxis(1) .or.zP<1) .and. gP==1) then
                  zP=1
                  rVec%z = grid(gP)%yAxis(1)
                  vHat%z = -vHat%z
                  
               end if
               
               if ( (rVec%z <= grid(gP)%zAxis(1)-grid(gP)%geoCorrZ &
                    & .or.zP<1) .and. gP>1) then
                  
                  xP = enPacket%xP(grid(gP)%motherP)
                  yP =  enPacket%yP(grid(gP)%motherP)
                  zP =  enPacket%zP(grid(gP)%motherP)
                  gP = grid(gP)%motherP
                  
               end if
               
               if ( (rVec%z >=  grid(gP)%zAxis(grid(gP)%nz) .or. zP>grid(gP)%nz) &
                    & .and. gP==1) then
                  
                  zP = grid(gP)%nz
                  rVec%z = grid(gP)%zAxis(grid(gP)%nz)
                  vHat%z = -vHat%z
                  
               end if

               if ((rVec%z >=  grid(gP)%zAxis(grid(gP)%nz)+grid(gP)%geoCorrZ &
                    & .or. zP>grid(gP)%nz) .and. gP>1) then
                  
                  xP = enPacket%xP(grid(gP)%motherP)
                  yP =  enPacket%yP(grid(gP)%motherP)
                  zP =  enPacket%zP(grid(gP)%motherP)
                  gP = grid(gP)%motherP                      
                  
               end if
               
               if (lgReturn) then            

                  ! the packet escapes without further interaction
                  
                  idirT = int(acos(enPacket%direction%z)/dTheta)+1
                  if (idirT>totangleBinsTheta) then
                     idirT=totangleBinsTheta
!                  print*, '! energyPacketRun: idir>totanglebins - error corrected', &
!                       & idir, totanglebins, enPacket%direction, dtheta
                  end if
                  if (idirT<1 .or. idirT>totAngleBinsTheta) then
                     print*, '! energyPacketRun: error in theta direction cosine assignment',&
                          &  idirT, enPacket, dTheta, totAngleBinsTheta
                     stop
                  end if
               

                  if (enPacket%direction%x<1.e-35) then
                     idirP = 0
                  else
                     idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)             
                  end if
                  if (idirP<0) idirP=totAngleBinsPhi+idirP
                  idirP=idirP+1
                  
                  if (idirP>totangleBinsPhi) then
                     idirP=totangleBinsPhi
!                  print*, '! energyPacketRun: idir>totanglebins - error corrected', &
!                       & idir, totanglebins, enPacket%direction, dtheta
                  end if
                  if (idirP<1 .or. idirP>totAngleBinsPhi) then
                     print*, '! energyPacketRun: error in phi direction cosine assignment',&
                          &  idirP, enPacket, dPhi, totAngleBinsPhi
                     stop
                  end if


                  if (nAngleBins>0) then
                     if (viewPointPtheta(idirT) == viewPointPphi(idirP).or. &
                          & (viewPointTheta(viewPointPphi(idirP))==viewPointTheta(viewPointPtheta(idirT))) .or. & 
                          & (viewPointPhi(viewPointPtheta(idirT))==viewPointPhi(viewPointPphi(idirP))) ) then
                        grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,viewPointPtheta(idirT)) = &
                             &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                             & enPacket%nuP,viewPointPtheta(idirT)) +  deltaE(iStar)
                        if (viewPointPtheta(idirT)/=0) grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                             enPacket%nuP,0) = &
                             & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                             & enPacket%nuP,0) +  deltaE(iStar)
                        
                     else
                        grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                             enPacket%nuP,0) = &
                             & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                             & enPacket%nuP,0) +  deltaE(iStar)
                        
                     end if

                  else

                     grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                          enPacket%nuP,0) = &
                          & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                          & enPacket%nuP,0) +  deltaE(iStar)      
                  end if

                  return
               end if
               
            end if

                ! check if the path is still within the ionized region 
                radius = 1.e10*sqrt((rVec%x/1.e10)*(rVec%x/1.e10) + &
                     &                     (rVec%y/1.e10)*(rVec%y/1.e10) + &
                     &                     (rVec%z/1.e10)*(rVec%z/1.e10))

                if (.not.lgPlaneIonization) then
                     
                   if ( (abs(rVec%x) >= grid(gP)%xAxis(grid(gP)%nx)+grid(gP)%geoCorrX) .or.&
                        &(abs(rVec%y) >= grid(gP)%yAxis(grid(gP)%ny)+grid(gP)%geoCorrY) .or.&
                        &(abs(rVec%z) >= grid(gP)%zAxis(grid(gP)%nz)+grid(gP)%geoCorrZ) .or. &
                        & xP>grid(gP)%nx .or. yP>grid(gP)%ny .or. zP>grid(gP)%nz ) then

                      if ((gP==1) .or.  (radius >= R_out .and. R_out >= 0.)) then

                         if (xP > grid(gP)%nx) xP = grid(gP)%nx
                         if (yP > grid(gP)%ny) yP = grid(gP)%ny
                         if (zP > grid(gP)%nz) zP = grid(gP)%nz

                         ! the energy packet escapes



                         ! the packet escapes without further interaction
                  
                         idirT = int(acos(enPacket%direction%z)/dTheta)+1
                         if (idirT>totangleBinsTheta) then
                            idirT=totangleBinsTheta
!                  print*, '! energyPacketRun: idir>totanglebins - error corrected', &
!                       & idir, totanglebins, enPacket%direction, dtheta
                         end if
                         if (idirT<1 .or. idirT>totAngleBinsTheta) then
                            print*, '! energyPacketRun: error in theta direction cosine assignment',&
                                 &  idirT, enPacket, dTheta, totAngleBinsTheta
                            stop
                         end if
                         
                         if (enPacket%direction%x<1.e-35) then
                            idirP = 0
                         else
                            idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)             
                         end if
                         if (idirP<0) idirP=totAngleBinsPhi+idirP
                         idirP=idirP+1
      
                         if (idirP>totangleBinsPhi) then
                            idirP=totangleBinsPhi
!                  print*, '! energyPacketRun: idir>totanglebins - error corrected', &
!                       & idir, totanglebins, enPacket%direction, dtheta
                         end if
      
                         if (idirP<1 .or. idirP>totAngleBinsPhi) then
                            print*, '! energyPacketRun: error in phi direction cosine assignment',&
                                 &  idirP, enPacket, dPhi, totAngleBinsPhi
                            stop
                         end if

                         if (nAngleBins>0) then
                            if (viewPointPtheta(idirT) == viewPointPphi(idirP).or. &
                                 & (viewPointTheta(viewPointPphi(idirP))==viewPointTheta(viewPointPtheta(idirT))) .or. & 
                                 & (viewPointPhi(viewPointPtheta(idirT))==viewPointPhi(viewPointPphi(idirP))) ) then
                               grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,& 
                                    & viewPointPtheta(idirT)) = &
                                    &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                    & enPacket%nuP,viewPointPtheta(idirT)) +  deltaE(iStar)
                               if (viewPointPtheta(idirT)/=0) grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                                    enPacket%nuP,0) = &
                                    & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                    & enPacket%nuP,0) +  deltaE(iStar)
                            else
                               grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                                    enPacket%nuP,0) = &
                                    & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                    & enPacket%nuP,0) +  deltaE(iStar)
                            end if
                         else

                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                                 enPacket%nuP,0) = &
                                 & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,0) +  deltaE(iStar)
                         end if

                         !b2.005
                         return

                      else if (gP>1) then

                         xP = enPacket%xP(grid(gP)%motherP)
                         yP = enPacket%yP(grid(gP)%motherP)
                         zP = enPacket%zP(grid(gP)%motherP)
                         gP = grid(gP)%motherP

                      else
                         print*, '! pathSegment: insanity occured - invalid gP - ', gP
                         stop
                      end if
                   end if
                end if

                if (lgSymmetricXYZ .and. gP == 1) then
                   if (lgPlaneIonization) then
                      print*, '! pathSegment: lgSymmetric and lgPlaneionization flags both raised'
                      stop
                   end if

                   if ( rVec%x <= grid(gP)%xAxis(1) .or. xP<1) then
                      if (vHat%x<0.) vHat%x = -vHat%x 
                      xP = 1
                      rVec%x = grid(gP)%xAxis(1)
                   end if
                   if ( rVec%y <= grid(gP)%yAxis(1) .or. yP<1) then
                      if (vHat%y<0.) vHat%y = -vHat%y 
                      yP = 1
                      rVec%y = grid(gP)%yAxis(1)
                   end if
                   if ( rVec%z <= grid(gP)%zAxis(1) .or. zP<1) then
                      if (vHat%z<0.) vHat%z = -vHat%z 
                      zP=1
                      rVec%z = grid(gP)%zAxis(1)
                   end if

                end if

                if (gP>1) then
                   if ( ( (rVec%x <= grid(gP)%xAxis(1)-grid(gP)%geoCorrX &
                        &.or. xP<1) .and. vHat%x <=0.) .or. & 
                        & ( (rVec%y <= grid(gP)%yAxis(1)-grid(gP)%geoCorrY &
                        & .or. yP<1) .and. vHat%y <=0.) .or. & 
                        & ( (rVec%z <= grid(gP)%zAxis(1)-grid(gP)%geoCorrZ &
                        &  .or. zP<1) .and. vHat%z <=0.) .or. & 
                        & ( (rVec%x >= grid(gP)%xAxis(grid(gP)%nx)+grid(gP)%geoCorrX &
                        &.or. xP>grid(gP)%nx) .and. vHat%x >=0.) .or. & 
                        & ( (rVec%y >= grid(gP)%yAxis(grid(gP)%ny)+grid(gP)%geoCorrY &
                        & .or. yP>grid(gP)%ny) .and. vHat%y >=0.) .or. & 
                        & ( (rVec%z >= grid(gP)%zAxis(grid(gP)%nz)+grid(gP)%geoCorrZ &
                        &  .or. zP>grid(gP)%nz) .and. vHat%z >=0.) ) then
                      
                      ! go back to mother grid
                      xP = enPacket%xP(grid(gP)%motherP)
                      yP =  enPacket%yP(grid(gP)%motherP)
                      zP =  enPacket%zP(grid(gP)%motherP)
                      gP = grid(gP)%motherP

                   end if
                      
                end if

             end if
          end do ! safelimit loop


          if (i>= safeLimit) then
             if (.not.lgPlaneIonization) then
                print*, '! pathSegment: [warning] packet trajectory has exceeded&
                     &  maximum number of events', safeLimit, gP, xP,yP,zP, grid(gP)%active(xP,yP,zP), & 
                     & rvec, vhat, enPacket, iphot
             end if
             return
          end if

          enPacket%xP(gP) = xP
          enPacket%yP(gP) = yP
          enPacket%zP(gP) = zP
          
          ! the energy packet has beenid absorbed - reemit a new packet from this position
          call energyPacketRun(packetType, rVec, enPacket%xP(1:nGrids), enPacket%yP(1:nGrids), &
               & enPacket%zP(1:nGrids), gP)
          return
          
        end subroutine pathSegment


      end subroutine energyPacketDriver

         
    end module photon_mod


