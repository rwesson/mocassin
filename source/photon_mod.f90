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


    type(vector), parameter :: origin=vector(0.,0.,0.)  ! origin of the cartesian grid axes

    integer     , parameter :: safeLim = 10000          ! safety limit for the loops
    integer                        :: totalEscaped

    contains

    subroutine energyPacketDriver(iStar, n, grid, gpLoc, cellLoc)
        implicit none

        integer, intent(in)            :: n           ! number of energy packets
        integer, intent(in)            :: iStar       ! central star index

        integer, intent(in), optional &
             & :: gpLoc                               ! local grid (only used for extra diffuse sources)
        integer, intent(inout), optional &
             & :: cellLoc(3)                          ! local cell (only used for extra diffuse sources)

        type(grid_type), dimension(:), intent(inout) :: grid        ! the grid(s)

        type(vector)                   :: posDiff     ! initial position vector for diff ext
        type(vector)                   :: posVector   ! initial position vector for dust emi

        ! local variables
        integer                        :: gPIn        !
        integer                        :: igp         ! 1= mother 2 =sub
        integer                        :: i           ! counter
        integer                        :: iCell       ! cell counter
        integer                        :: igrid,ix,iy,iz ! location indeces
        integer                        :: ierr        ! allocation error status
        integer                        :: iPhot       ! counter
        integer                        :: seedSize    ! pseudo random number generator seed
        integer, dimension(2)          :: inX,inY,inZ ! initial position indeces
        integer, allocatable               :: seed(:)     ! seed array
        integer                        :: msec        ! millisecs of the sec
        integer                        :: dt(8)       ! date and time values
        integer                        :: trapped
        integer                        :: reRun

        character(len=7)               :: chTypeD     ! character type for driver
        character(len=7)               :: chTypeIn    ! character type

        type(vector)                   :: positionIn  ! position of packet absorption


        if (iStar == 0) then
           deltaE(0) = grid(gpLoc)%LdiffuseLoc(grid(gpLoc)%active(cellLoc(1),cellLoc(2),cellLoc(3)))/NphotonsDiffuseLoc
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

        if (allocated(seed)) deallocate(seed)

        Qphot = 0.

        trapped = 0

        do iPhot = 1, n

           if (iStar>=1) then

              chTypeD = "stellar"
              if (starIndeces(iStar,4)==1) then
                 igp = 1
              else if (starIndeces(iStar,4)>1) then
                 igp = 2
              else
                 print*, "! energyPacketDriver: insane grid pointer, starIndeces(iStar,4)", starIndeces(iStar,4)
                 stop
              end if

              inX=-1
              inY=-1
              inZ=-1

              inX(igp) = starIndeces(iStar,1)
              inY(igp) = starIndeces(iStar,2)
              inZ(igp) = starIndeces(iStar,3)


              chTypeIn = chTypeD
              positionIn = starPosition(iStar)
              gPIn = starIndeces(iStar,4)
              reRun = 0
              do i = 1, recursionLimit
                 call energyPacketRun(chTypeIn, positionIn, inX, inY, &
                      &inZ, gPIn, reRun)
                 if (rerun == 0) exit

              end do
              if (i>=recursionLimit) trapped = trapped+1

           else if (iStar==0) then

              chTypeD = "diffExt"
              inX=-1
              inY=-1
              inZ=-1

              if (gpLoc==1) then
                 igp = 1
              else if (gpLoc>1) then
                 igp = 2
              else
                 print*, "! energyPacketDriver: insane grid pointer"
                 stop
              end if

              inX(igp) = cellLoc(1)
              inY(igp) = cellLoc(2)
              inZ(igp) = cellLoc(3)

              posDiff%x = grid(gpLoc)%xAxis(cellLoc(1))
              posDiff%y = grid(gpLoc)%yAxis(cellLoc(2))
              posDiff%z = grid(gpLoc)%zAxis(cellLoc(2))

              chTypeIn = chTypeD
              positionIn = posDiff
              gPIn = gPLoc
              reRun = 0
              do i = 1, recursionLimit
                 call energyPacketRun(chTypeIn, positionIn, inX, inY, &
                      &inZ, gPIn, reRun)
                 if (rerun == 0) exit

              end do
              if (i>=recursionLimit) trapped = trapped+1

           else

              print*, '! energyPacketDriver: insanity in iStar value'
              stop

           end if
        end do

        if (iStar>0) then
          if (taskid==0) then
            print*, 'Star: ', iStar
            print*, 'Qphot = ', Qphot
          endif
        end if


        if (lgDust.and.convPercent>=resLinesTransfer .and.&
             & .not.lgResLinesFirst&
             & .and. (.not.nIterateMC==1)) then

           print*, "! energyPacketDriver: starting resonance line packets transfer"


           iCell = 0
           do igrid = 1, nGrids

              if (igrid==1) then
                 igp = 1
              else if (igrid>1) then
                 igp = 2
              else
                 print*, "! energyPacketDriver: insane grid pointer"
                 stop
              end if


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
                                inX(igp)=ix
                                inY(igp)=iy
                                inZ(igp)=iz

                                if (igrid>1) then
                                   ! check location on mother grid
                                   call locate(grid(grid(igrid)%motherP)%xAxis, &
                                        & posVector%x,inX(1))
                                   if (posVector%x > (grid(grid(igrid)%motherP)%xAxis(inX(1))+&
                                        & grid(grid(igrid)%motherP)%xAxis(inX(1)+1) )/2.) &
                                        & inX(1) = inX(1)+1
                                   call locate(grid(grid(igrid)%motherP)%yAxis, &
                                        & posVector%y,inY(1))
                                   if (posVector%y > (grid(grid(igrid)%motherP)%yAxis(inY(1))+&
                                        & grid(grid(igrid)%motherP)%yAxis(inY(1)+1) )/2.) &
                                        & inY(1) = inY(1)+1
                                   call locate(grid(grid(igrid)%motherP)%zAxis, &
                                        & posVector%z,inZ(1))
                                   if (posVector%z > (grid(grid(igrid)%motherP)%zAxis(inZ(1))+&
                                        & grid(grid(igrid)%motherP)%zAxis(inZ(1)+1) )/2.) &
                                        & inZ(1) = inZ(1)+1
                                end if

                                chTypeIn = chTypeD
                                positionIn = posVector
                                gPIn = igrid
                                reRun = 0

                                do i = 1, recursionLimit
                                   call energyPacketRun(chTypeIn, positionIn, inX, inY, &
                                        &inZ, gPIn, reRun)
                                   if (rerun == 0) exit

                                end do
                                if (i>=recursionLimit) trapped = trapped+1

                             end do
                          end if
                       end if

                    end do
                 end do
              end do

           end do

           print*, "! energyPacketDriver: ending resonance line packets transfer"

        end if


        if (iStar>0 .and. taskid==0) print*, 'Qphot = ', Qphot
        if (lgTalk) print*, 'Packets trapped = ', trapped, taskid

        ! evaluate Jste and Jdif
        ! NOTE : deltaE is in units of [E36 erg/s] however we also need a factor of
        ! 1.e-45 from the calculations of the volume of the cell hence these
        ! two factors cancel each other out giving units of [E-9erg/s] so we need to
        ! multiply by 1.E-9
        ! NOTE : Jste and Jdif calculated this way are in units of
        ! [erg sec^-1 cm^-2] -> no Hz^-1 as they are summed over separate bins (see
        ! Lucy A&A (1999)

        if(iStar>0. .and. taskid==0) then
           print*, 'Lstar', Lstar(iStar)
        end if


        contains


          subroutine energyPacketRun(chType, position, xP, yP, zP, gP, rR)
            implicit none

            type(vector),intent(inout) :: position         ! the position of the photon

            integer, dimension(2), intent(inout)    :: xP, yP, &
                 & zP                                            ! cartesian axes indeces
                                                                 ! 1= mother; 2=sub


            integer, intent(inout) :: rR               ! rerun?
            integer, intent(inout) :: gP               ! grid index
            integer                          :: igpr             ! grid pointer 1= mother 2=sub
            integer                          :: difSourceL(3)    ! cell indeces

            integer                          :: idirP, idirT     ! direction cosines

            type(photon_packet)              :: enPacket         ! the energu packet

            character(len=7), intent(inout)     :: chType           ! stellar or diffuse?

!            countRecursive = countRecursive+1
!            if (countRecursive > recursionLimit) then
!               trapped = trapped+1
!               return
!            end if

            rR = 0

            if (gP==1) then
               igpr = 1
            else if (gP>1) then
               igpr = 2
            else
               print*,  "! energyPacketRun: insane grid index"
               stop
            end if


            ! create a new photon packet
            select case (chType)

            ! if the energy packet is stellar
            case ("stellar")
               ! check for errors in the sources position
               if( position /= starPosition(iStar) ) then
                  print*, "! energyPacketRun: stellar energy packet must&
                       & start at the stellar position"
                  stop
               end if

               ! create the packet
               enPacket = newPhotonPacket(chType,position, xp, yp, zp, gp, noCellLoc)

            ! if the photon is from an extra source of diffuse radiation
             case ("diffExt")

                difSourceL(1) = xP(igpr)
                difSourceL(2) = yP(igpr)
                difSourceL(3) = zP(igpr)


                enPacket = newPhotonPacket(chType, position, xp, yp, zp, gP, difSourceL)

            ! if the photon is diffuse

            case ("diffuse")

                ! create the packet

                enPacket = newPhotonPacket(chType=chType, position=position, xP=xP, yP=yP, &
                     & zP=zP, gP=gP, difSource=noCellLoc)

            case ("dustEmi")

                ! crate the packet
                enPacket = newPhotonPacket(chType=chType, position=position, xP=xP, yP=yP, &
                     &zP=zP, gP=gP, difSource = noCellLoc)

            end select

            if (.not.lgDust .and. enPacket%nu < ionEdge(1) .and. .not.enPacket%lgLine) then

               ! the packet escapes without further interaction
               if (lgSymmetricXYZ) then
                  idirT = int(acos(abs(enPacket%direction%z))/dTheta)+1
               else
                  idirT = int(acos(enPacket%direction%z)/dTheta)+1
               end if

               if (idirT>totangleBinsTheta) then
                  idirT=totangleBinsTheta
               end if
               if (idirT<1 .or. idirT>totAngleBinsTheta) then
                  print*, '! energyPacketRun: error in theta direction cosine assignment',&
                       &  idirT, enPacket, dTheta, totAngleBinsTheta
                  stop
               end if


               if (abs(enPacket%direction%x)<1.e-35) then
                  idirP = 0
               else
                  if (lgSymmetricXYZ) then
                     idirP = int(atan(abs(enPacket%direction%y)/abs(enPacket%direction%x))/dPhi)
                  else
                     idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                  end if
               end if
               if (idirP<0) idirP=totAngleBinsPhi+idirP

               idirP=idirP+1

               if (idirP>totangleBinsPhi) then
                  idirP=totangleBinsPhi
!                  print*, '! energyPacketRun: idir>totanglebins - error corrected', &
!                       & idir, totanglebins, enPacket%direction, dtheta
               end if

               if (idirP<1 .or. idirP>totAngleBinsPhi) then
                  print*, '! energyPacketRun: error in phi direction cosine assignment -1',&
                       &  idirP, enPacket, dPhi, totAngleBinsPhi
                  stop
               end if

               if (nAngleBins>0) then
                  if (viewPointPtheta(idirT) > 0 .and. viewPointPhi(viewPointPtheta(idirT)) < 0) then

                     grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,viewPointPtheta(idirT)) = &
                          &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                          & enPacket%nuP,viewPointPtheta(idirT)) + deltaE(iStar)
                     grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                          & enPacket%nuP,0) = &
                          & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                          & enPacket%nuP,0) +  deltaE(iStar)

                  elseif (viewPointPtheta(idirT) == viewPointPphi(idirP).or. &
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

                  if (enPacket%origin(1) == 0) then
                     print*, '! energyPacketRun: enPacket%origin(1) ==0'
                     stop
                  end if
                  if (enPacket%origin(2) < 0) then
                     print*, '! energyPacketRun: enPacket%origin(2) < 0'
                     stop
                  end if

                  grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                       & enPacket%nuP,0) = &
                       & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                       & enPacket%nuP,0) +  deltaE(iStar)

               end if

               return
            end if

            ! if the new packet is capable of re-ionizing or we have dust
            if (.not.enPacket%lgLine) then

                ! compute the next segment of trajectory
                call pathSegment(enPacket)

                return

            else ! if the packet is a line packet
                ! add to respective line packet bin

               if (lgDebug) &
                    & grid(gP)%linePackets(grid(gP)%active(enPacket%xP(igpr), &
                    & enPacket%yP(igpr), enPacket%zP(igpr)), enPacket%nuP) = &
                    & grid(gP)%linePackets(grid(gP)%active(enPacket%xP(igpr), &
                    & enPacket%yP(igpr), enPacket%zP(igpr)), enPacket%nuP) +deltaE(iStar)


            end if

        end subroutine energyPacketRun


        ! this function initializes a photon packet
        function initPhotonPacket(nuP,  position, direction, lgLine, lgStellar, xP, yP, zP, gP, lgHg)
            implicit none

            type(photon_packet)      :: initPhotonPacket  ! the photon packet

            real                     :: random            ! random number

            integer, intent(in)      :: nuP               ! the frequency of the photon packet
            integer, intent(in),dimension(2) :: xP, yP, &
                 & zP                                     ! indeces of position on the x, y and z axes
            integer, intent(in)      :: gP                ! grid index
            integer                  :: igpi              ! grid pointer 1=mother, 2=sub


            logical, intent(in)      :: lgLine, lgStellar,&
                 & lgHG ! line, stellar packet Heyney-Greenstein?

            type(vector), intent(in) :: position          ! the position at which the photon
                                                          ! packet is created
            type(vector), intent(in) :: direction
            ! local variables

            integer                  :: irepeat           ! counter


            initPhotonPacket%direction = direction

            if (.not. (direction%x >= 0. .or. direction%x < 0)) then
               print*, '! initPhotonPacket: [0] insane direction%x [direction%x]'
               print*, direction%x
               stop
            end if


            initPhotonPacket%position = position

            initPhotonPacket%iG  = gP

            if (gP==1) then
               igpi=1
            else if (gp>1) then
               igpi=2
            else
               print*, "! initPhotonPacket: insane gridp pointer"
               stop
            end if

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

            initPhotonPacket%xP  = xP
            initPhotonPacket%yP  = yP
            initPhotonPacket%zP  = zP

            if (.not.lgHG .or. lgIsotropic .or. initPhotonPacket%lgStellar) then
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
                  if (initPhotonPacket%position%x>grid(gP)%xAxis(grid(gP)%nx))&
                       & initPhotonPacket%position%x=grid(gP)%xAxis(grid(gP)%nx)

                  call locate(grid(gP)%xAxis, initPhotonPacket%position%x,&
                       & initPhotonPacket%xP(igpi))
                  if (initPhotonPacket%xP(igpi) < grid(gP)%nx) then
                     if (initPhotonPacket%position%x >= &
                          &(grid(gP)%xAxis(initPhotonPacket%xP(igpi))+&
                          & grid(gP)%xAxis(initPhotonPacket%xP(igpi)+1))/2.) &
                          & initPhotonPacket%xP(igpi) = initPhotonPacket%xP(igpi)+1
                  end if

                  ! y-direction
                  initPhotonPacket%position%y = 0.
                  initPhotonPacket%yP(igpi) = 1

                  ! z-direction
                  call random_number(random)
                  random = 1. - random
                  initPhotonPacket%position%z = -&
                       & (grid(gP)%zAxis(2)-grid(gP)%zAxis(1))/2. &
                       & + random*( &
                       & (grid(gP)%zAxis(2)-grid(gP)%zAxis(1))/2.+&
                       & (grid(gP)%zAxis(grid(gP)%nz)-&
                       & grid(gP)%zAxis(grid(gP)%nz-1))/2.+&
                       & grid(gP)%zAxis(grid(gP)%nz))
                  if (initPhotonPacket%position%z<grid(gP)%zAxis(1)) &
                       & initPhotonPacket%position%z=grid(gP)%zAxis(1)
                  if (initPhotonPacket%position%z>&
                       & grid(gP)%zAxis(grid(gP)%nz)) &
                       & initPhotonPacket%position%z=&
                       & grid(gP)%zAxis(grid(gP)%nz)

                  call locate(grid(gP)%zAxis, &
                       & initPhotonPacket%position%z, initPhotonPacket%zP(igpi))
                  if (initPhotonPacket%zP(igpi) < grid(gP)%nz) then
                     if (initPhotonPacket%position%z >= &
                          & (grid(gP)%xAxis(initPhotonPacket%zP(igpi))+&
                          & grid(gP)%zAxis(initPhotonPacket%zP(igpi)+1))&
                          & /2.) initPhotonPacket%zP(igpi) =&
                          & initPhotonPacket%zP(igpi)+1
                  end if

                  if (initPhotonPacket%xP(igpi)<1) &
                       & initPhotonPacket%xP(igpi)=1
                  if (initPhotonPacket%zP(igpi)<1) &
                       & initPhotonPacket%zP(igpi)=1

                  ! direction is parallel to y-axis direction
                  initPhotonPacket%direction%x = 0.
                  initPhotonPacket%direction%y = 1.
                  initPhotonPacket%direction%z = 0.

                  if (initPhotonPacket%xP(igpi) >  &
                       & grid(gP)%xAxis(grid(gP)%nx) .or. &
                       & initPhotonPacket%zP(igpi) >  &
                       & grid(gP)%zAxis(grid(gP)%nz)) then
                     print*, "! initPhotonPacket: insanity in &
                          & planeIonisation init"
                     print*, igpi, initPhotonPacket%xP(igpi),  &
                          & grid(gP)%xAxis(grid(gP)%nx), &
                          & initPhotonPacket%zP(igpi), &
                          & grid(gP)%zAxis(grid(gP)%nz),  random, &
                          & initPhotonPacket%position%z

                     stop
                  end if

                  planeIonDistribution(initPhotonPacket%xP(igpi),&
                       & initPhotonPacket%zP(igpi)) = &
                       & planeIonDistribution(initPhotonPacket%xP(igpi),&
                       & initPhotonPacket%zP(igpi)) + 1

               else

                  do irepeat = 1, 1000000
                     ! get a random direction
                     initPhotonPacket%direction = randomUnitVector()

                     if (.not. (initPhotonPacket%direction%x >= 0. .or. initPhotonPacket%direction%x < 0)) then
                        print*, '! initPhotonPacket: [1] insane initPhotonPacket%direction%x [initPhotonPacket%direction%x]'
                        print*, initPhotonPacket%direction%x
                        stop
                     end if

                     if (initPhotonPacket%direction%x/=0. .and. &
                          & initPhotonPacket%direction%y/=0. .and. &
                          & initPhotonPacket%direction%z/=0.) exit
                  end do

               end if

               if ((lgSymmetricXYZ) .and. initPhotonPacket%lgStellar &
                    & .and. .not.lgMultistars) then
                  if (initPhotonPacket%direction%x<0.) &
                       & initPhotonPacket%direction%x = -initPhotonPacket%direction%x
                  if (initPhotonPacket%direction%y<0.) &
                       & initPhotonPacket%direction%y = -initPhotonPacket%direction%y
                  if (initPhotonPacket%direction%z<0.) &
                       & initPhotonPacket%direction%z = -initPhotonPacket%direction%z
               end if

            end if

            initPhotonPacket%origin(1) = gP
            initPhotonPacket%origin(2) = grid(gP)%active(initPhotonPacket%xP(igpi),&
                 & initPhotonPacket%yP(igpi), initPhotonPacket%zP(igpi))


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

        ! this subroutine determines the frequency of a newly created photon packet
        ! according to the given probability density
        ! does not use bisection to locate nu on array
        subroutine getNu2(probDen, nuP)

            real, dimension(:), intent(in) :: probDen    ! probability density function

            real                           :: random     ! random number

            integer, intent(out)           :: nuP        ! frequency index of the new

            integer                        :: isearch,i  !

            ! get a random number
            call random_number(random)

            do i = 1, 10000
               if (random==0 .or. random==1. .or. random == 0.9999999) then
                  call random_number(random)
               else
                  exit
               end if
            end do
            if (i>=10000) then
               print*, '! getNu2: problem with random number generator', random, i
               stop
            end if

            ! see what frequency random corresponds to
            nuP=1
            do isearch = 1, nbins
               if (random>=probDen(isearch)) then
                  nuP=isearch
               else
                  exit
               end if
            end do

            if (nuP<nbins-1) then
               nuP=nuP+1
            end if

            if (nuP>=nbins) then
               print*, 'random: ', random
               print*, 'probDen: ', probDen
            end if

          end subroutine getNu2


        ! this function creates a new photon packet
        function newPhotonPacket(chType, position, xP, yP, zP, gP, difSource)

            type(photon_packet)                :: newPhotonPacket! the photon packet to be created

            character(len=7), intent(in)       :: chType         ! stellar or diffuse?

            type(vector), intent(in), optional :: position       ! the position of the photon
                                                                 ! packet
            ! local variables
            type(vector)                       :: positionLoc    ! the position of the photon
                                                                 ! packet

            integer                            :: nuP            ! the frequency index of the photon packet
            integer, dimension(2)         :: orX,orY,orZ    ! dummy
            integer, dimension(2), intent(in) :: xP, yP, &
                 & zP                                            ! cartesian axes indeces
            integer, intent(in)      :: difSource(3)  ! grid and cell indeces
            integer, intent(inout)   :: gP
            integer                            :: igpn           ! grid pointe 1=motehr, 2=sub

            real                               :: random         ! random number

            if (gP==1) then
               igpn = 1
            else if (gp>1) then
               igpn = 2
            else
               print*,  "! newPhotonPacket: insane grid pointer"
               stop
            end if

            select case (chType)

            ! if the photon is stellar
            case ("stellar")

                ! check for errors in the sources position
               if( position /= starPosition(iStar) ) then
                  print*, "! newPhotonPacket: stellar photon packet must&
                       & start at the stellar position"
                  stop
               end if


!                gP = starIndeces(iStar,4)

                if (starIndeces(iStar,4) == 1) then
                   igpn = 1
                else if (starIndeces(iStar,4) > 1) then
                   igpn = 2
                else
                   print*,  "! newPhotonPacket: insane grid pointer -star position- "
                   stop
                end if

                ! determine the frequency of the newly created photon packet
                call getNu2(inSpectrumProbDen(iStar,1:nbins), nuP)

                if (nuP>=nbins) then
                   print*, "! newPhotonPacket: insanity occurred in stellar photon &
                        &nuP assignment (nuP,xP,yP,zP,activeP)", nuP, xP(igpn),yP(igpn),zP(igpn), &
                        & grid(starIndeces(iStar,4))%active(xP(igpn),yP(igpn),zP(igpn))
                   print*, "inSpectrumProbDen: ",iStar,inSpectrumProbDen(iStar,:), nuP
                   stop
                end if

                if (nuP < 1) then
                    print*, "! newPhotonPacket: insanity occurred in stellar photon &
&                               nuP assignment"
                    stop
                end if

                ! initialize the new photon packet
                orX(igpn) = starIndeces(iStar,1)
                orY(igpn) = starIndeces(iStar,2)
                orZ(igpn) = starIndeces(iStar,3)


                if (grid(starIndeces(iStar,4))%active(orX(igpn), orY(igpn), orZ(igpn)) < 0.) then
                   print*, "! newPhotonPacket: new packet cannot be emitted from re-mapped cell -1-"
                   print*, "chType, nuP, starPosition(iStar), .false., .true., orX,orY,orZ, gp"
                   print*, chType, nuP, starPosition(iStar), .false., .true., orX,orY,orZ, gp
                   stop
                end if

                newPhotonPacket = initPhotonPacket(nuP, &
                     &starPosition(iStar), nullUnitVector,.false., .true., &
                     & orX,orY,orZ, &
                     & starIndeces(iStar,4), .false.)


                if (newPhotonPacket%nu>1.) then
                   Qphot = Qphot + deltaE(iStar)/(2.1799153e-11*newPhotonPacket%nu)
                end if

                ! if the photon is from an extra diffuse source
             case ("diffExt")

                call getNu2(inSpectrumProbDen(0,1:nbins), nuP)

                if (nuP>=nbins) then
                   print*, "! newPhotonPacket: insanity occurred in extra diffuse photon &
                        & nuP assignment (nuP,gp,activeP)", nuP, gp
                   print*, "difSpectrumProbDen: ", inSpectrumProbDen(0,:)
                   stop
                end if

                if (nuP < 1) then
                   print*, "! newPhotonPacket: insanity occurred in extra diffuse photon &
                        & nuP assignment (nuP,gp,activeP)", nuP, gp,grid(gP)%active(xP(igpn),yP(igpn),zP(igpn))
                   print*, "difSpectrumProbDen: ", inSpectrumProbDen(0,:)
                   stop
                end if

                positionLoc%x = grid(gP)%xAxis(difSource(1))
                positionLoc%y = grid(gP)%yAxis(difSource(2))
                positionLoc%z = grid(gP)%zAxis(difSource(3))

                ! initialize the new photon packet
                orX(igPn) = difSource(1)
                orY(igPn) = difSource(2)
                orZ(igPn) = difSource(3)

                ! initialize the new photon packet
                if (grid(gP)%active(orX(igpn), orY(igpn), orZ(igpn)) < 0.) then
                   print*, "! newPhotonPacket: new packet cannot be emitted from re-mapped cell -3-"
                   print*, "chType, nuP, starPosition(iStar), .false., .true., orX,orY,orZ, gp"
                   print*, chType, nuP, starPosition(iStar), .false., .true., orX,orY,orZ, gp
                   stop
                end if
!                newPhotonPacket = initPhotonPacket(nuP, positionLoc, .false., .false., orX,&
!                     & orY, orZ, gP)

                newPhotonPacket = initPhotonPacket(nuP, positionLoc, &
                     & nullUnitVector,.false., .false., orX,&
                     & orY, orZ, gP, .false.)


            ! if the photon is diffuse
            case ("diffuse")

                ! check that gas is present in the grid
                if (.not.lgGas) then
                   print*, "! newPhotonPacket: diffuse packet cannot be created in a no gas grid"
                   stop
                end if

                ! check that the position is not inside the inner region
                if (grid(gP)%active(xP(igPn),yP(igPn),zP(igPn))<= 0) then
                    print*, "! newPhotonPacket: position of the new diffuse &
                         & photon packet is inside the inner region", xP(igPn),yP(igPn),zP(igPn),gP
                    stop
                end if

                ! decide whether continuum or line photon
                call random_number(random)

                random = 1.-random

                if (random <= grid(gP)%totalLines(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)))) then
                   ! line photon
                   ! line photons escape so don't care which one it is unless debugging
                   if (lgDebug) then
                      call getNu2( grid(gP)%linePDF(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)),:), nuP )

                      if (nuP < 1) then
                         print*, "! newPhotonPacket: insanity occurred in line photon &
                              & nuP assignment"
                         stop
                      end if

                   else

                      nuP = 0

                   end if

                   ! initialize the new photon packet
                   if (grid(gP)%active(xP(igpn), yp(igpn), zp(igpn)) < 0.) then
                      print*, "! newPhotonPacket: new packet cannot be emitted from re-mapped cell -4-"
                      print*, "chType, nuP, starPosition(iStar), .false., .true., xp,yp,zp, gp"
                      print*, chType, nuP, starPosition(iStar), .false., .true., xp,yp,zp, gp
                      stop
                   end if


                   newPhotonPacket = initPhotonPacket(nuP, position, &
                        & nullUnitVector, .true., .false., xP, yP, zP, gP, .false.)
!                   newPhotonPacket = initPhotonPacket(nuP, position, .true., .false., xP, yP, zP, gP)
                else
                    ! continuum photon

                    call getNu2(grid(gP)%recPDF(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)),1:nbins), nuP)

                    if (nuP>=nbins) then
                       print*, "! newPhotonPacket: insanity occurred in diffuse photon &
                       & nuP assignment (nuP,xP,yP,zP,activeP)", nuP, xP(igPn),yP(igPn),zP(igPn),&
                       &  grid(gP)%active(xP(igPn),yP(igPn),zP(igPn))
                       print*, "recPDF: ", grid(gP)%recPDF(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)),1:nbins)
                       stop
                    end if

                    if (nuP < 1) then
                        print*, "! newPhotonPacket: insanity occurred in diffuse photon &
                             & nuP assignment", nuP, xP(igPn),yP(igPn),zP(igPn), &
                             & grid(gP)%active(xP(igPn),yP(igPn),zP(igPn))
                       print*, "recPDF: ", grid(gP)%recPDF(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)),:)
                        stop
                    end if

                    ! initialize the new photon packet
                   if (grid(gP)%active(xP(igpn), yp(igpn), zp(igpn)) < 0.) then
                      print*, "! newPhotonPacket: new packet cannot be emitted from re-mapped cell -5-"
                      print*, "chType, nuP, starPosition(iStar), .false., .true., xp,yp,zp, gp"
                      print*, chType, nuP, starPosition(iStar), .false., .true., xp,yp,zp, gp
                      stop
                   end if

!                    newPhotonPacket = initPhotonPacket(nuP, position, .false., .false., xP, yP, zP, gP)


                    newPhotonPacket = initPhotonPacket(nuP, position, nullUnitVector,&
                         & .false., .false., xP, yP, zP, gP, .false.)

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

               ! check that the position is not inside the inner region
               if (grid(gP)%active(xP(igPn),yP(igPn),zP(igPn))<= 0) then
                  print*, "! newPhotonPacket: position of the new dust emitted &
                       &photon packet is inside the inner region", xP(igPn),yP(igPn),zP(igPn)
                  stop
               end if

               call getNu2(grid(gP)%dustPDF(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)),1:nbins), nuP)

               if (nuP>=nbins) then
                   print*, "! newPhotonPacket: insanity occurred in dust emitted photon &
                       &nuP assignment (iphot, nuP,xP(gP),yP(gP),zP(gP),activeP)", iphot, &
                       & nuP, xP(igPn),yP(igPn),zP(igPn), &
                       & grid(gP)%active(xP(igPn),yP(igPn),zP(igPn))
                   print*, "dustPDF: ", grid(gP)%dustPDF(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)),1:nbins)
                   print*, "grain T: ", grid(gP)%Tdust(:,0,grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)))
                   stop
                end if

               if (nuP < 1) then
                  print*, "! newPhotonPacket: insanity occurred in dust emitted photon &
                       &nuP assignment", nuP,xP(igPn),yP(igPn),zP(igPn), grid(gP)%active(xP(igPn),yP(igPn),zP(igPn))
                  print*, "dustPDF: ", grid(gP)%dustPDF(grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)),1:nbins)
                  print*, "grain T: ", grid(gP)%Tdust(:, 0, grid(gP)%active(xP(igPn),yP(igPn),zP(igPn)))
                  stop
               end if

               ! initialize the new photon packet
                   if (grid(gP)%active(xP(igpn), yp(igpn), zp(igpn)) < 0.) then
                      print*, "! newPhotonPacket: new packet cannot be emitted from re-mapped cell -6-"
                      print*, "chType, nuP, starPosition(iStar), .false., .true., xp,yp,zp, gp"
                      print*, chType, nuP, starPosition(iStar), .false., .true., xp,yp,zp, gp
                      stop
                   end if

!               newPhotonPacket = initPhotonPacket(nuP, position, .false., .false., xP, yP, zP, gP)

               newPhotonPacket = initPhotonPacket(nuP, position, nullUnitVector,&
                    & .false., .false., xP, yP, zP, gP, .false.)



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
          real                            :: dx, dy, dz
          real                            :: dSx, dSy, dSz
                                                      ! distances from x,y and z wall
          real                            :: dS       ! distance from nearest wall
          real                            :: dV       ! lume of this cell
          real                            :: passProb ! prob of passing the next segment
          real                            :: probSca  ! prob that the packet scatters
          real                            :: radius   ! radius
          real                            :: random   ! random number
          real                            :: tauCell  ! local tau

          integer                         :: iierr, ihg
          integer                         :: idirT,idirP ! direction cosine counters
          integer                         :: i, j, nS ! counter
          integer                         :: xP,yP,zP ! cartesian axes indeces
          integer                         :: gP       ! grid index
          integer                         :: igpp     ! grid index 1=mother 2=sub
          integer                         :: safeLimit =1000
                                                      ! safe limit for the loop

          character(len=7)                :: packetType ! stellar, diffuse, dustEmitted?

          logical                         :: lgScattered ! is the packet scattering with dust?
          logical                         :: lgReturn


          if (enPacket%iG == 1) then
             igpp = 1
          else if (enPacket%iG>1) then
             igpp = 2
          else
             print*, "! pathSegment: insane grid index"
             stop
          end if

          ! check that the input position is not outside the grid
          if ( (enPacket%iG <= 0).or.(enPacket%iG > nGrids) ) then
             print*, "! pathSegment: starting position not in any defined grids"
             print*, 'iPhot                 ', iPhot
             print*, 'enPacket%nuP          ', enPacket%nuP
             print*, 'enPacket%iG           ', enPacket%iG
             print*, 'enPacket%xP           ', enPacket%xP
             print*, 'enPacket%yP           ', enPacket%yP
             print*, 'enPacket%zP           ', enPacket%zP
             print*, 'enPacket%lgStellar    ', enPacket%lgStellar
             print*, 'enPacket%lgLine       ', enPacket%lgLine
             print*, 'enPacket%position     ', enPacket%position
             print*, 'enPacket%direction    ', enPacket%direction
             stop
          else if ( (enPacket%xP(igpp) <= 0).or.&
               &(enPacket%xP(igpp) > grid(enPacket%iG)%nx) ) then
             print*, "! pathSegment: starting position in x is outside the grid",&
                  & enPacket
             stop
          else if ( (enPacket%yP(igpp) <= 0).or. &
               & (enPacket%yP(igpp) > grid(enPacket%iG)%ny) ) then
             print*, "! pathSegment: starting position in y is outside the grid",&
                  & enPacket
             stop
          else if ( (enPacket%zP(igpp) <= 0).or.&
               & (enPacket%zP(igpp) > grid(enPacket%iG)%nz) ) then
             print*, "! pathSegment: starting position in z is outside the grid",&
                  & enPacket
             stop
          end if

          ! define vHat and rVec
          rVec = enPacket%position
          vHat = enPacket%direction


          if (.not. (rVec%x >= 0. .or. rVec%x < 0)) then
             print*, '! pathSegment: [0] insane rVec%x [rVec%x, dlLoc, vHat%x]'
             print*, rVec%x, dlLoc, vHat%x
             stop
          end if
          if (.not. (vHat%x >= 0. .or. vHat%x < 0)) then
             print*, '! pathSegment: [0] insane vHat%x [rVec%x, dlLoc, vHat%x]'
             print*, rVec%x, dlLoc, vHat%x
             stop
          end if


          ! initialize xP, yP,zP
          xP = enPacket%xP(igpp)
          yP = enPacket%yP(igpp)
          zP = enPacket%zP(igpp)
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
                print*, " ! pathSegment: multiple grids are not allowed in a 1D simulation"
                stop
             end if
          end if

          ! initialize optical depth
          absTau = 0.

          ! get a random number
          call random_number(random)

          ! calculate the probability
          passProb = -log(1.-random)

          ! speed up photons that my be trapped
          if (lgPlaneIonization) then
             safeLimit=5000
          else
             safeLimit=500000
!             safeLimit=500
          end if

          do i = 1, safeLimit

             do j = 1, safeLimit

                if (xP > grid(gP)%nx .or. xP < 1 .or. &
                     & yP > grid(gP)%ny .or. yP < 1 .or. &
                     & zP > grid(gP)%nz .or. zP < 1 ) then
                   print*, "! pathSegment: insanity [gp,xp,yp,zp,j,i]", &
                        & gp, xp, yp, zp, j, i
                   stop
                end if

                if (grid(gP)%active(xP,yP,zP)<0) then

                   ! packet is entering a subgrid
                   enPacket%xP(1) = xP
                   enPacket%yP(1) = yP
                   enPacket%zP(1) = zP


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


                   enPacket%iG = gP
                   igpp = 2

                end if ! line 1224


                ! find distances from all walls

                if (lgSymmetricXYZ) then
                   if ( rVec%x <= grid(1)%xAxis(1) ) then
                      if (vHat%x<0.) vHat%x = -vHat%x
                      rVec%x = grid(1)%xAxis(1)
                   end if
                   if ( rVec%y <= grid(1)%yAxis(1) ) then
                      if (vHat%y<0.) vHat%y = -vHat%y
                      rVec%y = grid(1)%yAxis(1)
                   end if
                   if ( rVec%z <= grid(1)%zAxis(1) ) then
                      if (vHat%z<0.) vHat%z = -vHat%z
                      rVec%z = grid(1)%zAxis(1)
                   end if
                end if

                if (vHat%x>1.e-10) then
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
                else if (vHat%x<-1.e-10) then
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
                else
!                   dSx = grid(gP)%xAxis(grid(gP)%nx)
                   dSx =1.e35
                end if

                if (.not. (dSx >= 0. .or. dSx <0.)) then
                   print*, '! pathSegment: insane dSx [dSx, xP, grid(gP)%xAxis(xP), rVec%x, vHat%x]'
                   print*, dSx, xP, grid(gP)%xAxis(xP), rVec%x, vHat%x
                   stop
                end if

                if (.not.lg1D) then
                   if (vHat%y>1.e-10) then
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
                   else if (vHat%y<-1.e-10) then
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
                   else
!                      dSy = grid(gP)%yAxis(grid(gP)%ny)
                      dSy = 1.e35
                   end if

                   if (.not. (dSy >= 0. .or. dSy <0.)) then
                      print*, '! pathSegment: insane dSy [dSy, yP, grid(gP)%yAxis(yP), rVec%y, vHat%y]'
                      print*, dSy, yP, grid(gP)%yAxis(yP), rVec%y, vHat%y
                      stop
                   end if


                   if (vHat%z>1.e-10) then
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
                   else if (vHat%z<-1.e-10) then
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
                   else
!                      dSz = grid(gP)%zAxis(grid(gP)%nz)
                      dSz = 1.e35
                   end if

                   if (.not. (dSz >= 0. .or. dSz <0.)) then
                      print*, '! pathSegment: insane dSz [dSz, zP, grid(gP)%zAxis(zP), rVec%z, vHat%z]'
                      print*, dSz, zP, grid(gP)%zAxis(zP), rVec%z, vHat%z
                      stop
                   end if

                   if (xP > grid(gP)%nx .or. xP < 1 .or. &
                        & yP > grid(gP)%ny .or. yP < 1 .or. &
                        & zP > grid(gP)%nz .or. zP < 1 ) then
                      print*, "! pathSegment: insanity -2- [gp,xp,yp,zp]", &
                           & gp, xp, yp, zp
                      stop
                   end if

                end if

                if (grid(gP)%active(xP,yP,zP)>=0) exit
             end do

             ! cater for cells on cell wall
             if ( abs(dSx)<1.e-10 ) dSx = grid(gP)%xAxis(grid(gP)%nx)
             if ( abs(dSy)<1.e-10 ) dSy = grid(gP)%yAxis(grid(gP)%ny)
             if ( abs(dSz)<1.e-10 ) dSz = grid(gP)%zAxis(grid(gP)%nz)

             ! find the nearest wall
             dSx = abs(dSx)
             dSy = abs(dSy)
             dSz = abs(dSz)

             if (dSx<=0.) then
                print*, '! pathSegment: [warning] dSx <= 0.',dSx
                print*, 'grid(gP)%xAxis ', grid(gP)%xAxis
                print*, 'gP,xP,grid(gP)%xAxis(xP), rVec%x, vHat%x'
                print*, gP,xP,grid(gP)%xAxis(xP), rVec%x, vHat%x
                dS = min(dSy, dSz)
             else if (dSy<=0.) then
                print*, '! pathSegment: [warning] dSy <= 0.', dSy
                print*, 'grid(gP)%yAxis ', grid(gP)%yAxis
                print*, 'gP,yP,grid(gP)%yAxis(yP), rVec%y, vHat%y'
                print*, gP,yP,grid(gP)%yAxis(yP), rVec%y, vHat%y
                dS = min(dSx, dSz)
             else if (dSz<=0.) then
                print*, '! pathSegment: [warning] dSz <= 0.', dSz
                print*, 'grid(gP)%zAxis ', grid(gP)%zAxis
                print*, 'gP,zP,grid(gP)%zAxis(zP), rVec%z, vHat%z'
                print*, gP,zP,grid(gP)%zAxis(zP), rVec%z, vHat%z
                dS = min(dSx, dSy)
             else
                dS = min(dSx,dSy)
                dS = min(dS, dSz)
             end if

             ! this should now never ever happen
             if (dS <= 0.) then
                print*, 'pathSegment: dS <= 0', dSx, dSy, dSz
                print*, gP, rVec
                stop
             end if

             ! calculate the optical depth to the next cell wall
             tauCell = dS*grid(gP)%opacity(grid(gP)%active(xP,yP,zP), enPacket%nuP)


             ! find the volume of this cell
!             dV = getVolumeLoc(grid(gP), xP,yP,zP)


             if (lg1D) then
                if (nGrids>1) then
                   print*, '! getVolumeLoc: 1D option and multiple grids options are not compatible'
                   stop
                end if

                if (xP == 1) then

                   dV = 4.*Pi* ( (grid(gP)%xAxis(xP+1)/1.e15)**3)/3.


                else if ( xP==grid(gP)%nx) then

                   dV = Pi* ( (3.*(grid(gP)%xAxis(xP)/1.e15)-(grid(gP)%xAxis(xP-1)/1.e15))**3 - &
                        & ((grid(gP)%xAxis(xP)/1.e15)+(grid(gP)%xAxis(xP-1)/1.e15))**3 ) / 6.

                else

                   dV = Pi* ( ((grid(gP)%xAxis(xP+1)/1.e15)+(grid(gP)%xAxis(xP)/1.e15))**3 - &
                        & ((grid(gP)%xAxis(xP-1)/1.e15)+(grid(gP)%xAxis(xP)/1.e15))**3 ) / 6.

                end if

                dV = dV/8.

             else

                if ( (xP>1) .and. (xP<grid(gP)%nx) ) then

                   dx = abs(grid(gP)%xAxis(xP+1)-grid(gP)%xAxis(xP-1))/2.
                else if ( xP==1 ) then
                   if (lgSymmetricXYZ) then
                      dx = abs(grid(gP)%xAxis(xP+1)-grid(gP)%xAxis(xP))/2.
                   else
                      dx = abs(grid(gP)%xAxis(xP+1)-grid(gP)%xAxis(xP))
                   end if
                else if ( xP==grid(gP)%nx ) then
                   dx = abs(grid(gP)%xAxis(xP)  -grid(gP)%xAxis(xP-1))
                end if

                if ( (yP>1) .and. (yP<grid(gP)%ny) ) then
                   dy = abs(grid(gP)%yAxis(yP+1)-grid(gP)%yAxis(yP-1))/2.
                else if ( yP==1 ) then
                   if (lgSymmetricXYZ) then
                      dy = abs(grid(gP)%yAxis(yP+1)-grid(gP)%yAxis(yP))/2.
                   else
                      dy = abs(grid(gP)%yAxis(yP+1)-grid(gP)%yAxis(yP))
                   end if
                else if ( yP==grid(gP)%ny ) then
                   dy = abs(grid(gP)%yAxis(yP)  -grid(gP)%yAxis(yP-1))
                end if

                if ( (zP>1) .and. (zP<grid(gP)%nz) ) then
                   dz = abs(grid(gP)%zAxis(zP+1)-grid(gP)%zAxis(zP-1))/2.
                 else if ( zP==1 ) then
                   if (lgSymmetricXYZ) then
                      dz = abs(grid(gP)%zAxis(zP+1)-grid(gP)%zAxis(zP))/2.
                   else
                      dz = abs(grid(gP)%zAxis(zP+1)-grid(gP)%zAxis(zP))
                   end if
                else if ( zP==grid(gP)%nz ) then
                   dz = abs(grid(gP)%zAxis(zP)-grid(gP)%zAxis(zP-1))
                end if

                dx = dx/1.e15
                dy = dy/1.e15
                dz = dz/1.e15


                ! calculate the volume
                dV = dx*dy*dz

             end if

             ! check if the packet interacts within this cell
             if ((absTau+tauCell > passProb) .and. (grid(gP)%active(xP,yP,zP)>0)) then
! .and. &
!                  & grid(gP)%opacity(grid(gP)%active(xP,yP,zP), enPacket%nuP) > 1.e-10) then

                ! packet interacts

                ! calculate where within this cell the packet is absorbed
                dlLoc = (passProb-absTau)/grid(gP)%opacity(grid(gP)%active(xP,yP,zP), enPacket%nuP)

                ! update packet's position
                rVec = rVec + dlLoc*vHat

                if (.not. (rVec%x >= 0. .or. rVec%x < 0)) then
                   print*, '! pathSegment: insane rVec%x [rVec%x, dlLoc, vHat%x]'
                   print*, rVec%x, dlLoc, vHat%x
                   stop
                end if

                if (.not. (rVec%y >= 0. .or. rVec%y < 0)) then
                   print*, '! pathSegment: insane rVec%y [rVec%y, dlLoc, vHat%y]'
                   print*, rVec%y, dlLoc, vHat%y
                   stop
                end if

                if (.not. (rVec%z >= 0. .or. rVec%z < 0)) then
                   print*, '! pathSegment: insane rVec%z [rVec%z, dlLoc, vHat%z]'
                   print*, rVec%z, dlLoc, vHat%z
                   stop
                end if

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
                if ( sqrt( (rvec%x/1.e10)**2 + (rvec%y/1.e10)**2 + (rvec%z/1.e10)**2)*1.e10 >= R_out &
                     & .and. R_out > 0.) then


                   ! the packet escapes without further interaction
                   if (lgSymmetricXYZ) then
                      idirT = int(acos(abs(enPacket%direction%z))/dTheta)+1
                   else
                      idirT = int(acos(enPacket%direction%z)/dTheta)+1
                   end if

                   if (idirT>totangleBinsTheta) then
                      idirT=totangleBinsTheta
                   end if
                   if (idirT<1 .or. idirT>totAngleBinsTheta) then
                      print*, '! pathSegment: error in theta direction cosine assignment',&
                           &  idirT, enPacket, dTheta, totAngleBinsTheta
                      stop
                   end if

!                   if (enPacket%direction%x<1.e-35) then
!                      idirP = 0
!                   else
!                      idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
!                   end if
!                   if (idirP<0) idirP=totAngleBinsPhi+idirP


                   if (abs(enPacket%direction%x)<1.e-35) then
                      idirP = 0
                   else
                      if (lgSymmetricXYZ) then
                         idirP = int(atan(abs(enPacket%direction%y)/abs(enPacket%direction%x))/dPhi)
                      else
                         idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                      end if
                   end if
                   if (idirP<0) idirP=totAngleBinsPhi+idirP

                   idirP=idirP+1

                   if (idirP>totangleBinsPhi) then
                      idirP=totangleBinsPhi
                   end if

                   if (idirP<1 .or. idirP>totAngleBinsPhi) then
                      print*, '! pathSegment: error in Phi direction cosine assignment -2',&
                           &  idirP, enPacket, dPhi, totAngleBinsPhi
                      stop
                   end if

!print*, 'a', nanglebins
!print*, enPacket%lgStellar, enPacket%direction
!print*, idirT, viewPointPTheta(idirT), viewPointTheta(viewPointPTheta(idirT))
!print*, idirP, viewPointPPhi(idirP), viewPointPhi(viewPointPPhi(idirP)), viewPointPhi(viewPointPtheta(idirT))

                   if (nAngleBins>0) then
                      if (viewPointPtheta(idirT) > 0 .and. viewPointPhi(viewPointPtheta(idirT)) < 0) then
                         grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,viewPointPtheta(idirT)) = &
                              &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              & enPacket%nuP,viewPointPtheta(idirT)) + deltaE(iStar)
                         grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              & enPacket%nuP,0) = &
                              & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              & enPacket%nuP,0) +  deltaE(iStar)

                      elseif (viewPointPtheta(idirT) == viewPointPphi(idirP).or. &
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
                      print*, '! pathSegment: insanity occurred at the scattering/absorption&
                           & decision stage.'
                      print*,'gP,xP,yP,zP,nuP'
                      print*,gP,xP,yP,zP,nuP
                      print*,'scaOpac,opacity,Ndust'
                      print*,&
 &grid(gP)%scaOpac(grid(gP)%active(xP,yP,zP),enPacket%nuP),&
 &grid(gP)%opacity(grid(gP)%active(xP,yP,zP),enPacket%nuP),&
 &grid(gP)%Ndust(grid(gP)%active(xP,yP,zP))
                      print*,'random number=',random,'probSca=',probSca
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

                      if (lgMultiDustChemistry) then
                         do nS = 1, nSpeciesPart(grid(gP)%dustAbunIndex(grid(gP)%active(xP,yP,zP)))
                            if (grainabun(grid(gP)%dustAbunIndex(grid(gP)%active(xP,yP,zP)),nS)>0. &
                                 &.and. grid(gP)%Tdust(nS, 0, &
                                 & grid(gP)%active(xP,yP,zP))<TdustSublime(dustComPoint(&
                                 &grid(gP)%dustAbunIndex(grid(gP)%active(xP,yP,zP)))-1+nS)) exit
                         end do
                         if (nS>nSpeciesPart(grid(gP)%dustAbunIndex(grid(gP)%active(xP,yP,zP)))) then
                            print*, "! pathSegment: packet scatters with dust at position where all &
                                 &grains have sublimed -1-."
                            print*, xP,yP,zP, grid(gP)%active(xP,yP,zP), tauCell, absTau, passProb
                            stop
                         end if
                      else
                         do nS = 1, nSpeciesPart(1)
                            if (grainabun(1,nS)>0. &
                                 &.and. grid(gP)%Tdust(nS, 0, &
                                 & grid(gP)%active(xP,yP,zP))<TdustSublime(dustComPoint(&
                                 &1)-1+nS)) exit
                         end do
                         if (nS>nSpeciesPart(1)) then
                            print*, "! pathSegment: packet scatters with dust at position where all &
                                 &grains have sublimed -2-."
                            print*, xP,yP,zP, grid(gP)%active(xP,yP,zP), tauCell, absTau, passProb
                            stop
                         end if
                      end if

                      ! packet is scattered by the grain
                      ! calculate new direction

                      enPacket%xP(igpp) = xP
                      enPacket%yP(igpp) = yP
                      enPacket%zP(igpp) = zP


                      if (grid(gP)%active(enPacket%xp(igpp), enPacket%yp(igpp), enPacket%zp(igpp)) < 0.) then
                         print*, "! pathSegment: new packet cannot be emitted from re-mapped cell -1-"
                         print*, "nuP, starPosition(iStar), .false., .true., xp,yp,zp, gp"
                         print*, nuP, starPosition(iStar), .false., .true.,  xp,yp,zp, gp
                         stop
                      end if

                      enPacket = initPhotonPacket(enPacket%nuP, rVec, enPacket%direction, .false., .false., enPacket%xP(1:2), &
                           & enPacket%yP(1:2), enPacket%zP(1:2), gP, .true.)

                      if (.not.lgIsotropic .and. .not.enPacket%lgStellar) then
                         do ihg = 1,10
                            call hg(enPacket,iierr)
                            if(iierr==0) exit
                         end do
                         if (ihg >=10) then
                            print*, '! pathSegment [warning]: hg called ten times [enPacket]', enPacket
                         end if
                      end if

                      vHat%x = enPacket%direction%x
                      vHat%y = enPacket%direction%y
                      vHat%z = enPacket%direction%z


                      if (.not. (enPacket%direction%x >= 0. .or. enPacket%direction%x < 0)) then
                         print*, '! pathSegment: [1] insane direction%x [direction%x]'
                         print*, enPacket%direction%x
                         stop
                      end if


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
                      print*, "! pathSegment: Insanity occurred - no gas present when no dust interaction"
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

                   if (enPacket%xP(1) <= 0 .or. &
                        & enPacket%yP(1) <= 0 .or. &
                        & enPacket%zP(1) <= 0 .or. &
                        & enPacket%xP(1) > grid(grid(gP)%motherP)%nx .or. &
                        & enPacket%yP(1) > grid(grid(gP)%motherP)%ny .or. &
                        & enPacket%zP(1) > grid(grid(gP)%motherP)%nz) then

                      ! locate where we are at on the mother grid
                      call locate(grid(grid(gP)%motherP)%xAxis, rvec%x, enPacket%xP(1))
                      if ( enPacket%xP(1)< grid(grid(gP)%motherP)%nx) then
                         if ( rvec%x > ( grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)) + &
                              & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1)/2.) ) &
                              & enPacket%xP(1) = enPacket%xP(1)+1
                      end if

                      call locate( grid(grid(gP)%motherP)%yAxis, rvec%y, enPacket%yP(1))
                      if ( enPacket%yP(1)< grid(grid(gP)%motherP)%ny) then
                         if ( rvec%y > ( grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)) + &
                              & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1)/2.) ) &
                              & enPacket%yP(1) = enPacket%yP(1)+1
                      end if

                      call locate(grid(grid(gP)%motherP)%zAxis, rvec%z, enPacket%zP(1))
                      if ( enPacket%zP(1)< grid(grid(gP)%motherP)%nz) then
                         if ( rvec%z > ( grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)) + &
                              & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1)/2.) ) &
                              & enPacket%zP(1) = enPacket%zP(1)+1
                      end if

                   else

                      if (vHat%x>0.) then
                         if ( enPacket%xP(1) < grid(grid(gP)%motherP)%nx ) then
                            if ( rVec%x > (grid(grid(gP)%motherP)%xAxis(enPacket%xP(1))+&
                                 & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)+1))/2. ) then
                               enPacket%xP(1) = enPacket%xP(1)+1
                            end if
                         else
                            if ( rVec%x > grid(grid(gP)%motherP)%xAxis(enPacket%xP(1))) then
!                            print*, '! pathSegment: insanity occurred at mother grid transfer (x axis +)', &
!                                 & rVec%x, gP, grid(gP)%motherP
!                            stop
                            end if
                         end if
                      else
                         if ( enPacket%xP(1) > 1 ) then
                            if ( rVec%x <= (grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)-1)+&
                                 & grid(grid(gP)%motherP)%xAxis(enPacket%xP(1)))/2. ) then
                               enPacket%xP(1) = enPacket%xP(1)-1
                            end if
                         else
                            if (rVec%x < grid(grid(gP)%motherP)%xAxis(1)) then
!                            print*, '! pathSegment: insanity occurred at mother grid transfer (x axis-)',&
!                                 & rVec%x, gP, grid(gP)%motherP
!                            stop
                            end if
                         end if
                      end if
                      if (vHat%y>0.) then
                         if (  enPacket%yP(1) < grid(grid(gP)%motherP)%ny ) then
                            if ( rVec%y > (grid(grid(gP)%motherP)%yAxis( enPacket%yP(1))+&
                                 & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)+1))/2. ) then
                               enPacket%yP(1) =  enPacket%yP(1)+1
                            end if
                         else
                            if ( rVec%y > grid(grid(gP)%motherP)%yAxis( enPacket%yP(1))) then
!                            print*, '! pathSegment: insanity occurred at mother grid transfer (y axis +)',&
!                                 & rVec%y, gP, grid(gP)%motherP
!                            stop
                            end if
                         end if
                      else
                         if (  enPacket%yP(1) > 1 ) then
                            if ( rVec%y <= (grid(grid(gP)%motherP)%yAxis( enPacket%yP(1)-1)+&
                                 & grid(grid(gP)%motherP)%yAxis(enPacket%yP(1)))/2. ) then
                               enPacket%yP(1) =  enPacket%yP(1)-1
                            end if
                         else
                            if (rVec%y < grid(grid(gP)%motherP)%yAxis(1)) then
!                            print*, '! pathSegment: insanity occurred at mother grid transfer (y axis -)', &
!                                 & rVec%y, gP, grid(gP)%motherP
!                            stop
                            end if
                         end if
                      end if
                      if (vHat%z>0.) then
                         if (  enPacket%zP(1) < grid(grid(gP)%motherP)%nz ) then
                            if ( rVec%z > (grid(grid(gP)%motherP)%zAxis( enPacket%zP(1))+&
                                 & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)+1))/2. ) then
                               enPacket%zP(1) =  enPacket%zP(1)+1
                            end if
                         else
                            if ( rVec%z > grid(grid(gP)%motherP)%zAxis( enPacket%zP(1))) then
!                            print*, '! pathSegment: insanity occurred at mother grid transfer (z axis +)', &
!                                 & rVec%z, gP, grid(gP)%motherP
!                            stop
                            end if
                         end if
                      else
                         if (  enPacket%zP(1) > 1 ) then
                            if ( rVec%z <= (grid(grid(gP)%motherP)%zAxis( enPacket%zP(1)-1)+&
                                 & grid(grid(gP)%motherP)%zAxis(enPacket%zP(1)))/2. ) then
                               enPacket%zP(1) =  enPacket%zP(1)-1
                            end if
                         else
                            if (rVec%z < grid(grid(gP)%motherP)%zAxis(1)) then
!                        print*, '! pathSegment: insanity occurred at mother grid transfer (z axis -)', &
!                             & rVec%z, gP, grid(gP)%motherP
!                        stop
                            end if
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
                      print*, '! pathSegment: insanity occurred in dS assignement &
                           & [dS,dSx,dSy,dSz,vHat]', dS,dSx,dSy,dSz,vHat
                   end if
                else
                   radius = 1.e10*sqrt((rVec%x/1.e10)*(rVec%x/1.e10) + &
                        & (rVec%y/1.e10)*(rVec%y/1.e10) + &
                        & (rVec%z/1.e10)*(rVec%z/1.e10))
                   call locate(grid(gP)%xAxis, radius , xP)

                end if

                ! be 6/6/06
                if(.not.lgPlaneIonization.and..not.lgSymmetricXYZ) then
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
                         print*, '! pathSegment: insanity occurred - invalid gP', gP
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
                         print*, '! pathSegment: insanity occurred - invalid gP', gP
                         stop
                      end if

                   end if

                   if ( (rVec%x <= grid(gP)%xAxis(1)-grid(gP)%geoCorrX .or. xP<1) .and. gP==1) then
                      xP=1
                      lgReturn=.true.

                   end if


                   if ( (rVec%x <= grid(gP)%xAxis(1)-grid(gP)%geoCorrX .or. xP<1) &
                    & .and. gP>1) then

                      xP = enPacket%xP(grid(gP)%motherP)
                      yP =  enPacket%yP(grid(gP)%motherP)
                      zP =  enPacket%zP(grid(gP)%motherP)
                      gP = grid(gP)%motherP

                   end if


                   if ( (rVec%x >=  grid(gP)%xAxis(grid(gP)%nx)+grid(gP)%geoCorrX &
                        & .or. xP>grid(gP)%nx) .and. gP==1 )then
                      xP = grid(gP)%nx
                      lgReturn=.true.

                   end if

                   if ( (rVec%x >=  grid(gP)%xAxis(grid(gP)%nx)+grid(gP)%geoCorrX&
                        & .or. xP>grid(gP)%nx) .and.  gP>1) then

                      xP = enPacket%xP(grid(gP)%motherP)
                      yP =  enPacket%yP(grid(gP)%motherP)
                      zP =  enPacket%zP(grid(gP)%motherP)
                      gP = grid(gP)%motherP

                   end if

                   if ( (rVec%z <= grid(gP)%zAxis(1)-grid(gP)%geoCorrZ .or.zP<1) &
                        & .and. gP==1) then
                      zP=1
                      lgReturn=.true.

                   end if

                   if ( (rVec%z <= grid(gP)%zAxis(1)-grid(gP)%geoCorrZ &
                        & .or.zP<1) .and. gP>1) then

                      xP = enPacket%xP(grid(gP)%motherP)
                      yP =  enPacket%yP(grid(gP)%motherP)
                      zP =  enPacket%zP(grid(gP)%motherP)
                      gP = grid(gP)%motherP

                   end if

                   if ( (rVec%z >=  grid(gP)%zAxis(grid(gP)%nz)+grid(gP)%geoCorrZ &
                        & .or. zP>grid(gP)%nz) &
                        & .and. gP==1) then

                      zP = grid(gP)%nz
                      lgReturn=.true.

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
                      ! the packet escapes without further interaction
                      if (lgSymmetricXYZ) then
                         idirT = int(acos(abs(enPacket%direction%z))/dTheta)+1
                      else
                         idirT = int(acos(enPacket%direction%z)/dTheta)+1
                      end if

!                      idirT = int(acos(enPacket%direction%z)/dTheta)+1
                      if (idirT>totangleBinsTheta) then
                         idirT=totangleBinsTheta
                      end if
                      if (idirT<1 .or. idirT>totAngleBinsTheta) then
                         print*, '! pathSegment: error in theta direction cosine assignment',&
                              &  idirT, enPacket, dTheta, totAngleBinsTheta
                         stop
                      end if


!                      if (enPacket%direction%x<1.e-35) then
!                         idirP = 0
!                      else
!                         idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
!                      end if
!                      if (idirP<0) idirP=totAngleBinsPhi+idirP

                      if (abs(enPacket%direction%x)<1.e-35) then
                         idirP = 0
                      else
                         if (lgSymmetricXYZ) then
                            idirP = int(atan(abs(enPacket%direction%y)/abs(enPacket%direction%x))/dPhi)
                         else
                            idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                         end if
                      end if
                      if (idirP<0) idirP=totAngleBinsPhi+idirP

                      idirP=idirP+1

                      if (idirP>totangleBinsPhi) then
                         idirP=totangleBinsPhi
                      end if
                      if (idirP<1 .or. idirP>totAngleBinsPhi) then
                         print*, '! pathSegment: error in phi direction cosine assignment -3',&
                              &  idirP, enPacket, dPhi, totAngleBinsPhi
                         stop
                      end if


!print*, 'b''c', nanglebins
!print*, enPacket%lgStellar, enPacket%direction
!print*, idirT, viewPointPTheta(idirT), viewPointTheta(viewPointPTheta(idirT))
!print*, idirP, viewPointPPhi(idirP), viewPointPhi(viewPointPPhi(idirP)), viewPointPhi(viewPointPtheta(idirT))

                      if (nAngleBins>0) then
                         if (viewPointPtheta(idirT) > 0 .and. viewPointPhi(viewPointPtheta(idirT)) < 0) then

                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,viewPointPtheta(idirT)) = &
                                 &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,viewPointPtheta(idirT)) + deltaE(iStar)
                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,0) = &
                                 & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,0) +  deltaE(iStar)

                         elseif (viewPointPtheta(idirT) == viewPointPphi(idirP).or. &
                              & (viewPointTheta(viewPointPphi(idirP))==viewPointTheta(viewPointPtheta(idirT))) .or. &
                              & (viewPointPhi(viewPointPtheta(idirT))==viewPointPhi(viewPointPphi(idirP))) ) then
                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,&
                                 & viewPointPtheta(idirT)) =grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2),&
                                 & enPacket%nuP,viewPointPtheta(idirT)) +  deltaE(iStar)

                            if (viewPointPtheta(idirT)/=0) grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 &enPacket%nuP,0) = &
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
                              enPacket%nuP,0) = &
                              & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              & enPacket%nuP,0) +  deltaE(iStar)
                      end if

                      return
                   end if

                end if

                ! end be 6/6/06


                if(lgPlaneIonization) then
                   lgReturn=.false.

                   if ( rVec%y <= grid(gP)%yAxis(1)-grid(gP)%geoCorrY .or. yP<1) then

                      ! the energy packet escapes this grid
                      if (gP==1) then
                         yP=1
                         lgReturn=.true.
                      else if (gP>1) then
                         xP = enPacket%xP(1)
                         yP = enPacket%yP(1)
                         zP = enPacket%zP(1)
                         gP = 1
                         igpp = 1
                      else
                         print*, '! pathSegment: insanity occurred - invalid gP', gP
                         stop
                      end if

                   end if

                   if (rVec%y > grid(gP)%yAxis(grid(gP)%ny)+grid(gP)%geoCorrY .or. yP>grid(gP)%ny) then

                      if (gP==1) then
                         ! the energy packet escapes
                         yP = grid(gP)%ny
                         lgReturn=.true.
                      else if (gP>1) then
                         xP = enPacket%xP(1)
                         yP = enPacket%yP(1)
                         zP = enPacket%zP(1)
                         gP = 1
                         igpp = 1
                      else
                         print*, '! pathSegment: insanity occurred - invalid gP', gP
                         stop
                      end if

                   end if

                   if ( (rVec%x <= grid(1)%xAxis(1) .or. xP<1) ) then
                      xP=1
                      rVec%x = grid(gP)%xAxis(1)
                      vHat%x = -vHat%x

                   end if

                   if ( (rVec%x <= grid(gP)%xAxis(1)-grid(gP)%geoCorrX .or. xP<1) &
                        & .and. gP>1) then

                      xP = enPacket%xP(1)
                      yP = enPacket%yP(1)
                      zP = enPacket%zP(1)
                      gP = 1
                      igpp = 1

                   end if

                   if ( (rVec%x >=  grid(1)%xAxis(grid(gP)%nx) &
                        & .or. xP>grid(gP)%nx)  )then
                      xP = grid(gP)%nx
                      rVec%x = grid(gP)%xAxis(grid(gP)%nx)
                      vHat%x = -vHat%x

                   end if

                   if ( (rVec%x >=  grid(gP)%xAxis(grid(gP)%nx)+grid(gP)%geoCorrX&
                        & .or. xP>grid(gP)%nx) .and.  gP>1) then

                      xP = enPacket%xP(1)
                      yP = enPacket%yP(1)
                      zP = enPacket%zP(1)
                      gP = 1
                      igpp = 1
                   end if

                   if ( (rVec%z <= grid(1)%zAxis(1) .or.zP<1) ) then
                      zP=1
                      rVec%z = grid(gP)%yAxis(1)
                      vHat%z = -vHat%z

                   end if

                   if ( (rVec%z <= grid(gP)%zAxis(1)-grid(gP)%geoCorrZ &
                        & .or.zP<1) .and. gP>1) then

                      xP = enPacket%xP(1)
                      yP = enPacket%yP(1)
                      zP = enPacket%zP(1)
                      gP = 1
                      igpp = 1

                   end if

                   if ( (rVec%z >=  grid(1)%zAxis(grid(gP)%nz) .or. zP>grid(gP)%nz) &
                        & ) then

                      zP = grid(gP)%nz
                      rVec%z = grid(gP)%zAxis(grid(gP)%nz)
                      vHat%z = -vHat%z

                   end if

                   if ((rVec%z >=  grid(gP)%zAxis(grid(gP)%nz)+grid(gP)%geoCorrZ &
                        & .or. zP>grid(gP)%nz) .and. gP>1) then

                      xP = enPacket%xP(1)
                      yP = enPacket%yP(1)
                      zP = enPacket%zP(1)
                      gP = 1
                      igpp = 1

                   end if


                   if (lgReturn) then

                      ! the packet escapes without further interaction
                      ! the packet escapes without further interaction
                      if (lgSymmetricXYZ) then
                         idirT = int(acos(abs(enPacket%direction%z))/dTheta)+1
                      else
                         idirT = int(acos(enPacket%direction%z)/dTheta)+1
                      end if

!                      idirT = int(acos(enPacket%direction%z)/dTheta)+1
                      if (idirT>totangleBinsTheta) then
                         idirT=totangleBinsTheta
                      end if
                      if (idirT<1 .or. idirT>totAngleBinsTheta) then
                         print*, '! pathSegment: error in theta direction cosine assignment',&
                              &  idirT, enPacket, dTheta, totAngleBinsTheta
                         stop
                      end if



                      if (abs(enPacket%direction%x)<1.e-35) then
                         idirP = 0
                      else
                         if (lgSymmetricXYZ) then
                            idirP = int(atan(abs(enPacket%direction%y)/abs(enPacket%direction%x))/dPhi)
                         else
                            idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                         end if
                      end if
                      if (idirP<0) idirP=totAngleBinsPhi+idirP

!                      if (enPacket%direction%x<1.e-35) then
!                         idirP = 0
!                      else
!                         idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
!                      end if
!                      if (idirP<0) idirP=totAngleBinsPhi+idirP

                      idirP=idirP+1

                      if (idirP>totangleBinsPhi) then
                         idirP=totangleBinsPhi
                      end if
                      if (idirP<1 .or. idirP>totAngleBinsPhi) then
                         print*, '! pathSegment: error in phi direction cosine assignment -4',&
                              &  idirP, enPacket, dPhi, totAngleBinsPhi
                         stop
                      end if


!print*, 'c', nanglebins
!print*, enPacket%lgStellar, enPacket%direction
!print*, idirT, viewPointPTheta(idirT), viewPointTheta(viewPointPTheta(idirT))
!print*, idirP, viewPointPPhi(idirP), viewPointPhi(viewPointPPhi(idirP)), viewPointPhi(viewPointPtheta(idirT))

                      if (nAngleBins>0) then
                         if (viewPointPtheta(idirT) > 0 .and. viewPointPhi(viewPointPtheta(idirT)) < 0) then

                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,viewPointPtheta(idirT)) = &
                                 &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,viewPointPtheta(idirT)) + deltaE(iStar)
                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,0) = &
                                 & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,0) +  deltaE(iStar)

                         elseif (viewPointPtheta(idirT) == viewPointPphi(idirP).or. &
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

                ! check if the path is still within the simulation region
                radius = 1.e10*sqrt((rVec%x/1.e10)*(rVec%x/1.e10) + &
                     &                     (rVec%y/1.e10)*(rVec%y/1.e10) + &
                     &                     (rVec%z/1.e10)*(rVec%z/1.e10))

                if (.not.lgPlaneIonization) then

                   if ( (.not.lgSymmetricXYZ .and. (rVec%x<=grid(1)%xAxis(1)-grid(1)%geoCorrX .or.&
                        & rVec%y<=grid(1)%yAxis(1)-grid(1)%geoCorrY .or. rVec%z<=grid(1)%zAxis(1)-&
                        & grid(1)%geoCorrZ)) .or.&
                        & (rVec%x >= grid(gP)%xAxis(grid(gP)%nx)+grid(gP)%geoCorrX) .or.&
                        &(rVec%y >= grid(gP)%yAxis(grid(gP)%ny)+grid(gP)%geoCorrY) .or.&
                        &(rVec%z >= grid(gP)%zAxis(grid(gP)%nz)+grid(gP)%geoCorrZ) .or. &
                        & xP>grid(gP)%nx .or. yP>grid(gP)%ny .or. zP>grid(gP)%nz ) then

                  if (gP==1) then

                     if (enPacket%xP(1) > grid(1)%nx) xP = grid(1)%nx
                     if (enPacket%yP(1) > grid(1)%ny) yP = grid(1)%ny
                     if (enPacket%zP(1) > grid(1)%nz) zP = grid(1)%nz
                     if (enPacket%xP(1) < 1) xP = 1
                     if (enPacket%yP(1) < 1) yP = 1
                     if (enPacket%zP(1) < 1) zP = 1


                     ! the energy packet escapes

                     ! the packet escapes without further interaction
                     ! the packet escapes without further interaction
                     if (lgSymmetricXYZ) then
                        idirT = int(acos(abs(enPacket%direction%z))/dTheta)+1
                     else
                        idirT = int(acos(enPacket%direction%z)/dTheta)+1
                     end if

!                     idirT = int(acos(enPacket%direction%z)/dTheta)+1
                     if (idirT>totangleBinsTheta) then
                        idirT=totangleBinsTheta
                     end if
                     if (idirT<1 .or. idirT>totAngleBinsTheta) then
                        print*, '! pathSegment: error in theta direction cosine assignment',&
                             &  idirT, enPacket, dTheta, totAngleBinsTheta
                        stop
                     end if

!                     if (enPacket%direction%x<1.e-35) then
!                        idirP = 0
!                     else
!                        idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
!                     end if
!                     if (idirP<0) idirP=totAngleBinsPhi+idirP


                     if (abs(enPacket%direction%x)<1.e-35) then
                        idirP = 0
                     else
                        if (lgSymmetricXYZ) then
                           idirP = int(atan(abs(enPacket%direction%y)/abs(enPacket%direction%x))/dPhi)
                        else
                           idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                        end if
                     end if
                     if (idirP<0) idirP=totAngleBinsPhi+idirP

                     idirP=idirP+1

                     if (idirP>totangleBinsPhi) then
                        idirP=totangleBinsPhi
                     end if

                     if (idirP<1 .or. idirP>totAngleBinsPhi) then
                        print*, '! pathSegment: error in phi direction cosine assignment -5',&
                             &  idirP, enPacket, dPhi, totAngleBinsPhi
                        stop
                     end if

!print*, 'd', nanglebins
!print*, enPacket%lgStellar, enPacket%direction
!print*, idirT, viewPointPTheta(idirT), viewPointTheta(viewPointPTheta(idirT))
!print*, idirP, viewPointPPhi(idirP), viewPointPhi(viewPointPPhi(idirP)), viewPointPhi(viewPointPtheta(idirT))


                     if (nAngleBins>0) then
                        if (viewPointPtheta(idirT) > 0 .and. viewPointPhi(viewPointPtheta(idirT)) < 0) then

                           grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,viewPointPtheta(idirT)) = &
                                &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                & enPacket%nuP,viewPointPtheta(idirT)) + deltaE(iStar)
                           grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                & enPacket%nuP,0) = &
                                & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                & enPacket%nuP,0) +  deltaE(iStar)
!print*, 'd1', enPacket%origin(1), enPacket%origin(2)

                        elseif (viewPointPtheta(idirT) == viewPointPphi(idirP).or. &
                             & (viewPointTheta(viewPointPphi(idirP))==viewPointTheta(viewPointPtheta(idirT))) .or. &
                             & (viewPointPhi(viewPointPtheta(idirT))==viewPointPhi(viewPointPphi(idirP))) ) then
                           grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,&
                                & viewPointPtheta(idirT)) = &
                                &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                & enPacket%nuP,viewPointPtheta(idirT)) +  deltaE(iStar)
!print*, 'd2'
                           if (viewPointPtheta(idirT)/=0) grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                enPacket%nuP,0) = &
                                & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                & enPacket%nuP,0) +  deltaE(iStar)
                        else

!print*, 'd3'
                           grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                enPacket%nuP,0) = &
                                & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                & enPacket%nuP,0) +  deltaE(iStar)
                        end if
                     else
!print*, 'd4'
                        grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                             enPacket%nuP,0) = &
                             & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                             & enPacket%nuP,0) +  deltaE(iStar)

                     end if

                     !b2.005
                     return

                  else if (gP>1) then

                     xP = enPacket%xP(1)
                     yP = enPacket%yP(1)
                     zP = enPacket%zP(1)
                     gP = 1
                     igpp = 1


                     if (gP/=1) then
                        print*, '! pathSegment: nested multigrids still not implemented'
                        stop
                     end if

                     if ( (radius >= R_out .and. R_out >= 0.) .or.&
                          & (rVec%x >= grid(1)%xAxis(grid(1)%nx)+grid(1)%geoCorrX) .or.&
                          &(rVec%y >= grid(1)%yAxis(grid(1)%ny)+grid(1)%geoCorrY) .or.&
                          &(rVec%z >= grid(1)%zAxis(grid(1)%nz)+grid(1)%geoCorrZ) .or. &
                          & (.not.lgSymmetricXYZ .and.  (rVec%x<=grid(1)%xAxis(1)-grid(1)%geoCorrX .or.&
                          & rVec%y<=grid(1)%yAxis(1)-grid(1)%geoCorrY .or. rVec%z<=grid(1)%zAxis(1)-&
                          & grid(1)%geoCorrZ))) then


                        if (xP > grid(gP)%nx) xP = grid(gP)%nx
                        if (yP > grid(gP)%ny) yP = grid(gP)%ny
                        if (zP > grid(gP)%nz) zP = grid(gP)%nz

                        if (xP < 1) xP = 1
                        if (yP < 1) yP = 1
                        if (zP < 1) zP = 1

                        ! the energy packet escapes

                        ! the packet escapes without further interaction
                        ! the packet escapes without further interaction
                        if (lgSymmetricXYZ) then
                           idirT = int(acos(abs(enPacket%direction%z))/dTheta)+1
                        else
                           idirT = int(acos(enPacket%direction%z)/dTheta)+1
                        end if

!                        idirT = int(acos(enPacket%direction%z)/dTheta)+1
                        if (idirT>totangleBinsTheta) then
                           idirT=totangleBinsTheta
                        end if
                        if (idirT<1 .or. idirT>totAngleBinsTheta) then
                           print*, '! pathSegment: error in theta direction cosine assignment',&
                                &  idirT, enPacket, dTheta, totAngleBinsTheta
                           stop
                        end if

!                        if (enPacket%direction%x<1.e-35) then
!                           idirP = 0
!                        else
!                           idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
!                        end if
!                        if (idirP<0) idirP=totAngleBinsPhi+idirP


                        if (abs(enPacket%direction%x)<1.e-35) then
                           idirP = 0
                        else
                           if (lgSymmetricXYZ) then
                              idirP = int(atan(abs(enPacket%direction%y)/abs(enPacket%direction%x))/dPhi)
                           else
                              idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
                           end if
                        end if
                        if (idirP<0) idirP=totAngleBinsPhi+idirP


                        idirP=idirP+1

                        if (idirP>totangleBinsPhi) then
                           idirP=totangleBinsPhi
                        end if

                        if (idirP<1 .or. idirP>totAngleBinsPhi) then
                           print*, '! pathSegment: error in phi direction cosine assignment -6',&
                                &  idirP, enPacket, dPhi, totAngleBinsPhi
                           stop
                        end if

                        if (nAngleBins>0) then
                           if (viewPointPtheta(idirT) > 0 .and. viewPointPhi(viewPointPtheta(idirT)) < 0) then

                              grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,viewPointPtheta(idirT)) = &
                                   &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,viewPointPtheta(idirT)) + deltaE(iStar)
                              grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,0) = &
                                   & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,0) +  deltaE(iStar)

                           elseif (viewPointPtheta(idirT) == viewPointPphi(idirP).or. &
                                &(viewPointTheta(viewPointPphi(idirP))==viewPointTheta(viewPointPtheta(idirT))).or. &
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

                     end if
                  else
                     print*, '! pathSegment: insanity occurred - invalid gP - ', gP
                     stop
                  end if

               end if

!            if (lgSymmetricXYZ .and. gP == 1) then
               if (lgSymmetricXYZ ) then
                  if (lgPlaneIonization) then
                     print*, '! pathSegment: lgSymmetric and lgPlaneionization flags both raised'
                     stop
                  end if

                  if ( rVec%x <= grid(1)%xAxis(1) .or. (gP==1 .and. xP<1)) then
                     if (vHat%x<0.) vHat%x = -vHat%x
                     enPacket%xP(1) = 1
                     xP = 1
                     rVec%x = grid(gP)%xAxis(1)
                  end if
                  if ( rVec%y <= grid(1)%yAxis(1) .or. (gP==1 .and. yP<1)) then
                     if (vHat%y<0.) vHat%y = -vHat%y
                     enPacket%yP(1)=1
                     yP = 1
                     rVec%y = grid(gP)%yAxis(1)
                  end if
                  if ( rVec%z <= grid(1)%zAxis(1) .or. (gP==1 .and. zP<1)) then
                     if (vHat%z<0.) vHat%z = -vHat%z
                     enPacket%zP(1) = 1
                     zP=1
                     rVec%z = grid(1)%zAxis(1)
                  end if

               end if

            end if

            if (gP>1) then
               if ( ( (rVec%x <= grid(gP)%xAxis(1) &
                    &.or. xP<1) .and. vHat%x <=0.) .or. &
                    & ( (rVec%y <= grid(gP)%yAxis(1) &
                    & .or. yP<1) .and. vHat%y <=0.) .or. &
                    & ( (rVec%z <= grid(gP)%zAxis(1) &
                    &  .or. zP<1) .and. vHat%z <=0.) .or. &
                    & ( (rVec%x >= grid(gP)%xAxis(grid(gP)%nx) &
                    &.or. xP>grid(gP)%nx) .and. vHat%x >=0.) .or. &
                    & ( (rVec%y >= grid(gP)%yAxis(grid(gP)%ny) &
                    & .or. yP>grid(gP)%ny) .and. vHat%y >=0.) .or. &
                    & ( (rVec%z >= grid(gP)%zAxis(grid(gP)%nz) &
                    &  .or. zP>grid(gP)%nz) .and. vHat%z >=0.) ) then

                  ! go back to mother grid
                  xP = enPacket%xP(1)
                  yP = enPacket%yP(1)
                  zP = enPacket%zP(1)
                  gP = 1
                  igpp = 1

               end if

            end if




         end if

         if (.not. lgPlaneIonization .and. gP==1 .and. (xP > grid(gP)%nx  &
              & .or. yP > grid(gP)%ny .or. zP > grid(gP)%nz) ) then

            ! the energy packet escapes

            ! the packet escapes without further interaction

            ! the packet escapes without further interaction
            if (lgSymmetricXYZ) then
               idirT = int(acos(abs(enPacket%direction%z))/dTheta)+1
            else
               idirT = int(acos(enPacket%direction%z)/dTheta)+1
            end if

!            idirT = int(acos(enPacket%direction%z)/dTheta)+1
            if (idirT>totangleBinsTheta) then
               idirT=totangleBinsTheta
            end if
            if (idirT<1 .or. idirT>totAngleBinsTheta) then
               print*, '! pathSegment: error in theta direction cosine assignment',&
                    &  idirT, enPacket, dTheta, totAngleBinsTheta
               stop
            end if

!            if (enPacket%direction%x<1.e-35) then
!               idirP = 0
!            else
!               idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
!            end if
!            if (idirP<0) idirP=totAngleBinsPhi+idirP

            if (abs(enPacket%direction%x)<1.e-35) then
               idirP = 0
            else
               if (lgSymmetricXYZ) then
                  idirP = int(atan(abs(enPacket%direction%y)/abs(enPacket%direction%x))/dPhi)
               else
                  idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
               end if
            end if
            if (idirP<0) idirP=totAngleBinsPhi+idirP

            idirP=idirP+1

            if (idirP>totangleBinsPhi) then
               idirP=totangleBinsPhi
            end if

            if (idirP<1 .or. idirP>totAngleBinsPhi) then
               print*, '! pathSegment: error in phi direction cosine assignment -7',&
                    &  idirP, enPacket, dPhi, totAngleBinsPhi
               stop
            end if


!print*, 'f', nanglebins
!print*, enPacket%lgStellar, enPacket%direction
!print*, idirT, viewPointPTheta(idirT), viewPointTheta(viewPointPTheta(idirT))
!print*, idirP, viewPointPPhi(idirP), viewPointPhi(viewPointPPhi(idirP)), viewPointPhi(viewPointPtheta(idirT))

            if (nAngleBins>0) then
               if (viewPointPtheta(idirT) > 0 .and. viewPointPhi(viewPointPtheta(idirT)) < 0) then

                  grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,viewPointPtheta(idirT)) = &
                       &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                       & enPacket%nuP,viewPointPtheta(idirT)) + deltaE(iStar)
                  grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                       & enPacket%nuP,0) = &
                       & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                       & enPacket%nuP,0) +  deltaE(iStar)

               elseif (viewPointPtheta(idirT) == viewPointPphi(idirP).or. &
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

      if (gP==1) then
         igpp = 1
      else if (gP>1) then
         igpp = 2
      else
         print*, "! pathSegment: insane grid index "
         stop
      end if

      enPacket%xP(igpp) = xP
      enPacket%yP(igpp) = yP
      enPacket%zP(igpp) = zP

      ! the energy packet has beenid absorbed - reemit a new packet from this position

      chTypeIn = packetType
      positionIn = rVec
      inX =  enPacket%xP(1:2)
      inY =  enPacket%yP(1:2)
      inZ =  enPacket%zP(1:2)
      gPIn = gP
      reRun = 1
      return

    end subroutine pathSegment


 subroutine hg(inpacket,ierr)

   implicit none

   type(photon_packet), intent(inout) :: inpacket

   real :: sint,cost,sinp,cosp,phi
   real :: random, random0, vin(3), vout(3), s, denom
   real :: hgg, znorm

   integer :: ierr

   ierr = 0

! BEKS 2010
   vin(1) = inpacket%direction%x
   vin(2) = inpacket%direction%y
   vin(3) = inpacket%direction%z

   hgg = gSca(inpacket%nuP)

   ! henyey-greenstein http://robertoreif.com/Documents/Chapter3.pdf
   call random_number(random0) ! this is the theta dependence
   s=2.*random0-1.
   if (hgg.ge.0.0001) then
      cost=0.5/hgg*(1.+hgg**2-((1.-hgg**2)/(1.+hgg*s))**2)
   else
      cost=s !+1.5*hgg*(1.-s**2)-2*hgg**2*s*(1.-s**2)
   endif
   if (cost .ge. 1.0) then
      cost=1.0
      sint=0.
   elseif (cost .lt. -1.0) then
      cost=-1.0
      sint=0.
   else
      sint=sqrt(1.-cost**2)
   endif

   call random_number(random) ! this is the phi dependence
   phi = twoPi*random
   cosp=cos(phi)
   sinp=sin(phi)
   denom=sqrt(1.-vin(3)**2)

   if (denom.gt.0.001) then
      vout(1)=sint/denom*(vin(1)*vin(3)*cosp-vin(2)*sinp) + vin(1)*cost
      vout(2)=sint/denom*(vin(2)*vin(3)*cosp+vin(1)*sinp) + vin(2)*cost
      vout(3)=-sint*cosp*denom+vin(3)*cost
   else
      vout(1)=sint*cosp
      vout(2)=sint*sinp
      if (vin(3).ge.0.) then
         vout(3)=cost
      else
         vout(3)=-cost
      endif
   endif



   if ( ( (abs(vout(1)) <= 1.) .and. (abs(vout(2)) <= 1.) .and. (abs(vout(3)) <= 1.) )&
        & .and. (vout(1) >= 0. .or. vout(1) < 0.) .and. &
        & (vout(2) >= 0. .or. vout(2) < 0.) .and. &
        & (vout(3) >= 0. .or. vout(3) < 0.)) then

      inpacket%direction%x =  vout(1)
      inpacket%direction%y =  vout(2)
      inpacket%direction%z =  vout(3)


      ierr = 0
   elseif ((abs(vout(1))>=1.).or.(abs(vout(2))>=1.).or.(abs(vout(3))>=1.)&
        & .and. (vout(1) >= 0. .or. vout(1) < 0.) .and. &
        & (vout(2) >= 0. .or. vout(2) < 0.) .and. &
        & (vout(3) >= 0. .or. vout(3) < 0.)) then
      znorm=sqrt(vout(1)**2 + vout(2)**2 + vout(3)**2)
      vout=vout/znorm
   else
      ierr = 1
      print*,'FAIL',random0,random,s
      print*,' ',vin
      print*,' ',cost,sint
      print*,' ',cosp,sinp,denom
      print*,' ',vout

      vout(1) = vin(1)
      vout(2) = vin(2)
      vout(3) = vin(3)

!      print*, vout
!      print*, vin
!      print*, random0, random
!      print*, hgg
!      print*, inpacket%nuP
!      print*, sint, cost, sinp, cosp

   end if

   return




! OLD
!   bmu= ((1.+g2)-((1.-g2)/(1.-hgg+2.*hgg*random))**2)/(2.*hgg)
!ßpathSe
!   cosb2=bmu**2
!   b=cosb2-1.
!
!   if(abs(bmu) > 1.) then
!      if(bmu>1.) then
!         bmu=1.
!         cosb2=1.
!         b=0.
!      else
!         bmu=-1.
!         cosb2=1.
!         b=0.
!      end if
!   end if
!   sinbt=sqrt(1.-cosb2)
!   call random_number(random)
!
!   ri1=twoPi*random
!
!   if(ri1>pi) then
!      ri3=twoPi-ri1
!      cosi3=cos(ri3)
!      sini3=sin(ri3)
!      sin2i3=2.*sini3*cosi3
!      cos2i3=2.*cosi3*cosi3-1.
!
!
!      if(abs(bmu)==1.) return
!
!      cost=costp*bmu+sintp*sinbt*cosi3
!
!      if(abs(cost)<1.) then
!         sint=abs(sqrt(1.-cost*cost))
!         sini2=sini3*sintp/sint
!         bott=sint*sinbt
!         cosi2=costp/bott-cost*bmu/bott
!      else
!         sint=0.
!         sini2=0.
!         if(cost.ge.1.)  cosi2=-1.
!         if(cost.le.-1.) cosi2=1.
!      end if
!
!      cosdph=-cosi2*cosi3+sini2*sini3*bmu
!      if(abs(cosdph)>1.) then
!         if(cosdph>1.) then
!            cosdph=1.
!         else
!            cosdph=-1.
!         end if
!      end if
!      phi=phip+acos(cosdph)
!      if(phi>twoPi) phi=phi-twoPi
!      if(phi<0.)    phi=phi+twoPi
!
!      sin2i2=2.*sini2*cosi2
!      cos2i2=2.*cosi2*cosi2-1.
!      sin2=sin2i2*sin2i3
!      cos2=cos2i2*cos2i3
!      sin2cos1=sin2i2*cos2i3
!      cos2sin1=cos2i2*sin2i3
!
!   else
!
!      cosi1=cos(ri1)
!      sini1=sin(ri1)
!      sin2i1=2.*sini1*cosi1
!      cos2i1=2.*cosi1*cosi1-1.
!
!      if(abs(bmu)==1.) return
!
!      cost=costp*bmu+sintp*sinbt*cosi1
!
!      if(abs(cost)<1.) then
!         sint=abs(sqrt(1.-cost*cost))
!         sini2=sini1*sintp/sint
!         bott=sint*sinbt
!         cosi2=costp/bott-cost*bmu/bott
!      else
!         sint=0.
!         sini2=0.
!         if(cost.ge.1.)  cosi2=-1.
!         if(cost.le.-1.) cosi2=1.
!      end if
!
!      cosdph=-cosi1*cosi2+sini1*sini2*bmu
!      if(abs(cosdph)>1.) then
!         if(cosdph>1.) then
!            cosdph=1.
!         else
!            cosdph=-1.
!         end if
!      end if
!      phi=phip-acos(cosdph)
!      if(phi>twoPi) phi=phi-twoPi
!      if(phi<0.)    phi=phi+twoPi
!
!      sin2i2=2.*sini2*cosi2
!      cos2i2=2.*cosi2*cosi2-1.
!      sin2=sin2i2*sin2i1
!      cos2=cos2i2*cos2i1
!      sin2cos1=sin2i2*cos2i1
!      cos2sin1=cos2i2*sin2i1
!
!   end if
!
!
!   cosp=cos(phi)
!   sinp=sin(phi)
!
!   nxp=sint*cosp
!   nyp=sint*sinp
!   nzp=cost
!
!
!   inpacket%direction%x =  nxp
!   inpacket%direction%y =  nyp
!   inpacket%direction%z =  nzp

 end subroutine hg

end subroutine energyPacketDriver


 end module photon_mod



