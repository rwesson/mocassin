! Copyright (C) 2007 Barbara Ercolano 
! 
! Version 3.00
module fluorescence_mod
    
    use common_mod
    use constants_mod
    use continuum_mod
    use grid_mod
    use interpolation_mod
    use pathIntegration_mod
    use vector_mod


    type(vector), parameter :: origin=vector(0.,0.,0.)  ! origin of the cartesian grid axes

    integer     , parameter :: safeLim = 10000          ! safety limit for the loops

    contains      
        
      recursive subroutine fluorescencePacketRun(grid, fType, rvec, xP, yP, zP, gP)
        implicit none
        
        
        real :: number 
        real, save :: ionPhot = 0.       

        type(vector), intent(inout)                  :: rVec        ! 

        integer, dimension(2), intent(inout)              :: xP, yP, &
             & zP                                            ! cartesian axes indeces 
        ! 1= mother; 2=sub
        
        integer, intent(inout)           :: gP               ! grid index
        integer                          :: igpr             ! grid pointer 1= mother 2=sub
        integer                          :: difSourceL(3)    ! cell indeces
        integer                          :: err              ! allocation error status
        integer                          :: i, j             ! counters  
        integer                          :: idirP, idirT     ! direction cosines

        type(grid_type), dimension(:), intent(inout) :: grid        ! the grid(s)
        
        type(photon_packet)              :: fluoPacket         ! the energu packet

        character(len=4), intent(inout)                   :: fType ! what resonance line?                
        
        
        if (fType == 'diffuse') then
           ! the energy packet has beenid absorbed - reemit a new packet from this position
           call energyPacketRun(ftype, rVec, enPacket%xP(1:2), enPacket%yP(1:2), &
                & enPacket%zP(1:2), gP)
           return
        end if

        if (gP==1) then
           igpr = 1
        else if (gP>1) then
           igpr = 2
        else
           print*,  "! fluorescencePacketRun: insane grid index"
           stop
        end if
           
        difSourceL(1) = xP(igpr)
        difSourceL(2) = yP(igpr)
        difSourceL(3) = zP(igpr)                

        fluoPacket = newFluorescencePacket(chType=fType, gP=gP, difSource=difSourceL, rvec=rvec)

        ! compute the next segment of trajectory
        call fluoPathSegment(fluoPacket)
        
        return


      contains

      
        ! this function initializes a photon packet
        function initFluorescencePacket(nuP,  position, lgLine, lgStellar, xP, yP, zP, gP)
          implicit none
          
          
          real                     :: random            ! random number
          type(vector), intent(in) :: position          ! the position at which the photon
          
          integer, intent(in)      :: nuP               ! the frequency of the photon packet
          integer, intent(in),dimension(2) :: xP, yP, &
               & zP                                     ! indeces of position on the x, y and z axes            
          integer, intent(in)      :: gP                ! grid index
          integer                  :: igpi              ! grid pointer 1=mother, 2=sub
          integer                  :: i, irepeat        ! counter
          
          type(photon_packet)      :: initFluorescencePacket  ! the photon packet
          
          logical, intent(in)      :: lgLine, lgStellar ! line, stellar packet?
          


          initFluorescencePacket%position = position
          
          initFluorescencePacket%iG  = gP
            
          if (gP==1) then
             igpi=1 
          else if (gp>1) then
             igpi=2
          else
             print*, "! initFluorescencePacket: insane gridp pointer"
             stop
          end if
           
          initFluorescencePacket%nuP      = nuP       
           
          initFluorescencePacket%lgStellar = lgStellar

          ! check if photon packen is line or continuum photon
          if ( lgLine ) then
             ! line photon
             initFluorescencePacket%nu       = 0.
             initFluorescencePacket%lgLine   = .true.
          else
             ! continuum photon
             initFluorescencePacket%nu       = nuArray(nuP)
             initFluorescencePacket%lgLine   = .false.          
          end if

          initFluorescencePacket%xP  = xP
          initFluorescencePacket%yP  = yP
          initFluorescencePacket%zP  = zP


          ! cater for plane parallel ionization case
          if (initFluorescencePacket%lgStellar .and. lgPlaneIonization) then
               
             ! get position
             
             ! x-direction
             call random_number(random)
             random = 1. - random
             initFluorescencePacket%position%x = &
                  & -(grid(gP)%xAxis(2)-grid(gP)%xAxis(1))/2. + random*( &
                  & (grid(gP)%xAxis(2)-grid(gP)%xAxis(1))/2.+&
                  & (grid(gP)%xAxis(grid(gP)%nx)-grid(gP)%xAxis(grid(gP)%nx-1))/2.+&
                  & grid(gP)%xAxis(grid(gP)%nx))
             if (initFluorescencePacket%position%x<grid(gP)%xAxis(1)) &
                  & initFluorescencePacket%position%x=grid(gP)%xAxis(1)
             if (initFluorescencePacket%position%x>grid(gP)%xAxis(grid(gP)%nx)) & 
                  initFluorescencePacket%position%x=grid(gP)%xAxis(grid(gP)%nx)
             
             call locate(grid(gP)%xAxis, initFluorescencePacket%position%x, initFluorescencePacket%xP(igpi))
             if (initFluorescencePacket%xP(igpi) < grid(gP)%nx) then
                if (initFluorescencePacket%position%x >= (grid(gP)%xAxis(initFluorescencePacket%xP(igpi))+&
                     & grid(gP)%xAxis(initFluorescencePacket%xP(igpi)+1))/2.) &
                     & initFluorescencePacket%xP(igpi) = initFluorescencePacket%xP(igpi)+1
             end if

             ! y-direction
             initFluorescencePacket%position%y = 0.
             initFluorescencePacket%yP(igpi) = 1
             
             ! z-direction
             call random_number(random)
             random = 1. - random
             initFluorescencePacket%position%z = -(grid(gP)%zAxis(2)-grid(gP)%zAxis(1))/2. + random*( &
                  & (grid(gP)%zAxis(2)-grid(gP)%zAxis(1))/2.+&
                  & (grid(gP)%zAxis(grid(gP)%nz)-grid(gP)%zAxis(grid(gP)%nz-1))/2.+&
                  & grid(gP)%zAxis(grid(gP)%nz))
             if (initFluorescencePacket%position%z<grid(gP)%zAxis(1)) & 
                  & initFluorescencePacket%position%z=grid(gP)%zAxis(1)
             if (initFluorescencePacket%position%z>grid(gP)%zAxis(grid(gP)%nz)) & 
                  & initFluorescencePacket%position%z=grid(gP)%zAxis(grid(gP)%nz)
             
             call locate(grid(gP)%zAxis, initFluorescencePacket%position%z, initFluorescencePacket%zP(igpi))
             if (initFluorescencePacket%zP(igpi) < grid(gP)%nz) then               
                if (initFluorescencePacket%position%z >= (grid(gP)%xAxis(initFluorescencePacket%zP(igpi))+&
                     & grid(gP)%zAxis(initFluorescencePacket%zP(igpi)+1))/2.) initFluorescencePacket%zP(igpi) =& 
                     & initFluorescencePacket%zP(igpi)+1
             end if

             if (initFluorescencePacket%xP(igpi)<1) initFluorescencePacket%xP(igpi)=1             
             if (initFluorescencePacket%zP(igpi)<1) initFluorescencePacket%zP(igpi)=1
             
             ! direction is parallel to y-axis direction
             initFluorescencePacket%direction%x = 0.
             initFluorescencePacket%direction%y = 1.
             initFluorescencePacket%direction%z = 0.
             
             if (initFluorescencePacket%xP(igpi) >  grid(gP)%xAxis(grid(gP)%nx) .or. &
                  & initFluorescencePacket%zP(igpi) >  grid(gP)%zAxis(grid(gP)%nz)) then
                print*, "! initFluorescencePacket: insanity in planeIonisation init"
                print*, igpi, initFluorescencePacket%xP(igpi),  grid(gP)%xAxis(grid(gP)%nx), &
                     & initFluorescencePacket%zP(igpi), grid(gP)%zAxis(grid(gP)%nz),  &
                     &random, initFluorescencePacket%position%z
                
                stop
             end if
             
             planeIonDistribution(initFluorescencePacket%xP(igpi),initFluorescencePacket%zP(igpi)) = &
                  & planeIonDistribution(initFluorescencePacket%xP(igpi),initFluorescencePacket%zP(igpi)) + 1
             
          else

             do irepeat = 1, 1000000
                ! get a random direction
                initFluorescencePacket%direction = randomUnitVector()
                if (initFluorescencePacket%direction%x/=0. .and. & 
                     & initFluorescencePacket%direction%y/=0. .and. & 
                     & initFluorescencePacket%direction%z/=0.) exit
             end do
          end if

          initFluorescencePacket%origin(1) = gP
          initFluorescencePacket%origin(2) = grid(gP)%active(initFluorescencePacket%xP(igpi),&
               & initFluorescencePacket%yP(igpi), initFluorescencePacket%zP(igpi))                           
          
        end function initFluorescencePacket

          
        ! this function creates a new fluorescence packet
        function newFluorescencePacket(chType, gP, difSource, rvec)
          
            type(photon_packet)                :: newFluorescencePacket! the photon packet to be created
            type(vector), intent(in)           :: rVec           ! 
            real                               :: random         ! random number
            
            type(vector)                       :: positionLoc    ! the position of the photon
                                                                 ! packet
            integer                            :: nuP            ! the frequency index of the photon packet
            integer, dimension(2)              :: orX,orY,orZ    ! dummy
            integer, intent(in)                :: difSource(3)   ! grid and cell indeces
            integer, intent(inout)             :: gP
            integer                            :: igpn           ! grid pointe 1=motehr, 2=sub

            character(len=4), intent(in)       :: chType         ! what line?




            if (gP==1) then
               igpn = 1
            else if (gp>1) then
               igpn = 2
            else
               print*,  "! newFluorescencePacket: insane grid pointer"
               stop
            end if

            positionLoc%x = rVec%x
            positionLoc%y = rVec%y
            positionLoc%z = rVec%z

            ! initialize the new fluorescence packet
            orX(igPn) = difSource(1)
            orY(igPn) = difSource(2)
            orZ(igPn) = difSource(3)

            ! initialize the new fluorescence packet
            if (grid(gP)%active(orX(igpn), orY(igpn), orZ(igpn)) < 0.) then
               print*, "! newFluorescencePacket: new packet cannot be emitted from re-mapped cell -3-"
               print*, "chType, nuP,  .false., .true., orX,orY,orZ, gp"
               print*, chType, nuP,  .false., .true., orX,orY,orZ, gp
               stop
            end if
            
            select case (chType)
            ! if the photon is stellar
            case ("FeKa")

               newFluorescencePacket = initFluorescencePacket(FeKaColdP, positionLoc, .false., .false., orX,&
                    & orY, orZ, gP)

            case ("FeL1")

               newFluorescencePacket = initFluorescencePacket(FeL1P, positionLoc, .false., .false., orX,&
                    & orY, orZ, gP)
            case ("FeL2")

               newFluorescencePacket = initFluorescencePacket(FeL2P, positionLoc, .false., .false., orX,&
                    & orY, orZ, gP)

            case ("CKa")

               newFluorescencePacket = initFluorescencePacket(CKaP, positionLoc, .false., .false., orX,&
                    & orY, orZ, gP)

            case ("NKa")

               newFluorescencePacket = initFluorescencePacket(NKaP, positionLoc, .false., .false., orX,&
                    & orY, orZ, gP)

            case ("OKa")

               newFluorescencePacket = initFluorescencePacket(OKaP, positionLoc, .false., .false., orX,&
                    & orY, orZ, gP)

            case ("NeKa")

               newFluorescencePacket = initFluorescencePacket(NeKaP, positionLoc, .false., .false., orX,&
                    & orY, orZ, gP)

            case ("MgKa")

               newFluorescencePacket = initFluorescencePacket(MgKaP, positionLoc, .false., .false., orX,&
                    & orY, orZ, gP)

            case ("AlKa")

               newFluorescencePacket = initFluorescencePacket(AlKaP, positionLoc, .false., .false., orX,&
                    & orY, orZ, gP)

            case ("SiKa")

               newFluorescencePacket = initFluorescencePacket(SiKaP, positionLoc, .false., .false., orX,&
                    & orY, orZ, gP)

            case ("SKa")

               newFluorescencePacket = initFluorescencePacket(SKaP, positionLoc, .false., .false., orX,&
                    & orY, orZ, gP)

            case ("ArKa")

               newFluorescencePacket = initFluorescencePacket(ArKaP, positionLoc, .false., .false., orX,&
                    & orY, orZ, gP)

            case ("CaKa")

               newFluorescencePacket = initFluorescencePacket(CaKaP, positionLoc, .false., .false., orX,&
                    & orY, orZ, gP)

            case ("NiKa")

               newFluorescencePacket = initFluorescencePacket(NiKaP, positionLoc, .false., .false., orX,&
                    & orY, orZ, gP)

            ! if the photon packet type is wrong or missing
            case default
        
                print*, "! newFluorescencePacket: unsupported fluorescent transition.", chType
                stop
                
             end select
             
           end function newFluorescencePacket
            
           subroutine fluoPathSegment(enPacket)
             implicit none


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

             integer                         :: idirT,idirP ! direction cosine counters
             integer                         :: i, j, nS ! counter
             integer                         :: xP,yP,zP ! cartesian axes indeces
             integer                         :: gP       ! grid index
             integer                         :: igpp     ! grid index 1=mother 2=sub
             integer                         :: safeLimit =1000 ! safe limit for the loop

             type(photon_packet), intent(inout) :: enPacket ! the energy packet

             character(len=7)                :: packetType ! what line?

             logical                         :: lgScattered ! is the packet scattering with dust?
             logical                         :: lgReturn

             if (enPacket%iG == 1) then
                igpp = 1
             else if (enPacket%iG>1) then
                igpp = 2
             else 
                print*, "! fluoPathSegment: insane grid index"
                stop
             end if

             ! check that the input position is not outside the grid
             if ( (enPacket%iG <= 0).or.(enPacket%iG > nGrids) ) then   
                print*, "! fluoPathSegment: starting position not in any defined gridhses",&
                     & enPacket
                stop
             else if ( (enPacket%xP(igpp) <= 0).or.&
                  &(enPacket%xP(igpp) > grid(enPacket%iG)%nx) ) then
                print*, "! fluoPathSegment: starting position in x is outside the grid",&
                     & enPacket
                stop
             else if ( (enPacket%yP(igpp) <= 0).or. & 
                  & (enPacket%yP(igpp) > grid(enPacket%iG)%ny) ) then
                print*, "! fluoPathSegment: starting position in y is outside the grid",&
                     & enPacket
                stop
             else if ( (enPacket%zP(igpp) <= 0).or.& 
                  & (enPacket%zP(igpp) > grid(enPacket%iG)%nz) ) then   
                print*, "! fluoPathSegment: starting position in z is outside the grid",&
                     & enPacket
                stop          
             end if

             ! define vHat and rVec
             rVec = enPacket%position
             vHat = enPacket%direction
             
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
                   print*, " ! fluorescencePacketRun: multiple grids are not allowed in a 1D simulation"
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
!             safeLimit=500000
                safeLimit=1000
             end if

             do i = 1, safeLimit

                do j = 1, safeLimit

                   if (xP > grid(gP)%nx .or. xP < 1 .or. &
                        & yP > grid(gP)%ny .or. yP < 1 .or. &
                        & zP > grid(gP)%nz .or. zP < 1 ) then
                      print*, "! fluoPathSegment: insanity [gp,xp,yp,zp,j,i]", & 
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
                   
                   end if

                   enPacket%iG = gP
                   igpp = 2

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
                      
                      if (xP > grid(gP)%nx .or. xP < 1 .or. &
                           & yP > grid(gP)%ny .or. yP < 1 .or. &
                           & zP > grid(gP)%nz .or. zP < 1 ) then
                         print*, "! fluoPathSegment: insanity -2- [gp,xp,yp,zp]", & 
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
                   print*, '! fluoPathSegment: [warning] dSx <= 0.',dSx
                   print*, 'grid(gP)%xAxis ', grid(gP)%xAxis
                   print*, 'gP,xP,grid(gP)%xAxis(xP), rVec%x, vHat%x'
                   print*, gP,xP,grid(gP)%xAxis(xP), rVec%x, vHat%x
                   dS = amin1(dSy, dSz)
                else if (dSy<=0.) then
                   print*, '! fluoPathSegment: [warning] dSy <= 0.', dSy
                   print*, 'grid(gP)%yAxis ', grid(gP)%yAxis
                   print*, 'gP,yP,grid(gP)%yAxis(yP), rVec%y, vHat%y'
                   print*, gP,yP,grid(gP)%yAxis(yP), rVec%y, vHat%y
                   dS = amin1(dSx, dSz)
                else if (dSz<=0.) then
                   print*, '! fluoPathSegment: [warning] dSz <= 0.', dSz
                   print*, 'grid(gP)%zAxis ', grid(gP)%zAxis
                   print*, 'gP,zP,grid(gP)%zAxis(zP), rVec%z, vHat%z'
                   print*, gP,zP,grid(gP)%zAxis(zP), rVec%z, vHat%z
                   dS = amin1(dSx, dSy)
                else
                   dS = min(dSx,dSy)
                   dS = min(dS, dSz)
                end if

                ! this should now never ever happen
                if (dS <= 0.) then
                   print*, 'fluoPathSegment: dS <= 0', dSx, dSy, dSz
                   print*, gP, rVec
                   stop
                end if

                ! calculate the optical depth to the next cell wall 
                tauCell = dS*grid(gP)%opacity(grid(gP)%active(xP,yP,zP), enPacket%nuP)

                if (lg1D) then
                   if (nGrids>1) then
                      print*, '! getVolumeLoc: 1D option and multiple grids options are not compatible'
                      stop
                   end if

                   if (xP == 1) then              

                      dV = 4.*Pi* ( (grid(gP)%xAxis(xP+1)/1.e15)**3.)/3.


                   else if ( xP==grid(gP)%nx) then
                   
                      dV = Pi* ( (3.*(grid(gP)%xAxis(xP)/1.e15)-(grid(gP)%xAxis(xP-1)/1.e15))**3. - &
                           & ((grid(gP)%xAxis(xP)/1.e15)+(grid(gP)%xAxis(xP-1)/1.e15))**3. ) / 6.

                   else 

                      dV = Pi* ( ((grid(gP)%xAxis(xP+1)/1.e15)+(grid(gP)%xAxis(xP)/1.e15))**3. - &
                           & ((grid(gP)%xAxis(xP-1)/1.e15)+(grid(gP)%xAxis(xP)/1.e15))**3. ) / 6.
                      
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

                   ! packet interacts
                
                   ! calculate where within this cell the packet interacts
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
                           grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dlLoc*deltaEUsed / dV
                   else ! if the energy packet is diffuse
                      if (lgDebug) then
                         grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                              & grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dlLoc*deltaEUsed / dV
                      else
                         grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                              & grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dlLoc*deltaEUsed / dV
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
                         print*, '! fluorescencePacketRun: error in theta direction cosine assignment',&
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
                         print*, '! fluorescencePacketRun: error in Phi direction cosine assignment',&
                              &  idirP, enPacket, dPhi, totAngleBinsPhi
                         stop
                      end if
                      
                      if (nAngleBins>0) then
                         if (viewPointPtheta(idirT) == viewPointPphi(idirP).or. &
                              & (viewPointTheta(viewPointPphi(idirP))==viewPointTheta(viewPointPtheta(idirT))) .or. & 
                              & (viewPointPhi(viewPointPtheta(idirT))==viewPointPhi(viewPointPphi(idirP))) ) then
                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,viewPointPtheta(idirT)) = &
                                 &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                            if (viewPointPtheta(idirT)/=0) grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                                 & enPacket%nuP,0) = &
                                 & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,0) +  deltaEUsed
                         else
                            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                                 & enPacket%nuP,0) = &
                                 & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                 & enPacket%nuP,0) +  deltaEUsed
                         end if
                      else
                         
                         grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                              & enPacket%nuP,0) = &
                              & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                              & enPacket%nuP,0) +  deltaEUsed
                         
                      end if

                      return
                   end if
                

                   ! check if the packet is absorbed or scattered 
                   ! NOTE : we are ignoring dust
                   probSca = KNsigmaT(enPacket%nuP)*&
                        &grid(gP)%Ne(grid(gP)%active(xp,yp,zp))/&
                        & (grid(gp)%opacity(grid(gp)%active(xp,yp,zp),&
                        &enPacket%nuP))
                   
                   call random_number(random)
                   
                   random = 1.-random
                   
                   if (random > probSca) then
                      lgScattered = .false.         
                   else if (random <= probSca) then
                      lgScattered = .true.         
                   else
                      print*, '! fluoPathSegment: insanity occured and scattering/absorption &
                           & decision stage.'
                      stop
                   end if

                   if (.not. lgScattered) then

                      absInt = absInt + 1.                            
                      exit
                      
                   else
                      
                      scaInt = scaInt + 1.                                                 
                      ! packet is compton scattered                         
                      ! calculate new direction
                      ! for now assume scattering is isotropic, when phase
                      ! function is introduced the following must be changed                         
                         
                      enPacket%xP(igpp) = xP
                      enPacket%yP(igpp) = yP
                      enPacket%zP(igpp) = zP            
                      
                      if (grid(gP)%active(enPacket%xp(igpp), enPacket%yp(igpp), enPacket%zp(igpp)) < 0.) then
                         print*, "! fluoPathSegment: new packet cannot be emitted from re-mapped cell -1-"
                         print*, "nuP, .false., .true., xp,yp,zp, gp"
                         print*, nuP, .false., .true.,  xp,yp,zp, gp
                         stop
                      end if

                      call comptonScatter(enPacket)
                                            
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

                   ! the packet is not absorbed within this cell
                   ! add contribution of the packet to the radiation field
                   
                   if (enPacket%lgStellar) then
                      grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                           grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dS*deltaEUsed / dV
                   else ! if the energy packet is diffuse
                      if (lgDebug) then
                         grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                              & grid(gP)%Jdif(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dS*deltaEUsed / dV
                      else
                         grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) = &
                              & grid(gP)%Jste(grid(gP)%active(xP,yP,zP),enPacket%nuP) + dS*deltaEUsed / dV
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
                           & enPacket%zP(1) <= 0 ) then
                         
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
!                            print*, '! fluoPathSegment: insanity occured at mother grid transfer (x axis +)', & 
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
!                            print*, '! fluoPathSegment: insanity occured at mother grid transfer (x axis-)',&  
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
!                            print*, '! fluoPathSegment: insanity occured at mother grid transfer (y axis +)',&
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
!                            print*, '! fluoPathSegment: insanity occured at mother grid transfer (y axis -)', & 
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
!                            print*, '! fluoPathSegment: insanity occured at mother grid transfer (z axis +)', &
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
!                        print*, '! fluoPathSegment: insanity occured at mother grid transfer (z axis -)', &
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
                         print*, '! fluoPathSegment: insanity occured in dS assignement &
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
                            print*, '! fluoPathSegment: insanity occured - invalid gP', gP
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
                            print*, '! fluoPathSegment: insanity occured - invalid gP', gP
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
                        idirT = int(acos(enPacket%direction%z)/dTheta)+1
                        if (idirT>totangleBinsTheta) then
                           idirT=totangleBinsTheta
                        end if
                        if (idirT<1 .or. idirT>totAngleBinsTheta) then
                           print*, '! fluorescencePacketRun: error in theta direction cosine assignment',&
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
                           print*, '! fluorescencePacketRun: error in phi direction cosine assignment',&
                                &  idirP, enPacket, dPhi, totAngleBinsPhi
                           stop
                        end if
                        
                        
                        if (nAngleBins>0) then
                           if (viewPointPtheta(idirT) == viewPointPphi(idirP).or. &
                                & (viewPointTheta(viewPointPphi(idirP))==viewPointTheta(viewPointPtheta(idirT))) .or. &
                                & (viewPointPhi(viewPointPtheta(idirT))==viewPointPhi(viewPointPphi(idirP))) ) then
                              grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,& 
                                   & viewPointPtheta(idirT)) =grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2),& 
                                   & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                              
                              if (viewPointPtheta(idirT)/=0) grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   &enPacket%nuP,0) = &
                                   & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,0) +  deltaEUsed
                              
                           else
                              grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,0) = &
                                   & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,0) +  deltaEUsed
                              
                           end if
                           
                        else
                           
                           grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                enPacket%nuP,0) = &
                                & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                & enPacket%nuP,0) +  deltaEUsed
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
                           print*, '! fluoPathSegment: insanity occured - invalid gP', gP
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
                           print*, '! fluoPathSegment: insanity occured - invalid gP', gP
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
                                   & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                              if (viewPointPtheta(idirT)/=0) grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                                   enPacket%nuP,0) = &
                                   & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,0) +  deltaEUsed
                              
                           else
                              grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                                enPacket%nuP,0) = &
                                & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                & enPacket%nuP,0) +  deltaEUsed
                              
                           end if
                           
                        else

                           grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                                enPacket%nuP,0) = &
                                & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                & enPacket%nuP,0) +  deltaEUsed   
                           
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
                                      & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                                 if (viewPointPtheta(idirT)/=0) grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                                      enPacket%nuP,0) = &
                                      & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                      & enPacket%nuP,0) +  deltaEUsed
                              else
                                 grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                                      enPacket%nuP,0) = &
                                      & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                      & enPacket%nuP,0) +  deltaEUsed
                              end if
                           else
                              
                              grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                                   enPacket%nuP,0) = &
                                   & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                   & enPacket%nuP,0) +  deltaEUsed
                              
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
                              print*, '! fluoPathSegment: nested multigrids still not implemented'
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
                                 print*, '! energyPacketRun: error in phi direction cosine assignment',&
                                      &  idirP, enPacket, dPhi, totAngleBinsPhi
                                 stop
                              end if
                              
                              if (nAngleBins>0) then
                                 if (viewPointPtheta(idirT) == viewPointPphi(idirP).or. &
                                      &(viewPointTheta(viewPointPphi(idirP))==viewPointTheta(viewPointPtheta(idirT))).or. & 
                                      & (viewPointPhi(viewPointPtheta(idirT))==viewPointPhi(viewPointPphi(idirP))) ) then
                                    grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,& 
                                         & viewPointPtheta(idirT)) = &
                                         &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                         & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                                    if (viewPointPtheta(idirT)/=0) grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                                         enPacket%nuP,0) = &
                                         & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                         & enPacket%nuP,0) +  deltaEUsed
                                 else
                                    grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                                         enPacket%nuP,0) = &
                                         & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                         & enPacket%nuP,0) +  deltaEUsed
                                 end if
                              else
                                 
                                 grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                                      enPacket%nuP,0) = &
                                      & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                                      & enPacket%nuP,0) +  deltaEUsed
                                 
                              end if

                              !b2.005
                              return

                           end if
                        else
                           print*, '! fluoPathSegment: insanity occured - invalid gP - ', gP
                           stop
                        end if
                        
                     end if

                     if (lgSymmetricXYZ ) then
                         if (lgPlaneIonization) then
                            print*, '! fluoPathSegment: lgSymmetric and lgPlaneionization flags both raised'
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
                   & enPacket%nuP,viewPointPtheta(idirT)) +  deltaEUsed
                 if (viewPointPtheta(idirT)/=0) grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                      enPacket%nuP,0) = &
                      & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                      & enPacket%nuP,0) +  deltaEUsed
              else
                 grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                      enPacket%nuP,0) = &
                      & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                      & enPacket%nuP,0) +  deltaEUsed
              end if
           else
              
              grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                   enPacket%nuP,0) = &
                   & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                   & enPacket%nuP,0) +  deltaEUsed
              
           end if

           !b2.005
            return
            
            
         end if
         
      end do ! safelimit loop
      
      if (i>= safeLimit) then
         if (.not.lgPlaneIonization) then
            print*, '! fluoPathSegment: [warning] packet trajectory has exceeded&
                 &  maximum number of events', safeLimit, gP, xP,yP,zP, grid(gP)%active(xP,yP,zP), & 
              & rvec, vhat, enPacket
         end if
         return
         
      end if
     
      if (gP==1) then
         igpp = 1
      else if (gP>1) then
         igpp = 2 
      else
         print*, "! fluoPathSegment: insane grid index "
         stop
      end if
      
      enPacket%xP(igpp) = xP
      enPacket%yP(igpp) = yP
      enPacket%zP(igpp) = zP   
      
      packetType = 'diffuse'
      
      ! the energy packet has beenid absorbed - reemit a new packet from this position
      call fluorescencePacketRun(grid, packetType, rvec, enPacket%xP(1:2), enPacket%yP(1:2), &
           & enPacket%zP(1:2), gP)

      return
   
 end subroutine fluoPathSegment

 ! determine energy and direction for compton redistribution of fluorescent lines
 subroutine comptonScatter(fPacket)
   implicit none

   
   real                              :: newNu,newTheta     ! new energy and theta
   real                              :: random             ! random number
   real                              :: u,v,w,t            ! direction units
   type(photon_packet),intent(inout) :: fPacket            ! the fluorescent packet

   integer                           :: newNuP,newThetaP   ! pointer to new energy and theta

   ! calculate new direction
   call getNu2(KNsigmaArray(fPacket%nuP,1:180),newThetaP)     
   newTheta = real(newThetaP)*Pi/180.

   ! calculate u,v,w (Harries & Howarth, 1997)
   call random_number(random)
   w = 2.*random - 1.
   t = sqrt(1-w*w)
   u = t*cos(newTheta)
   v = t*sin(newTheta)

   fPacket%direction%x = u
   fPacket%direction%y = v
   fPacket%direction%z = w      

   ! calculate new energy
   newNu = fPacket%nu*PcompArray(fPacket%nuP,newThetaP) 

   ! find this in the nuArray 
   call locate(nuArray,newNu,newNuP)

   fPacket%Nu  = newNu     
   fPacket%NuP = newNuP

 end subroutine comptonScatter
 
 ! this subroutine determines the index in the PDFs
 ! same as getNu2
 subroutine getNu2(probDen, nuP)
   
   real, dimension(:), intent(in) :: probDen    ! probability density function
   
   real                           :: random     ! random number
   
   integer, intent(out)           :: nuP        ! frequency index of the new
   
   integer                        :: isearch,i  !  
   
   ! get a random number
   call random_number(random)
   
   do i = 1, 10000
      if (random==0 .or. random==1.) then
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
   
   if (nuP>=nbins) then
      print*, 'random: ', random
      print*, 'probDen: ', probDen
   end if
   
 end subroutine getNu2
 
end subroutine fluorescencePacketRun

        
end module fluorescence_mod


