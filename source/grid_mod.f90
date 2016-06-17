! Copyright (C) 2005 Barbara Ercolano 
!
! Version 2.02
module grid_mod

    use common_mod             ! common variablesesca
    use composition_mod        ! cemical abundances 
    use constants_mod          ! physical constants
    use continuum_mod          ! ionising field
    use elements_mod           ! hydrogen data 
    use interpolation_mod      ! interpolation maths
    use pathIntegration_mod    ! path integration
    use set_input_mod          ! model inputs
    use vector_mod             ! vectyor maths
    use xSec_mod               ! x Section data


    contains

    ! initCartesianGrid initializes a cartesian grid
    subroutine initCartesianGrid(grid, nx, ny, nz)
        implicit none

        integer, intent(in) :: nx, ny, nz           ! x, y and z arrays sizes

        type(grid_type),intent(inout) :: grid       ! grid


        ! local variables

        integer :: err, ios                         ! allocation error status
        integer :: ii
        integer :: i, j, iCount, nuCount, elem, ion ! counters
        integer :: nradio
        integer :: g0,g1
        integer :: nEdges  
        integer :: nElec
        integer :: outshell

        logical, save :: lgFirst=.true.
        logical       :: lgAssigned

        integer, parameter :: maxLim = 10000
        integer, parameter :: nSeries = 17

        real, dimension(450) :: ordered
        real, dimension(17)  :: seriesEdge

        real                 :: nuMinArray, nuMaxArray        
        real                 :: dradio
        real                 :: nuStepSizeLoc

        real, pointer        :: nuTemp(:)

        print*, "in initCartesianGrid"
        
        nullUnitVector%x = 1.
        nullUnitVector%y = 0.
        nullUnitVector%z = 0.


        elemLabel = (/'*H','He','Li','Be','*B','*C','*N','*O','*F','Ne','Na','Mg',&
             &'Al','Si','*P','*S','Cl','Ar','*K','Ca','Sc','Ti','*V','Cr','Mn','Fe',&
             &'Co','Ni','Cu','Zn'/)

       ! set up atomic weight array
       aWeight = (/1.0080, 4.0026, 6.941, 9.0122, 10.811, 12.0111, 14.0067, 15.9994, &
            & 18.9984, 20.179, 22.9898, 24.305, 26.9815, 28.086, 30.9738, 32.06, 35.453, &
            & 39.948, 39.102, 40.08, 44.956, 47.90, 50.9414, 51.996, 54.9380, 55.847, 58.9332, &
            & 58.71, 63.546, 65.37 /)

        ! only calculate the frequency/energy grid and assign the abundances if this is the
        ! first time that 
        ! this procedure is called - the grid will be constant throughout the execution

        allocate(grid%elemAbun(nAbComponents, nElements), stat = err)
        if (err /= 0) then 
           print*, "! initCartesianGrid: can't allocate grid memory"
            stop
        end if
        grid%elemAbun=0.

        if (lgFirst) then  

           if (lgGas) then
              ! set chemical abundances according to the grid%composition
              ! variable  given 
              allocate(forbiddenLines(nElements,nstages, nForLevels,nForLevels), stat=err)
              if (err /= 0) then
                 print*, "! emissionDriver: can't allocate array memory"
                 stop
              end if
              allocate(forbiddenLinesLarge(nForLevelsLarge,nForLevelsLarge), stat=err)
              if (err /= 0) then
                 print*, "! emissionDriver: can't allocate array memory"
                 stop
              end if

              call setComposition(grid)

           else

              lgElementOn=.false.
              elementXref=0
              nElementsUsed=0

           end if
 
            ! calculate the frequency grid

            ! allocate just enough space for the energy array
            allocate(nuArray(1:nbins), stat = err)
            if (err /= 0) then
                print*, "! grid: can't allocate grid memory"
                stop
            end if
 
            ! allocate just enough space for the widFlx array
            allocate(widFlx(1:nbins), stat = err)
            if (err /= 0) then
                print*, "! grid: can't allocate grid memory"
                stop
            end if

            allocate(continuum(1:nbins), stat = err)
            if (err /= 0) then
                print*, "! grid: can't allocate continuum emission array memory"
                stop
            end if
            
            if (lgGas) then

               call phInit() ! set up photoionization data
               
               seriesEdge = (/0.0069, 0.0083, 0.01, 0.0123, 0.0156, 0.0204, 0.0278, 0.04, &
                    & 0.0625, 0.11117, &
                    & 0.11610, 0.12248, 0.13732, 0.24763, 0.24994, 0.26630, 0.29189/)

               nEdges = 1
               
               ! find the ionization edges included in our frequency range
               do elem = 1, nElements  ! begin element loop
                  do ion = 1, min(elem, nStages-1) ! begin ion loop
                     if(.not.lgElementOn(elem)) exit
                         
                     if (elem > 2) then

                        ! find the number of electrons in this ion
                        nElec = elem - ion + 1

                        ! find the outer shell number and statistical weights 
                        call getOuterShell(elem, nELec, outShell, g0, g1)
                         
                        ! get threshold energy
                        ionEdge(nEdges) = ph1(1,elem, nElec, outShell)/RydtoEv
                     else if (elem == 1) then ! HI
            
                        ionEdge(nEdges) = 0.999434
                        
                     else if ( (elem == 2) .and. (ion == 1) ) then ! HeI
 
                        ionEdge(nEdges) = 1.80804
                        
                     else if ( (elem == 2) .and. (ion == 2) ) then ! HeII

                        ionEdge(nEdges) = 4.
                     
                     end if

                     if (ionEdge(nEdges) <= nuMax) nEdges = nEdges + 1


                  end do
               end do

               nEdges = nEdges -1

               call sortUp(ionEdge(1:nEdges))

            end if


            if ( (lgDust) .and. (.not.lgGas) ) then 
               close(72)
               open (unit= 72,  action="read", file=trim(home)//"dustData/nuDustRyd.dat", status = "old", position = "rewind", &
                    & iostat = ios)
               if (ios /= 0) then
                  print*, "! initCartesianGrid: can't open dust nu grid file - ",trim(home),"dustData/nuDustRyd.dat"
                  stop
               end if

               do i = 1, 10000000
                  if (i<=nbins+1) then
                     read(unit=72,fmt=*,iostat=ios) nuArray(i)
                     if (nuArray(i)>nuMax) exit ! nuMax reached
                     if (ios < 0) exit ! end of file reached
                  else
                     print*, "! initCartesianGrid: nbins is smaller that the number of &
                          & frequency points in dustData/nuDustRyd.dat file - enlarge nbins"
                     stop
                  end if
               end do
               
               nbins = i-1
               print*, "! initCartesianGrid: nbins reset to ", nbins

               close(72)
               
               allocate (nuTemp(1:nbins))
               nuTemp = nuArray(1:nbins)
               if (associated(nuArray)) deallocate(nuArray)
               allocate (nuArray(1:nbins))
               nuArray = nuTemp
               if (associated(nuTemp)) deallocate(nuTemp)

            else if (lgGas .and. (.not.lgDust)) then

               ! first count how many edge, thresholds etc
               if (nuMin<radio4p9GHz) then
                  nuArray(1)     = radio4p9GHz
                  nuCount = 2
               else
                  nuCount = 1
               end if

               nuMinArray = nuMin
               nuMaxArray = nuMax
               ! H series edges 
               do i = 1, nSeries 
                  nuArray(nuCount) = seriesEdge(i)
                  nuArray(nuCount+1) = seriesEdge(i)- 0.0003
                  nuArray(nuCount+2) = seriesEdge(i)+ 0.0003
                  if ( nuArray(nuCount) < nuMinArray ) nuMinArray= nuArray(nuCount)- 0.0006
                  if ( nuArray(nuCount) > nuMaxArray ) then
                     nuMaxArray= nuArray(nuCount)
                     print*, 'initCartesianGrid [warning]: H series - nuMaxArray increased to ', nuMaxArray, i
                  end if
                  nuCount = nuCount+3

               end do
               ! IP Thresholds
               do i = 1, nEdges
                  if (ionEdge(i) < nuMaxArray) then
                     nuArray(nuCount) = ionEdge(i)
                     nuArray(nuCount+1) = ionEdge(i)- 0.0003
                     nuArray(nuCount+2) = ionEdge(i)+ 0.0003
                     if ( nuArray(nuCount) < nuMinArray ) nuMinArray= nuArray(nuCount)- 0.0006
                     nuCount = nuCount+3
                  end if
               end do

               ! build the log energy mesh
               iCount = nbins-nuCount+1
               nuStepSize = (log10(nuMaxArray)-log10(nuMinArray))/(iCount-1)
               nuArray(nuCount) = nuMinArray
               do i = nuCount+1, nbins                  
                  nuArray(i) = 10.**(log10(nuArray(i-1))+nuStepSize)
               enddo

               ! now sort in ascending order
               call sortUp(nuArray)

            else if (lgDust .and. lgGas) then

               ! first count how many edge, thresholds etc
               if (nuMin<radio4p9GHz) then
                  nuArray(1)     = radio4p9GHz
                  nuCount = 2
               else
                  nuCount = 1
               end if

               nuMinArray = nuMin
               nuMaxArray = nuMax
               ! H series edges 
               do i = 1, nSeries 
                  nuArray(nuCount) = seriesEdge(i)
                  nuArray(nuCount+1) = seriesEdge(i)- 0.0003
                  nuArray(nuCount+2) = seriesEdge(i)+ 0.0003
                  if ( nuArray(nuCount) < nuMinArray ) nuMinArray= nuArray(nuCount)- 0.0006
                  if ( nuArray(nuCount) > nuMaxArray ) then
                     nuMaxArray= nuArray(nuCount)
                     print*, 'initCartesianGrid [warning]: H series - nuMaxArray increased to ', nuMaxArray, i
                  end if
                  nuCount = nuCount+3

               end do
               ! IP Thresholds
               do i = 1, nEdges
                  if (ionEdge(i) < nuMaxArray) then
                     nuArray(nuCount) = ionEdge(i)
                     nuArray(nuCount+1) = ionEdge(i)- 0.0003
                     nuArray(nuCount+2) = ionEdge(i)+ 0.0003
                     if ( nuArray(nuCount) < nuMinArray ) nuMinArray= nuArray(nuCount)- 0.0006
                     nuCount = nuCount+3
                  end if
               end do
               
               ! dust data points
               close(72)
               open (unit= 72,  action="read", file=trim(home)//"dustData/nuDustRyd.dat", status = "old", position = "rewind", &
                    & iostat = ios)
               if (ios /= 0) then
                  print*, "! initCartesianGrid: can't open dust nu grid file - ",trim(home),"nuDustGrid.dat"
                  stop
               end if

               do i = 1, 10000000
                  if (i<=nbins) then
                     read(unit=72,fmt=*,iostat=ios) nuArray(nuCount)
                     nuCount = nuCount+1
                     if (ios < 0) exit ! end of file reached
                     if (nuArray(i)>= seriesEdge(1)) exit
                  else
                     print*, "! initCartesianGrid: nbins is smaller that the number of &
                          & frequency points in nuDustGrid.dat file - enlarge nbins"
                     stop
                  end if
               end do

               close(72)

               ! build the log energy mesh
               iCount = nbins-nuCount+1
               nuStepSize = (log10(nuMaxArray)-log10(nuMinArray))/(iCount-1)
               nuArray(nuCount) = nuMinArray
               do i = nuCount+1, nbins                  
                  nuArray(i) = 10.**(log10(nuArray(i-1))+nuStepSize)
               enddo

               ! now sort in ascending order
               call sortUp(nuArray)               

            end if

            widFlx(1) = nuArray(2)-nuArray(1)
            do i = 2, nbins-1
               widFlx(i) = (nuArray(i+1)-nuArray(i-1))/2.
               print*, i, nuArray(i), widFlx(i)
            end do
            widFlx(nbins) =  nuArray(nbins)-nuArray(nbins-1)
            
             ! set the 4.9 GHz pointer
            if (nuArray(1) <= radio4p9GHz) then
               call locate(nuArray,radio4p9GHz,radio4p9GHzP)
               if (radio4p9GHz > (nuArray(radio4p9GHzP)+nuArray(radio4p9GHzP+1))/2.) &
                    & radio4p9GHzP =radio4p9GHzP+1
            end if
            
            if (widFlx(1)<=0. .or. widFlx(nbins)<=0. ) then
               print*, " ! initCartesianGrid: null or negative frequency bin [1, nbins]", widFlx(1), widFlx(nbins)
               stop
            end if

            lgFirst = .false.

            ! initialize the continuum gamma coeffs
            call initGammaCont()

         end if

        ! allocate active cells pointers
        allocate(grid%active(1:nx, 1:ny, 1:nz), stat = err)
        if (err /= 0) then
            print*, "Can't allocate grid memory, active"
            stop
        end if    

        ! allocate active cells volume
!        if (lgEcho) then
           allocate(grid%echoVol(1:nx, 1:ny, 1:nz), stat = err)
           if (err /= 0) then
              print*, "Can't allocate grid memory, echoVol"
              stop
           end if
!        endif

        ! allocate axes

        allocate(grid%xAxis(1:nx), stat = err)
        if (err /= 0) then
            print*, "Can't allocate grid memory, xAxis"
            stop
        end if 

        allocate(grid%yAxis(1:ny), stat = err)
        if (err /= 0) then
            print*, "Can't allocate grid memory, yAxis"
            stop
        end if 

        allocate(grid%zAxis(1:nz), stat = err)
        if (err /= 0) then
            print*, "Can't allocate grid memory, zAxis"
            stop
        end if 

        if (lgGas) then
           allocate(grid%abFileIndex(1:nx,1:ny,1:nz), stat = err)
           if (err /= 0) then
              print*, "Can't allocate grid memory, abFileIndex "
              stop
           end if
           
           grid%abFileIndex = 1
        end if

        grid%nx = nx
        grid%ny = ny
        grid%nz = nz

        ! initialize the arrays with zeros
        grid%active = 0
        grid%xAxis = 0.
        grid%yAxis = 0.
        grid%zAxis = 0.

        !new
         dTheta = Pi/totAngleBinsTheta
!         dTheta = Pi/totAngleBins
 
         do i = 1, nAngleBins 
            if (viewPointPhi(i)<0.) then 
               totAngleBinsPhi=1
               print*, '! initCartesianGrid : [warning] phi-dependance in viewing angle turned off'
               viewPointPhi=-1.
               exit
            end if
         end do

         ! new
         dPhi = twoPi/totAngleBinsPhi

!         allocate(viewPointP(0:totAngleBins), stat = err)
!         if (err /= 0) then
!            print*, "Can't allocate grid memory, viewAngleP "
!            stop
!         end if        


         allocate(viewPointPTheta(0:totAngleBinsTheta), stat = err)
         if (err /= 0) then
            print*, "! initCartesianGrid: Can't allocate grid memory, viewAnglePTheta "
            stop
         end if
         allocate(viewPointPPhi(0:totAngleBinsPhi), stat = err)
         if (err /= 0) then
            print*, "! initCartesianGrid: Can't allocate grid memory, viewAnglePPhi "
            stop
         end if
         
         viewPointPtheta = 0
         viewPointPphi = 0

!         ii = 1
!         do i = 1, nAngleBins         
!            viewPointP(int(viewPoint(i)/dTheta)+1) = ii
!            ii = ii+1
!         end do

         do i = 1, nAngleBins
            if (viewPointTheta(i) > Pi/2. .and. lgSymmetricXYZ) then
               print*, '! initCartesianGrid: the inclination theta required is not available for symmetricXYZ models (theta > Pi/2)'
               stop
            end if
            viewPointPtheta(int(viewPointTheta(i)/dTheta)+1) = i            
            viewPointPphi(int(viewPointPhi(i)/dPhi)+1) = i
         end do
         
!         print*, 'viewpointptheta'
!         print*, viewpointptheta
!         print*, ' '
!         print*, 'viewpointpphi'
!         print*, viewpointpphi

         print*, 'dTheta : ', dTheta
         print*, 'dPhi : ', dPhi

         print*, "out initCartesianGrid"

    end subroutine initCartesianGrid

    
    ! fillGrid sets up a 3d grid ith a bipolar geometry. this 
    ! procedure is only implemented for a cartesian geometry at present.
    ! 
    ! Version 1.0 and later: this subroutine now sets up any type of grids
    subroutine fillGrid(grid)
        implicit none

        include 'mpif.h'

        type(grid_type), dimension(:),intent(inout) :: grid

        ! local variables
        integer :: i, j, k, l, m, n, elem, ion, iG, jG    ! counters
        integer :: ios                                    ! I/O status
        integer :: err                                    ! memory allocation status
        integer :: ngridsloc
        integer, parameter :: max_points = 10000          ! safety limit

!        character(len=30)            :: in_file

        real :: radius     ! sqrt(x^2 + y^2 + z^2)

        ! pointer to the array of 2nd derivatives of the interpolating function
        ! calculated by spline for use into splint

        integer :: iCount, nuCount              ! counters
        integer :: g0,g1
        integer :: nEdges  
        integer :: nElec
        integer :: outshell
        integer :: totCells
        integer :: totCellsLoc=0

        integer, parameter :: maxLim = 10000
        integer, parameter :: nSeries = 17
        
        real                 :: nuStepSizeLoc

        print*, "in fillGrid"

        if (lg1D) lgSymmetricXYZ = .false.

        if (lgPlaneIonization) then
           allocate(planeIonDistribution(grid(1)%nx,grid(1)%nz), stat = err)
           if (err /= 0) then
              print*, "! fillGrid: can't allocate dl memory"
              stop
           end if
           planeIonDistribution = 0
        end if

        if (.not.lgDfile) then
                  
           if (lgSymmetricXYZ ) then

              ! create the grid axes, forcing it to have grid points at the
              ! centre
              do i = 1, grid(1)%nx
                 grid(1)%xAxis(i) = real(i-1)/real(grid(1)%nx-1) 
                 grid(1)%xAxis(i) = grid(1)%xAxis(i) * Rnx
              end do

              do i = 1, grid(1)%ny
                 grid(1)%yAxis(i) = real(i-1)/real(grid(1)%ny-1) 
                 grid(1)%yAxis(i) = grid(1)%yAxis(i) * Rny
              end do

              do i = 1, grid(1)%nz
                 grid(1)%zAxis(i) = real(i-1)/real(grid(1)%nz-1) 
                 grid(1)%zAxis(i) = grid(1)%zAxis(i) * Rnz
              end do

           else if (lg1D) then

              do i = 1, grid(1)%nx
                 grid(1)%xAxis(i) = real(i-1)/real(grid(1)%nx-1) 
                 grid(1)%xAxis(i) = grid(1)%xAxis(i) * Rnx
              end do
           
              grid(1)%yAxis = 0.
              grid(1)%zAxis = 0.
           
           else ! not lgSymmetricXYZ           
           
              if (mod(grid(1)%nx,2) == 0.) then
                 print*, "! fillGrid: the automatic grid option &
                      & requires odd integer nx if not symmetric"
                 stop
              end if

              if (mod(grid(1)%ny,2) == 0.) then
                 print*, "! fillGrid: the automatic grid option &
                      & requires odd integer ny if not symmetric"
                 stop
              end if

              if (mod(grid(1)%nz,2) == 0.) then
                 print*, "! fillGrid: the automatic grid option &
                      & requires odd integer nz if not symmetric"
                 stop
              end if
               
              ! create the grid axes, forcing it to have grid points at the
              ! centre
              do i = 1, grid(1)%nx
                 grid(1)%xAxis(i) = 2.*real(i-1)/real(grid(1)%nx-1) - 1.
                 grid(1)%xAxis(i) = grid(1)%xAxis(i) * Rnx
              end do
 
 
              do i = 1, grid(1)%ny
                 grid(1)%yAxis(i) = 2.*real(i-1)/real(grid(1)%ny-1) - 1.
                 grid(1)%yAxis(i) = grid(1)%yAxis(i) * Rny
              end do
           
              do i = 1, grid(1)%nz
                 grid(1)%zAxis(i) = 2.*real(i-1)/real(grid(1)%nz-1) - 1.
                grid(1)%zAxis(i) = grid(1)%zAxis(i) * Rnz
              end do

           end if

        end if

        if (lg1D) lgSymmetricXYZ = .true.
        

        
        if (lgNeutral) then
           
           ! set up density distribution and initial grid properties

           call setMotherGrid(grid(1))

           if (taskid == 0) then
              print*, "Mother Grid:"
              print*, "xAxis: ", grid(1)%xAxis
              print*, "yAxis: ", grid(1)%yAxis
              print*, "zAxis: ", grid(1)%zAxis
           end if

           ! set the subGrids
           if (nGrids>1) call setSubGrids(grid(1:nGrids))

           totCells = 0       
           totcellsloc=0
           if (emittingGrid>0) then
              nGridsloc = emittingGrid
           else
              nGridsloc = ngrids
           end if

           do i = 1, nGrids
              totCells = totCells+grid(i)%nCells
           end do
           do i = 1, nGridsloc
              totCellsloc = totCellsloc+grid(i)%nCells
           end do

           if (taskid==0) print*, '! fillGrid: total number of active cells over all grids: ', totCells
              
           if (lgDust) then
              if (lgSymmetricXYZ) totalDustMass=totalDustMass*8.
              print*, 'Total dust mass [1.e45 g]: ', totalDustMass
              print*, 'Total dust mass [Msol]: ', totalDustMass*5.028e11                
           end if           

           if (nPhotonsDiffuse > 0 .and. nPhotonsDiffuse < totCellsloc) then
              print*, '! fillGrid: total number of active cells is larger than the &
                   &total number of packets to be used in the simulation please enlarge nPhotonsDiffuse,&
                   &[nPhotonsDiffuse, totcellsloc]', nPhotonsDiffuse, totcellsloc
              stop
           end if

           if (totCellsloc<=0) then
              print*, '! setSubGrids: totCellsloc <=0'
              stop
           end if

           nPhotonsDiffuseLoc = nPhotonsDiffuse/totCellsloc
           if (taskid==0) print*, '! setSubGrids: number of diffuse packets per active cell: ', nPhotonsDiffuseLoc

        else

           print*, '! fillGrid: only neutral option availale &
                & in this version. plese re-run with neutral keyword'
           stop 
        end if
        
        ! print out Grid
        if (taskid == 0) then
           print*, "Mother Grid:"
           print*, "xAxis: ", grid(1)%xAxis
           print*, "yAxis: ", grid(1)%yAxis
           print*, "zAxis: ", grid(1)%zAxis
           do iG = 2, nGrids
              print*, "Grid : ", iG
              print*, "xAxis: ", grid(iG)%xAxis
              print*, "yAxis: ", grid(iG)%yAxis
              print*, "zAxis: ", grid(iG)%zAxis
           end do
        end if

        if (lgPlaneIonization) then
           R_out = 1.e10*( (grid(1)%xAxis(grid(1)%nx)/1.e10)*&
                & (grid(1)%yAxis(grid(1)%ny)/1.e10)*&
                & (grid(1)%zAxis(grid(1)%nz)/1.e10) )
           print*, "! fillGrid [warning]: Plane Ionization model -> R_out &
                & reset to maximum extension", R_out
        end if

        ! locate the origin of the axes
        call locate(grid(1)%xAxis, 0., iOrigin)
        call locate(grid(1)%yAxis, 0., jOrigin)
        call locate(grid(1)%zAxis, 0., kOrigin)

        if (taskid == 0) print*, 'Origin at mother grid cell:  ' , iOrigin, jOrigin, kOrigin

        ! allocate dl
        allocate(dl(nGrids), stat = err)
        if (err /= 0) then
           print*, "! fillGrid: can't allocate dl memory"
           stop
        end if
        dl = 0.        
        
        do iG = 1, nGrids
           ! find geometric corrections
           grid(iG)%geoCorrX = (grid(iG)%xAxis(grid(iG)%nx) - grid(iG)%xAxis(grid(iG)%nx-1))/2.
           if (.not. lg1D) then
              grid(iG)%geoCorrY = (grid(iG)%yAxis(grid(iG)%ny) - grid(iG)%yAxis(grid(iG)%ny-1))/2.
              grid(iG)%geoCorrZ = (grid(iG)%zAxis(grid(iG)%nz) - grid(iG)%zAxis(grid(iG)%nz-1))/2.
           else
              grid(iG)%geoCorrY = 0.
              grid(iG)%geoCorrZ = 0.
           end if
 
           if (taskid==0) print*, "Geometric grid corrections for grid ", &
                & iG, ' : ', grid(iG)%geoCorrX, grid(iG)%geoCorrY, grid(iG)%geoCorrZ
           
           ! find linear increment
           dl(iG) =  abs(grid(iG)%xAxis(2) - grid(iG)%xAxis(1))
           do i = 2, grid(iG)%nx-1
              dl(iG) = min(dl(iG), abs(grid(iG)%xAxis(i+1)-grid(iG)%xAxis(i)) )
           end do
           do i = 1, grid(iG)%ny-1
              dl(iG) = min(dl(iG), abs(grid(iG)%yAxis(i+1)-grid(iG)%yAxis(i)) )
           end do
           do i = 1, grid(iG)%nz-1
              dl(iG) = min(dl(iG), abs(grid(iG)%zAxis(i+1)-grid(iG)%zAxis(i)) )
           end do
           dl(iG) = dl(iG)/50.                                                                                                 
        end do

        if (nGrids>1) then

           do iG = 1, nGrids
              
              do i = 1, grid(iG)%nx
                 do j = 1, grid(iG)%ny
                    do k = 1, grid(iG)%nz
                 
                       do jG = 2, nGrids
                          if (jG /= iG) then
                             if (lgSymmetricXYZ) then                               

                                if ( ( &
                                     & (grid(iG)%xAxis(i) > grid(jG)%xAxis(1)) .or.&
                                     & (grid(iG)%xAxis(i) >= grid(jG)%xAxis(1) .and. grid(jG)%xAxis(1)==0.) & 
                                     & ) .and. &
                                     & grid(iG)%xAxis(i)<grid(jG)%xAxis(grid(jG)%nx) .and.&
                                     & ( & 
                                     & (grid(iG)%yAxis(j) > grid(jG)%yAxis(1)) .or. & 
                                     & (grid(iG)%yAxis(j) >= grid(jG)%yAxis(1) .and. grid(jG)%yAxis(1)==0. ) & 
                                     & ) .and. &
                                     & grid(iG)%yAxis(j)<grid(jG)%yAxis(grid(jG)%ny) .and. &
                                     & ( & 
                                     & (grid(iG)%zAxis(k) > grid(jG)%zAxis(1)) .or.& 
                                     & (grid(iG)%zAxis(k) >= grid(jG)%zAxis(1) .and. grid(jG)%zAxis(1)==0. ) & 
                                     & ) .and. &
                                     & grid(iG)%zAxis(k)<grid(jG)%zAxis(grid(jG)%nz) ) then
                                   grid(iG)%active(i,j,k) = -jG

                                end if



                          else
                                if ( (grid(iG)%xAxis(i) > grid(jG)%xAxis(1) .and. &
                                     &grid(iG)%xAxis(i)<grid(jG)%xAxis(grid(jG)%nx)) .and.&
                                     & (grid(iG)%yAxis(j)>grid(jG)%yAxis(1) .and. & 
                                     & grid(iG)%yAxis(j)<grid(jG)%yAxis(grid(jG)%ny)) .and. &
                                     & (grid(iG)%zAxis(k)>grid(jG)%zAxis(1) .and. & 
                                     & grid(iG)%zAxis(k)<grid(jG)%zAxis(grid(jG)%nz)) ) then
                                   grid(iG)%active(i,j,k) = -jG
                                end if

                             end if

                          end if
                       end do

                    end do
                 end do
              end do

           end do

        end if


        print*, "out fillGrid"
        
      end subroutine fillGrid


      subroutine setMotherGrid(grid)
        implicit none

        type(grid_type), intent(inout) :: grid      ! the grid

        ! local variables
        real                           :: denominator   ! denominator
        real                           :: dV           ! volume element
        real                           :: expFactor    ! exp factor in density law calculations
        real                           :: expTerm      ! exp term
        real                           :: factor       ! conputation factor
        real                           :: gasCell      ! mass of gas in current cell
        real                           :: H0in         ! estimated H0 at the inner radius for regionI
        real, pointer                  :: MdMg(:,:,:)  ! Md/Mg
        real                           :: MhMg         ! mass oh hydrogen over mass of gas
        real                           :: norm, scale  ! normalisation and scaling for meanField
        real                           :: radius       ! distance from the origin
        real                           :: random       ! random nmumber 
        real                           :: readReal     ! real number reader
        real                           :: surfIn       ! surface at inner radius [e36 cm^2]
        real                           :: totalMass    ! total ionized mass 
        real                           :: totalVolume  ! total active volume
        real                      :: echoVolume, vol   ! just echo volume

        real, dimension(nElements) :: aWeight
        real, parameter :: amu = 1.66053e-24 ! [g]

        real, pointer                  :: HdenTemp(:,:,:) ! temporary Hden
        real, pointer                  :: NdustTemp(:,:,:) ! temporary dust number density arra
        real, pointer                  :: dustAbunIndexTemp(:,:,:) ! temporary dust abundance index array
        real, pointer                  :: twoDscaleJTemp(:)
        


        integer                        :: icomp 
        integer                        :: i,j,k        ! counters
        integer                        :: index        ! general index
        integer                        :: iOrigin,&    ! indeces of the origin of the grid
             & jOrigin, kOrigin
        integer                        :: ios, err     ! I/O and allocation error status
        integer                        :: elem, ion    ! counters
        integer                        :: nspec, ai    ! counters
        integer                        :: nsp          ! pointer
        integer                        :: nu0P         ! 
        integer                        :: RinP         ! pointer to the inner radius intercept  

        integer                        :: yTop, xPmap  ! 2D indeces

                                                       ! with one of the axes
        character(len=40)              :: keyword      ! character string readers

        print*, 'in setMotherGrid'

        ! this is the mother grid
        grid%motherP = 0

        if (lgGas) then

           ! allocate space for HdenTemp
           allocate(HdenTemp(1:grid%nx, 1:grid%ny, 1:grid%nz), stat = err)
           if (err /= 0) then
              print*, "! setMotherGrid: can't allocate grid memory"
              stop
           end if
           HdenTemp = 0.        
        
           if (lgDfile) then
              open (unit= 77,  action="read", file=densityFile, status = "old", position = "rewind", &
                   & iostat = ios)
              if (ios /= 0) then
                 print*, "! setMotherGrid: can't open density file"
                 stop
              end if
              read(77,*)keyword
              if (keyword .ne. "#") backspace 77
           end if

        end if


        if (lg2D) then
           yTop = 1
        else
           yTop = grid%ny
        end if                


        grid%active = 1
        if (lgGas) then
           do i = 1, grid%nx
              do j = 1, yTop
                 do k = 1, grid%nz
                    
                    ! set density 
                    if (lgDlaw) then
                       
                       ! calculate radius
                       if (lg1D) then
                          radius = grid%xAxis(i)
                       else
                          radius = 1.e10*sqrt( (grid%xAxis(i)/1.e10)*(grid%xAxis(i)/1.e10) + &
                                 &                                        (grid%yAxis(j)/1.e10)*(grid%yAxis(j)/1.e10) + &
                                 &                                        (grid%zAxis(k)/1.e10)*(grid%zAxis(k)/1.e10) ) 
                       end if
                         
                         ! edit the following to use a different density law

                         ! N(r) = N0 * (r/r1)^n * exp[-0.5*(r/r1)^2]
                         ! Clegg, Harrington, Barlow and Walsh 1987, ApJ, 314, 551
                      
                         expFactor = exp(-0.5*(radius/densityLaw(1))*(radius/densityLaw(1)))
                        
                         
                         HdenTemp(i,j,k) = densityLaw(3)*expFactor*(radius/densityLaw(1))**densityLaw(2)
                              
                      
                      else if (lgHdenConstant) then                     

                         HdenTemp(i,j,k) = Hdensity
                         
                      else if (lgDfile) then
                         
                         if (.not.lgMultiChemistry) then

                            read(77, *) grid%xAxis(i), grid%yAxis(j), grid%zAxis(k), HdenTemp(i,j,k)
                         
                         else

                            read(77, *) grid%xAxis(i), grid%yAxis(j), grid%zAxis(k), &
                              & HdenTemp(i,j,k), grid%abFileIndex(i,j,k)
                            
                         end if
                      
                      else 
 
                         print*, "! setMotherGrid: no density distribution was set"
                         stop
                         
                      end if
                        
                      if (fillingFactor<1.) then

                         call random_number(random)
                         
                         if (random > fillingFactor) HdenTemp(i,j,k) = 0.

                      end if

                                           
                      if (grid%active(i,j,k) <= 0 ) HdenTemp(i,j,k) = 0.
                        
                   end do
                end do
             end do
          end if ! lgGas

          if (lgGas) close(77)

          ! set up dust data
          if (lgDust) then
              allocate(NdustTemp(1:grid%nx,1:grid%ny,1:grid%nz), stat = err)
              if (err /= 0) then
                 print*, "! setMotherGrid: can't allocate NdustTemp memory"
                 stop
              end if

              NdustTemp = 0.              

              if (lgMultiDustChemistry) then
                 allocate(dustAbunIndexTemp(1:grid%nx,1:grid%ny,1:grid%nz), stat = err)
                 if (err /= 0) then
                    print*, "! setMotherGrid: can't allocate dustAbunIndexTemp memory"
                    stop
                 end if
                 dustAbunIndexTemp = 0.              
              end if

              ! set grains mass density [g/cm^3]
!              allocate(dustComPoint(nDustComponents))
!              dustComPoint = 0
!              dustComPoint(1) = 1

!              nSpecies = 0
!              nSpeciesMax = 0
!print*, 'heer', ndustcomponents
!              do icomp = 1, nDustComponents
!print*, icomp
!                 close(13)
!                 open(file =   dustSpeciesFile(icomp), action="read",unit=13, &
!                      &position="rewind",status="old", iostat = ios)
!                 if (ios /= 0 ) then
!                    print*, "! setMotherGrid: can't open file ", dustSpeciesFile(icomp)
!                    stop
!                 end if
!                 read(13, *) nSpeciesPart(icomp)
!print*, nspeciespart(icomp)
!                 close(13)
!                 nSpecies = nSpecies+nSpeciesPart(icomp)
!print*, nspecies
!                 if (nSpeciesMax < nSpeciesPart(icomp)) nSpeciesMax = nSpeciesPart(icomp)
!print*, nspeciesmax
!              end do
!
!              allocate(rho(1:nSpecies), stat = err)
!              if (err /= 0) then
!                 print*, "! setMotherGrid: can't allocate rho memory"
!                 stop
!              end if
!              rho=0.
!              allocate(grainVn(1:nSpecies), stat = err)
!              if (err /= 0) then
!                 print*, "! setMotherGrid: can't allocate grainVn memory"
!                 stop
!              end if
!              grainVn=0.
!              allocate(MsurfAtom(1:nSpecies), stat = err)
!              if (err /= 0) then
!                 print*, "! setMotherGrid: can't allocate surfAtom memory"
!                 stop
!              end if
!              MsurfAtom=0
!              
!              do icomp = 1, nDustComponents
!                 if (icomp > 1) dustComPoint(icomp) = dustComPoint(icomp-1)+nSpeciesPart(icomp)
!                 close(13)
!                 open(file =   dustSpeciesFile(icomp), action="read",unit=13, &
!                      & position="rewind",status="old", iostat = ios)
!                 if (ios /= 0 ) then
!                    print*, "! setMotherGrid: can't open file ", dustSpeciesFile(icomp)
!                    stop
!                 end if
!                 read(13, *) nSpeciesPart(icomp)              
!
!                 do i = 1, nSpeciesPart(icomp)
!                    read(13,*) extFile
!                    close(14)
!                    open(file=extFile,unit=14,  action="read", position="rewind",status="old", iostat = ios)
!                    if (ios /= 0 ) then
!                       print*, "! setMotherGrid: can't open file ", extFile
!                       stop
!                    end if
!                    read(14,*) readChar
!                    read(14,*) readChar, readReal, rho(dustComPoint(icomp)+i-1), grainVn(dustComPoint(icomp)+i-1), &
!                         &MsurfAtom(dustComPoint(icomp)+i-1)
!                    close(14)
!                 end do
!                 close(13)
!              end do
!
              if (lgMdMg .or. lgMdMh) then

                 allocate(MdMg(1:grid%nx,1:grid%ny,1:grid%nz), stat = err)
                 if (err /= 0) then
                    print*, "! setMotherGrid: can't allocate MdMg memory"
                    stop
                 end if
                 MdMg = 0.
                 
                 if (lgDustConstant) then
                    MdMg = MdMgValue
                 else
                    close(20)
                    open(unit=20, file=MdMgFile,  action="read", position="rewind",status="old", iostat = ios)
                    if (ios /= 0 ) then
                       print*, "! setMotherGrid: can't open MdMgFile file ", MdMgFile
                       stop
                    end if
                    read(20,*)keyword
                    if (keyword .ne. "#") backspace 20
                    
                    if (lgMultiDustChemistry) then
                       do i = 1, grid%nx
                          do j = 1, yTop
                             do k = 1, grid%nz
                                read(20, *) index, index, index, MdMg(i,j,k), dustAbunIndexTemp(i,j,k)
                             end do
                          end do
                       end do
                    else                    
                       do i = 1, grid%nx
                          do j = 1, yTop
                             do k = 1, grid%nz
                                read(20, *) index, index, index, MdMg(i,j,k)
                             end do
                          end do
                       end do
                    end if
                    close(20)
                 end if

              else

                 ! Ndust was directly defined by  the user
                 if (lgDustConstant) then
                    NdustTemp = NdustValue
                 else
                    close(20)
                    open(unit=20, file=NdustFile,  action="read", position="rewind",status="old", iostat = ios)
                    if (ios /= 0 ) then
                       print*, "! setMotherGrid: can't open NdustFile file ", NdustFile
                       stop
                    end if
                    read(20,*)keyword
                    if (keyword .ne. "#") backspace 20

                    if (lgMultiDustChemistry) then
                       do i = 1, grid%nx
                          do j = 1, yTop
                             do k = 1, grid%nz
                                read(20, *) grid%xAxis(i), grid%yAxis(j), &
                                     & grid%zAxis(k), NdustTemp(i,j,k), dustAbunIndexTemp(i,j,k)
                             end do
                          end do
                       end do
                    else
                       print*,grid%nx,yTop,grid%nz
                       do i = 1, grid%nx
                          do j = 1, yTop
                             do k = 1, grid%nz
                                read(20, *) grid%xAxis(i), grid%yAxis(j), grid%zAxis(k), NdustTemp(i,j,k)
!                                write(6,'(3i4,4es11.3)')i,j,k,grid%xAxis(i), grid%yAxis(j), grid%zAxis(k), NdustTemp(i,j,k)
                             end do
                          end do
                       end do
                    end if
                    close(20)
                 end if

              end if              

           end if
           
           if (lg2D) grid%yAxis = grid%xAxis             


          ! set active cells pointers
          grid%nCells = 0
          do i = 1, grid%nx
             do j = 1, yTop
                do k = 1, grid%nz

                   ! calculate radius
                   if (lg1D) then
                      radius = grid%xAxis(i)
                   else
                      radius = 1.e10*sqrt( (grid%xAxis(i)/1.e10)*(grid%xAxis(i)/1.e10) + &
&                                        (grid%yAxis(j)/1.e10)*(grid%yAxis(j)/1.e10) + &
&                                        (grid%zAxis(k)/1.e10)*(grid%zAxis(k)/1.e10) ) 
                   end if
  

                   ! check if this grid point is  valid nebular point
                   if (.not. lgPlaneIonization) then

                      if (radius < R_in) then
                         grid%active(i,j,k) = 0
                      else if (radius > R_out .and. R_out>0.) then
                         grid%active(i,j,k) = 0
                      end if

                   end if

                   if (grid%active(i,j,k) <= 0 .and. lgDust ) then
                      NdustTemp(i,j,k) = 0.
                      if (lgMultiDustChemistry) dustAbunIndexTemp(i,j,k) = 0.
                   end if

                   if (grid%active(i,j,k) <= 0 .and. lgGas ) HdenTemp(i,j,k) = 0.

                   if (lgDust .and. lgGas) then
                      if (HdenTemp(i,j,k) > 0. .or. NdustTemp(i,j,k)>0.) then
                         grid%nCells = grid%nCells + 1
                         grid%active(i,j,k) = grid%nCells
                      else
                         grid%active(i,j,k) = 0
                         HdenTemp(i,j,k) = 0.
                         NdustTemp(i,j,k) = 0.
                         if (lgMultiDustChemistry) dustAbunIndexTemp(i,j,k) = 0.
                      end if
                   else if ( lgDust .and. (.not.lgGas) ) then
                      if (NdustTemp(i,j,k)>0.) then
                         grid%nCells = grid%nCells + 1
                         grid%active(i,j,k) = grid%nCells
                      else
                         grid%active(i,j,k) = 0
                         NdustTemp(i,j,k) = 0.
                         if (lgMultiDustChemistry) dustAbunIndexTemp(i,j,k) = 0.
                      end if
                   else if ( (.not.lgDust) .and. lgGas) then 
                      if (HdenTemp(i,j,k) > 0.) then
                         grid%nCells = grid%nCells + 1
                         grid%active(i,j,k) = grid%nCells
                      else
                         grid%active(i,j,k) = 0
                         HdenTemp(i,j,k) = 0.
                      end if
                   else

                      print*, '! setMotherGrid: no gas and no dust? The grid is empty.'
                      stop
                   end if

                end do
             end do
          end do


          allocate(TwoDscaleJtemp(grid%nCells))
          TwoDscaleJtemp = 1.


          if (lg2D) then
             do i = 1, grid%nx
                do j = 2, grid%ny
                   do k = 1, grid%nz
                      radius = 1.e10*sqrt( (grid%xAxis(i)/1.e10)*&
                           &(grid%xAxis(i)/1.e10) + &
                           &(grid%yAxis(j)/1.e10)*(grid%yAxis(j)/1.e10) ) 

                      call locate(grid%xAxis, radius, xPmap)
                      if (xPmap < grid%nx) then
                         if (radius >= (grid%xAxis(xPmap)+grid%xAxis(xPmap+1))/2.) &
                              & xPmap = xPmap+1
                      end if
                      grid%active(i,j,k) = grid%active(xPmap, 1, k)
                      
                      if (grid%active(xPmap,1,k)>0) &
                           & TwoDScaleJtemp(grid%active(xPmap,1,k)) = &
                           & TwoDScaleJtemp(grid%active(xPmap,1,k))+1.
                      
                   end do
                end do
             end do

             grid%nCells = 0
             do i = 1,  grid%nx
                do k = 1,  grid%nz
                   if (grid%active(i,1,k) > 0) grid%nCells = grid%nCells +1                   

                end do
             end do

             allocate(TwoDscaleJ(grid%nCells))
             do i = 1, grid%nCells
                TwoDscaleJ(i) = TwoDscaleJtemp(i)
             end do
             deallocate(TwoDscaleJtemp)

          end if



          print*, '! setMotherGrid: active cells :', grid%nCells



          ! allocate grid arrays
          if (lgGas .and. lgDust .and. lgPhotoelectric) then
             allocate(grid%JPEots(1:grid%nCells, 1:nbins), stat = err)
             if (err /= 0) then
                print*, "! setMotherGrid: can't allocate JPEots memory"
                stop
             end if
             grid%JPEots=0
          end if

          if (lgGas) then

             allocate(grid%Hden(0:grid%nCells), stat = err)
             if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory"
                stop
             end if

             allocate(grid%recPDF(0:grid%nCells, 1:nbins), stat = err)
             if (err /= 0) then
                print*, "Can't allocate grid memory, 8"
                stop
             end if

             allocate(grid%totalLines(0:grid%nCells), stat = err)
             if (err /= 0) then
                print*, "Can't allocate grid memory, 10"
                stop
             end if

             allocate(grid%ionDen(0:grid%nCells, 1:nElementsUsed, 1:nstages), stat = err)
             if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory,ionDen"  
                stop
             end if
             allocate(ionDenUsed(1:nElementsUsed, 1:nstages), stat = err)
             if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory,ionDen"  
                stop
             end if
             allocate(grid%Ne(0:grid%nCells), stat = err)
             if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory"
                stop
             end if
             allocate(grid%Te(0:grid%nCells), stat = err)
             if (err /= 0) then
                print*, "! setMotherGrid:can't allocate grid memory"
                stop
             end if

             grid%Hden = 0.
             grid%Ne = 0.
             grid%Te = 0.        
             grid%ionDen = 0.
             grid%recPDF = 0.
             grid%totalLines = 0.
             
          end if

          if (Ldiffuse>0.) then
             allocate(grid%LdiffuseLoc(0:grid%nCells), stat = err)
             if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory : LdiffuseLoc "
                stop
             end if
             grid%LdiffuseLoc=0.
          end if


          allocate(grid%opacity(0:grid%nCells, 1:nbins), stat = err)
          if (err /= 0) then
             print*, "! setMotherGrid: can't allocate grid memory : opacity "
             stop
          end if
          allocate(grid%Jste(0:grid%nCells, 1:nbins), stat = err)
          if (err /= 0) then
             print*, "! setMotherGrid: can't allocate grid memory : Jste"
             stop
          end if

          allocate(grid%escapedPackets(0:grid%nCells, 0:nbins,0:nAngleBins), stat = err)
          if (err /= 0) then
             print*, "! setMotherGrid: can't allocate grid memory : Jste"
             stop
          end if

          if (lgEquivalentTau) then

             allocate(SEDnoExt(1:nbins), stat = err)
             if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory : SEDnoExt"
                stop
             end if
             allocate(equivalentTau(1:nbins), stat = err)
             if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory : equivalentTau"
                stop
             end if

          end if
          if (lgDust) then

             allocate(grid%Ndust(0:grid%nCells), stat = err)
             if (err /= 0) then
                print*, "! grid: can't allocate Ndust memory"
                stop
             end if
             grid%Ndust=0.
             ! be 2.02.44
             allocate(grid%dustAbunIndex(0:grid%nCells), stat = err)
             if (err /= 0) then
                print*, "! grid: can't allocate dustAbunIndex memory"
                stop
             end if
             grid%dustAbunIndex=0.
             ! be 2.02.44 end
             
             if (.not.lgGas) then
                allocate(grid%dustPDF(0:grid%nCells, 1:nbins), stat = err)
                if (err /= 0) then
                   print*, "! grid: can't allocate dustPDF memory"
                   stop
                end if
                grid%dustPDF = 0.
             end if

          end if
!          print* ,'a'

!BS10          if (lgDebug) then
             allocate(grid%Jdif(0:grid%nCells, 1:nbins), stat = err)
             if (err /= 0) then
                print*, "! grid: can't allocate grid memory"
                stop
             end if

             grid%Jdif = 0. 

             ! allocate pointers depending on nLines
             if (nLines > 0) then

                allocate(grid%linePackets(0:grid%nCells, 1:nLines), stat = err)
                if (err /= 0) then
                   print*, "! initCartesianGrid: can't allocate grid%linePackets memory"
                   stop
                end if

                allocate(grid%linePDF(0:grid%nCells, 1:nLines), stat = err)
                if (err /= 0) then
                   print*, "! initCartesianGrid: can't allocate grid%linePDF memory"
                   stop
                end if
                
                grid%linePackets = 0.
                grid%linePDF     = 0.
                
             end if
             
!BS10          end if

          allocate(grid%lgConverged(0:grid%nCells), stat = err)
          if (err /= 0) then
             print*, "Can't allocate memory to lgConverged array"
             stop
          end if

          allocate(grid%lgBlack(0:grid%nCells), stat = err)
          if (err /= 0) then
             print*, "Can't allocate memory to lgBlack array"
             stop
          end if


          if (lgNeInput) then
             allocate(grid%NeInput(0:grid%nCells), stat = err)
             if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory, grid%NeInput"
                stop
             end if
             grid%NeInput = 0.
          end if


          grid%opacity = 0.        
          grid%Jste = 0.        
          grid%lgConverged = 0
          grid%lgBlack = 0           

          do i = 1, grid%nx
             do j = 1, yTop
                do k = 1, grid%nz
                   if (grid%active(i,j,k)>0) then

                      if (lgGas) then
                         grid%Hden(grid%active(i,j,k)) = HdenTemp(i,j,k)
                         grid%Te(grid%active(i,j,k)) = TeStart
                      end if
                      if (lgDust) then
                         grid%Ndust(grid%active(i,j,k)) = NdustTemp(i,j,k)
                         if (lgMultiDustChemistry) &
                              & grid%dustAbunIndex(grid%active(i,j,k)) = &
                              &dustAbunIndexTemp(i,j,k)
                      end if
                   end if
                end do
             end do
          end do
!          print*, 'b'
          if (lgNeInput) then 
             grid%NeInput = grid%Hden
             ! 1.11 is just a starting guess for the ionization 
             ! factor
             grid%Hden = grid%Hden/1.11
          end if

          if (lgGas) then

             H0in  = 1.e-5

             do i = 1, grid%nx
                do j = 1, yTop
                   do k = 1, grid%nz
                      
                      if (grid%active(i,j,k)>0 ) then

                         ! calculate ionDen for H0
                         if (lgElementOn(1)) then
                            grid%ionDen(grid%active(i,j,k),elementXref(1),1) = H0in                                              
                            grid%ionDen(grid%active(i,j,k),elementXref(1),2) = &
                                 & 1.-grid%ionDen(grid%active(i,j,k),elementXref(1),1)
                         end if
                         if (lgElementOn(2)) then
                            grid%ionDen(grid%active(i,j,k),elementXref(2),1) = &
                                 & grid%ionDen(grid%active(i,j,k),(1),1)
                            grid%ionDen(grid%active(i,j,k),elementXref(2),2) = &
                                 & (1.-grid%ionDen(grid%active(i,j,k),elementXref(2),1))
                            grid%ionDen(grid%active(i,j,k),elementXref(2),3) = 0.
                         end if

                         ! initialize Ne
                         grid%Ne(grid%active(i,j,k)) =  grid%Hden(grid%active(i,j,k))
                         
                         ! initialize all heavy ions (careful that the sum over all ionisation 
                         ! stages of a given atom doesn't exceed 1.)
                         do elem = 3, nElements
                            do ion = 1, min(elem+1,nstages)
                               if (lgElementOn(elem)) then
                                  if (ion == 1) then
                                     grid%ionDen(grid%active(i,j,k),elementXref(elem),ion) &
                                          & = grid%ionDen(grid%active(i,j,k),1,1)
                                  else
                                     grid%ionDen(grid%active(i,j,k),elementXref(elem),ion) = 0.
                                  end if
                               end if
                            end do
                         end do
                      end if     ! active condition
                   
                   end do
                end do
             end do
             
             ! deallocate temp array
             if(associated(HdenTemp)) deallocate(HdenTemp)

          end if ! lgGas

           if (lgDust) then
              if(associated(NdustTemp)) deallocate(NdustTemp)
              if(lgMultiDustChemistry .and. associated(dustAbunIndexTemp)) deallocate(dustAbunIndexTemp)
           end if

           ! set up atomic weight array
           aWeight = (/1.0080, 4.0026, 6.941, 9.0122, 10.811, 12.0111, 14.0067, 15.9994, &
                & 18.9984, 20.179, 22.9898, 24.305, 26.9815, 28.086, 30.9738, 32.06, 35.453, &
                & 39.948, 39.102, 40.08, 44.956, 47.90, 50.9414, 51.996, 54.9380, 55.847, 58.9332, &
                & 58.71, 63.546, 65.37 /)

           totalDustMass = 0.
           totalMass = 0.
           totalVolume = 0.
           echoVolume = 0.

           do i = 1, grid%nx
              do j = 1, yTop
                 do k = 1, grid%nz
                    grid%echoVol(i,j,k)=0. ! initialize
                    if (grid%active(i,j,k)>0) then

                       dV = getVolume(grid,i,j,k)

                       if (lgGas) then
                          gasCell = 0.
                          do elem = 1, nElements
                             gasCell = gasCell + grid%elemAbun(grid%abFileIndex(i,j,k),elem)*&
                                  & aWeight(elem)*amu
                             totalMass = totalMass + &
                                  & grid%Hden(grid%active(i,j,k))*dV*grid%elemAbun(grid%abFileIndex(i,j,k),elem)*&
                                  & aWeight(elem)*amu
                          end do
                       end if

                       totalVolume = totalVolume + dV
!
! echoes only
!
                       if (lgEcho) then
                          grid%echoVol(i,j,k)=vEcho(grid,i,j,k,echot1,echot2,vol)
                          echoVolume = echoVolume + vol 
                       endif

                       if (lgDust .and. (lgMdMg.or.lgMdMh) ) then

                          
                          if (lgMdMh) then
                             MhMg=0.
                             do elem = 1, nElements
                                ! transform to MdMg
                                MhMg = MhMg+grid%elemAbun(grid%abFileIndex(i,j,k),elem)*&
                                     & aWeight(elem)
                             end do
                             MhMg = 1./MhMg                             
                             MdMg(i,j,k) = MdMg(i,j,k)*MhMg
                          end if

                          if (.not.lgGas) then
                             print*, '! setMotherGrid: Mass to dust ratio (MdMg) cannot be used in a pure dust (noGas)&
                                  & simulation. Ndust must be used instead.'
                             stop
                          end if

                          grid%Ndust(grid%active(i,j,k)) = gasCell*MdMg(i,j,k)*grid%Hden(grid%active(i,j,k))

                          denominator = 0.

                          if (lgMultiDustChemistry) then
                             nsp = grid%dustAbunIndex(grid%active(i,j,k))
                          else
                             nsp = 1
                          end if
!print*, nsp, 'here!'
                          do nspec = 1, nSpeciesPart(nsp)
                             do ai = 1, nSizes
                                denominator = denominator + &
                                     & (1.3333*Pi*( (grainRadius(ai)*1.e-4)**3)*&
                                     & rho(dustComPoint(nsp)+nspec-1)*grainWeight(ai)*&
                                     & grainAbun(nsp, nspec))
                             end do
                          end do
                          grid%Ndust(grid%active(i,j,k)) = grid%Ndust(grid%active(i,j,k))/&
                               & (denominator)
                       end if

                       ! calculate total dust mass
                       if (lgDust) then
                          if (lgMultiDustChemistry) then
                             nsp = grid%dustAbunIndex(grid%active(i,j,k))
                          else
                             nsp = 1
                          end if

                          do ai = 1, nsizes
                             do nspec = 1, nspeciesPart(nsp)
                                totalDustMass = totalDustMass + &
                                     &(1.3333*Pi*((grainRadius(ai)*1.e-4)**3)*&
                                     & rho(dustComPoint(nsp)-1+nspec)*grainWeight(ai)*&
                                     & grainAbun(nsp,nspec))*grid%Ndust(grid%active(i,j,k))*dV

                             end do
                          end do

                       end if

                    end if
                 end do
              end do
           end do

           if (lgDust .and. lginputDustMass) then
           do i = 1, grid%nx
              do j = 1, yTop
                 do k = 1, grid%nz
             grid%Ndust(grid%active(i,j,k)) = grid%Ndust(grid%active(i,j,k)) * inputDustMass / totalDustMass
                 end do
              end do
           end do
           totalDustMass = inputDustMass
           end if


           if(associated(MdMg)) deallocate(MdMg)

           if (taskid == 0) then

              print*, 'Mothergrid :'
              if (lgGas) then
                 print*, 'Total gas mass of ionized region by mass [1.e45 g]: ', totalMass
              end if
              print*, 'Total volume of the active region [e45 cm^3]: ', totalVolume
              if (lgEcho) then 
                 print*, 'Total volume of the echo [e45 cm^3]: ', echoVolume*846732407.," or ",echoVolume," ly^3"
                 open(unit=99, status='unknown', position='rewind', file='output/echo.out', action="write",iostat=ios)
                 write(99,*)'Total volume of the active region [e45 cm^3]: ', totalVolume
                 write(99,*)'Total volume of the echo [e45 cm^3]: ', echoVolume*846732407.," or ",echoVolume," ly^3"
                 if (echoVolume .eq. 0.) then
                    print*,'No dust in echo region. Stopping'
                    write(99,*)'No dust in echo region. Stopping'
                    stop
                 endif
                 close(99)
              endif

           end if
              
           ! if we are using a plane parallel ionization then we must find the luminosity 
           ! of the ionizing plane from the input meanField
           if (lgPlaneIonization) then

              print*, 'Flux above ', nu0, ' is ', meanField

              if (nu0 > 0.) then
                 call locate(nuArray, nu0, nu0P) 
                 if (nu0P >= nbins .or. nu0P <1) then
                    print*, "! setMotherGrid: insanity in nu0P", nu0P, nuArray(i), nuArray(nbins)
                    stop
                 end if
                 norm = 0.
                 do i = nu0P, nbins
                    norm = norm+inSpectrumPhot(i)*widflx(i)
                 end do
                 scale  = meanField/norm
                 norm = 0.
                 do i = 1, nbins
                    norm = norm+inSpectrumErg(i)*widFlx(i)
                 end do
                 meanField = norm*scale
              end if

              print*, 'Flux bolometric is ', meanField              

              Lstar(1) = (meanField/1.e36)*grid%xAxis(grid%nx)*grid%zAxis(grid%nz)
              deltaE(1) = Lstar(1)/nPhotons(1)
           end if
                 
           if (taskid == 0) then
              print*, 'Total ionizing flux :', Lstar(1)
              print*, 'deltaE :', deltaE(1)
           end if

           print*, 'out setMotherGrid'

         end subroutine setMotherGrid

         subroutine setSubGrids(grid)
           implicit none
           
           type(grid_type), dimension(:),intent(inout) :: grid        ! the grid

           

           ! local variables
           real, pointer                  :: a(:)         ! grain radii
           real                           :: x,y,z        ! x,y,z 
           real                           :: delta        ! displacement for realigning subgrids to mothergrid
           real                           :: denominator   ! denominator
           real                           :: denfac       ! enhancement factor
           real                           :: dV           ! volume element
           real                           :: expFactor    ! exp factor in density law calculations
           real                           :: expTerm      ! exp term
           real                           :: factor       ! conputation factor
           real                           :: gasCell      ! mass of gas in current cell
           real                           :: H0in         ! estimated H0 at the inner radius for regionI
           real                           :: MhMg         ! hdrogen to gas mass ratio
           real                           :: normWeight   ! normalization const for grain size distribution
           real                           :: radius       ! distance from the origin
           real                           :: random       ! random nmumber 
           real                           :: surfIn       ! surface at inner radius [e36 cm^2]
           real                           :: totalMass    ! total ionized mass 
           real                           :: totalVolume  ! total active volume
           real                      :: echoVolume, vol   ! just echo volume
           real, pointer                  :: weight(:)    ! weight of grain
           
           
           real, pointer                  :: HdenTemp(:,:,:) ! temporary Hden
           real, pointer                  :: NdustTemp(:,:,:) ! temporary dust number density array
           real, pointer                  :: dustAbunIndexTemp(:,:,:) ! temporary dust number density array
           
           integer                        :: edgeP        ! subgrid edge pointer on mothergrid
           integer                        :: i,j,k,iG,ai  ! counters
           integer                        :: ix,iy,iz     ! counters          
           integer                        :: index        ! general index
           integer                        :: iOrigin,&    ! indeces of the origin of the grid jOrigin,&    
                & kOrigin
           integer                        :: ios, err     ! I/O and allocation error status
           integer                        :: elem, ion    ! counters
           integer                        :: nspec, nsize ! counters
           integer                        :: nsp          ! pointer
           integer                        :: NgrainSize   ! number of grain sizes
           integer                        :: RinP         ! pointer to the inner radius intercept  
                                                       ! with one of the axes
           character(len=50)              :: dFileRead    ! subgrid dFile  reader

           ! open grid info file
           close(71)
           open (unit= 71, file=gridList,  action="read", status = "old", position = "rewind", &
                & iostat = ios)
           if (ios /= 0) then
              print*, "! setsubGrids: can't open grid list file", gridList
              stop
           end if

           do iG = 2, nGrids

              ! set the elemental abundances 
              grid(iG)%elemAbun = grid(1)%elemAbun


              read(71, *) grid(iG)%motherP ,  grid(iG)%nx,  grid(iG)%ny,  grid(iG)%nz, dFileRead, denfac

              read(71, *) grid(iG)%xAxis(1), grid(iG)%xAxis(grid(iG)%nx), &
                   & grid(iG)%yAxis(1), grid(iG)%yAxis(grid(iG)%ny), &
                   &grid(iG)%zAxis(1), grid(iG)%zAxis(grid(iG)%nz)

              ! realign boundaries with mother grid cell half-way boundaries
              if (.false.) then

              call locate(grid(grid(iG)%motherP)%xAxis, grid(iG)%xAxis(1), edgeP)

              if (edgeP>=grid(grid(iG)%motherP)%nx) then
                 print*, 'setSubGrids : edgeP for lower x boundary equals to nx of the mother grid'
                 stop
              end if
              if (edgeP == 1 .and. grid(iG)%xAxis(1) > 1.001*grid(grid(iG)%motherP)%xAxis(1) ) then

                 delta = (grid(grid(iG)%motherP)%xAxis(edgeP)+&
                      & grid(grid(iG)%motherP)%xAxis(edgeP+1))/2. -&
                      & grid(iG)%xAxis(1)

                 if (taskid==0) then 
                    print*, 'setSubGrids : Lower x-edge displaced by ', delta, ' for subGrid ', iG
                 end if
                 
                 grid(iG)%xAxis(1) = &
                      & (grid(grid(iG)%motherP)%xAxis(edgeP)+grid(grid(iG)%motherP)%xAxis(edgeP+1))/2.

              else
                 
                 grid(iG)%xAxis(1) = &
                      & (grid(grid(iG)%motherP)%xAxis(edgeP))

              end if

              call locate(grid(grid(iG)%motherP)%xAxis, grid(iG)%xAxis(grid(iG)%nx), edgeP)

              if (edgeP==grid(grid(iG)%motherP)%nx) then

                 delta = grid(grid(iG)%motherP)%xAxis(edgeP)-&
                   & grid(iG)%xAxis(grid(iG)%nx)
              else

                 delta = (grid(grid(iG)%motherP)%xAxis(edgeP)+grid(grid(iG)%motherP)%xAxis(edgeP+1))/2. -&
                      & grid(iG)%xAxis(grid(iG)%nx)
              end if

              if (taskid==0) then 
                 print*, 'setSubGrids : Upper x-edge displaced by ', delta, ' for subGrid ', iG
              end if

              if (edgeP==grid(grid(iG)%motherP)%nx) then
                 grid(iG)%xAxis(grid(iG)%nx) = grid(grid(iG)%motherP)%xAxis(edgeP)
              else
                 grid(iG)%xAxis(grid(iG)%nx) = (grid(grid(iG)%motherP)%xAxis(edgeP)+&
                      &grid(grid(iG)%motherP)%xAxis(edgeP+1))/2.
              end if

              if (taskid==0) then 
                 print*, 'setSubGrids : x1,x2,iG: ',  grid(iG)%xAxis(1),  grid(iG)%xAxis(grid(iG)%nx), iG
              end if

              call locate(grid(grid(iG)%motherP)%yAxis, grid(iG)%yAxis(1), edgeP)

              if (edgeP>=grid(grid(iG)%motherP)%ny) then
                 print*, 'setSubGrids : edgeP for lower y boundary equals to ny of the mother grid'
                 stop
              end if

              if (edgeP == 1 .and. grid(iG)%yAxis(1) > 1.001*grid(grid(iG)%motherP)%yAxis(1) ) then

                 delta = (grid(grid(iG)%motherP)%yAxis(edgeP)+&
                      & grid(grid(iG)%motherP)%yAxis(edgeP+1))/2. -&
                      & grid(iG)%yAxis(1)
 
                 if (taskid==0) then 
                    print*, 'setSubGrids : Lower y-edge displaced by ', delta, ' for subGrid ', iG
                 end if
                 
                 grid(iG)%yAxis(1) =  (grid(grid(iG)%motherP)%yAxis(edgeP)+&
                         & grid(grid(iG)%motherP)%yAxis(edgeP+1))/2.


              else
                 
                 grid(iG)%yAxis(1) = &
                      & (grid(grid(iG)%motherP)%yAxis(edgeP))

              end if

              call locate(grid(grid(iG)%motherP)%yAxis, grid(iG)%yAxis(grid(iG)%ny), edgeP)

              if (edgeP==grid(grid(iG)%motherP)%ny) then

                 delta = grid(grid(iG)%motherP)%yAxis(edgeP)-&
                   & grid(iG)%yAxis(grid(iG)%ny)
              else
                 delta = (grid(grid(iG)%motherP)%zAxis(edgeP)+& 
                      & grid(grid(iG)%motherP)%zAxis(edgeP+1))/2. -&
                      & grid(iG)%zAxis(grid(iG)%nz)
              end if

              if (taskid==0) then 
                 print*, 'setSubGrids : Upper y-edge displaced by ', delta, ' for subGrid ', iG
              end if

              if (edgeP==grid(grid(iG)%motherP)%ny) then
                 grid(iG)%yAxis(grid(iG)%ny) = grid(grid(iG)%motherP)%yAxis(edgeP)
              else
                 grid(iG)%yAxis(grid(iG)%ny) =  (grid(grid(iG)%motherP)%yAxis(edgeP)+&
                      &grid(grid(iG)%motherP)%yAxis(edgeP+1))/2.
              end if

              if (taskid==0) then 
                 print*, 'setSubGrids : y1,y2,iG: ',  grid(iG)%yAxis(1),  grid(iG)%yAxis(grid(iG)%ny)
              end if


              call locate(grid(grid(iG)%motherP)%zAxis, grid(iG)%zAxis(1), edgeP)

              if (edgeP>=grid(grid(iG)%motherP)%nz) then
                 print*, 'setSubGrids : edgeP for lower z boundary equals to nz of the mother grid'
                 print*, grid(grid(iG)%motherP)%zAxis, grid(iG)%zAxis(1), edgeP
                 stop
              end if

              if (edgeP == 1 .and. grid(iG)%zAxis(1) > 1.001*grid(grid(iG)%motherP)%zAxis(1) ) then

                 delta = (grid(grid(iG)%motherP)%zAxis(edgeP)+&
                      & grid(grid(iG)%motherP)%zAxis(edgeP+1))/2. -&
                      & grid(iG)%zAxis(1)

                 if (taskid==0) then 
                    print*, 'setSubGrids : Lower z-edge displaced by ', delta, ' for subGrid ', iG
                 end if


                 grid(iG)%zAxis(1) =  (grid(grid(iG)%motherP)%zAxis(edgeP)+& 
                      & grid(grid(iG)%motherP)%zAxis(edgeP+1))/2.

              else
                 
                 grid(iG)%zAxis(1) = &
                      & (grid(grid(iG)%motherP)%zAxis(edgeP))

              end if

              call locate(grid(grid(iG)%motherP)%zAxis, grid(iG)%zAxis(grid(iG)%nz), edgeP)

              if (edgeP==grid(grid(iG)%motherP)%nz) then
                 delta = grid(grid(iG)%motherP)%zAxis(edgeP) -&
                      & grid(iG)%zAxis(grid(iG)%nz)
              else
                 delta = (grid(grid(iG)%motherP)%zAxis(edgeP)+& 
                      & grid(grid(iG)%motherP)%zAxis(edgeP+1))/2. -&
                      & grid(iG)%zAxis(grid(iG)%nz)
              end if

              if (taskid==0) then 
                 print*, 'setSubGrids : Upper z-edge displaced by ', delta, ' for subGrid ', iG
              end if

              if (edgeP==grid(grid(iG)%motherP)%nz) then
                 grid(iG)%zAxis(grid(iG)%nz) = grid(grid(iG)%motherP)%zAxis(edgeP) 
              else
                 grid(iG)%zAxis(grid(iG)%nz) = (grid(grid(iG)%motherP)%zAxis(edgeP)+& 
                      & grid(grid(iG)%motherP)%zAxis(edgeP+1))/2.
              end if

              if (taskid==0) then 
                 print*, 'setSubGrids : z1,z2,iG: ',  grid(iG)%zAxis(1),  grid(iG)%zAxis(grid(iG)%nz)
              end if

           end if

              close(72)
              open (unit= 72, file=dFileRead,  action="read", status = "old", position = "rewind", &
                   & iostat = ios)
              if (ios /= 0) then
                 print*, "! setsubGrids: can't open subgrid density file", dFileRead
                 stop
              end if

              if (lgDust) then
                 allocate(NdustTemp(1:grid(iG)%nx,1:grid(iG)%ny,1:grid(iG)%nz), stat = err)
                 if (err /= 0) then
                    print*, "! setSubGrids: can't allocate NdustTemp memory for subGrid", iG
                    stop
                 end if
                 NdustTemp = 0.
                 if (lgMultiDustChemistry) then
                    allocate(dustAbunIndexTemp(1:grid(iG)%nx,1:grid(iG)%ny,1:grid(iG)%nz), stat = err)
                    if (err /= 0) then
                       print*, "! setSubGrids: can't allocate dustAbunIndexTemp memory for subGrid", iG
                       stop
                    end if
                    dustAbunIndexTemp = 0.
                 end if
              end if

              if (lgEcho) then ! allocate active cells volume
                 allocate(grid(iG)%echoVol(1:grid(iG)%nx, 1:grid(iG)%ny, 1:grid(iG)%nz), stat = err)
                 if (err /= 0) then
                    print*, "Can't allocate grid memory, echoVol"
                    stop
                 end if
              endif

              ! allocate space for HdenTemp 
              if (lgGas) then
                 allocate(HdenTemp(1:grid(iG)%nx, 1:grid(iG)%ny, 1:grid(iG)%nz), stat = err)
                 if (err /= 0) then
                    print*, "! setSubGrid: can't allocate sub grid memory: HdenTemp", iG
                    stop
                 end if
                 HdenTemp = 0.        
              end if

              grid(iG)%active = 1                 
              grid(iG)%nCells = 0
              do ix = 1, grid(iG)%nx
                 do iy = 1, grid(iG)%ny
                    do iz = 1, grid(iG)%nz
                       
                       if (lg1D) then
                          print*, "! setSubGrids: No 1D option with multiple grids!!"
                          stop
                       end if
                       
                       if (.not.lgMultiChemistry) then

                          if (lgMultiDustChemistry) then

                             if (lgGas .and. lgDust) then
                                read(72, *) x,y,z, HdenTemp(ix,iy,iz), &
                                     & NdustTemp(ix,iy,iz), dustAbunIndexTemp(ix,iy,iz)
!print*, x,y,z, HdenTemp(ix,iy,iz), NdustTemp(ix,iy,iz), dustAbunIndexTemp(ix,iy,iz)

                             else if (lgGas .and. .not. lgDust) then
                                read(72, *) x,y,z, HdenTemp(ix,iy,iz)
                             else if (.not.lgGas .and. lgDust) then
                                read(72, *) x,y,z, NdustTemp(ix,iy,iz), dustAbunIndexTemp(ix,iy,iz)
                             else
                                print*, "! setSubGrids: No gas or dust with multiple grids!!"
                                stop
                             end if

                          else

                             if (lgGas .and. lgDust) then
                                read(72, *) x,y,z, HdenTemp(ix,iy,iz), NdustTemp(ix,iy,iz)
                             else if (lgGas .and. .not. lgDust) then
                                read(72, *) x,y,z, HdenTemp(ix,iy,iz)
                             else if (.not.lgGas .and. lgDust) then
                                read(72, *) x,y,z, NdustTemp(ix,iy,iz)
                             else
                                print*, "! setSubGrids: No gas or dust with multiple grids!!"
                                stop
                             end if

                          end if
                       else
                          
                          if (lgMultiDustChemistry) then
                             
                             if (lgGas .and. lgDust) then
                                read(72, *) x,y,z, HdenTemp(ix,iy,iz), NdustTemp(ix,iy,iz),&
                                     &grid(iG)%abFileIndex(ix,iy,iz), dustAbunIndexTemp(ix,iy,iz)
                             else if (.not.lgGas .and. lgDust) then
                                print*, '! setSubGrids: gas must be present when the multichemistry option is used'
                                stop
                             else
                                print*, "! setSubGrids: No gas or dust with multiple grids and multichemistry!!"
                                stop
                             end if
                          else
                             if (lgGas .and. lgDust) then
                                read(72, *) x,y,z, HdenTemp(ix,iy,iz), NdustTemp(ix,iy,iz),grid(iG)%abFileIndex(ix,iy,iz)
                             else if (lgGas .and. .not. lgDust) then
                                read(72, *) x,y,z, HdenTemp(ix,iy,iz),grid(iG)%abFileIndex(ix,iy,iz)
                             else if (.not.lgGas .and. lgDust) then
                                print*, '! setSubGrids: gas must be present when the multichemistry option is used'
                                stop
                             else
                                print*, "! setSubGrids: No gas or dust with multiple grids and multichemistry!!"
                                stop
                             end if

                          end if
                       end if

                       x = grid(iG)%xAxis(1)+x*(grid(iG)%xAxis(grid(iG)%nx)-grid(iG)%xAxis(1))
                       y = grid(iG)%yAxis(1)+y*(grid(iG)%yAxis(grid(iG)%ny)-grid(iG)%yAxis(1))
                       z = grid(iG)%zAxis(1)+z*(grid(iG)%zAxis(grid(iG)%nz)-grid(iG)%zAxis(1))

                       if (ix == grid(iG)%nx) then
                          if ( abs(x-grid(iG)%xAxis(ix)) >= abs(grid(iG)%xAxis(ix)-grid(iG)%xAxis(ix-1))  ) then
                             print*, "! setSubGrids: insanity occured in setting xAxis for subGrid ", iG,dFileRead,&
                                  & grid(iG)%xAxis
                             stop
                          end if
                       end if
                       if (iy == grid(iG)%ny) then
                          if ( abs(y-grid(iG)%yAxis(iy)) >= abs(grid(iG)%yAxis(iy)-grid(iG)%yAxis(iy-1))  ) then
                             print*, "! setSubGrids: insanity occured in setting yAxis for subGrid ", iG,dFileRead,&
                                  & grid(iG)%yAxis 
                             stop
                          end if
                       end if
                       if (iz == grid(iG)%nz) then
                          if ( abs(z-grid(iG)%zAxis(iz)) >= abs(grid(iG)%zAxis(iz)-grid(iG)%zAxis(iz-1))  ) then
                             print*, "! setSubGrids: insanity occured in setting zAxis for subGrid ", iG,dFileRead,&
                                  & grid(iG)%zAxis 
                             stop
                          end if
                       end if

                       grid(iG)%xAxis(ix) = x                         
                       grid(iG)%yAxis(iy) = y                         
                       grid(iG)%zAxis(iz) = z

                       radius = 1.e10*sqrt( (grid(iG)%xAxis(ix)/1.e10)*(grid(iG)%xAxis(ix)/1.e10) + &
                            & (grid(iG)%yAxis(iy)/1.e10)*(grid(iG)%yAxis(iy)/1.e10) + &
                            & (grid(iG)%zAxis(iz)/1.e10)*(grid(iG)%zAxis(iz)/1.e10) ) 


                       ! check if this grid point is  valid nebular point
                       if (.not. lgPlaneIonization) then
                          
                          if (radius < R_in) then
                             grid(iG)%active(ix,iy,iz) = 0
                          else if (radius > R_out .and. R_out>0.) then
                             grid(iG)%active(ix,iy,iz) = 0
                          end if

                       end if
                       
                       if (grid(iG)%active(ix,iy,iz) <= 0 .and. lgDust ) then
                          NdustTemp(ix,iy,iz) = 0.
                          if (lgMultiDustChemistry) dustAbunIndexTemp(ix,iy,iz) = 0.
                       end if
                       if (grid(iG)%active(ix,iy,iz) <= 0 .and. lgGas ) HdenTemp(ix,iy,iz) = 0.

                       if (lgDust .and. lgGas) then
                          if (HdenTemp(ix,iy,iz) > 0. .or. NdustTemp(ix,iy,iz)>0.) then
                             grid(iG)%nCells =  grid(iG)%nCells + 1
                             grid(iG)%active(ix,iy,iz) =  grid(iG)%nCells
                          else
                             grid(iG)%active(ix,iy,iz) = 0
                             HdenTemp(ix,iy,iz) = 0.
                             NdustTemp(ix,iy,iz) = 0.
                             if (lgMultiDustChemistry) dustAbunIndexTemp(ix,iy,iz) = 0.
                          end if
                       else if ( lgDust .and. (.not.lgGas) ) then
                          if (NdustTemp(ix,iy,iz)>0.) then
                             grid(iG)%nCells =  grid(iG)%nCells + 1
                             grid(iG)%active(ix,iy,iz) =  grid(iG)%nCells
                          else
                             grid(iG)%active(ix,iy,iz) = 0
                             NdustTemp(ix,iy,iz) = 0.
                             if (lgMultiDustChemistry) dustAbunIndexTemp(ix,iy,iz) = 0.
                          end if
                       else if ( (.not.lgDust) .and. lgGas) then                            
                          if (HdenTemp(ix,iy,iz) > 0.) then
                             grid(iG)%nCells =  grid(iG)%nCells + 1
                             grid(iG)%active(ix,iy,iz) =  grid(iG)%nCells
                          else
                             grid(iG)%active(ix,iy,iz) = 0
                             HdenTemp(ix,iy,iz) = 0.
                          end if
                       else
                          
                          print*, '! setSubGrids: no gas and no dust? The grid is empty.'
                          stop
                       end if
                    end do
                 end do
              end do
              
              close(72)

              ! allocate grid arrays
              
              if (lgGas) then

                 allocate(grid(iG)%Hden(0:grid(iG)%nCells), stat = err)
                 if (err /= 0) then
                    print*, "! setSubGrid: can't allocate grid memory"
                    stop
                 end if
                 
                 allocate(grid(iG)%recPDF(0:grid(iG)%nCells, 1:nbins), stat = err)
                 if (err /= 0) then
                    print*, "Can't allocate grid memory, 8"
                    stop
                 end if
                 
                 allocate(grid(iG)%totalLines(0:grid(iG)%nCells), stat = err)
                 if (err /= 0) then
                    print*, "Can't allocate grid memory, 10"
                    stop
                 end if
                 
                 allocate(grid(iG)%ionDen(0:grid(iG)%nCells, 1:nElementsUsed, 1:nstages), stat = err)
                 if (err /= 0) then
                    print*, "! setSubGrid: can't allocate grid memory,ionDen"  
                    stop
                 end if
                 allocate(ionDenUsed(1:nElementsUsed, 1:nstages), stat = err)
                 if (err /= 0) then
                    print*, "! setSubGrid: can't allocate grid memory,ionDen"  
                    stop
                 end if
                 allocate(grid(iG)%Ne(0:grid(iG)%nCells), stat = err)
                 if (err /= 0) then
                    print*, "! setSubGrid: can't allocate grid memory"
                    stop
                 end if
                 allocate(grid(iG)%Te(0:grid(iG)%nCells), stat = err)
                 if (err /= 0) then
                    print*, "! setSubGrid:can't allocate grid memory"
                    stop
                 end if
                 
                 grid(iG)%Hden = 0.
                 grid(iG)%Ne = 0.
                 grid(iG)%Te = 0.        
                 grid(iG)%ionDen = 0.
                 grid(iG)%recPDF = 0.
                 grid(iG)%totalLines = 0.
                 
              end if

              if (Ldiffuse>0.) then
                 allocate(grid(iG)%LdiffuseLoc(0:grid(iG)%nCells), stat = err)
                 if (err /= 0) then
                    print*, "! setSubGrid: can't allocate grid memory : LdiffuseLoc "
                    stop
                 end if
              grid(iG)%LdiffuseLoc=0.
              end if              


              allocate(grid(iG)%opacity(0:grid(iG)%nCells, 1:nbins), stat = err)
              if (err /= 0) then
                 print*, "! setSubGrid: can't allocate grid memory : opacity "
                 stop
              end if
              allocate(grid(iG)%Jste(0:grid(iG)%nCells, 1:nbins), stat = err)
              if (err /= 0) then
                 print*, "! setSubGrid: can't allocate grid memory : Jste"
                 stop
              end if
              
              allocate(grid(iG)%escapedPackets(0:grid(iG)%nCells, 0:nbins,0:nAngleBins), stat = err)
              if (err /= 0) then
                 print*, "! setSubGrid: can't allocate grid memory : escapedPackets"
                 stop
              end if
              
              if (lgDust) then
                 allocate(grid(iG)%Ndust(0:grid(iG)%nCells), stat = err)
                 if (err /= 0) then
                    print*, "! grid: can't allocate Ndust memory"
                    stop
                 end if
                 grid(iG)%Ndust=0.

                 allocate(grid(iG)%dustAbunIndex(0:grid(iG)%nCells), stat = err)
                 if (err /= 0) then
                    print*, "! grid: can't allocate dustAbunIndex memory"
                    stop
                 end if
                 grid(iG)%dustAbunIndex=0.

                 if (.not.lgGas) then
                    allocate(grid(iG)%dustPDF(0:grid(iG)%nCells, 1:nbins), stat = err)
                    if (err /= 0) then
                       print*, "! grid: can't allocate dustPDF memory"
                       stop
                    end if                
                    grid(iG)%dustPDF = 0.
                 end if

              end if

!BS10              if (lgDebug) then
                 allocate(grid(iG)%Jdif(0:grid(iG)%nCells, 1:nbins), stat = err)
                 if (err /= 0) then
                    print*, "! grid: can't allocate grid memory"
                    stop
                 end if
                 
                 grid(iG)%Jdif = 0. 
                 
                 ! allocate pointers depending on nLines
                 if (nLines > 0) then
                    
                    allocate(grid(iG)%linePackets(0:grid(iG)%nCells, 1:nLines), stat = err)
                    if (err /= 0) then
                       print*, "! initCartesianGrid: can't allocate grid(iG)%linePackets memory"
                       stop
                    end if

                    allocate(grid(iG)%linePDF(0:grid(iG)%nCells, 1:nLines), stat = err)
                    if (err /= 0) then
                       print*, "! initCartesianGrid: can't allocate grid(iG)%linePDF memory"
                       stop
                    end if
                    
                    grid(iG)%linePackets = 0.
                    grid(iG)%linePDF     = 0.
                    
                 end if
             
!BS10              end if

              allocate(grid(iG)%lgConverged(0:grid(iG)%nCells), stat = err)
              if (err /= 0) then
                 print*, "Can't allocate memory to lgConverged array"
                 stop
              end if
              
              allocate(grid(iG)%lgBlack(0:grid(iG)%nCells), stat = err)
              if (err /= 0) then
                 print*, "Can't allocate memory to lgBlack array"
                 stop
              end if


              if (lgNeInput) then
                 allocate(grid(iG)%NeInput(0:grid(iG)%nCells), stat = err)
                 if (err /= 0) then
                    print*, "! setSubGrid: can't allocate grid memory, grid(iG)%NeInput"
                    stop
                 end if
                 grid(iG)%NeInput = 0.
              end if

              grid(iG)%opacity = 0.        
              grid(iG)%Jste = 0.        
              grid(iG)%lgConverged = 0
              grid(iG)%lgBlack = 0           
              
              do ix = 1, grid(iG)%nx
                 do iy = 1, grid(iG)%ny
                    do iz = 1, grid(iG)%nz
                       if (grid(iG)%active(ix,iy,iz)>0) then
                          
                          if (lgGas) then
                             grid(iG)%Hden(grid(iG)%active(ix,iy,iz)) = HdenTemp(ix,iy,iz)*denfac
                             grid(iG)%Te(grid(iG)%active(ix,iy,iz)) = TeStart
                          end if
                          if (lgDust) then
                             grid(iG)%Ndust(grid(iG)%active(ix,iy,iz)) = NdustTemp(ix,iy,iz)*denfac
                             if (lgMultiDustChemistry) grid(iG)%dustAbunIndex(grid(iG)%active(ix,iy,iz)) = &
                                  & dustAbunIndexTemp(ix,iy,iz)*denfac                          
                          end if
                       end if
                    end do
                 end do
              end do
              
              if (lgNeInput) then 
                 grid(iG)%NeInput = grid(iG)%Hden
                 ! 1.11 is just a starting guess for the ionization 
                 ! factor
                 grid(iG)%Hden = grid(iG)%Hden/1.11
              end if

              if (lgGas) then

                 H0in  = 1.e-5
                 
                 do ix = 1, grid(iG)%nx
                    do iy = 1, grid(iG)%ny
                       do iz = 1, grid(iG)%nz
                          
                          if (grid(iG)%active(ix,iy,iz)>0 ) then
                             
                             ! calculate ionDen for H0
                             if (lgElementOn(1)) then
                                grid(iG)%ionDen(grid(iG)%active(ix,iy,iz),elementXref(1),1) = H0in    
                                grid(iG)%ionDen(grid(iG)%active(ix,iy,iz),elementXref(1),2) = &
                                     & 1.-grid(iG)%ionDen(grid(iG)%active(ix,iy,iz),elementXref(1),1)
                             end if
                             if (lgElementOn(2)) then
                                grid(iG)%ionDen(grid(iG)%active(ix,iy,iz),elementXref(2),1) = &
                                     & grid(iG)%ionDen(grid(iG)%active(ix,iy,iz),(1),1)
                                grid(iG)%ionDen(grid(iG)%active(ix,iy,iz),elementXref(2),2) = &
                                 & (1.-grid(iG)%ionDen(grid(iG)%active(ix,iy,iz),elementXref(2),1))
                                grid(iG)%ionDen(grid(iG)%active(ix,iy,iz),elementXref(2),3) = 0.
                             end if

                             ! initialize Ne
                             grid(iG)%Ne(grid(iG)%active(ix,iy,iz)) =  grid(iG)%Hden(grid(iG)%active(ix,iy,iz))
                         
                             ! initialize all heavy ions (careful that the sum over all ionisation 
                             ! stages of a given atom doesn't exceed 1.)
                             do elem = 3, nElements
                                do ion = 1, min(elem+1,nstages)
                                   if (lgElementOn(elem)) then
                                      if (ion == 1) then
                                         grid(iG)%ionDen(grid(iG)%active(ix,iy,iz),elementXref(elem),ion) &
                                              & = grid(iG)%ionDen(grid(iG)%active(ix,iy,iz),1,1)
                                      else
                                         grid(iG)%ionDen(grid(iG)%active(ix,iy,iz),elementXref(elem),ion) = 0.
                                      end if
                                   end if
                                end do
                             end do
                          end if     ! active condition
                      
                       end do
                    end do
                 end do
                 
                 ! deallocate temp array
                 if(associated(HdenTemp)) deallocate(HdenTemp)
                 
              end if ! lgGas

              if (lgDust) then
                 if(associated(NdustTemp)) deallocate(NdustTemp)
                 if (lgMultiDustChemistry .and. associated(dustAbunIndexTemp)) deallocate(dustAbunIndexTemp)
              end if
              
              totalMass = 0.
              totalVolume = 0.
              echoVolume = 0.

              do ix = 1, grid(iG)%nx
                 do iy = 1, grid(iG)%ny
                    do iz = 1, grid(iG)%nz
                       if (grid(iG)%active(ix,iy,iz)>0) then

                          dV = getVolume(grid(iG),ix,iy,iz)

                          if (lgGas) then
                             gasCell = 0.
                             do elem = 1, nElements
                                gasCell = gasCell + &
                                     & grid(iG)%elemAbun(grid(iG)%abFileIndex(ix,iy,iz),elem)*&
                                     & aWeight(elem)*amu
                                totalMass = totalMass + &
                                     & grid(iG)%Hden(grid(iG)%active(ix,iy,iz))*dV*&
                                     & grid(iG)%elemAbun(grid(iG)%abFileIndex(ix,iy,iz),elem)*&
                                     & aWeight(elem)*amu
                             end do
                          end if
                          
                          totalVolume = totalVolume + dV
!
! echoes only
!
                          if (lgEcho) then
                             grid(iG)%echoVol(i,j,k)=vEcho(grid(iG),i,j,k,echot1,echot2,vol)
                             echoVolume = echoVolume + vol 
                          endif

                          
                          if (lgDust .and. lgMdMg .or. lgMdMh) then
                             
                             if (lgMdMh) then
                                MhMg=0.
                                do elem = 1, nElements
                                   ! transform to MdMg
                                   MhMg = MhMg+grid(iG)%elemAbun(grid(iG)%abFileIndex(ix,iy,iz),elem)*&
                                        & aWeight(elem)
                                end do
                                MhMg = 1./MhMg                             
                                grid(iG)%Ndust(grid(iG)%active(ix,iy,iz)) = grid(iG)%Ndust(grid(iG)%active(ix,iy,iz))*MhMg
                             end if

                             if (.not.lgGas) then
                                print*, '! setSubGrid: Mass to dust ratio (MdMg) cannot be used in a pure dust (noGas)&
                                     & simulation. Ndust must be used instead.'
                                stop
                             end if
                             
                             grid(iG)%Ndust(grid(iG)%active(ix,iy,iz)) = gasCell*grid(iG)%Ndust(grid(iG)%active(ix,iy,iz)) &
                                 & *grid(iG)%Hden(grid(iG)%active(ix,iy,iz))

                             denominator = 0.
                             
                             if (lgMultiDustChemistry) then
                                nsp = grid(iG)%dustAbunIndex(grid(iG)%active(ix,iy,iz))
                             else
                                nsp = 1
                             end if

                             do nspec = 1, nSpeciesPart(nsp)
                                do ai = 1, nSizes
                                   denominator = denominator + &
                                        & (1.3333*Pi*( (grainRadius(ai)*1.e-4)**3)*&
                                        & rho(dustComPoint(nsp)-1+nspec)&
                                        &*grainWeight(ai)*grainAbun(nsp,nspec))
                                end do
                             end do
                             grid(iG)%Ndust(grid(iG)%active(ix,iy,iz)) = grid(iG)%Ndust(grid(iG)%active(ix,iy,iz))/&
                                  & (denominator)
                             
                          end if

                          ! calculate total dust mass
                          if (lgDust) then

                             if (lgMultiDustChemistry) then
                                nsp = grid(iG)%dustAbunIndex(grid(iG)%active(ix,iy,iz))
                             else
                                nsp = 1
                             end if

                             do ai = 1, nsizes
                                do nspec = 1, nspeciesPart(nsp)
                                   totalDustMass = totalDustMass + &
                                        &(1.3333*Pi*((grainRadius(ai)*1.e-4)**3)*&
                                        & rho(dustComPoint(nsp)-1+nspec)*&
                                        & grainWeight(ai)*grainAbun(nsp,nspec)&
                                        & *grid(iG)%Ndust(grid(iG)%active(ix,iy,iz))*dV)
                                end do
                             end do

                          end if

                       end if
                    end do
                 end do
              end do

           if (lgDust .and. lginputDustMass) then
           do ix = 1, grid(iG)%nx
              do iy = 1, grid(iG)%ny
                 do iz = 1, grid(iG)%nz
        grid(iG)%Ndust(grid(iG)%active(ix,iy,iz)) = grid(iG)%Ndust(grid(iG)%active(ix,iy,iz)) * inputDustMass / totalDustMass
                 end do
              end do
           end do
           totalDustMass = inputDustMass
           end if

              if (taskid == 0) then
              
                 print*, 'Sub Grid: ', iG                    
                 if (lgGas) &
                    &print*, 'Total gas mass of ionized region by mass [1.e45 g]: ', totalMass


                 print*, 'Total volume of the active region [e45 cm^3]: ', totalVolume
                 if (lgEcho) print*, 'Total volume of the echo region [e45 cm^3]: ', echoVolume
                 print*, '! Number of active cells :', grid(iG)%nCells
              end if           

           end do

         end subroutine setSubGrids         

         subroutine freeGrid(grid)
 
           type(grid_type), intent(inout) :: grid
           print*, "in freeGrid"

           if (lgDust) then
              if(associated(grid%absOpac)) deallocate(grid%absOpac)
              if(associated(grid%scaOpac)) deallocate(grid%scaOpac)
              if(associated(grid%Ndust)) deallocate(grid%Ndust)
              if(associated(grid%Tdust)) deallocate(grid%Tdust)
              if (lgMultiChemistry .and. associated(grid%dustAbunIndex)) &
                   &  deallocate(grid%dustAbunIndex)
              if (.not.lgGas) then
                 if(associated(grid%dustPDF)) deallocate(grid%dustPDF)
              end if
           end if
           if (associated(grid%active)) deallocate(grid%active)
           if (associated(grid%lgConverged)) deallocate(grid%lgConverged)
           if (associated(grid%lgBlack)) deallocate(grid%lgBlack)     
           if (lgGas) then
              if (associated(grid%abFileIndex)) deallocate(grid%abFileIndex)           
              if (associated(grid%Te)) deallocate(grid%Te)
              if (associated(grid%Ne)) deallocate(grid%Ne)
              if (associated(grid%Hden)) deallocate(grid%Hden) 
              if (associated(grid%ionDen)) deallocate(grid%ionDen)
              if (associated(ionDenUsed)) deallocate(ionDenUsed)
              if (associated(grid%recPDF)) deallocate(grid%recPDF)
              if (associated(grid%totalLines)) deallocate(grid%totalLines)
           end if
           if (associated(grid%opacity)) deallocate(grid%opacity)
           if (associated(grid%Jste)) deallocate(grid%Jste)
!BS10           if (lgDebug) then
              if (associated(grid%Jdif)) deallocate(grid%Jdif)
              if (associated(grid%linePackets)) deallocate(grid%linePackets)
              if (associated(grid%linePDF)) deallocate(grid%linePDF)
!BS10           end if
           if (lgNeInput) then
               if (associated(grid%NeInput)) deallocate(grid%NeInput)
           end if
           if (associated(grid%xAxis)) deallocate(grid%xAxis)
           if (associated(grid%yAxis)) deallocate(grid%yAxis)
           if (associated(grid%zAxis)) deallocate(grid%zAxis)

           print*, "out freeGrid"
         end subroutine freeGrid

         subroutine writeGrid(grid)
           implicit none

           type(grid_type), dimension(:), intent(in) :: grid                ! grid
           
           ! local variables
           integer                     :: cellP               ! cell pointer
           integer                     :: elem                ! element counter
           integer                     :: ion                 ! ion counter
           integer                     :: ios                 ! I/O error status
           integer                     :: i,j,k,iG,ai         ! counters
           integer                     :: yTop                ! 2D index

           print* , 'in writeGrid'

           close(21)
           open(unit=21, file="output/grid0.out",  action="write", position="rewind",status="unknown", iostat = ios)
           if (ios /= 0 ) then
              print*, "! writeGrid: can't open file for writing - grid1.out"
              stop
           end if

           if (lgGas) then
              ! open files for writing
              close(20)
              open(unit=20, file="output/grid1.out", action="write",position="rewind",status="unknown", iostat = ios)
              if (ios /= 0 ) then
                 print*, "! writeGrid: can't open file for writing - grid1.out"
                 stop
              end if
              close(30)
              open(unit=30, file="output/grid2.out", action="write",position="rewind",status="unknown", iostat = ios)   
              if (ios /= 0 ) then
                 print*, "! writeGrid: can't open file for writing - grid2.out"
                 stop
              end if
           end if
           if (lgDust) then
              close(50)
              open(unit=50, file="output/dustGrid.out", action="write",position="rewind",status="unknown", iostat = ios)
              if (ios /= 0 ) then
                 print*, "! writeGrid: can't open file for writing - dustGrid.out"
                 stop
              end if
           end if

           do iG = 1, nGrids

              ! write nx, ny, nz, to file
              write(21, *) nGrids
              write(21, *) grid(iG)%nx, grid(iG)%ny, grid(iG)%nz, grid(iG)%nCells, grid(iG)%motherP, R_out

              ! write x, y, and z axis to the file
              do i = 1, grid(iG)%nx
                 write(21,*) grid(iG)%xAxis(i)
              end do
              do i = 1, grid(iG)%ny
                 write(21,*) grid(iG)%yAxis(i)
              end do
              do i = 1, grid(iG)%nz
                 write(21,*) grid(iG)%zAxis(i)
              end do

              if (iG>1 .or. (.not. lg2D)) then
                 yTop = grid(iG)%ny
              else if (iG == 1 .and. lg2D) then
                 yTop = 1
              end if


              ! write the rest of the grid to files
              do i = 1, grid(iG)%nx
                 do j = 1, yTop
                 do k = 1, grid(iG)%nz

                    cellP = grid(iG)%active(i,j,k)

                    if (cellP<0) cellP=0
                    
                    write(21, *) grid(iG)%active(i,j,k), grid(iG)%lgConverged(cellP), &
                         & grid(iG)%lgBlack(cellP)


                    if (lgGas) then
                       if (lgMultiChemistry) then
                          write(20, *) grid(iG)%Te(cellP), grid(iG)%Ne(cellP), &
                               & grid(iG)%Hden(cellP), grid(iG)%abFileIndex(i,j,k)
                       else
                          write(20, *) grid(iG)%Te(cellP), grid(iG)%Ne(cellP), &
                               & grid(iG)%Hden(cellP)
                       end if
                    
                       do elem = 1, nElements
                          if (lgElementOn(elem)) then

                             write(30, *) (grid(iG)%ionDen(cellP,elementXref(elem),ion), ion =1, &
                                  & min(elem+1,nstages))
                          end if
                       end do
                    end if
                    if (lgDust) then
                       
                       if (lgMultiChemistry) then                          
                          write(50, *) grid(iG)%Ndust(cellP), grid(iG)%dustAbunIndex(cellP)               
                       else
                          write(50, *) grid(iG)%Ndust(cellP)
                       end if
                       do ai = 0, nSizes
                          write(50, *) (grid(iG)%Tdust(elem,ai,cellP),' ', elem=0,nSpeciesMax)
                       end do

                    end if

                 end do
              end do
           end do

        end do

        ! close files
        close(21)
        if (lgGas) then
           close(20)
           close(30)
        end if
        if (lgDust) then
           write(50,*) ' '
           write(50,*) 'Total dust mass [1.e45 g]: ', totalDustMass
           write(50,*) 'Total dust mass [Msol]: ', totalDustMass*5.028e11
           close(50)
        end if

        ! stellar parameters
        close(42) 
        open(unit=42, file="output/photoSource.out", action="write",position="rewind",status="unknown", iostat = ios)   
        if (ios /= 0 ) then
            print*, "! writeGrid: can't open file for writing - photoSource.out"
            stop
        end if

        
        write(42, *) nStars, ' number of photon sources'
        do i = 1, nStars
           write(42, *) "'",trim(contShapeIn(i)),"'", TStellar(i), LStar(i), nPhotons(i), &
                & starPosition(i)%x/grid(1)%xAxis(grid(1)%nx),&
                & starPosition(i)%y/grid(1)%yAxis(grid(1)%ny),&
                & starPosition(i)%z/grid(1)%zAxis(grid(1)%nz),&
                & trim(spID(i)), tStep(i)
           if (contShapeIn(i)=='powerlaw') write(42,*) pwlIndex
        end do

        write(42, *) '(contShape, T_eff[K], L_* [E36 erg/s], nPackets, (x,y,z) position, spID, tstep)'
        close(42)

        ! general simulation parameters
        close(40)
        open(unit=40, file="output/grid3.out", action="write",position="rewind",status="unknown", iostat = ios)   
        if (ios /= 0 ) then
            print*, "! writeGrid: can't open file for writing - grid3.out"
            stop
        end if

        write(40, *) nGrids
        write(40, *) convWriteGrid, ' convWriteGrid'
        write(40, *) lgAutoPackets, convIncPercent, nPhotIncrease, maxPhotons, ' lgAutoPackets'
        write(40, *) lgSymmetricXYZ, ' lgSymmetricXYZ'
        write(40, *) lgTalk, ' lgTalk'
        write(40, *) lg1D, ' lg1D'
        write(40, *) nbins, ' nbins'
        write(40, *) nuStepSize, ' nuStepSize'
        write(40, *) nuMax,' nuMax'
        write(40, *) nuMin, ' nuMin'
        write(40, *) R_in, ' R_in'
        write(40, *) XHIlimit, ' XHIlimit'
        write(40, *) maxIterateMC, minConvergence, ' maxIterateMC'
        write(40, *) lgDebug, ' lgDebug'
        write(40, *) lgPlaneIonization, ' lgPlaneIonization'
        write(40, *) nAbComponents, ' nAbComponents'
        do i=1,nAbComponents
           write(40,*) '"',abundanceFile(i),'"'
        end do
        write(40, *) lgOutput, ' lgOutput'
        write(40, *) dxSlit,dySlit,' dxSlit,dySlit'
        write(40, *) lgDust, lgDustConstant, ' lgDust, lgDustConstant'
        write(40, *) lgMultiDustChemistry, nDustCOmponents, ' lgMultiDustChemistry, nDustComponents'
        if (lgDust) then
           do i = 1, nDustComponents
              write(40, *) '"',trim(dustSpeciesFile(i)),'"', ' dustFile'
           end do
        else
           write(40, *) 'none', ' dustFile'
        end if
        write(40, *) '"',trim(dustFile(2)),'"', ' dustFile'
        write(40, *) lgGas, ' lgGas'
        write(40, *) lgRecombination, ' lgRecombination'
        if (lgDust) then
           write(40, *) nSpeciesMax, nSpecies, nSizes, ' nSpeciesMax, nSpecies, nSizes'
           write(40, *) (nSpeciesPart(i), i=1, nDustComponents) , ' Partial nspecies'
        else
           write(40, *) 1, 1, 1, ' nSpeciesMax, nSpecies, nSizes'
           write(40, *) 1, ' Partial nspecies'
        end if
        write(40, *) resLinesTransfer, 'resLinesTransfer'
        write(40, *) lgDustScattering, 'lgDustScattering'
        write(40, *) nAngleBins, ' nAngleBins'
        if (nanglebins>0) then
           write(40, *) viewPointTheta, ' inclination theta'
           write(40, *) viewPointPhi, ' inclination theta'
        end if
        write(40, *) contCube(1),contCube(2), ' continuumCube'
        write(40, *) lgPhotoelectric, ' lgPhotoelectric'
        write(40, *) lgTraceHeating, ' lgTraceHeating'
        write(40, *) Ldiffuse, ' Ldiffuse'
        write(40, *) Tdiffuse, ' Tdiffuse'
        write(40, *) shapeDiffuse, ' shapeDiffuse'
        write(40, *) nPhotonsDiffuse, 'nPhotonsDiffuse'
        write(40, *) emittingGrid, ' emittingGrid'
        write(40, *) nstages, ' emittingGrid'
        write(40, *) lgMultistars, ' lgMultiStars'
        write(40,*)  lg2D, ' 2D geometry?'
        write(40,*) "'",trim(home),"' home()"
        write(40,*) lgEcho, echot1, echot2, echoTemp," Echo on/off"
        write(40,*) lgNosource," NoSourceSED"
        ! close file
        close(40)
     
        print*, 'out writeGrid'
      end subroutine writeGrid


      ! this function returns the volume of a cell in [e45 cm^3]
      function getVolume(grid,xP, yP, zP)
        implicit none

        type(grid_type),intent(in) :: grid              ! the grid
        
        integer, intent(in)        :: xP, yP, zP        ! cell indeces  

        real                       :: getVolume         ! volume of the cell [e45 cm^3]

        ! local variables
         
        real                       :: dx, &             ! cartesian axes increments
&                                      dy, &             ! in [cm] 
&                                      dz                ! 

        if (lg1D) then
           if (nGrids>1) then
              print*, '! getVolume: 1D option and multiple grids options are not compatible'
              stop
           end if

           if (xP == 1) then              

              getVolume = 4.*Pi* ( (grid%xAxis(xP+1)/1.e15)**3.)/3.


           else if ( xP==grid%nx) then
 
              getVolume = Pi* ( (3.*(grid%xAxis(xP)/1.e15)-(grid%xAxis(xP-1)/1.e15))**3. - &
                   & ((grid%xAxis(xP)/1.e15)+(grid%xAxis(xP-1)/1.e15))**3. ) / 6.

           else 

              getVolume = Pi* ( ((grid%xAxis(xP+1)/1.e15)+(grid%xAxis(xP)/1.e15))**3. - &
                   & ((grid%xAxis(xP-1)/1.e15)+(grid%xAxis(xP)/1.e15))**3. ) / 6.

           end if

           getVolume = getVolume/8.

        else

           if ( (xP>1) .and. (xP<grid%nx) ) then
              dx = abs(grid%xAxis(xP+1)-grid%xAxis(xP-1))/2.
           else if ( xP==1 ) then
              if (lgSymmetricXYZ) then
                 dx = abs(grid%xAxis(xP+1)-grid%xAxis(xP))/2.
              else 
                 dx = abs(grid%xAxis(xP+1)-grid%xAxis(xP))
              end if
           else if ( xP==grid%nx ) then
              dx = abs(grid%xAxis(xP)  -grid%xAxis(xP-1))
           end if
        
           if ( (yP>1) .and. (yP<grid%ny) ) then
              dy = abs(grid%yAxis(yP+1)-grid%yAxis(yP-1))/2.
           else if ( yP==1 ) then
              if (lgSymmetricXYZ) then
                 dy = abs(grid%yAxis(yP+1)-grid%yAxis(yP))/2.
              else
                dy = abs(grid%yAxis(yP+1)-grid%yAxis(yP))
             end if
          else if ( yP==grid%ny ) then
             dy = abs(grid%yAxis(yP)  -grid%yAxis(yP-1))
          end if

          if ( (zP>1) .and. (zP<grid%nz) ) then    
             dz = abs(grid%zAxis(zP+1)-grid%zAxis(zP-1))/2.    
          else if ( zP==1 ) then    
             if (lgSymmetricXYZ) then
                dz = abs(grid%zAxis(zP+1)-grid%zAxis(zP))/2.
             else
                dz = abs(grid%zAxis(zP+1)-grid%zAxis(zP))
             end if
          else if ( zP==grid%nz ) then    
             dz = abs(grid%zAxis(zP)-grid%zAxis(zP-1))
          end if

          dx = dx/1.e15
          dy = dy/1.e15
          dz = dz/1.e15
      

          ! calculate the volume
          getVolume = dx*dy*dz


       end if

    end function getVolume

    subroutine resetGrid(grid)
      implicit none

      type(grid_type), intent(inout) :: grid(maxGrids)  ! the 3d grids

      ! local variables

      logical, save :: lgfirst = .true.
      integer                        :: cellP ! cell pointer
      integer                        :: err,ios   ! I/O error status
      integer                        :: i,j,k ! counters
      integer                        :: elem,&! 
&                                       ion,i1 ! counters
      real                           :: p0,p00,p1,p2,p3,p4,p5,p6,p7, ind
      real, pointer                  :: p(:)                                        
      integer :: iCount, nuCount, iG, ai      ! counters
      integer :: g0,g1
      integer :: nEdges  
      integer :: nElec
      integer :: outshell
      integer :: totCells, totcellsloc
      integer :: yTop, xPmap
        

      integer, parameter :: maxLim = 10000
      integer, parameter :: nSeries = 17
      
      real, dimension(450) :: ordered
      real, dimension(17)  :: seriesEdge
      real                 :: radius

      real                 :: nuStepSizeLoc

      print*, 'resetGrid in'
      

      ! read stellar parameters
      close(72) 
      open(unit=72, file="output/photoSource.out",action="read", position="rewind",status="old", iostat = ios)   
      if (ios /= 0 ) then
         print*, "! writeGrid: can't open file for reading - photoSource.out"
         stop
      end if

        
       read(72, *) nStars
       print*, nStars, ' photon sources'
       print*, '(contShape, T_eff[K], L_* [E36 erg/s], nPackets, (x,y,z) position [cm])'
       allocate(TStellar(nStars))
       allocate(LStar(nStars))
       allocate(nPhotons(nStars))
       allocate(starPosition(nStars))
       allocate(contShape(nStars))
       allocate(contShapeIn(nStars))
       allocate(spID(nStars))
       allocate(tStep(nStars))
       allocate(deltaE(0:nStars))

       do i = 1, nStars
          read(72, *) contShape(i), TStellar(i), LStar(i), nPhotons(i), starPosition(i)%x,starPosition(i)%y,&
               &starPosition(i)%z,spID(i), tStep(i)
          contShapeIn(i)=contShape(i)
          if (contShape(i)=='powerlaw') read(72,*) pwlIndex
          print*, i, contShape(i), TStellar(i), LStar(i), nPhotons(i), starPosition(i)%x,starPosition(i)%y,&
               &starPosition(i)%z, spID(i), tStep(i)
          print*, pwlIndex
          deltaE(i) = Lstar(i)/nPhotons(i)
          print*, 'deltaE', deltaE(i), i 
       end do
       close(72)

      ! read in file containing general simulation parameters
      close(77)
      open(unit=77, file='output/grid3.out', action="read",position='rewind',  &
&          status='old', iostat = err)
      if (err /= 0) then
         print*, "! resetMotherGrid: error opening file grid3.out"
         stop
      end if

      read(77, *) nGrids     
      read(77, *) convWriteGrid
      read(77, *) lgAutoPackets, convIncPercent, nPhotIncrease, maxPhotons
      read(77, *) lgSymmetricXYZ
      read(77, *) lgTalk
      read(77, *) lg1D
      read(77, *) nbins
      read(77, *) nuStepSize
      read(77, *) nuMax
      read(77, *) nuMin
      read(77, *) R_in
      read(77, *) XHIlimit
      read(77, *) maxIterateMC, minConvergence
      read(77, *) lgDebug
      read(77, *) lgPlaneIonization
      read(77, *) nAbComponents
      allocate(abundanceFile(nAbComponents))
      if (nAbComponents>1) then
         lgMultiChemistry = .true.
      else
         lgMultiChemistry = .false.
      end if
      do i = 1, nAbComponents
         read(77, *) abundanceFile(i)
      end do                  
      read(77, *) lgOutput
      read(77, *) dxSlit, dySlit
      read(77, *) lgDust, lgDustConstant
      read(77, *) lgMultiDustChemistry, nDustComponents
      allocate(dustSpeciesFile(1:nDustComponents))
      do i = 1, nDustComponents
         read(77, *) dustSpeciesFile(i)
      end do
      read(77, *) dustFile(2)
      read(77, *) lgGas
      read(77, *) lgRecombination
      read(77, *) nSpeciesMax,nSpecies, nSizes
      allocate(nSpeciesPart(1:nDustComponents))
      read(77, *) (nSpeciesPart(i), i = 1, nDustComponents)
      read(77, *) resLinesTransfer
      read(77, *) lgDustScattering
      read(77, *) nAngleBins
      if (nanglebins>0) then
         allocate(viewPointTheta(0:nAngleBins), stat=err)
         if (err /= 0) then
            print*, '! readInput: allocation error for viewPoint pointer'
            stop
         end if
         viewPointTheta = 0.
         allocate(viewPointPhi(0:nAngleBins), stat=err)
         if (err /= 0) then
            print*, '! readInput: allocation error for viewPoint pointer'
            stop
         end if
         viewPointPhi = 0.
         read(77, *) (viewPointTheta(i), i = 1, nAngleBins)
         read(77, *) (viewPointPhi(i), i = 1, nAngleBins)
      end if
      read(77, *) contCube(1),contCube(2)
      read(77, *) lgPhotoelectric
      read(77, *) lgTraceHeating
      read(77, *) Ldiffuse
      read(77, *) Tdiffuse
      read(77, *) shapeDiffuse
      read(77, *) nPhotonsDiffuse
      if (Ldiffuse>0. .and. nPhotonsDiffuse>0 ) then
         deltaE(0) = Ldiffuse/nPhotonsDiffuse
      end if
      read(77, *) emittingGrid
      read(77, *) nstages
      read(77, *) lgMultistars
      read(77, *) lg2D
      read(77, *) home
      read(77, *) lgEcho, echot1, echot2, echoTemp
      read(77,*) lgNosource

      if (taskid == 0) then
         print*,  nGrids,'nGrids'
         print*,  convWriteGrid, ' convWriteGrid'
         print*,  lgAutoPackets, convIncPercent, nPhotIncrease, maxPhotons, &
              & ' lgAutoPackets, convIncPercent, nPhotIncrease, maxPhotons'
         print*,  lgSymmetricXYZ, ' lgSymmetricXYZ'
         print*,  lgTalk, ' lgTalk'
         print*,  lg1D, ' lg1D'
         print*,  nbins, ' nbins'
         print*,  nuStepSize, ' nuStepSize.'
         print*,  nuMax, ' nuMax'
         print*,  nuMin, ' nuMin'
         print*,  R_in, ' R_in'
         print*,  XHIlimit, ' XHIlimit'
         print*,  maxIterateMC, minConvergence, ' maxIterateMC, minConvergence'
         print*,  lgDebug,  ' lgDebug'
         print*,  lgPlaneIonization, ' lgPlaneIonization'
         print*,  nAbComponents, ' nAbComponents'
         print*,  lgMultiChemistry, ' ',abundanceFile, ' lgMultiChemistry, abundanceFile'
         print*,  lgOutput, ' lgOutput'
         print*,  dxSlit, dySlit, ' dxSlit, dySlit'
         print*,  lgDust, lgDustConstant, ' lgDust, lgDustConstant'
         print*,  lgMultiDustChemistry, nDustComponents, ' lgMultiDustChemistry, nDustComponents'
         do i =1 , nDustComponents
            print*,  dustSpeciesFile(i), ' dustFiles'
         end do
         print*,  dustFile(2)
         print*,  lgGas, ' lgGas'
         print*,  lgRecombination, ' lgRecombination'
         print*,  nSPeciesMax, nSpecies, nSizes, ' nSPeciesMax, nSpecies, nSizes'
         print*,  (nSpeciesPart(i), i = 1, nDustComponents)
         print*,  resLinesTransfer, ' resLinesTransfer'
         print*,  lgDustScattering, ' lgDustScattering'
         print*,  nAngleBins, ' nAngleBins'
         print*,  contCube(1),contCube(2), 'continuumCube'
         print*,  lgPhotoelectric, ' lgPhotoelectric'
         print*,  lgTraceHeating, ' lgTraceHeating'
         print*,  Ldiffuse, ' Ldiffuse'
         print*,  Tdiffuse, ' Tdiffuse'
         print*,  shapeDiffuse, ' shapeDiffuse'
         print*,  nPhotonsDiffuse, ' nPhotonsDiffuse'
         print*,  emittingGrid, ' emittingGrid'
         print*,  nstages, ' nstages'
         print*,  lgMultistars, ' lgMultiStars'
         print*,  lg2D, ' lg2D'
         print*,  "'",trim(home),"' home"
         print*,  lgEcho, echot1, echot2, echoTemp
         print*,  lgNosource," NoSourceSED"
      end if
      close(77)


      allocate(p(nstages))
      p=0.
      allocate(lgDataAvailable(3:nElements, nstages))

      close(89)
      open(unit=89, file='output/grid0.out',  action="read",position='rewind',  &
           &          status='old', iostat = err)
      if (err /= 0) then
         print*, "! resetGrid: error opening file grid0.out"
         stop
      end if

      if (lgGas) then
         close(78)
         open(unit=78, file='output/grid1.out',  action="read",position='rewind',  &
              &          status='old', iostat = err)
         if (err /= 0) then
            print*, "! resetGrid: error opening file grid1.out"
            stop
         end if
         ! open the grid2.out file for later
         close(79)
         open(unit=79, file='output/grid2.out', action="read", position='rewind',  &
              &          status='old', iostat = err)
         if (err /= 0) then
            print*, "! resetGrid: error opening file grid2.out"
            stop
         end if                  
      end if
      if (lgDust) then
         close(88)
         open(unit=88, file='output/dustGrid.out', action="read", position='rewind',  &
              &          status='old', iostat = err)
         if (err /= 0) then
            print*, "! resetGrid: error opening file dustGrid.out"
            stop
         end if
      end if

      allocate(dl(1:nGrids), stat = err)
      if (err /= 0) then
         print*, "! resetGrid: can't allocate grid memory"
         stop
      end if

      totCells = 0
      totCellsloc = 0
      do iG = 1, nGrids
         read(89, *) nGrids
         read(89, *) nxIn(iG), nyIn(iG), nzIn(iG), grid(iG)%nCells, grid(iG)%motherP, R_out

         print*, nxIn(iG), nyIn(iG), nzIn(iG), grid(iG)%nCells, grid(iG)%motherP, R_out

         ! initialize cartesian grid
         call initCartesianGrid(grid(iG),nxIn(iG), nyIn(iG), nzIn(iG)) 
         
         
         if (iG>1) then
            grid(iG)%elemAbun = grid(1)%elemAbun
         end if

         if (lgPlaneIonization .and. iG==1) then
            allocate(planeIonDistribution(grid(iG)%nx,grid(iG)%nz), stat = err)
            if (err /= 0) then
               print*, "! resetMotherGrid: can't allocate dl memory"
               stop
            end if
            planeIonDistribution = 0
         end if


         ! Allocate grid arrays       
         if (lgGas) then
            allocate(grid(iG)%Hden(0:grid(iG)%nCells), stat = err)
            if (err /= 0) then
               print*, "! resetGrid: can't allocate grid memory"
               stop
            end if
            
            allocate(grid(iG)%Ne(0:grid(iG)%nCells), stat = err)
            if (err /= 0) then
               print*, "! resetGrid: can't allocate grid memory"
               stop
            end if
         
            allocate(grid(iG)%Te(0:grid(iG)%nCells), stat = err)
            if (err /= 0) then
               print*, "! resetGrid:can't allocate grid memory"
               stop
            end if

            allocate(grid(iG)%ionDen(0:grid(iG)%nCells, 1:nElementsUsed, 1:nstages), stat = err)
            if (err /= 0) then
               print*, "! resetGrid: can't allocate grid memory,ionDen"  
               stop
            end if

            allocate(ionDenUsed(1:nElementsUsed, 1:nstages), stat = err)
            if (err /= 0) then
               print*, "! resetGrid: can't allocate grid memory,ionDen"  
               stop
            end if
            
            allocate(grid(iG)%recPDF(0:grid(iG)%nCells, 1:nbins), stat = err)
            if (err /= 0) then
               print*, "Can't allocate grid memory, 8"
               stop
            end if
            
            allocate(grid(iG)%totalLines(0:grid(iG)%nCells), stat = err)
            if (err /= 0) then
               print*, "Can't allocate grid memory, 10"
               stop
            end if
            
            grid(iG)%Hden = 0.
            grid(iG)%Ne = 0.
            grid(iG)%Te = 0.        
            grid(iG)%ionDen = 0.
            grid(iG)%recPDF = 0.
            grid(iG)%totalLines = 0.
            
         end if

         if (Ldiffuse>0.) then
            allocate(grid(iG)%LdiffuseLoc(0:grid(iG)%nCells), stat = err)
            if (err /= 0) then
               print*, "! resetGrid: can't allocate grid memory : LdiffuseLoc "
               stop
            end if
            grid(iG)%LdiffuseLoc=0.
         end if

         totCells = totCells + grid(iG)%nCells
         if (emittingGrid>0 .and. iG<=emittingGrid) then
            totCellsLoc = totCellsLoc + grid(iG)%nCells
         elseif (emittingGrid==0) then
            totCellsLoc = totCellsLoc + grid(iG)%nCells
         end if
            


         allocate(grid(iG)%opacity(0:grid(iG)%nCells, 1:nbins), stat = err)
         if (err /= 0) then
            print*, "! resetGrid: can't allocate grid memory,opacity "
            stop
         end if
        
         allocate(grid(iG)%Jste(0:grid(iG)%nCells, 1:nbins), stat = err)
         if (err /= 0) then
            print*, "! resetGrid: can't allocate Jste grid memory"
            stop
         end if
         
!BS10         if (lgDebug) then
            allocate(grid(iG)%Jdif(0:grid(iG)%nCells, 1:nbins), stat = err)
            if (err /= 0) then
               print*, "! grid: can't allocate grid memory"
               stop
            end if
            
            ! allocate pointers depending on nLines
            if (nLines > 0) then
               
               allocate(grid(iG)%linePackets(0:grid(iG)%nCells, 1:nLines), stat = err)
               if (err /= 0) then
                  print*, "! resetGrid: can't allocate grid(iG)%linePackets memory"
                  stop
               end if
               
               allocate(grid(iG)%linePDF(0:grid(iG)%nCells, 1:nLines), stat = err)
               if (err /= 0) then
                  print*, "! resetGrid: can't allocate grid(iG)%linePDF memory"
                  stop
               end if

               grid(iG)%linePackets = 0.
               grid(iG)%linePDF     = 0.

            end if
            
            grid(iG)%Jdif = 0. 

!BS10         end if

         allocate(grid(iG)%lgConverged(0:grid(iG)%nCells), stat = err)
         if (err /= 0) then
            print*, "Can't allocate memory to lgConverged array"
            stop
         end if

         allocate(grid(iG)%lgBlack(0:grid(iG)%nCells), stat = err)
         if (err /= 0) then
            print*, "Can't allocate memory to lgBlack array"
            stop
         end if
         
         allocate(grid(iG)%escapedPackets(0:grid(iG)%nCells, 0:nbins,0:nAngleBins), stat = err)
         if (err /= 0) then
            print*, "! resetGrid: can't allocate grid memory : Jste"
            stop
         end if

         if (lgDust) then
            allocate(grid(iG)%Ndust(0:grid(iG)%nCells), stat = err)
            if (err /= 0) then
               print*, "Can't allocate grid memory, Ndust"
               stop
            end if
            allocate(grid(iG)%Tdust(0:nSpeciesMax, 0:nSizes, 0:grid(iG)%nCells), stat = err)
            if (err /= 0) then
               print*, "Can't allocate grid memory, Tdust"
               stop
            end if
            grid(iG)%Ndust = 0.
            grid(iG)%Tdust = 0.

            if (.not.lgGas) then
               
               allocate(grid(iG)%dustPDF(0:grid(iG)%nCells, 1:nbins), stat = err)
               if (err /= 0) then
                  print*, "Can't allocate grid memory, dustPDF"
                  stop
               end if
               grid(iG)%dustPDF=0.
            end if
            
         end if
         
         grid(iG)%opacity = 0.        
         grid(iG)%Jste = 0.        
         grid(iG)%lgConverged = 0
         grid(iG)%lgBlack = 0



         ! axis points
         do i = 1, grid(iG)%nx
            read(89, *) grid(iG)%xAxis(i)
         end do
         do i = 1, grid(iG)%ny
            read(89, *) grid(iG)%yAxis(i)
         end do
         do i = 1, grid(iG)%nz
            read(89, *) grid(iG)%zAxis(i)
         end do

         if (lg2D) then
            yTop = 1
         else
            yTop = grid(iG)%ny
         end if

         ! read the rest of the files into grid 
         do i = 1, grid(iG)%nx
            do j = 1, yTOp
               do k = 1, grid(iG)%nz
                  
                  read(89, *) cellP, p0,p00
!print*, cellP
                  grid(iG)%active(i,j,k) = cellP
                  if (cellP<0) cellP = 0

                  grid(iG)%lgConverged(cellP)=p0
                  grid(iG)%lgBlack(cellP)=p00
                  
                  if (lgGas) then
                  
                     i1=1
                     if (lgMultiChemistry) then
                        read(78, *) p1,p2,p3,i1
                     else
                        read(78, *) p1,p2,p3
                     end if

                     grid(iG)%Te(cellP)=p1
                     grid(iG)%Ne(cellP)=p2
                     grid(iG)%Hden(cellP)=p3
                     grid(iG)%abFileIndex(i,j,k)=i1

                     do elem = 1, nElements
                        if (lgElementOn(elem)) then
                           
                           read(79, *) (p(ion), ion =1,min(elem+1,nstages))
                           
                           grid(iG)%ionDen(cellP,elementXref(elem),:)=p

                        end if
                     end do
                  end if
                  
                  if (lgDust) then
                     if (lgMultiDustChemistry) then
                        read(88, *) p00, ind
                        grid(iG)%dustAbunIndex(cellP) = ind
                     else
                        read(88, *) p00
                     end if
                     
                     do ai = 0, nSizes
                        read(88, *) (grid(iG)%Tdust(elem,ai,cellP), elem = 0,nSpeciesMax)
                     end do
                     grid(iG)%Ndust(cellP) = p00
                  end if

               end do
            end do
         end do



         if (lg2D .and. iG==1) then
            allocate(TwoDscaleJ(grid(iG)%nCells))
            TwoDscaleJ = 1.

            do i = 1, grid(ig)%nx
               do j = 2, grid(ig)%ny
                  do k = 1, grid(ig)%nz
                     radius = 1.e10*sqrt( (grid(ig)%xAxis(i)/1.e10)*&
                          &(grid(ig)%xAxis(i)/1.e10) + &
                          &(grid(ig)%yAxis(j)/1.e10)*(grid(ig)%yAxis(j)/1.e10) ) 
                      
                     call locate(grid(ig)%xAxis, radius, xPmap)
                     if (xPmap < grid(ig)%nx) then
                        if (radius >= (grid(ig)%xAxis(xPmap)+grid(ig)%xAxis(xPmap+1))/2.) &
                             & xPmap = xPmap+1
                     end if
                     grid(ig)%active(i,j,k) = grid(ig)%active(xPmap, 1, k)
                     
                     if (grid(ig)%active(xPmap,1,k)>0) &
                          & TwoDScaleJ(grid(ig)%active(xPmap,1,k)) = &
                          & TwoDScaleJ(grid(ig)%active(xPmap,1,k))+1.
                      
                   end do
                end do
             end do             
          end if

 


         ! find geometric corrections
         grid(iG)%geoCorrX = (grid(iG)%xAxis(grid(iG)%nx) - grid(iG)%xAxis(grid(iG)%nx-1))/2.
         if (.not. lg1D) then
            grid(iG)%geoCorrY = (grid(iG)%yAxis(grid(iG)%ny) - grid(iG)%yAxis(grid(iG)%ny-1))/2.
            grid(iG)%geoCorrZ = (grid(iG)%zAxis(grid(iG)%nz) - grid(iG)%zAxis(grid(iG)%nz-1))/2.
         else
            grid(iG)%geoCorrY = 0.
            grid(iG)%geoCorrZ = 0.
         end if
         if (taskid == 0) print*, "Geometric grid corrections at grid: ", iG, &
           & grid(iG)%geoCorrX, grid(iG)%geoCorrY, grid(iG)%geoCorrZ

         ! find linear increment
         dl(iG) =  abs(grid(iG)%xAxis(2) - grid(iG)%xAxis(1))
         do i = 2, grid(iG)%nx-1
            dl(iG) = min(dl(iG), abs(grid(iG)%xAxis(i+1)-grid(iG)%xAxis(i)) )
         end do
         do i = 1, grid(iG)%ny-1
            dl(iG) = min(dl(iG), abs(grid(iG)%yAxis(i+1)-grid(iG)%yAxis(i)) )
         end do
         do i = 1, grid(iG)%nz-1
            dl(iG) = min(dl(iG), abs(grid(iG)%zAxis(i+1)-grid(iG)%zAxis(i)) )
         end do
         dl(iG) = dl(iG)/50.                                                              
         
         
      end do ! closes nGrids loop
         
      if (emittingGrid>0) then
         if (totCellsloc<=0) then
            print*, '! resetGrid: insanity totCellsloc'
            stop
         end if
         nPhotonsDiffuseLoc = nPhotonsDiffuse/totCellsLoc
      else
         if (totCells<=0) then
            print*, '! resetGrid: insanity totCellsloc'
            stop
         end if
         nPhotonsDiffuseLoc = nPhotonsDiffuse/totCells
      end if

      ! close files
      close(89)
      
      if (lgGas) then
         close(78)
         close(79)
      end if

      if (lgDust) close(88)

      ! locate the origin of the axes
      call locate(grid(1)%xAxis, 0., iOrigin)
      call locate(grid(1)%yAxis, 0., jOrigin)
      call locate(grid(1)%zAxis, 0., kOrigin)

      if (taskid == 0) print*, 'Mothergrid origin at cell:  ' , iOrigin, jOrigin, kOrigin

      if (associated(p)) deallocate(p)

    end subroutine resetGrid      

    subroutine setStarPosition(xA,yA,zA,grid)
      implicit none
      
      type(grid_type), intent(inout) :: grid(maxGrids)  ! the 3d grids

      real, dimension(:) :: xA,yA,zA
      
      Integer :: i, xP,yP,zP, nxA,nyA,nzA

      nxA = size(xA)
      nyA = size(yA)
      nzA = size(zA)

      allocate(starIndeces(nStars,4))

      do i = 1, nStars

         starPosition(i)%x = starPosition(i)%x*xA(nxA)
         starPosition(i)%y = starPosition(i)%y*yA(nyA)
         starPosition(i)%z = starPosition(i)%z*zA(nzA)

         call locate(xA, starPosition(i)%x, xP)
         if (xP<nxA) then
            if (starPosition(i)%x > & 
                 & (xA(xP)+xA(xP+1))/2.) &
                 xP=xP+1
         end if

         call locate(yA, starPosition(i)%y, yP)
         if (yP<nyA) then         
            if (starPosition(i)%y > &
                 & (yA(yP)+yA(yP+1))/2.) &
                 yP=yP+1
         end if

         call locate(zA, starPosition(i)%z, zP)
         if (zP<nzA) then   
            if (starPosition(i)%z > & 
                 & (zA(zP)+zA(zP+1))/2.) &
                 zP=zP+1
         end if
         
         
         if (grid(1)%active(xp,yp,zp)>=0) then
            starIndeces(i,1) = xP
            starIndeces(i,2) = yP
            starIndeces(i,3) = zP
            starIndeces(i,4) = 1

         else

            starIndeces(i,4) = abs(grid(1)%active(xp,yp,zp))

            nxA = grid(starIndeces(i,4))%nx
            nyA = grid(starIndeces(i,4))%ny
            nzA = grid(starIndeces(i,4))%nz

            call locate(grid(starIndeces(i,4))%xAxis, starPosition(i)%x, xP)
            if (xP<nxA) then
               if (starPosition(i)%x > & 
                    & (grid(starIndeces(i,4))%xAxis(xP)+grid(starIndeces(i,4))%xAxis(xP+1))/2.) &
                    xP=xP+1
            end if
            
            call locate(grid(starIndeces(i,4))%yAxis, starPosition(i)%y, yP)
            if (yP<nyA) then         
               if (starPosition(i)%y > &
                    & (grid(starIndeces(i,4))%yAxis(yP)+grid(starIndeces(i,4))%yAxis(yP+1))/2.) &
                    yP=yP+1
            end if

            call locate(grid(starIndeces(i,4))%zAxis, starPosition(i)%z, zP)
            if (zP<nzA) then   
               if (starPosition(i)%z > & 
                    & (grid(starIndeces(i,4))%yAxis(zP)+grid(starIndeces(i,4))%yAxis(zP+1))/2.) &
                    zP=zP+1
            end if

            starIndeces(i,1) = xP
            starIndeces(i,2) = yP
            starIndeces(i,3) = zP

         end if

      end do         

    end subroutine setStarPosition

!
! =Light Echo routines, BEKS 2010==========================================
!
      function vEcho(grid,xP,yP,zP,t1i,t2i,volume)
! return the fraction of the grid cell containing an echo
      implicit none
      type(grid_type),intent(in) :: grid              ! the grid
      integer, intent(in)        :: xP, yP, zP        ! cell indeces  
      double precision :: x,y,z,dx,dy,dz,t1,t2,vfrac,vol
      real             :: t1i,t2i,vEcho,odx,ody,odz,volume
      double precision :: h,rho,ddx,ddy,xx,yy
      double precision :: dfac=9.4605284e17,tfac=365.25
      integer :: i,j
      integer :: n=25 ! number of points in grid cell over which to integrate
      
! N.B. work in yrs & lt-yrs, they are easier...

      x=dble(grid%xAxis(xP))/dfac
      y=dble(grid%yAxis(yP))/dfac
      z=dble(grid%zAxis(zP))/dfac

      t1=dble(t1i)/tfac
      t2=dble(t2i)/tfac

      if ( (xP>1) .and. (xP<grid%nx) ) then
         odx = abs(grid%xAxis(xP+1)-grid%xAxis(xP-1))/2.
      else if ( xP==1 ) then
         if (lgSymmetricXYZ) then
            odx = abs(grid%xAxis(xP+1)-grid%xAxis(xP))/2.
         else 
            odx = abs(grid%xAxis(xP+1)-grid%xAxis(xP))
         end if
      else if ( xP==grid%nx ) then
         odx = abs(grid%xAxis(xP)  -grid%xAxis(xP-1))
      end if
      
      if ( (yP>1) .and. (yP<grid%ny) ) then
         ody = abs(grid%yAxis(yP+1)-grid%yAxis(yP-1))/2.
      else if ( yP==1 ) then
         if (lgSymmetricXYZ) then
            ody = abs(grid%yAxis(yP+1)-grid%yAxis(yP))/2.
         else
            ody = abs(grid%yAxis(yP+1)-grid%yAxis(yP))
         end if
      else if ( yP==grid%ny ) then
         ody = abs(grid%yAxis(yP)  -grid%yAxis(yP-1))
      end if
      
      if ( (zP>1) .and. (zP<grid%nz) ) then    
         odz = abs(grid%zAxis(zP+1)-grid%zAxis(zP-1))/2.    
      else if ( zP==1 ) then    
         if (lgSymmetricXYZ) then
            odz = abs(grid%zAxis(zP+1)-grid%zAxis(zP))/2.
         else
            odz = abs(grid%zAxis(zP+1)-grid%zAxis(zP))
         end if
      else if ( zP==grid%nz ) then    
         odz = abs(grid%zAxis(zP)-grid%zAxis(zP-1))
      end if

      dx=dble(odx)/dfac
      dy=dble(ody)/dfac
      dz=dble(odz)/dfac
      
      vEcho=0.
      vol=0.0

! integrate all echo space inside torus

      ddx=dx/dble(n)
      ddy=dy/dble(n)
!      write(6,'(a2,8f9.3)')'s',x,y,z,dx,dy,dz,ddx,ddy
      do i=1,n
         do j=1,n
            xx=x-dx/2.+ddx/2.+dble(i-1)*ddx ! test at center of each subcell
            yy=y-dy/2.+ddy/2.+dble(j-1)*ddy
            rho=sqrt(xx**2+yy**2)
            h=max((min(z+dz,zee(rho,t1))-max(z,zee(rho,t2))),0.)
            vol=vol+h
!            if (vol.gt.0.) write(6,'(2f8.4,x,4f8.4,x,f8.4,1pg12.4)')&
!                 &xx,yy,zee(rho,t2),z,z+dz,zee(rho,t1),h,vol
         enddo
      enddo
      volume=real(vol*ddx*ddy)
      vEcho=real(vol*ddx*ddy/(dx*dy*dz)) ! ratio of volume in echo to grid cell
!      if (vEcho.gt.0.) write(6,'(a2,4f9.3)')'a',x,y,z,vEcho

      return
    end function vEcho
!
    function zee(rho,t)
      double precision :: zee,rho,t
      double precision :: cl=2.998e10

      zee=rho**2/(2.*t)-t/2.
      return 
    end function zee


end module grid_mod




























