! Copyright (C) 2005 Barbara Ercolano 
!  
! MoCaSSiN = MOnte CArlo SImulationS of Nebulae 
! this is the plotting driver
!  
! Version 2.02
program MoCaSSiNplot
    use common_mod
    use constants_mod
    use emission_mod
    use grid_mod
    use iteration_mod
    use output_mod
    use set_input_mod
    use xSec_mod

    implicit none

    include 'mpif.h'

    type(grid_type) :: grid3D(maxGrids)         ! the 3D Cartesian  grid

    type(plot_type) :: plot           ! the plot

    real(kind=8),pointer    :: flinePlot(:,:,:,:)
    real(kind=8),pointer    :: flinePlotLarge(:,:)
    real, pointer   :: image(:,:,:,:)    

    real, pointer   :: linePDFTemp(:,:)
    real, pointer   :: opacityTemp(:,:)  
    real, pointer   :: recPDFTemp(:,:)
    real, pointer   :: totalLinesTemp(:)
    
    ! filter transmission stuff
    real            :: frequency(5000), frequencyTemp(5000),tranCoeff(5000),tranCoeffTemp(5000)
    Real, pointer   :: coeff(:)

    real            :: dV                     ! volume of local cell [e45 cm^3]
    real            :: freq1,freq2
    real            :: log10TeLoc

    character(len=30) :: filename       ! input file
    character(len=30) :: bandFile       ! band transmission coeff file

    integer         :: abFileIndexUsed
    integer         :: err            ! allocation error statu
    integer         :: elem, ion      ! element ion counters
    integer         :: i,j,k,n,freq,iG! counters
    integer         :: iCount         ! cell index
    integer         :: iLine          ! line counter
    integer         :: iup, ilow, l   ! energy levels pointerst
    integer         :: ios            ! I/O error status
    integer         :: nxMax,nyMax,nzMax ! 
    integer         :: maxCells       ! maxCells
    integer         :: plotNum        ! counter
    integer         :: size,load,rest ! mpi stuff
    integer         :: tranP          ! transmnission coeff index
    integer         :: tranMaxP       ! maximum tran coeff index 
    integer         :: contP          ! continuum index counter 
    integer         :: yP             ! 

    logical         :: lgContinuum    ! 

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, taskid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, numtasks, ierr)


    if (taskid == 0) then
        print*, "MOCASSIN 2005 plot Version 2.00"
        print*, " "
    end if 

    filename = 'input/plot.in'


    ! reset the 3D cartesian grid
    call resetGrid(grid3D)

    call setStarPosition(grid3D(1)%xAxis,grid3D(1)%yAxis,grid3D(1)%zAxis,grid3D)      

    ! prepare atomica data stuff
    call makeElements()

    maxCells = grid3D(1)%nCells
    do iG = 1, ngrids
       if (grid3D(iG)%nCElls>maxCElls) maxCElls = grid3D(iG)%nCElls
    end do       

    print*, nstages
    allocate(flinePlot(nElements,nstages, nForLevels,nForLevels), stat=err)
    if (err /= 0) then
       print*, "! emissionDriver: can't allocate array memory"
       stop
    end if
    flinePlot = 0.
    if (lgElementOn(26) .and. nstages>=2) then
       allocate(flinePlotLarge(nForLevelsLarge,nForLevelsLarge), stat=err)
       if (err /= 0) then
          print*, "! emissionDriver: can't allocate array memory"
          stop
       end if
       flinePlotLarge = 0.
    end if

    ! read the file into the plot variable
    plot = readPlot(filename)

    lgContinuum = .false.
    do plotNum = 1, plot%nPlots
       if (.not.plot%lgLine(plotNum)) then
          print*, " ! mocassinPlot: no continuum plots are available in this version. &
               &Please contact Barbara Ercolano (be@star.ucl.ac.uk)"
          stop
       end if
    end do

    if (plot%lgFilter) then
       open(unit=72, file=bandFile,  action="read",status='old', position='rewind',iostat=ios)
       if (ios /= 0) then
          print*, "! readPlot: can't open band file for reading: ", bandFile
          stop
       end if
  
       frequency = 0.
       tranCoeff = 0.

       tranMaxP = 0
       do i = 1, 5000
           
          read(unit=72,fmt=*,iostat=ios) frequency(i), tranCoeff(i)

          if (ios < 0) exit

          ! the following assumes that wavelengths are in band file are in angstroms       
          frequency(i) = 910.998/frequency(i)

          tranMaxP = tranMaxP + 1

       end do

       do i=1, tranMaxP
          frequencyTemp(i) = frequency(1+tranMaxP-i)
          tranCoeffTemp(i) = tranCoeff(1+tranMaxP-i)
       end do
       
       frequency = frequencyTemp
       tranCoeff = tranCoeffTemp

    end if

    if (taskid == 0) then
       open(unit=29, file='output/grid4.out', status='unknown',  action="write",position='rewind',iostat=ios)
       if (ios /= 0) then
          print*, "! readPlot: can't open file for writing: output/grid4.out"
          stop
       end if
    end if

    ! read in data file for higher level H transitions
    call hdatx
    ! read in data from benj skil smit 99
    call readHeIRecLines()


    do iG = 1, nGrids


       if (lg2D) then 
          yP = 1
       else
          yP = grid3D(iG)%ny
       end if

       do i = 1, grid3D(iG)%nx
          do j = 1, yP
             do k = 1, grid3D(iG)%nz

                ! find the volume of this cell 
                dV = getVolume(grid3D(iG), i,j,k)

                if (taskid ==0) write(29, *) grid3D(iG)%xAxis(i),grid3D(iG)%yAxis(j),grid3D(iG)%zAxis(k),dV

                ! check this cell is in the ionized region
                if (grid3D(iG)%active(i,j,k) >0) then
                   
                   ! find the physical properties of this cell
                   HdenUsed        = grid3D(iG)%Hden(grid3D(iG)%active(i,j,k))
                   ionDenUsed      = grid3D(iG)%ionDen(grid3D(iG)%active(i,j,k), :, :)
                   NeUsed          = grid3D(iG)%Ne(grid3D(iG)%active(i,j,k))
                   TeUsed          = grid3D(iG)%Te(grid3D(iG)%active(i,j,k))
                   abFileIndexUsed = grid3D(iG)%abFileIndex(i,j,k)
                   log10TeLoc      = log10(TeUsed)
                   
                   ! recalculate line emissivities at this cell 
                   
                   ! calculate the emission due to HI and HeI rec lines
                   call RecLinesEmission()

                   ! calculate the emission due to the heavy elements
                   ! forbidden lines
                   call forLines()

                   ! find the volume of this cell 
                   dV = getVolume(grid3D(iG), i,j,k)

                   do plotNum = 1, plot%nPlots

                      if (plot%lgLine(plotNum)) then
                      
                         ! initialise line number counter
                         iLine = 1


                         ! HI rec lines
                         do ilow = 2,8
                            do iup = ilow+1, 30
                               
                               if (plot%lineNumber(plotNum) == iLine) &
                                    & plot%intensity(iG, grid3D(iG)%active(i,j,k),plotNum) = &
                                    & HIRecLines(iup,ilow)*HdenUsed*dV
                               iLine = iLine+1
                               
                            end do
                         end do
                      
                         do l = 1, 34 
                            if (plot%lineNumber(plotNum) == iLine) &
                                 & plot%intensity(iG,grid3D(iG)%active(i,j,k),plotNum) = & 
                                 & HeIRecLines(l)*HdenUsed*dV

                            iLine = iLine+1

                         end do
                         
                         
                         ! HeII rec lines
                         do iup = 3, 30
                            do ilow = 2, min0(16, iup-1)
                               
                               
                               if (plot%lineNumber(plotNum) == iLine) &
                                    & plot%intensity(iG, grid3D(iG)%active(i,j,k),plotNum) = &
                                    & HeIIRecLines(iup,ilow)*HdenUsed*dV

                               iLine = iLine+1

                            end do
                         end do
                         
                         ! Heavy elements forbidden lines
                         do elem = 3, nElements
                            do ion = 1, min(elem+1, nstages)

                               if (.not.lgElementOn(elem)) exit
                               if (lgDataAvailable(elem,ion)) then
                                  
                                  if (elem == 26 .and. ion ==2) then

                                     do iup = 1, nforlevelsLarge
                                        do ilow = 1, nforlevelsLarge
                                                                                   
                                           if (plot%lineNumber(plotNum) == iLine) then 
                                              
!                                              print*, flinePlotLarge(iup,ilow), elem,ion,iup,ilow
                                              plot%intensity(iG, grid3D(iG)%active(i,j,k),plotNum) = &
                                                   & flinePlotLarge(iup,ilow)*HdenUsed*dV
                                              
                                           end if
                                           
                                           iLine = iLine+1

                                        end do
                                     end do

                                  else

                                     do iup = 1, nforlevels
                                        do ilow = 1, nforlevels
                                        
                                           
                                           if (plot%lineNumber(plotNum) == iLine) then 
                                              
!                                              print*, flinePlot(elem,ion,iup,ilow), elem,ion,iup,ilow
                                              plot%intensity(iG, grid3D(iG)%active(i,j,k),plotNum) = &
                                                   & flinePlot(elem,ion,iup,ilow)*HdenUsed*dV
                                              
                                           end if
                                           
                                           iLine = iLine+1

                                        end do
                                     end do
                                     
                                  end if
                               end if
                            end do
                         end do
                         
                         ! calculate the intensity per unit frequency 
!                         plot%intensity(iG, grid3D(iG)%active(i,j,k),plotNum) = &
!                              & plot%intensity(iG,grid3D(iG)%active(i,j,k),plotNum) / &
!                              & 2.1799153e-11*nuArray(plot%nuP(plotNum,1))
                         
                      else
                         
                         print*, " ! mocassinPlot: no continuum plots are available in this version. &
                              &Please contact Barbara Ercolano (be@star.ucl.ac.uk)"
                         stop
                         
                         do contP = 1, nbins
                            
                            if ( (contP >= plot%nuP(plotNum,1)) .and. (contP <= plot%nuP(plotNum,2)) &
                                 &  .and. (continuum(contP) >=1.e-20) ) then
                               
                               if (plot%lgFilter) then 
                                  
                                  call locate(frequency(1:tranMaxP), nuArray(contP), tranP)
                                  if ( nuArray(contP) >= (frequency(tranP) + frequency(tranP+1))/2. ) &
                                       tranP = tranP+1
                                  if( (tranP <= 0) .or. (tranP >= tranMaxP) ) coeff = 0.
                                  
                                  coeff(plotNum) = tranCoeff(tranP)
                                  
                                  ! continuum is in units of [e-40 erg/cm^3/s/Hz]; multiply by 1.e-15 
                                  ! to bring to the same units as the lines, also cRyd at the bottom (2.1799153e-11)
                                  ! so just multipy by 4.587e-5
                                  plot%intensity(iG,grid3D(iG)%active(i,j,k),plotNum) = &
                                       & plot%intensity(iG,grid3D(iG)%active(i,j,k),plotNum) + &
                                    & continuum(contP)*4.587e-5*HdenUsed*dV*coeff(plotNum)/nuArray(contP)
                               else
                                  ! continuum is in units of [e-40 erg/cm^3/s/Hz]; multiply by 1.e-15 
                                  ! to bring to the same units as the lines, also cRyd at the bottom (2.1799153e-11)
                                  ! so just multipy by 4.587e-5
                                  plot%intensity(iG,grid3D(iG)%active(i,j,k),plotNum) = &
                                       & plot%intensity(iG,grid3D(iG)%active(i,j,k),plotNum) + &
                                       & continuum(contP)*4.587e-5*HdenUsed*dV/nuArray(contP)
                               end if
                               
                            end if
                         end do
                         
                      end if
                   end do ! plots

                end if ! inner/outer cell condition
                
             end do
          end do
       end do
    end do

    if (taskid == 0) close(29)
       
    ! calculate intensities in units of [E30 erg/s/Hz]
    ! Note:  intensities are now in E-25 * E45 ergs/sec (the E-25 comes from the emission
    ! module calculations and the E45 comes from the volume calculations to avoid
    ! overflow. Hence intensities are in [E20 ergs/s], so we need to multiply by E-16
    ! to give units of  [E36 erg/s].
    plot%intensity = plot%intensity*1.e-16
    
    ! per sterradian
!    plot%intensity = plot%intensity/fourPi
    
    if (plot%lgFilter) then

       do plotNum = 1, plot%nPlots
          
          if (plot%lgLine(plotNum)) then
             call locate(frequency(1:tranMaxP), nuArray(plot%nuP(plotNum,1)), tranP)
             if ( nuArray(plot%nuP(plotNum,1)) >= (frequency(tranP) + frequency(tranP+1))/2. ) &
                  tranP = tranP+1
             if( (tranP <= 0) .or. (tranP >= tranMaxP) ) then 
                print*, 'mocassinPlot: [warning] frequency outside the band', nuArray(plot%nuP(plotNum,1))
                coeff = 0.
             end if
             
             coeff(plotNum) = tranCoeff(tranP)
             
          else
             coeff(plotNum) = 1.
          end if
          
       end do
       
    else
          
       coeff = 1.
       
    end if
       
    open(unit=28, file='output/plot.out', status='unknown',  action="write",position='rewind',iostat=ios)
    if (ios /= 0) then
       print*, "! readPlot: can't open file for writing: output/plot.out"
       stop
    end if
    
    iCount = 0
       
    ! allocate pointers depending on nLines
    nxMax=grid3D(1)%nx
    nyMax=grid3D(1)%ny
    nzMax=grid3D(1)%nz
    do iG=1,nGrids
       if (grid3D(iG)%nx > nxMax) nxMax=grid3D(iG)%nx
       if (grid3D(iG)%ny > nxMax) nyMax=grid3D(iG)%ny
       if (grid3D(iG)%nz > nxMax) nzMax=grid3D(iG)%nz
    end do

    allocate(image(1:nGrids,1:nxMax, 1:nyMax, 1:nzMax), stat = err)
    if (err /= 0) then
       print*, "! mocassinPlot: can't allocate image memory"
       stop
    end if

    image = 0.

    do iG=1, nGrids
       do i = 1, grid3D(iG)%nx
          do j = 1, yP
             do k = 1, grid3D(iG)%nz
                
                if (plot%lgFilter) then
                   do plotNum = 1, plot%nPlots
                      image(iG,i,j,k) = image(iG,i,j,k) + plot%intensity(iG,grid3D(iG)%active(i,j,k),plotNum)*coeff(plotNum)
                   end do
                   
                   iCount = iCount + 1 
                   write(28,*) iCount,' ',(plot%intensity(iG,grid3D(iG)%active(i,j,k),plotNum),' ', &
                        & plotNum = 1, plot%nPlots), image(iG,i,j,k)

                else

                   iCount = iCount + 1 
                   write(28,*) iCount,' ',(plot%intensity(iG,grid3D(iG)%active(i,j,k),plotNum),' ', &
                        & plotNum = 1, plot%nPlots)
                   
                end if
             
             end do
          end do
       end do
    end do
!******************************************************************************
       
    ! free all space allocated to the 3D grid
    do iG=1, nGrids
       call freeGrid(grid3D(iG))
    end do
    
    ! free all space allocated to the plot
    call freePlot(plot)
    
    if (associated(ionDenUsed)) deallocate(ionDenUsed)
    if (associated(fLinePlot)) deallocate(fLinePlot)
    

    call mpi_finalize(ierr)
    stop 'mpi done'
    

    
    contains

      function readPlot(filename)
        implicit none

        type(plot_type)           :: readPlot    ! plot

        character(len=15), intent(in) &
&                                 :: filename    ! plot parameters file


        ! local variables
        real                      :: freq        ! cont freq

        character(len=10)         :: lineORcont  ! line or continuum ?
        character(len=10)         :: anOrMC      ! analytical or MC?
        character(len=10)         :: filter      ! filter?


        integer                   :: err         ! allocation error status
        integer                   :: i           ! counters
        integer                   :: ios         ! I/O error status
        integer                   :: maxLim      ! loop limit
        integer                   :: code        ! line code number (0 if cont)
        

        close(77)
        open(unit=77, file=filename, status='old', position='rewind', action="read",iostat=ios)
        if (ios /= 0) then
           print*, "! readPlot: can't open file for reading", filename
           stop
        end if

        maxLim = 10000

        ! find nPlots
        readPlot%nPlots = 0

        read(unit=77,fmt=*,iostat=ios) filter

        select case(filter)
        case('band')
           readPlot%lgFilter = .true.
           backspace(77)
           read(unit=77,fmt=*,iostat=ios) filter, bandFile
        case('mono')
           readPlot%lgFilter = .false.
           backspace(77)
        case default
           print*, '! readPlot: band or mono must be specified for plots'
           stop
        end select
                    
        do i = 1, maxLim
           
        read(unit=77,fmt=*,iostat=ios) lineOrcont

           if (ios < 0) exit

           readPlot%nPlots = readPlot%nPlots + 1
           
        end do

        readPlot%nPlots = readPlot%nPlots - 1

        print*, 'Number of plots: ', readPlot%nPlots

        close(77)

        ! allocate pointers
        allocate(readPlot%lineNumber(readPlot%nPlots), stat=err)
        if (err /= 0) then
           print*, "! readPlot: can't allocate memory for lineNumber array"
           stop
        end if
        allocate(readPlot%nuP(readPlot%nPlots,2), stat=err)
        if (err /= 0) then
           print*, "! readPlot: can't allocate memory for nuP array"
           stop
        end if
        allocate(readPlot%lgLine(readPlot%nPlots), stat=err)
        if (err /= 0) then
           print*, "! readPlot: can't allocate memory for lgLine array"
           stop
        end if
        allocate(readPlot%intensity(1:nGrids,0:maxCells,readPlot%nPlots), stat=err)
        if (err /= 0) then
           print*, "! readPlot: can't allocate memory for intensity array"
           stop
        end if
        allocate(coeff(readPlot%nPlots), stat=err)
        if (err /= 0) then
           print*, "! readPlot: can't allocate memory for transmission coeff array"
           stop
        end if

        readPlot%lineNumber   = 0
        readPlot%nuP          = 0
        readPlot%lgLine       = .true.
        readPlot%intensity    = 0.
        freq                  = 0.
        code                  = 0

        open(unit=77, file=filename, status='old', position='rewind', action="read",iostat=ios)
        if (ios /= 0) then
           print*, "! readPlot: can't open file for reading", filename
           stop
        end if

        read(77,*) filter

        do i = 1, readPlot%nPlots

           read(77, *) lineORcont, code, freq1, freq2

           ! assume freq1 and freq2 are wavelengths in Angstroms
           freq1 = 910.998/freq1
           freq2 = 910.998/freq2

           ! decide whether line or continuum and assign to plot
           select case (lineORcont)
           case ('line')
              readPlot%lgLine(i)     = .true.
              readPlot%lineNumber(i) = code
              call locate(nuArray, freq1, readPlot%nuP(i,1))
              if (freq1 >= (nuArray(readPlot%nuP(i,1)) + nuArray(readPlot%nuP(i,1)+1))/2. ) &
                   readPlot%nuP(i,1) = readPlot%nuP(i,1)+1
              readPlot%nuP(i,2) = readPlot%nuP(i,1)
           case ('continuum')
              readPlot%lgLine(i)     = .false.
              readPlot%lineNumber(i) = code
              call locate(nuArray, freq2, readPlot%nuP(i,1))
              if (freq2 >= (nuArray(readPlot%nuP(i,1)) + nuArray(readPlot%nuP(i,1)+1))/2. ) &
                   readPlot%nuP(i,1) = readPlot%nuP(i,1)+1
              call locate(nuArray, freq1, readPlot%nuP(i,2))
              if (freq1 >= (nuArray(readPlot%nuP(i,2)) + nuArray(readPlot%nuP(i,2)+1))/2. ) &
                   readPlot%nuP(i,2) = readPlot%nuP(i,2)+1
           end select
        end do

        ! close file
        close(77)

        print*, 'out readPlot'

      end function readPlot

      subroutine freePlot(plot)
        implicit none

        type(plot_type), intent(inout) :: plot
        
        if (associated(plot%intensity))  deallocate(plot%intensity)
        if (associated(plot%lgLine))     deallocate(plot%lgLine)
        if (associated(plot%lineNumber)) deallocate(plot%lineNumber)
        if (associated(plot%nuP))        deallocate(plot%nuP)

      end subroutine freePlot

      

      ! this subroutine is the driver for the calculation of the emissivity
      ! from the heavy elements forbidden lines. 
      subroutine forLines()
        implicit none

        integer                     :: elem, ion ! counters
          
        ! re-initialize forbiddenLines
        fLinePlot = 0.
        
        do elem = 3, nElements
           do ion = 1, min(elem+1, nstages)
              if (.not.lgElementOn(elem)) exit
                
              if (lgDataAvailable(elem, ion)) then
                 
                 if (elem == 26 .and. ion ==2) then

                    if (ion<min(elem+1,nstages)) then
                       call equilibrium(dataFile(elem, ion), &
                            &ionDenUsed(elementXref(elem), ion+1)/&
                            &ionDenUsed(elementXref(elem),ion), &
                            & TeUsed, NeUsed, flinePlotLarge(:,:))
                    else
                       call equilibrium(dataFile(elem, ion), 0., &
                            & TeUsed, NeUsed, flinePlotLarge(:,:))
                    end if

                    flinePlotLarge(:, :) =flinePlotLarge(:,:)*&
                         &grid3D(iG)%elemAbun(abFileIndexUsed,elem)*&
                         & ionDenUsed(elementXref(elem), ion)
                 else

                    if (ion<min(elem+1,nstages)) then
                       call equilibrium(dataFile(elem, ion), &
                            &ionDenUsed(elementXref(elem), ion+1)/&
                            &ionDenUsed(elementXref(elem),ion), &
                            & TeUsed, NeUsed, flinePlot(elem,ion,:,:))
                    else
                       call equilibrium(dataFile(elem, ion), 0., &
                            & TeUsed, NeUsed, flinePlot(elem,ion,:,:))
                    end if

                    flinePlot(elem, ion, :, :) =flinePlot(elem,ion,:,:)*&
                         &grid3D(iG)%elemAbun(abFileIndexUsed,elem)*&
                         & ionDenUsed(elementXref(elem), ion)

                 end if
              end if
           end do
        end do


        ! scale the forbidden lines emissivity to give units of [10^-25 erg/s/Ngas]
        ! comment: the forbidden line emissivity is so far in units of cm^-1/s/Ngas
        !          the energy [erg] of unit wave number [cm^-1] is 1.9865e-16, hence 
        !          the right units are obtained by multiplying by 1.9865e-16*1e25 
        !          which is 1.9865*1e9 
        flinePlot = flinePlot*1.9865e9

        if (lgElementOn(26) .and. ion>=2) flinePlotLarge = flinePlotLarge*1.9865e9


      end subroutine forLines
    
    subroutine RecLinesEmission()
        implicit none

        ! local variables
        integer                    :: ios         ! I/O error status
        integer                    :: i           ! counters
        integer                    :: ilow,&      ! pointer to lower level
&                                      iup        ! pointer to upper level
        integer                    :: denint

        real                       :: A4471, A4922! HeI reference lines 
        real                       :: aFit        ! general interpolation coeff
        real                       :: C5876, C6678! collition exc. corrections
        real                       :: Hbeta       ! Hbeta emission
        real                       :: HeII4686    ! HeII 4686 emission
        real                       :: Lalpha      ! Lalpha emission
        real                       :: T4          ! TeUsed/10000.
        real                       :: x1,x2
        real                       :: hb          ! emissivity of H 4->2
        real                       :: fh          ! emissivity of H


        T4 = TeUsed / 10000.

        ! do hydrogenic ions first
        
        ! calculate Hbeta
        if (TeUsed > 26000.) then
           print*, "! recLineEmission: [warning] temperature exceeds 26000K - Hbeta &
                & calculations may be uncertain"
        end if
        Hbeta = 2530./(TeUsed**0.833) ! TeUsed < 26000K CASE 

        ! fits to Storey and Hummer MNRAS 272(1995)41
!        Hbeta = 10**(-0.870*log10Te + 3.57)
        Hbeta = Hbeta*NeUsed*ionDenUsed(elementXref(1),2)*grid3D(iG)%elemAbun(abFileIndexUsed,1)

        call hlinex(4,2,TeUsed,NeUsed,fh,2)
        hb = fh

        do ilow = 2, 8
           do iup = 30, ilow+1, -1
              call hlinex(iup,ilow,TeUsed,NeUsed,fh,2)              
              HIRecLines(iup, ilow) = (fh/hb)*Hbeta

           enddo
        enddo


        ! add contribution of Lyman alpha 
        ! fits to Storey and Hummer MNRAS 272(1995)41
        Lalpha = 10**(-0.897*log10TeLoc + 5.05) 

        HIRecLines(30, 8) =HIRecLines(30, 8) + grid3D(iG)%elemAbun(abFileIndexUsed,1)*&
             & ionDenUsed(elementXref(1),2)*&
             & NeUsed*Lalpha 
        

        ! read in HeII recombination lines [e-25 ergs*cm^3/s]
        ! (Storey and Hummer MNRAS 272(1995)41)
        close(95)
        open(unit = 95,  action="read", file = trim(home)//"data/r2b0100.dat", status = "old", position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! RecLinesEmission: can't open file: data/r2b0100.dat"
            stop
        end if
        do iup = 30, 3, -1
            read(95, fmt=*) (HeIIRecLines(iup, ilow), ilow = 2, min0(16, iup-1))
        end do

        close(95)

        ! calculate HeII 4686 [E-25 ergs*cm^3/s]
        HeII4686 = 10.**(-.997*log10(TeUsed)+5.16)
        HeII4686 = HeII4686*NeUsed*grid3D(iG)%elemAbun(abFileIndexUsed,2)*ionDenUsed(elementXref(2),3)

        ! calculate emission due to HeI recombination lines [e-25 ergs/s/cm^3]
        do iup = 30, 3, -1
            do ilow = 2, min0(16, iup-1)
                HeIIRecLines(iup, ilow) = HeIIRecLines(iup, ilow)*HeII4686
            end do
        end do    

        ! now do HeI
        
        if (NeUsed <= 100.) then
           denint=0
        elseif (NeUsed > 100. .and. NeUsed <= 1.e4) then
           denint=1
        elseif (NeUsed > 1.e4 .and. NeUsed < 1.e6) then
            denint=2
        elseif (NeUsed >= 1.e6) then
           denint=3
        end if

        ! data from Benjamin, Skillman and Smits ApJ514(1999)307 [e-25 ergs*cm^3/s]
        if (denint>0.and.denint<3) then
           do i =1, 34
              x1=HeIrecLineCoeff(i,denint,1)*(T4**(HeIrecLineCoeff(i,denint,2)))*exp(HeIrecLineCoeff(i,denint,3)/T4)

              x2=HeIrecLineCoeff(i,denint+1,1)*(T4**(HeIrecLineCoeff(i,denint+1,2)))*exp(HeIrecLineCoeff(i,denint+1,3)/T4)

              HeIRecLines(i) = x1+((x2-x1)*(NeUsed-100.**denint)/(100.**(denint+1)-100.**(denint)))


           end do
        elseif(denint==0) then
           do i = 1, 34
              HeIRecLines(i) = HeIrecLineCoeff(i,1,1)*(T4**(HeIrecLineCoeff(i,1,2)))*exp(HeIrecLineCoeff(i,1,3)/T4)
           end do
        elseif(denint==3) then
           do i = 1, 34
              HeIRecLines(i) = HeIrecLineCoeff(i,3,1)*(T4**(HeIrecLineCoeff(i,3,2)))*exp(HeIrecLineCoeff(i,3,3)/T4)
           end do
        end if
        HeIRecLines=HeIRecLines*NeUsed*grid3D(iG)%elemAbun(abFileIndexUsed,2)*ionDenUsed(elementXref(2),2)
    end subroutine RecLinesEmission

end program MoCaSSiNplot
   
    


