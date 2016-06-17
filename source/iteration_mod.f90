! Copyright (C) 2005 Barbara Ercolano 
!
! Version 2.02
module iteration_mod
    use common_mod
    use continuum_mod
    use emission_mod
    use ionization_mod
    use output_mod
    use photon_mod   
    use update_mod

    implicit none

    contains

    subroutine MCIterationDriver(grid)
       implicit none

       type(grid_type), intent(inout) :: grid(*) ! the 3D cartesian grid

       ! local variables
 
       integer :: i, j, k                     ! counters
       
       logical, save :: lgFirst = .true.      ! first time the procedure is called? 

       call iterateMC()

       contains

       recursive subroutine iterateMC()
           implicit none

           include 'mpif.h'

           ! local variables
           real, pointer :: budgetTemp(:,:)      ! temporary dust heating budget array
           real, pointer :: dustPDFTemp(:,:)     ! temporary dust emission PDF array
           real, pointer :: escapedPacketsTemp(:,:,:)!temporary escaped packets array
           real, pointer :: fEscapeResPhotonsTemp(:,:) ! temporary escape res line phot 
           real, pointer :: JDifTemp(:,:)        ! temporary diffuse field array
           real, pointer :: JSteTemp(:,:)        ! temporary stellar field array
           real, pointer :: linePacketsTemp(:,:) ! temporary line packets array
           real, pointer :: opacityTemp(:,:)     ! temporary opacities array
           real, pointer :: recPDFTemp(:,:)      ! temporary rec prob distribution function
           real, pointer :: linePDFTemp(:,:)     ! temporary line prob distribution function
           real, pointer :: totalLinesTemp(:)    ! temporary fraction of non-ionizing line phots
   
           real, pointer :: absTauTemp(:), &     ! temp absorption opacity
                & lambdaTemp(:)                  ! temp displacement
           real, pointer :: absTau(:), &         ! absorption opacity
                & lambda(:)                      ! displacement

            
           real, pointer          :: noHitPercent(:)    ! percentage of no Hit cells
           real, pointer          :: noIonBalPercent(:) ! percentage of cell where ion bal not conv
           real, pointer          :: noTeBalPercent(:)  ! percentage of cell where Te bal not conv
           real,save              :: totPercentOld   ! percentage of converged cells from prev iteration
           real                   :: photRatio       ! photons/cells ratio
           real                   :: radius          ! ditance in cm from star
           real                   :: tau             ! optical depth
           real                   :: totCells        ! total # of active cells 
           real                   :: totheatdust     ! total dust heating
           real          :: echod1, echod2            ! light echo distances  
           
           integer, pointer       :: planeIonDistributionTemp(:,:) 
           integer, pointer       :: resLinePacketsTemp(:) ! temporary array for extra packets
           integer                :: err             ! allocation error status
           integer                :: elem,freq,ion,nS! counters
           integer                :: ifreq, ian      ! counters
           integer                :: ios,iG          ! I/O error status           
           integer                :: load,rest       ! 
           integer                :: size            ! size for mpi
           integer                :: iCell           ! cell index including non-active
           integer                :: iStar           ! star index
           integer                :: ai              ! grain size counter
           integer                :: icontrib,icomp ! counters
           integer                :: imu             ! direction cosine            
           integer                :: cellLoc(3)      ! local cell counters
           integer                :: gpLoc           ! local grid counter
           integer                :: ii,jj,kk        ! counters
           integer                :: ngridloc
           integer                :: yTop

           integer                :: dcp, nsp

           allocate(noHitPercent(nGrids))
           allocate(noIonBalPercent(nGrids))
           allocate(noTeBalPercent(nGrids))           

           ! re-initialize MC estimators
           do iG = 1, nGrids
              grid(iG)%lgConverged(0:grid(iG)%nCells)    = 0
              grid(iG)%lgBlack(0:grid(iG)%nCells)        = 0
              if (lgGas) then
                 ! zero out PDF arrays
                 grid(iG)%recPDF(0:grid(iG)%nCells, 1:nbins) = 0.
                 if (lgDebug) grid(iG)%linePDF(0:grid(iG)%nCells, 1:nLines)    = 0.
                 grid(iG)%totalLines(0:grid(iG)%nCells) = 0.  
                 
                 ! zero out Balmer jump
                 BjumpTemp = 0.
                 Bjump     = 0.
              end if
           
              if (lgDust .and. .not.lgGas) then
                 ! zero out dust PDF arrays
                 grid(iG)%dustPDF(0:grid(iG)%nCells, 1:nbins)     = 0.
              end if
           end do

!*****************************************************************************

           iCell = 0

           do iG = 1, nGrids

              if (ig>1 .or. (.not.lg2D)) then
                 yTop = grid(iG)%ny
              else if (iG ==1 .and. lg2D) then
                 yTop = 1
              end if


              grid(iG)%opacity(0:grid(iG)%nCells, 1:nbins) = 0.

              if (lgGas) then
                 if (taskid==0) print*, '! iterateMC: ionizationDriver in', iG
                 ! calculate the opacities at every grid cell
                 do i = 1, grid(iG)%nx
                    do j = 1, yTop
                       do k = 1, grid(iG)%nz
                          iCell = iCell+1
                          if (mod(iCell-(taskid+1),numtasks)==0) &
                               & call ionizationDriver(grid(iG),i,j,k)
                       end do
                    end do
                 end do
                 if (taskid==0) print*, '! iterateMC: ionizationDriver out', iG

                 allocate(opacityTemp(0:grid(iG)%nCells, nbins), stat = err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array memory: opacityTemp ", iG
                    stop
                 end if
                 
                 opacityTemp = 0.
                 
                 size = (grid(iG)%nCells+1)*nbins
                 
                 call mpi_allreduce(grid(iG)%opacity, opacityTemp, size, &
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)
                 
                 
                 do i = 1, grid(iG)%nx
                    do j = 1, yTop
                       do k = 1, grid(iG)%nz
                          do freq = 1, nbins
                             
                             if (grid(iG)%active(i,j,k)>0) then
                                grid(iG)%opacity(grid(iG)%active(i,j,k),freq) = opacityTemp(grid(iG)%active(i,j,k),freq)
                             end if

                          end do
                       end do
                    end do
                 end do
                 if ( associated(opacityTemp) ) deallocate(opacityTemp)

              end if
           
              ! add dust contribution to total opacity
              if (taskid==0) print*, '! iterateMC: adding dust contribution to total opacity ',iG           
              if (lgDust) then
                 grid(iG)%scaOpac(0:grid(iG)%nCells, 1:nbins) = 0.
                 grid(iG)%absOpac(0:grid(iG)%nCells, 1:nbins) = 0.

                 if((nIterateMC==1 .and. lgEquivalentTau)) then
                    grid(iG)%scaOpac(0:grid(iG)%nCells, 1:nbins) = 0.
                    grid(iG)%absOpac(0:grid(iG)%nCells, 1:nbins) = 0.
                 else


                    if (ig>1 .or. (.not. lg2D)) then
                       yTop = grid(iG)%ny
                    else if (iG ==1 .and. lg2D) then
                       yTop = 1
                    end if

                    do i = 1, grid(iG)%nx
                       do j = 1, yTop
                          do k = 1, grid(iG)%nz
                             if (grid(iG)%active(i,j,k)>0) then
                                if (lgMultiDustChemistry) then 
                                   dcp = dustComPoint(grid(iG)%dustAbunIndex(grid(iG)%active(i,j,k)))
                                   nsp = grid(iG)%dustAbunIndex(grid(iG)%active(i,j,k))
                                else
                                   dcp = dustComPoint(1)
                                   nsp = 1
                                end if

                                do nS = 1, nSpeciesPart(nsp)

                                   do ai = 1, nSizes                                
                                      if (grid(iG)%Tdust(nS,ai,grid(iG)%active(i,j,k))<TdustSublime(dcp-1+nS)) then
                                         do freq = 1, nbins 
                                            grid(iG)%scaOpac(grid(iG)%active(i,j,k),freq) = &
                                                 & grid(iG)%scaOpac(grid(iG)%active(i,j,k),freq) + & 
                                                 & grainAbun(nsp,nS)*&
                                                 & grainWeight(ai)*grid(iG)%Ndust(grid(iG)%active(i,j,k))*&
                                                 & xSecArray(dustScaXsecP(nS+dcp-1,ai)+freq-1)
                                            grid(iG)%absOpac(grid(iG)%active(i,j,k),freq) = &
                                                 & grid(iG)%absOpac(grid(iG)%active(i,j,k),freq) + & 
                                                 & grainAbun(nsp,nS)&
                                                 & *grainWeight(ai)*grid(iG)%Ndust(grid(iG)%active(i,j,k))*&
                                                 & xSecArray(dustAbsXsecP(nS+dcp-1,ai)+freq-1)
                                         end do
                                      end if
                                   end do
                                end do
                                
                                do freq = 1, nbins
                                   grid(iG)%opacity(grid(iG)%active(i,j,k),freq) = &
                                        &grid(iG)%opacity(grid(iG)%active(i,j,k),freq) + &
                                        & (grid(iG)%scaOpac(grid(iG)%active(i,j,k),freq) + &
                                        &grid(iG)%absOpac(grid(iG)%active(i,j,k),freq))
                                end do
                             end if
                             
                          end do
                       end do
                    end do
                 end if
                 
              end if
              if (taskid==0) print*, '! iterateMC: dust contribution to total opacity added ',iG                            
              
           end do ! ngrids

           call mpi_barrier(mpi_comm_world, ierr)


!*****************************************************************************

           

           
           if (taskid==0 .and. lgWritePss) then
              open(unit=89, file='output/qHeatPss.out',  action="write",status='unknown', position='rewind', iostat=ios)
              if (ios /= 0) then
                 print*, "! iterationMC: can't open file for writing, output/qHeatPss.out"
                 stop
              end if
           end if
              
           ! check if min convergence was reached to carry out resonant line transfer
           if (lgDust .and. lgGas .and. convPercent>=resLinesTransfer .and. lgResLinesFirst &
                &.and. (.not.nIterateMC==1)) then

              call initResLines(grid(1:nGrids))
              
              do iG = 1, nGrids
                 allocate(fEscapeResPhotonsTemp(0:grid(iG)%nCells, 1:nResLines), stat &
                      &= err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array&
                         & memory:fEscapeResPhotonsTemp"
                    stop
                 end if
                 fEscapeResPhotonsTemp = 0.

                 allocate(resLinePacketsTemp(0:grid(iG)%nCells), stat &
                      &= err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array&
                         & memory:resLinePacketsTemp"
                    stop
                 end if
                 resLinePacketsTemp = 0.

              end do
           end if           

           ! set the diffuse PDFs at every grid cell
           do iG = 1, nGrids

              if (ig>1 .or. (.not.lg2D) ) then
                 yTop = grid(iG)%ny
              else if (iG ==1 .and. lg2D) then
                 yTop = 1
              end if


              if (taskid==0) print*, '! iterateMC: emissionDriver in',iG           
              iCell = 0
              do i = 1, grid(iG)%nx
                 do j = 1, yTop
                    do k = 1, grid(iG)%nz
                       iCell = iCell+1
                       if (mod(iCell-(taskid+1),numtasks)==0) &
                            & call emissionDriver(grid,i,j,k,iG)

                    end do
                 end do
              end do

              if (taskid==0) print*, '! iterateMC: emissionDriver out', iG           
              if (taskid==0 .and. lgWritePss) close(89)

              if (lgDust .and. .not.lgGas) then 
                 allocate(dustPDFTemp(0:grid(iG)%nCells, 1:nbins),&
                      & stat = err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array memory:&
                         & dustPDFTemp ", iG
                    stop
                 end if
                 dustPDFTemp = 0.
              end if

              if (lgGas) then
                 call mpi_allreduce(BjumpTemp, Bjump, 1, mpi_real&
                      &, mpi_sum, mpi_comm_world, ierr)
                 
                 ! Balmer jump is in [erg/s/A] since factor of e-40 from gas emissivity
                 ! and factor of e45 from volume calculations, so multiply by e5
                 if (lgSymmetricXYZ) Bjump = Bjump*8.
                 Bjump = Bjump*1.e5
                 
                 print*, "Balmer Jump: [erg/s/A] ", Bjump
                 
                 allocate(recPDFTemp(0:grid(iG)%nCells, 1:nbins),&
                      & stat = err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array memory:&
                         & recPDFTemp "
                    stop
                 end if

                 allocate(totalLinesTemp(0:grid(iG)%nCells), stat =&
                      & err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array memory:&
                         & opacityTemp "
                    stop
                 end if
                 
                 recPDFTemp = 0.
                 totalLinesTemp = 0.
                 
              end if

              if (lgDebug) then
                 allocate(linePDFTemp(0:grid(iG)%nCells, nLines),&
                      & stat = err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array memory:&
                         & opacityTemp "
                    stop
                 end if
                 linePDFTemp = 0.
              end if

              size =  (grid(iG)%nCells+1)*nbins
              

              
              if (lgGas) then
                 call mpi_allreduce(grid(iG)%recPDF, recPDFTemp, size, mpi_real&
                      &, mpi_sum, mpi_comm_world, ierr)
              
                 if (lgDebug) then
                    size =  (grid(iG)%nCells+1)*nLines

                    call mpi_allreduce(grid(iG)%linePDF, linePDFTemp, size,&
                         & mpi_real, mpi_sum, mpi_comm_world, ierr)
                 end if
              end if

              if (lgDust .and. .not.lgGas) then

                 call mpi_allreduce(grid(iG)%dustPDF, dustPDFTemp, size, mpi_real&
                      &, mpi_sum, mpi_comm_world, ierr)

                 do i = 0, grid(iG)%nCells
                    do freq = 1, nbins
                       grid(iG)%dustPDF(i,freq) = dustPDFTemp(i,freq)
!if (nIterateMC>1 .and. ig == 181 .and. i==9) then
!print*,i, grid(iG)%dustPDF(i,freq)
!end if
                    end do
                 end do

                 if (associated(dustPDFTemp)) deallocate(dustPDFTemp)
              end if

              call mpi_barrier(mpi_comm_world, ierr)

              size = (grid(iG)%nCells+1)

              if (lgGas) then
                 call mpi_allreduce(grid(iG)%totalLines, totalLinesTemp, size,&
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)
                 
                 do i = 1, grid(iG)%nx
                    do j = 1, yTop
                       do k = 1, grid(iG)%nz
                          if (grid(iG)%active(i,j,k)>0) then
                             grid(iG)%totalLines(grid(iG)%active(i,j,k)) = totalLinesTemp(grid(iG)%active(i,j,k))
                             do freq = 1, nbins
                                grid(iG)%recPDF(grid(iG)%active(i,j,k),freq) = recPDFTemp(grid(iG)%active(i,j,k) ,freq)
                             end do
                             if (lgDebug) then
                                do freq = 1, nLines
                                   grid(iG)%linePDF(grid(iG)%active(i,j,k),freq) = linePDFTemp(grid(iG)%active(i,j,k),freq)
                                end do
                             end if
                          end if
                       end do
                    end do
                 end do

                 call mpi_barrier(mpi_comm_world, ierr)
                 
                 if ( associated(totalLinesTemp) )&
                      & deallocate(totalLinesTemp)
                 if (lgDebug) then
                    if ( associated(linePDFTemp) ) deallocate(linePDFTemp)
                 end if
                 if ( associated(recPDFTemp) ) deallocate(recPDFTemp) 
              end if

              ! check if min convergence was reached to carry out resonant line transfer
              if (lgDust .and. lgGas .and. convPercent>=resLinesTransfer .and. lgResLinesFirst &
                   &.and. (.not.nIterateMC==1)) then

                 if (iG==nGrids) lgResLinesFirst = .false.
                 
                 size = (grid(iG)%nCells+1)*nResLines

                 call mpi_allreduce(grid(iG)%fEscapeResPhotons, fEscapeResPhotonsTemp, size, &
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)
                 
                 grid(iG)%fEscapeResPhotons(0:grid(iG)%nCells, 1:nResLines) = fEscapeResPhotonsTemp
                 
                 size = grid(iG)%nCells+1

                 if (associated(fEscapeResPhotonsTemp)) deallocate(fEscapeResPhotonsTemp)

                 call mpi_allreduce(grid(iG)%resLinePackets, resLinePacketsTemp, size, &
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)

                 grid(iG)%resLinePackets(0:grid(iG)%nCells) = resLinePacketsTemp

                 if (associated(resLinePacketsTemp)) deallocate(resLinePacketsTemp)
                 
                 call mpi_barrier(mpi_comm_world, ierr)

              end if


           end do ! nGrids loop

!**********************************************************************************************

           do iG = 1, nGrids

              grid(iG)%Jste(0:grid(iG)%nCells, 1:nbins)    = 0.
              if (lgDebug) then
                 grid(iG)%Jdif(0:grid(iG)%nCells, 1:nbins) = 0.
                 grid(iG)%linePackets(0:grid(iG)%nCells, 1:nLines) = 0.
              end if
           end do
          

           totalEscaped = 0.

           do iG = 1, nGrids
              grid(iG)%escapedPackets(0:grid(iG)%nCells, 0:nbins,0:nAngleBins) = 0.
           end do

           do iStar = 1, nStars
              if(taskid==0) print*, 'iterateMC: Starting transfer for ionising source ', iStar
              
              load = int(nPhotons(iStar)/numtasks)
              rest = mod(nPhotons(iStar), numtasks)           
           

              if (lgPlaneIonization) then
                 planeIonDistribution = 0
              end if


              ! send the photons through and evaluate the MC 
              ! estimators of Jste and Jdif at every grid cell
              if (taskid < rest) then
                 load = load+1
                 call energyPacketDriver(iStar,load, grid(1:nGrids))
              else
                 call energyPacketDriver(iStar,load, grid(1:nGrids))
              end if
           
              call mpi_barrier(mpi_comm_world, ierr)
           end do

           if (Ldiffuse>0.) then
              
              if (emittingGrid>0) then                 
                 ngridloc = emittingGrid
              else
                 ngridloc = ngrids
              end if

              do gpLoc = 1, nGridloc

                 if (gploc>1 .or. (.not. lg2D)) then
                    yTop = grid(gploc)%ny
                 else if (gploc ==1 .and. lg2D) then
                    yTop = 1
                 end if


                 if(taskid==0) print*, 'iterateMC: Starting transfer for diffuse source grid: ', gpLoc
                 do ii = 1,grid(gpLoc)%nx
                    do jj = 1,yTop
                       do kk = 1,grid(gpLoc)%nz
                 
                          if (grid(gpLoc)%active(ii,jj,kk)>0) then
                          
                             cellLoc(1)  = ii
                             cellLoc(2)  = jj
                             cellLoc(3)  = kk

                             load = int(nPhotonsDiffuseLoc/numtasks)
                             rest = mod(nPhotonsDiffuseLoc, numtasks)           

              
                             ! send the photons through and evaluate the MC 
                             ! estimators of Jste and Jdif at every grid cell
                             if (taskid < rest) then
                                load = load+1
                                call energyPacketDriver(iStar=0,n=load, grid=grid(1:nGrids), &
                                     & gpLoc=gpLoc, cellLoc=cellLoc)
                             else
                                call energyPacketDriver(iStar=0,n=load, grid=grid(1:nGrids), &
                                     & gpLoc=gpLoc, cellLoc=cellLoc)
                             end if
           
                             call mpi_barrier(mpi_comm_world, ierr)

                          end if

                       end do
                    end do
                 end do
              end do

           end if


           if (lgPlaneIonization) then

              allocate(planeIonDistributionTemp(grid(1)%nx, grid(1)%nz), stat = err)
              if (err /= 0) then
                 print*, "! iterateMC: can't allocate grid memory, planeIonDistributionTemp"
                 stop
              end if
              planeIonDistributionTemp = 0

              size = grid(1)%nx*grid(1)%nz

              call mpi_allreduce(planeIonDistribution, planeIonDistributionTemp, size, &
                   & mpi_integer, mpi_sum, mpi_comm_world, ierr)

              planeIonDistribution = planeIonDistributionTemp

              if (taskid ==0) then
                 open(file="output/planeIonDistribution.out",  action="write",unit=18, status="unknown")
                 do i = 1, grid(1)%nx
                    do k = 1, grid(1)%nz
                       write(18,*) i,k,planeIonDistribution(i,k)
                    end do
                 end do                 
                 close(18)
              end if

              if (associated(planeIonDistributionTemp)) deallocate(planeIonDistributionTemp)

           end if

           do iG = 1, nGrids

              if (ig>1 .or. (.not. lg2D)) then
                 yTop = grid(iG)%ny
              else if (iG ==1 .and. lg2D) then
                 yTop = 1
              end if

  

              allocate(escapedPacketsTemp(0:grid(iG)%nCells, 0:nbins, 0:nAngleBins), stat = err)
              if (err /= 0) then
                 print*, "! iterateMC: can't allocate grid memory: escapedPacketsTemp", iG
                 stop
              end if
           
              escapedPacketsTemp  = 0.

              allocate(JSteTemp(0:grid(iG)%nCells, nbins), stat = err)
              if (err /= 0) then
                 print*, "! iterateMC: can't allocate array memory: JsteTemp ", iG
                 stop
              end if
              JSteTemp           = 0.

              if (lgDebug) then
                 allocate(JDifTemp(0:grid(iG)%nCells, nbins), stat = err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array memory: JdifTemp "
                    stop
                 end if
                 allocate(linePacketsTemp(0:grid(iG)%nCells, nLines), stat = err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array memory: linePacketsTemp "
                    stop
                 end if
              
                 JDifTemp           = 0.              
                 linePacketsTemp    = 0.
              end if
          

              size =  (grid(iG)%nCells+1)*(1+nbins)*(nAngleBins+1)

              call mpi_allreduce(grid(iG)%escapedPackets, escapedPacketsTemp, size, &
                   & mpi_real, mpi_sum, mpi_comm_world, ierr)

              do i = 0, grid(iG)%nCells
                 do freq = 0, nbins                       
                    do imu = 0, nAngleBins                       
                       grid(iG)%escapedPackets(i, freq,imu) = escapedPacketsTemp(i, freq, imu)

                    end do
                 end do
              end do
                 

              if ( associated(escapedPacketsTemp) ) deallocate(escapedPacketsTemp)
              
!              if (taskid==0) call writeSED(grid)
!              if (taskid==0 .and. contCube(1)>0. .and. contCube(2)>0. ) &
!                   & call writeContCube(grid, contCube(1),contCube(2))    

              size =  (grid(iG)%nCells+1)*nbins

              if (lgDebug) then
                 call mpi_allreduce(grid(iG)%JDif, JDifTemp, size, &
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)
              end if

              call mpi_allreduce(grid(iG)%JSte, JSteTemp, size, &
                   & mpi_real, mpi_sum, mpi_comm_world, ierr)

              size =  (grid(iG)%nCells+1)*nLines

              if (lgDebug) then
                 call mpi_allreduce(grid(iG)%linePackets, linePacketsTemp, size, &
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)
              end if


              do i = 1, grid(iG)%nx
                 do j = 1, yTop
                    do k = 1, grid(iG)%nz
                       if (grid(iG)%active(i,j,k)>0) then
                          do freq = 1, nbins
                             if (lgDebug) then
                                grid(iG)%JDif(grid(iG)%active(i,j,k),freq) = &
                                     & JDifTemp(grid(iG)%active(i,j,k),freq)
                             end if
                             grid(iG)%JSte(grid(iG)%active(i,j,k),freq) = &
                                  & JSteTemp(grid(iG)%active(i,j,k),freq)
                             if (lg2D) grid(iG)%JSte(grid(iG)%active(i,j,k),freq) = &
                                  & grid(iG)%JSte(grid(iG)%active(i,j,k),freq)/&
                                  & TwoDscaleJ(grid(iG)%active(i,j,k))
                             
                          end do
                          if (lgDebug) then
                             do freq = 1, nLines
                                grid(iG)%linePackets(grid(iG)%active(i,j,k),freq) &
                                     & = linePacketsTemp(grid(iG)%active(i,j,k),freq)
                                if (lg2D) grid(iG)%JDif(grid(iG)%active(i,j,k),freq) = &
                                     & grid(iG)%JDif(grid(iG)%active(i,j,k),freq)/&
                                     & TwoDscaleJ(grid(iG)%active(i,j,k))
 
                             end do
                          end if
                       end if
                    end do
                 end do
              end do

              call mpi_barrier(mpi_comm_world, ierr)


              if (lgDebug) then
                 if ( associated(linePacketsTemp) )    deallocate(linePacketsTemp)
                 if ( associated(JDifTemp) )           deallocate(JDifTemp)
              end if
              if ( associated(JSteTemp) )           deallocate(JSteTemp)           


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
                    totalEscaped = totalEscaped/8.
                    if (lgDebug) grid(iG)%Jdif(i,:) = grid(iG)%Jdif(i,:)/8.
                 end if
                 
              end do

           end do


           if (lgDust .and. (nIterateMC>1 .or. .not.lgEquivalentTau)) then
              print*, "! iterateMC: [Interactions] : total -- abs -- sca: "
              print*, "! iterateMC: [Interactions] ", absInt+scaInt, " -- ", &
                   &  absInt*100./(absInt+scaInt),"% -- ", &
                   &  scaInt*100./(scaInt+absInt),"%"
           end if
           
           print*, " total Escaped Packets :",  totalEscaped              

           if (taskid==0) call writeSED(grid)                          
           if (taskid==0 .and. contCube(1)>0. .and. contCube(2)>0. ) & 
                & call writeContCube(grid, contCube(1),contCube(2))    

!*******************************************************************

           do iG = 1, nGrids
              
              if (ig>1 .or. (.not.lg2D)) then
                 yTop = grid(iG)%ny
              else if (iG ==1 .and. lg2D) then
                 yTop = 1
              end if


              allocate(lgConvergedTemp(0:grid(iG)%nCells), stat &
                   &= err)
              if (err /= 0) then
                 print*, "! iterateMC: can't allocate array&
                      & memory:lgConvergedTemp"
                 stop
              end if

              allocate(lgBlackTemp(0:grid(iG)%nCells), stat &
                   &= err)
              if (err /= 0) then
                 print*, "! iterateMC: can't allocate array&
                      & memory:lgBlackTemp  "
                 stop
              end if

              if (lgGas) then

                 allocate(NeTemp(0:grid(iG)%nCells), stat = err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array memory:&
                         & NeTemp ", iG
                    stop
                 end if
                 NeTemp           = 0.
                 allocate(TeTemp(0:grid(iG)%nCells), stat = err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array memory:&
                         & TeTemp ", iG
                    stop
                 end if
                 TeTemp           = 0.
                 allocate(ionDenTemp(0:grid(iG)%nCells, &
                      & nElementsUsed, nStages), stat = err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array memory: &
                         &ionDenTemp ", iG
                    stop
                 end if
                 ionDenTemp       = 0.
              
              end if
             
              if (lgDust) then
                 allocate(TdustTemp(0:nSpeciesMax,0:nSizes,0:grid(iG)%nCells), stat = err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array memory:&
                         &TdustTemp ", iG
                    stop
                 end if
                 TdustTemp = 0.
              end if


              grid(iG)%lgConverged(0:grid(iG)%nCells) = 0
              grid(iG)%lgBlack(0:grid(iG)%nCells)     = 0
              grid(iG)%noHit       = 0.
              grid(iG)%noIonBal    = 0.
              grid(iG)%noTeBal     = 0.
              noHitPercent(iG)     = 0.
              noIonBalPercent(iG)  = 0.
              noTeBalPercent(iG)   = 0.           
              lgConvergedTemp      = 0
              lgBlackTemp          = 0 

               if (lgTraceHeating.and.taskid==0) then
                 open(file="output/thermalBalance.out",  action="write",unit=57, status="unknown", iostat=ios)
                 if (ios /= 0) then
                    print*, "! iterationMC: can't open file for writing, output/thermalBalance.out"
                    stop
                 end if
              end if

              if(taskid==0) print*, 'iterateMC: updateCell in', iG
              iCell = 0
              do i = 1, grid(iG)%nx
                 do j = 1, yTop
                    do k = 1, grid(iG)%nz
                       iCell = iCell+1
                       if (mod(iCell-(taskid+1),numtasks)==0) &
                            & call updateCell(grid(iG),i,j,k)

                    end do
                 end do
              end do
              if(taskid==0) print*, 'iterateMC: updateCell out', iG

              if (lgTraceHeating.and.taskid==0) then
                 close (57)
              end if


              
              size = 1

              call mpi_allreduce(grid(iG)%noHit, noHitPercent(iG), size, &
                   & mpi_real, mpi_sum,mpi_comm_world, ierr)           

              if (lgGas) then
                 call mpi_allreduce(grid(iG)%noIonBal, noIonBalPercent(iG), size, &
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)           

                 call mpi_allreduce(grid(iG)%noTeBal, noTeBalPercent(iG), size, &
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)                                 
              end if

              size =  (grid(iG)%nCells+1)

              call mpi_allreduce(grid(iG)%lgConverged, lgConvergedTemp, size, &
                   & mpi_integer, mpi_sum, mpi_comm_world, ierr)
              call mpi_allreduce(grid(iG)%lgBlack, lgBlackTemp, size, &
                   & mpi_integer, mpi_sum, mpi_comm_world, ierr)


              if (lgGas) then
                 call mpi_allreduce(NeTemp, grid(iG)%Ne, size, &
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)
                 call mpi_allreduce(TeTemp, grid(iG)%Te, size, &
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)
           
           
                 size = (grid(iG)%nCells+1)*nElementsUsed*nStages
                 
                 call mpi_allreduce(ionDenTemp, grid(iG)%ionDen, size, &
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)

              end if

              if (lgDust) then
                 size =  (grid(iG)%nCells+1)*(nspeciesMax+1)*(nsizes+1)

                 call mpi_allreduce(TdustTemp,grid(iG)%Tdust,size, &
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)

                 if (lgGas .and. convPercent>=resLinesTransfer .and. (.not.lgResLinesFirst) &
                      & .and. .not.(nIterateMC==1)) then

                    allocate(budgetTemp(0:nAbComponents,0:nResLines), stat=err)
                    if (err /= 0) then
                       print*, "! iterateMC: can't allocate array memory:&
                            &budgetTemp "
                       stop
                    end if
                    budgetTemp=0.

                    size = (nAbcomponents+1)*(nResLines+1)
                    
                    call mpi_allreduce(dustHeatingBudget,budgetTemp,size, &
                         & mpi_real, mpi_sum, mpi_comm_world, ierr)

                    dustHeatingBudget = budgetTemp


                    deallocate(budgetTemp)
                    
                 end if

              end if

              do i = 1, grid(iG)%nx
                 do j = 1, yTop
                    do k = 1, grid(iG)%nz
                       if (grid(iG)%active(i,j,k)>0) then
                          grid(iG)%lgConverged(grid(iG)%active(i,j,k)) = &
                               & lgConvergedTemp(grid(iG)%active(i,j,k))           
                          grid(iG)%lgBlack(grid(iG)%active(i,j,k)) = &
                               & lgBlackTemp(grid(iG)%active(i,j,k))
                       end if
                    end do
                 end do
              end do


              call mpi_barrier(mpi_comm_world, ierr)
           
              if ( associated(lgConvergedTemp) )  deallocate(lgConvergedTemp)
              if ( associated(lgBlackTemp) )  deallocate(lgBlackTemp)
              if (lgGas) then
                 if ( associated(NeTemp) )           deallocate(NeTemp)
                 if ( associated(TeTemp) )           deallocate(TeTemp)
                 if ( associated(ionDenTemp) )       deallocate(ionDenTemp)
              end if
              if (lgDust) then
                 if ( associated(TdustTemp)) deallocate(TdustTemp)
              end if
              
           end do
           
!******************************************************************************

           if (nGrids>1) then
!              print*, " ! iterateMC: integratePathTau stuff still not implemented for multiple grids.... skipping"
           else
              call writeTauNu(grid)
           end if


           ! decide over final convergence of the model
           
           ! reinitialize convPercent and totCells           
           totPercent  = 0.

           ! calculate the percentage of converged cells
           do iG = 1, nGrids

              totCells    = 0.
              convPercent = 0.

              if (ig>1 .or. (.not.lg2D)) then
                 yTop = grid(iG)%ny
              else if (iG ==1 .and. lg2D) then
                 yTop = 1
              end if

              do i = 1, grid(iG)%nx
                 do j = 1, yTop
                    do k = 1, grid(iG)%nz
                       if (grid(iG)%active(i,j,k)>0)  then
                          if (.not.lgEcho) then 
                             convPercent = convPercent + grid(iG)%lgConverged(grid(iG)%active(i,j,k))
                             totCells    = totCells + 1.
                          else ! if light echo then only count echo cells!
                             if (grid(iG)%echoVol(i,j,k).gt.0.0) then
                                convPercent = convPercent + grid(iG)%lgConverged(grid(iG)%active(i,j,k))
                                totCells    = totCells + 1.
                             end if
                          endif
                       end if
                    end do
                 end do
              end do

              convPercent               = 100.*convPercent / totCells
              noHitPercent              = 100.*noHitPercent(iG) / totCells
              grid(iG)%noIonBal         = 100.*noIonBalPercent(iG) / totCells
              grid(iG)%noTeBal          = 100.*noTeBalPercent(iG) / totCells 
              totPercent = totPercent + convPercent*grid(iG)%nCells/100.

              if (taskid == 0) then
                 if (nIterateMC == 1) then
                    close(21)
                    open(unit=21, status='unknown', position='rewind', file='output/summary.out', iostat=ios)
                    if (ios /= 0) then
                       print*, "! iterationMC: can't open file for writing, summary.out -1"
                       stop
                    end if
                 else
                    close(21)
                    open(unit=21, status='unknown', position='append', file='output/summary.out', iostat=ios)
                    if (ios /= 0) then
                       print*, "! iterationMC: can't open file for writing, summary.out -2"
                       stop
                    end if
                 end if                
                 
                 print*, "! iterateMC:  Summary] Iteration ",nIterateMC,'; ', &
                      &int(convPercent),"% conveged cells in grid ", iG
                 print*, "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                      &int(noHitPercent(iG)),"% no hit cells in grid ", iG
                 print*, "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                      &int(grid(iG)%noIonBal),"% Ion balance not reached in grid ", iG
                 print*, "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                      &int(grid(iG)%noTeBal),"% Te balance not reached in grid", iG
                 write(21,*) "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                      &int(convPercent),"% conveged cells in grid ", iG
                 write(21,*) "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                      &int(noHitPercent(iG)),"% no hit cells in grid ", iG
                 write(21,*) "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                      &int(grid(iG)%noIonBal),"% Ion balance not reached in grid ", iG
                 write(21,*) "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                      &int(grid(iG)%noTeBal),"% Te balance not reached in grid", iG
                 if (lgGas .and. convPercent>=resLinesTransfer .and. .not.lgResLinesFirst .and. &
                      & (.not. nIterateMC==1) ) then
                    write(21,*) "! iterateMC: [Summary] Iteration ",nIterateMC,'; '
                    write(21,*) "Dust Budgets: "
                    do icomp = 0, nAbComponents
                       write(21,*) " Component ", icomp
                       totheatdust = 0.
                       do icontrib =0, nResLines
                          totheatdust = totheatdust+dustHeatingBudget(icomp,icontrib)
                       end do
                       do icontrib =0, nResLines
                          write(21,*) " Contribution", icontrib, dustHeatingBudget(icomp,icontrib)/&
                               & totheatdust
                       end do
                    end do
                 end if
              end if
              
           end do
           
           if (lgDust .and. convPercent>=resLinesTransfer .and. lgGas) dustHeatingBudget = 0.

           totCells = 0
           do iG =1,nGrids
              totCells = totCells+grid(iG)%nCells
           end do

           totPercent      = 100.*totPercent / totCells           

           if (taskid==0) then
              print*, "! iterateMC: [Summary] Iteration ",nIterateMC,'; Total:  ', &
                   &totPercent,"% converged cells over all grids"
              print*, "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                   &nPhotons, " energy packets used"
              write(21, *) "! iterateMC: [Summary] Iteration ",nIterateMC,'; Total:  ', &
                   &totPercent,"% converged cells over all grids"
              write(21,*) "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                   &nPhotons, " energy packets used"
              close(21)
           
              ! write grid to files for warm start
              
              if ( totPercent >= convWriteGrid) then

                 ! write old grid to disk
                 call writeGrid(grid(1:nGrids))

              end if
           
           end if
        
           nPhotonsTot = nPhotons(1)
           do iStar=1, nStars              
              nPhotonsTot = nPhotonsTot+nPhotons(iStar)
           end do

           if (nIterateMC > 1 .and. totPercent < 95. .and. lgAutoPackets & 
                & .and. nPhotonsTot < maxPhotons .and. totPercentOld > 0.) then

              if ( (totPercent-totPercentOld)/totPercentOld <= convIncPercent ) then
                 nPhotons = nPhotons*nPhotIncrease
                 deltaE   = deltaE/nPhotIncrease

                 if (taskid==0) &
                      & print*, "! iterateMC: [talk] Total number of energy packets &
                      &increased to ", nPhotons

                 if (Ldiffuse>0.) then
                    nPhotonsDiffuseLoc = nPhotonsDiffuseLoc*nPhotIncrease
                    if (taskid==0) &
                         & print*, "! iterateMC: [talk] number of diffuse energy packets &
                         &per cell increased to ", nPhotonsDiffuseLoc
                 end if              

              end if              
           end if


!           if (Ldiffuse>0. .and. nIterateMC > 1 .and. totPercent < 95. .and. lgAutoPackets & 
!                &  .and. totPercentOld > 0.) then
!
!              if ( (totPercent-totPercentOld)/totPercentOld <= convIncPercent ) then
!                 nPhotonsDiffuseLoc = nPhotonsDiffuseLoc*nPhotIncrease
!
!                 if (taskid==0) &
!                      & print*, "! iterateMC: [talk] number of diffuse energy packets &
!                      &per cell increased to ", nPhotonsDiffuseLoc
!              end if
!              
!           end if


           totPercentOld = totPercent

           if ( totPercent >= minConvergence ) then
              
              if (taskid==0) then
                 print*, "! iterateMC: [talk] convergence reached after ", &
                      &                           nIterateMC, " iterations. & Finishing up ... "
                 
                 ! output results at this itearation stage (every 3 iterations)

                 call writeGrid(grid(1:nGrids))
                 if (lgGas) call outputGas(grid(1:nGrids)) 
              end if
              
              call mpi_barrier(mpi_comm_world, ierr)

           else if (nIterateMC >= maxIterateMC) then

              if (taskid==0) then

                 print*, " ! iterateMC: maximum number of iterations reached. Finishing up ... ",&
                      & maxIterateMC

                 call writeGrid(grid(1:nGrids))
                 if (lgGas) call outputGas(grid(1:nGrids))
              end if

              call mpi_barrier(mpi_comm_world, ierr)
   
           else
              
              if (lgOutput .and. taskid == 0 ) then
                 ! output results at this itearation stage (every ? iterations)
                 if ( mod(nIterateMC, 1) == 0 ) then
                    if (lgGas) call outputGas(grid(1:nGrids)) 
                 end if
              end if
              
              call mpi_barrier(mpi_comm_world, ierr)
                 
              ! step up MC iterations counter
              nIterateMC = nIterateMC + 1
              
              if (lgTalk) print*, "! iterateMC: [talk] now starting iteration #",&
                   &                            nIterateMC

              ! start the next iteration

              call iterateMC()
              return
              
           end if

         end subroutine iterateMC

       end subroutine MCIterationDriver

 end module iteration_mod
