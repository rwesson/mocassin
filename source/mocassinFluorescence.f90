! Copyright (C) 2003 Barbara Ercolano
!  
! MoCaSSiNwarm = MOnte CArlo SImulationS of Nebulae 
! this is the warm-start driver of the simulation 
! (requires grid1.out,grid2.out,grid3.out files)
!  
! Version 2.02
program MoCaSSiNfluorescence
    use common_mod
    use constants_mod
    use dust_mod
    use grid_mod
    use iteration_mod
    use output_mod
    use set_input_mod
    use xSec_mod

    implicit none

    include 'mpif.h'

    type(grid_type) :: grid3D(maxGrids) ! the 3D Cartesian  grid

    integer         :: err              ! allocation error status
    integer         :: iGrid            ! 

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, taskid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, numtasks, ierr)

!    starttime = mpi_wtime()

    lgWarm = .true.
    lgFluorescence = .true.

    if (taskid == 0) then
        print*, "MOCASSIN 2006 fluorescence Version 3.00"
        print*, " "
    end if 

    if (taskid == 0) then
        print*, " "
    end if

    ! read input for fluorescence calculation file *N*
    call setFluorescenceInput()

    ! reset the 3D cartesian grid
    call resetGrid(grid3D)

    call setStarPosition(grid3D(1)%xAxis,grid3D(1)%yAxis,grid3D(1)%zAxis,grid3D)      

    ! initialize opacities x sections array
    call initXSecArray()
    ! set the ionzing continuum according to the contShape variable
    call setContinuum()

    do igrid = 1, ngrids
       allocate (grid3D(igrid)%fluorescenceCube(0:grid3D(igrid)%nCells, 1:nfluo, 0:nVP), stat=err) 
       if (err /= 0) then
          print*, '! setFluorescenceInput: allocation error for fluorescenceCube pointer'
          stop
       end if
       grid3D(igrid)%fluorescenceCube=0.
       allocate (grid3D(igrid)%Lfluorescence(0:grid3D(igrid)%nCells, 1:nfluo), stat=err) 
       if (err /= 0) then
          print*, '! setFluorescenceInput: allocation error for fluorescenceCube pointer'
          stop
       end if
       grid3D(igrid)%Lfluorescence = 0.

    end do

    ! prepare atomica data stuff
    call makeElements()

    if (taskid == 0) then
        print*, "Total number of grids (mother+sub-grids): ", nGrids
        print*, "Mother grid used:"
        print*, grid3D(1)%xAxis
        print*, grid3D(1)%yAxis
        print*, grid3D(1)%zAxis
    end if 

    ! if grains are included, calculate the dust opacity     
    if (lgDust) then
       do iGrid = 1, nGrids
          call dustDriver(grid3D(iGrid))
       end do
    end if

    call mpi_barrier(mpi_comm_world, ierr)
    
    ! start the Monte Carlo simulation
    call MCIterationDriver(grid3D)

    call mpi_barrier(mpi_comm_world, ierr)

    ! free all space allocated to the 3D grid
!    do iGrid=1, nGrids
!       call freeGrid(grid3D(iGrid))
!    end do

    call mpi_finalize(ierr)
    stop 'mpi done'
  
  end program MoCaSSiNfluorescence
   
    
