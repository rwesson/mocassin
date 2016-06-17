! Copyright (C) 2005 Barbara Ercolano 
!  
! MoCaSSiNoutput = MOnte CArlo SImulationS of Nebulae 
! this is the output-only driver of the simulation 
! (requires grid1.out,grid2.out,grid3.out files)
!  
! Version 2.02
program MoCaSSiNoutput
    use common_mod
    use constants_mod
    use grid_mod
    use iteration_mod
    use output_mod
    use set_input_mod
    use xSec_mod

    implicit none

    include 'mpif.h'

    type(grid_type) :: grid3D(maxGrids) ! the 3D Cartesian  grid

    integer         :: iGrid
    integer         :: err            ! allocation error status

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, taskid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, numtasks, ierr)

!    starttime = mpi_wtime()

    lgWarm = .true.
    lgRecombination = .true.

    if (taskid == 0) then
        print*, "MOCASSIN 2005 output Version 2.00"
        print*, "Creating output files from current grid*.out files "
        print*, " stored in the output/ directory"
        print*, " "
    end if 

    if (taskid == 0) then
        print*, " "
    end if

    ! reset the 3D cartesian grid
    call resetGrid(grid3D)

    ! initialize opacities x sections array
    call initXSecArray()

    ! set the ionzing continuum according to the contShape variable
    call setContinuum()

    ! prepare atomica data stuff
    call makeElements()
    
    if (taskid ==  0) then
       ! determine final statistics
       if (lgGas) then
          if (lg3DextinctionMap) then
             call outputGas(grid3D(1:nGrids), extMapFile) 
          else
             call outputGas(grid3D(1:nGrids)) 
          end if
       end if
    end if

    call mpi_barrier(mpi_comm_world, ierr)

    ! free all space allocated to the 3D grid
    do iGrid=1, nGrids
       call freeGrid(grid3D(iGrid))
    end do

!    endtime  = mpi_wtime()

    if (taskid == 0) then
!        print*, "time: ", endtime-starttime
    end if 

    call mpi_finalize(ierr)
    stop 'mpi done'


end program MoCaSSiNoutput
   
    
