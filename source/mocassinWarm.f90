! Copyright (C) 2003 Barbara Ercolano
!
! MoCaSSiNwarm = MOnte CArlo SImulationS of Nebulae
! this is the warm-start driver of the simulation
! (requires grid1.out,grid2.out,grid3.out files)
!
! Version 2.02
program MoCaSSiNwarm
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

    integer         :: iGrid            !

    real            :: etime
    real, dimension(2) :: tarray,timing        ! cputimer
    integer         :: nhours, nminutes, nseconds

    timing(1)=etime(tarray)/60.0


    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, taskid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, numtasks, ierr)

!    starttime = mpi_wtime()

    lgWarm = .true.

    if (taskid == 0) then
        print*, "MOCASSIN 2005 warm Version ",VERSION
        print*, " "
    end if

    if (taskid == 0) then
        print*, " "
    end if

    ! reset the 3D cartesian grid
    call resetGrid(grid3D)

    call setStarPosition(grid3D(1)%xAxis,grid3D(1)%yAxis,grid3D(1)%zAxis,grid3D)

    ! initialize opacities x sections array
    call initXSecArray()

    ! set the ionzing continuum according to the contShape variable
    call setContinuum()

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

    if (taskid ==  0 .and. .not.lgOutput) then
        ! determine final statistics
        if (lgGas) call outputGas(grid3D)
        call writeSED(grid3D)
        if (contCube(1)>0. .and. contCube(2)>0. ) &
             & call writeContCube(grid3D, contCube(1),contCube(2))
    end if

    call mpi_barrier(mpi_comm_world, ierr)

    ! free all space allocated to the 3D grid
    do iGrid=1, nGrids
       call freeGrid(grid3D(iGrid))
    end do

!    endtime  = mpi_wtime()
!
!    if (taskid == 0) then
!        print*, "time: ", endtime-starttime
!    end if

     timing(2) = etime(tarray)/60.0
     nhours = int((timing(2)-timing(1))/60.)
     nminutes = int(mod(timing(2)-timing(1),60.))
     nseconds = nint(60.*((timing(2)-timing(1))-real(nhours*60)-real(nminutes)))

     write(6,100)nhours,nminutes,nseconds
 100 format('total run time per processor ',1i3.2,':',1i2.2,':',1i2.2' (HMS)')


    call mpi_finalize(ierr)
    stop 'mpi done'


end program MoCaSSiNwarm


