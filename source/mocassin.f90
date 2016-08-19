! Copyright (C) 2005 Barbara Ercolano 
! Department of Physics and Astronomy
! University College London
! Gower Street
! London, WC1E 6BT
! UK
! +44-207-679-7146
! be@star.ucl.ac.uk
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version. This requires
! that any chnages or improvements made to the program should also be 
! made freely available.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

!  
! MoCaSSiN = MOnte CArlo SImulationS of Nebulae 
! Version 2.02
! this is the main driver of the simulation
!  
program MoCaSSiN
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

    type(grid_type) :: grid3D(maxGrids)       ! the 3D Cartesian  grid

    integer         :: iGrid, err             ! allocation error status

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, taskid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, numtasks, ierr)

    if (taskid == 0) then
        print*, "MOCASSIN 2005 Version 2.0"
        print*, " "
    end if 

    ! read the input parameters of the simulation
    call readInput()      

    if (taskid == 0) then
        print*, " "
    end if

    ! initialize the 3D cartesian grid
    
    do iGrid = 1, nGrids
       call initCartesianGrid(grid3D(iGrid), nxIn(iGrid), nyIn(iGrid), nzIn(iGrid))
    end do
    
    ! initialize opacities x sections array
    call initXSecArray()

    ! set the ionzing continuum according to the contShape variable
    call setContinuum()

    call fillGrid(grid3D(1:nGrids))

    if (taskid==0) then
       do iGrid = 1, nGrids
          print*, 'Grid : ', iGrid 
          print*, 'active cells: ', grid3D(iGrid)%nCells
       end do
    end if

    ! prepare atomica data stuff
    if (lgGas) call makeElements()
    
    print*, 'active elements: ', nElementsUsed

    ! if grains are included, calculate the dust opacity     
    if (lgDust) then
       do iGrid = 1, nGrids
          call dustDriver(grid3D(iGrid))
       end do
    end if

    call mpi_barrier(mpi_comm_world, ierr)

    ! start the Monte Carlo simulation
    call MCIterationDriver(grid3D(1:nGrids))

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

    call mpi_finalize(ierr)
    stop '! MoCaSSin: end simulation reached - clean exit -'

end program MoCaSSiN
   
    
