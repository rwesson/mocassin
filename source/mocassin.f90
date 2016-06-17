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

    real            :: test                   ! test
    integer         :: i, iGrid, err          ! allocation error status
    real            :: etime
    real, dimension(2) :: tarray,timing        ! cputimer 
    integer         :: nhours, nminutes, nseconds

    call cpu_time(timing(1))
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

    call setStarPosition(grid3D(1)%xAxis,grid3D(1)%yAxis,grid3D(1)%zAxis, grid3D(1:nGrids))

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
       if (taskid==0) print*, '! mocassin: calling dustDriver'
       do iGrid = 1, nGrids
          call dustDriver(grid3D(iGrid))
       end do
       if (taskid==0) print*, '! mocassin: dustDriver done'
    end if

    call mpi_barrier(mpi_comm_world, ierr)

    ! set local diffuse ionisation field
    if (Ldiffuse>0.) then
       if (taskid==0) print*, '! mocassin: calling setLdiffuse'
       call setLdiffuse(grid3D(1:nGrids))
       test=0.
       do iGrid = 1, nGrids
          do i = 1, grid3D(igrid)%ncells
             test = test+grid3d(igrid)%LdiffuseLoc(i)
          end do
       end do
       if (taskid==0) then
          print*, '! mocassin: setLdiffuse done, total Ldiffuse: ', test
       end if
    end if

    if (taskid==0) print*, '! mocassin: calling MCIterationDriver'    
    ! start the Monte Carlo simulation
    call MCIterationDriver(grid3D(1:nGrids))
    if (taskid==0) print*, '! mocassin: MCIterationDriver done'    

    call mpi_barrier(mpi_comm_world, ierr)

    ! free all space allocated to the 3D grid
    do iGrid=1, nGrids
       call freeGrid(grid3D(iGrid))
    end do

    call mpi_finalize(ierr)
    call cpu_time(timing(2))
    timing(1)=timing(1)/60.
    timing(2)=timing(2)/60.
    nhours = int((timing(2)-timing(1))/60.)
    nminutes = int(mod(timing(2)-timing(1),60.))
    nseconds = nint(60.*((timing(2)-timing(1))-real(nhours*60)-real(nminutes)))

    write(6,100)nhours,nminutes,nseconds
    100 format('total run time per processor ',1i3.2,':',1i2.2,':',1i2.2' (HMS)')
    stop '! MoCaSSin: end simulation reached - clean exit -'

end program MoCaSSiN
   
    
