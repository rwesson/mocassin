! Copyright (C) 2005 Barbara Ercolano 
!
! Version 2.02
module composition_mod
    use constants_mod
    use common_mod   

    contains

    subroutine setComposition(grid)
        implicit none

        type(grid_type) :: grid      ! grid

        ! local variables
        integer :: icount            ! counter
        integer :: nelem             ! counter
        integer :: ios               ! I/O error status     
        integer :: i,j,k             ! counters
        integer :: off               ! off everywhere? 

        ! open file containing composition data


        do i = 1, nAbComponents 

           close(80)        
           open(unit=80, action="read", file = abundanceFile(i), status="old", position="rewind", &
&                                                     iostat = ios)
           if (ios /= 0) then
              print*, "! setComposition: can't open file: ", abundanceFile(i)
              stop
           end if
           do j = 1, nElements
              read(80, *) grid%elemAbun(i, j)
           end do
           close(80)



        end do

        ! set the switches

        ! first initialise to true
        lgElementOn = .true.

        ! initialise nElementsUsed
        nElementsUsed = nElements
        icount = 1         
        
        ! check what elements can be switched off
        do nelem = 1, nElements
           off = 1
           do i = 1, nAbComponents
              if (grid%elemAbun(i,nelem) > 1.e-12) then 
                 off = 0
              end if
           end do
           
           
           if (off==0) then
              elementXref(nelem) = icount
              icount = icount+1
           else
              lgElementOn(nelem) = .false.
              nElementsUsed = nElementsUsed - 1
              elementXref(nelem) = 0            
           end if
           
        end do
        
    end subroutine setComposition
      
end module composition_mod

