! Copyright (C) 2005 Barbara Ercolano 
!
! Version 2.02
module pathIntegration_mod
    use common_mod
    use interpolation_mod
    use vector_mod

    contains

    ! This subroutine performs the integration of opacities from a position
    ! aVec to the edge of the nebula in the arbitrary direction uHat
    subroutine integratePathTau(grid, aVec, uHat, freq, nTau, absTau, lambda)
        implicit none

        type(grid_type), intent(in) :: grid(*)        ! grid

        type(vector),  intent(in)   :: uHat           ! direction vector
        type(vector),  intent(in)   :: aVec           ! starting position vector

        integer, intent(out)        :: nTau           ! actual size of opacity and length arrays 

        real, intent(in)            :: freq           ! the frequency [Hz]
        real, intent(out), &
             & dimension(maxTau)    :: absTau         ! absorption optical depth

        real, intent(out), &
             & dimension(maxTau)    :: lambda         ! distance array

        ! local variables

        type(vector)                :: rVec           ! position vector
        type(vector)                :: vHat           ! direction vector


        integer                     :: err            ! allocation error status
        integer                     :: freqP          ! frequency index
        integer                     :: i              ! counter
        integer                     :: xP, yP, zP     ! x, y, and z axis indeces
        integer, parameter          :: nl = 4         ! defines the size of the distance increment 


        real                        :: dlSmall        ! distance increment

        if (nGrids>1) then
           print*, 'integratePathTau: tau integration routine still not developed for multiple grids... returning '
           return
        end if

        ! locate the input frequency in the frequency array
        call locate(nuArray(1:nbins), freq, freqP)
        if (freqP<=0 .or. freqP>nbins) then
           print*, "! integratePathTau: out of bounds frequency ", freqP, freq, nuArray(1), nuArray(nbins)
           print*, nuArray
           stop
        end if


        ! locate the starting position in the cartesian grid
        call locate(grid(1)%xAxis, aVec%x, xP)
        if ( xP < grid(1)%nx ) then
            if ( aVec%x >= (grid(1)%xAxis(xP)+ &
                 & grid(1)%xAxis(xP+1))/2. ) xP = xP+1
        end if
        call locate(grid(1)%yAxis, aVec%y, yP)
        if ( yP < grid(1)%ny ) then
            if ( aVec%y >= (grid(1)%yAxis(yP)+ &                           
                 & grid(1)%yAxis(yP+1))/2. ) yP = yP+1 
        end if
        call locate(grid(1)%zAxis, aVec%z, zP)
        if ( zP < grid(1)%nz ) then
            if ( aVec%z >= (grid(1)%zAxis(zP)+ &                           
                 & grid(1)%zAxis(zP+1))/2. ) zP = zP+1 
        end if
 

        ! check that the input position is not outside the grid
        if ( (xP <= 0).or.(xP > grid(1)%nx) ) then
            print*, "! integratePathTau: starting position in x is outside the grid",xP,yP,zP
            stop
        else if ( (yP <= 0).or.(yP > grid(1)%ny) ) then
            print*, "! integratePathTau: starting position in y is outside the grid",xP,yP,zP
            stop
        else if ( (zP <= 0).or.(zP > grid(1)%nz) ) then   
            print*, "! integratePathTau: starting position in z is outside the grid",xP,yP,zP
            stop
        end if

        ! define dlSmall
        dlSmall = dl(1)

        ! initialize the optical depth arrays
        absTau = 0.
        lambda = 0.

        ! initialize nTau
        nTau = 1

        ! initialize direction vactor
        vHat = uHat        

        ! the first optical depth are all zero (as displacement is zero)
        ! these have already been set to zero in the initialization of the arrays
        ! so add the first displacements vector (dlSmall*vHat)
        rVec = aVec + dlSmall*vHat
        
        ! execute this loop until the edge of the grid is reached
        do i = 1, maxTau

           ! check if the path is still within the ionized region 
           if (rVec%x > grid(1)%xAxis(grid(1)%nx)) exit

           if (rVec%y > grid(1)%yAxis(grid(1)%ny)) exit
           
           if (rVec%z > grid(1)%zAxis(grid(1)%nz)) exit
           
           if ( sqrt( (rvec%x/1.e10)**2. + (rvec%y/1.e10)**2. + (rvec%z/1.e10)**2.)*1.e10 >= R_out &
                & .and. R_out > 0.) exit
           
           
           if (lgSymmetricXYZ) then
              if( rVec%x < grid(1)%xAxis(1) ) then
                 vHat%x = -vHat%x
                 rVec%x = -rVec%x
                 call locate(grid(1)%xAxis, rVec%x, xP)
              end if
              if( rVec%y < grid(1)%yAxis(1) ) then
                 vHat%y = -vHat%y
                 rVec%y = -rVec%y
                 call locate(grid(1)%yAxis, rVec%y, yP)
              end if
              if( rVec%z < grid(1)%zAxis(1) ) then
                 vHat%z = -vHat%z
                 rVec%z = -rVec%z
                 call locate(grid(1)%zAxis, rVec%z, zP)
              end if
           end if
             
           ! x-axis
           if (xP < grid(1)%nx) then
              if ( rVec%x > (grid(1)%xAxis(xP)+grid(1)%xAxis(xP+1))/2. ) then
                 xP = xP + 1
              end if
           else if (xP == grid(1)%nx) then
              if ( rVec%x > grid(1)%xAxis(xP) ) exit
           end if
           if (xP > 1) then
              if ( rVec%x < (grid(1)%xAxis(xP)+grid(1)%xAxis(xP-1))/2. ) then
                 xP = xP - 1
              end if
           end if

           ! y-axis
           if (yP < grid(1)%ny) then
              if ( rVec%y > (grid(1)%yAxis(yP)+grid(1)%yAxis(yP+1))/2. ) then
                 yP = yP + 1
              end if
           else if (yP == grid(1)%ny) then
              if (rVec%y > grid(1)%yAxis(yP)) exit
           end if
           if (yP > 1) then
              if ( rVec%y < (grid(1)%yAxis(yP)+grid(1)%yAxis(yP-1))/2. ) then
                 yP = yP - 1
              end if
           end if
           
           ! z-axis
            if (zP < grid(1)%nz) then
               if ( rVec%z > (grid(1)%zAxis(zP)+grid(1)%zAxis(zP+1))/2. ) then
                  zP = zP + 1
               end if
            else if (zP == grid(1)%nz) then
               if (rVec%z > grid(1)%zAxis(zP)) exit
            end if
            if (zP > 1) then
                if ( rVec%z < (grid(1)%zAxis(zP)+grid(1)%zAxis(zP-1))/2. ) then
                    zP = zP - 1
                end if
            end if
                
            if (.not.lgSymmetricXYZ) then
                if ((xP < 1) .or. (zP < 1) .or. (yP < 1)) exit
            end if

             ! check if there's any symmetries and if the path must be reflected
             if (lgSymmetricXYZ) then
                 if (xP<1) then
                     vHat%x=-vHat%x
                     rVec%x=-rVec%x
                     call locate(grid(1)%xAxis, rVec%x, xP)
                     if (xP < grid(1)%nx) then
                         if ( rVec%x > (grid(1)%xAxis(xP)+grid(1)%xAxis(xP+1))/2. ) xP = xP + 1
                     else if (xP == grid(1)%nx) then
                         if (rVec%x > grid(1)%xAxis(xP)) exit
                     end if
                 end if
                 if (yP<1) then
                     vHat%y=-vHat%y
                     rVec%y=-rVec%y
                     call locate(grid(1)%yAxis, rVec%y, yP)
                     if (yP < grid(1)%ny) then
                         if ( rVec%y > (grid(1)%yAxis(yP)+grid(1)%yAxis(yP+1))/2. ) yP = yP + 1
                     else if (yP == grid(1)%ny) then
                         if (rVec%y > grid(1)%yAxis(yP)) exit
                     end if 
                 end if
                 if (zP<1) then
                     vHat%z=-vHat%z
                     rVec%z=-rVec%z
                     call locate(grid(1)%zAxis, rVec%z, zP)
                     if (zP < grid(1)%nz) then
                         if ( rVec%z > (grid(1)%zAxis(zP)+grid(1)%zAxis(zP+1))/2. ) zP = zP + 1
                     else if (zP == grid(1)%nz) then
                         if (rVec%z > grid(1)%zAxis(zP)) exit
                     end if 
                 end if
             end if                                                                    

            ! increment nTau
            nTau = nTau + 1
 
            ! this should never happen (as the loop is stopped at maxTau anyway)
            ! just a sanity check
            if (nTau > maxTau) then 
                print*, "! integratePathTau: tau arrays are full", xP,yP,zP
                exit
            end if

            ! add the delta optical depths 
            absTau(nTau) = absTau(nTau-1) + grid(1)%opacity(grid(1)%active(xP,yP,zP),freqP)*dlSmall 
            lambda(nTau) = lambda(nTau-1) + dlSmall
            ! increment position by adding the displacements vector (dlSmall*vHat)
            rVec = rVec + dlSmall*vHat
            
        end do

    end subroutine integratePathTau

end module pathIntegration_mod       



