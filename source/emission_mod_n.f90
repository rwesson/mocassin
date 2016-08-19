! Copyright (C) 2003 Barbara Ercolano 
!
! Version 1.10
module emission_mod
    use constants_mod
    use common_mod
    use continuum_mod
    use grid_mod
    use xSec_mod

    type(resLine_type), pointer         :: resLine(:)  ! resonant lines

    real                                :: abFileUsed  ! abundance file index used
    real                                :: BjumpTemp   ! 
    real                                :: cellPUsed   !  cell index
    real                                :: HbetaProb   ! probability of Hbeta
    real                                :: log10Ne     ! log10(Ne) at this cell
    real                                :: log10Te     ! log10(Te) at this cell
    real                                :: sqrTeUsed   ! sqr(Te) at this cell

    ! units are [e-40erg/cm^3/s/Hz]
    double precision, pointer :: emissionHI(:)                     ! HI continuum emission coefficient 
    double precision, pointer :: emissionHeI(:)                    ! HeI continuum emission coefficient
    double precision, pointer :: emissionHeII(:)                   ! HeII continuum emission coefficient  

    ! units are [e-25erg/s/N_gas]
    double precision, dimension(nElements,nstages, nForLevels,nForLevels)::&
         & forbiddenLines                                          ! emissivity from heavies  rec lines
    double precision, dimension(3:15, 2:8) :: HIRecLines           ! emissivity from HI rec lines 
    double precision, dimension(9)         :: HeIRecLinesS         ! emissivity from HeI sing rec lines
    double precision, dimension(11)        :: HeIRecLinesT         ! emissivity from HeI trip rec lines
    double precision, dimension(3:30, 2:16):: HeIIRecLines         ! emissivity from HeII rec lines

    contains

    ! this is the subroutine to drive the calculation 
    ! of the emission spectrum at the given cell
    subroutine emissionDriver(grids, ix, iy, iz, ig)
        implicit none

        integer, intent(in) :: ix, iy, iz, ig  ! pointers to this cell
    

        type(grid_type), intent(inout) :: grids(*) ! the cartesian grid
        type(grid_type)                :: grid

        ! local variables

        integer :: err                       ! allocation error status
        integer :: i, n                      ! counters

        real    :: dV                        ! volume element        

        double precision, pointer :: gammaHI(:)          ! HI fb+ff emission coefficient [e-40 erg*cm^3/s/Hz]
        double precision, pointer :: gammaHeI(:)         ! HeI fb+ff emission coefficient [e-40erg*cm^3/s/Hz]
        double precision, pointer :: gammaHeII(:)        ! HeII fb+ff emission coefficient [e-40erg*cm^3/s/Hz]
        double precision, pointer :: gammaHeavies(:)     ! Heavy ions f-b emission coefficient [e-40 erg*cm^3/s/Hz]
        double precision, pointer :: ffCoeff1(:)         ! ff emission coefficient for Z=1 ions [e-40 erg*cm^3/s/Hz] 
        double precision, pointer :: ffCoeff2(:)         ! ff emission coefficient for Z=2 ions [e-40 erg*cm^3/s/Hz]
        double precision, pointer :: twoPhotHI(:)        ! HI 2photon emission coefficient [e-40 erg*cm^3/s/Hz]
        double precision, pointer :: twoPhotHeI(:)       ! HeI 2photon emission coefficient [e-40 erg*cm^3/s/Hz]  
        double precision, pointer :: twoPhotHeII(:)      ! HeII 2photon emission coefficient [e-40 erg*cm^3/s/Hz]


!print*, "in emissionDriver"

        grid = grids(iG)

        ! check whether this cell is outside the nebula
        if (grid%active(ix, iy, iz)<=0) return  


        ! set the dust emission PDF
        if (lgDust) call setDustPDF()

        if (.not.lgGas) return

        ! find the physical properties of this cell
        ionDenUsed= grid%ionDen(grid%active(ix, iy, iz), :, :)

        abFileUsed= grid%abFileIndex(ix,iy,iz)
        NeUsed    = grid%Ne(grid%active(ix, iy, iz))
        TeUsed    = grid%Te(grid%active(ix, iy, iz))
        sqrTeUsed = sqrt(TeUsed)

        log10Te = log10(TeUsed)
        log10Ne = log10(NeUsed)

        ! check the electron density is not too low
        if (NeUsed < 1.e-5) then
             NeUsed = 1.e-5
        end if

        ! allocate space for emissionHI, emissionHeI and emissionHeII
        allocate(emissionHI(nbins), stat = err)
        if (err /= 0) then
            print*, "! emissionDriver: can't allocate array memory"
            stop
        end if
        allocate(emissionHeI(nbins), stat = err)              
        if (err /= 0) then
            print*, "! emissionDriver: can't allocate array memory"
            stop
        end if
        allocate(emissionHeII(nbins), stat = err)
        if (err /= 0) then
            print*, "! emissionDriver: can't allocate array memory"
            stop
        end if
 
        ! allocate space for the f-b H, He and Heavy ions emission
        ! coefficients
        allocate(gammaHI(nbins), stat = err)
        if (err /= 0) then
            print*, "! emissionDriver: can't allocate array memory"
            stop
        end if
        allocate(gammaHeI(nbins), stat = err)
        if (err /= 0) then
            print*, "! emissionDriver: can't allocate array memory"
            stop
        end if
        allocate(gammaHeII(nbins), stat = err)
        if (err /= 0) then
            print*, "! emissionDriver: can't allocate array memory"
            stop
        end if
        allocate(gammaHeavies(nbins), stat = err)        
        if (err /= 0) then
            print*, "! emissionDriver: can't allocate array memory"
            stop
        end if
        allocate(twoPhotHI(nbins), stat = err)
        if (err /= 0) then
            print*, "! emissionDriver: can't allocate array memory"
            stop
        end if
        allocate(twoPhotHeI(nbins), stat = err)
        if (err /= 0) then
            print*, "! emissionDriver: can't allocate array memory"
            stop
        end if
        allocate(twoPhotHeII(nbins), stat = err)              
        if (err /= 0) then
            print*, "! emissionDriver: can't allocate array memory"
            stop
        end if
        allocate(ffCoeff1(nbins), stat = err)
        if (err /= 0) then
            print*, "! emissionDriver: can't allocate array memory"
            stop
        end if
        allocate(ffCoeff2(nbins), stat = err)
        if (err /= 0) then
            print*, "! emissionDriver: can't allocate array memory"
            stop
        end if


        ! zero out arrays
        continuum      = 0.
        emissionHI     = 0.
        emissionHeI    = 0.
        emissionHeII   = 0.
        gammaHI        = 0.
        gammaHeI       = 0.
        gammaHeII      = 0.
        gammaHeavies   = 0.
        twoPhotHI      = 0.
        twoPhotHeI     = 0.
        twoPhotHeII    = 0.
        ffCoeff1       = 0.
        ffCoeff2       = 0.

        ! calculates the H, He and heavy ions f-b and f-f emission
        ! coefficients
        call fb_ff()

        ! calculate two photon emission coefficients for HI, HeI and
        ! HeII
        call twoPhoton()

        ! calculate the total emission coefficients for H and He
        do i = 1, nbins
            emissionHI(i)  = (gammaHI(i) + twoPhotHI(i) ) * grid&
                 &%elemAbun(grid%abFileIndex(ix,iy,iz), 1) * ionDenUsed(elementXref(1),2)*NeUsed &
                 &+ gammaHeavies(i)*NeUsed
            emissionHeI(i) = (gammaHeI(i) + twoPhotHeI(i) ) * grid&
                 &%elemAbun(grid%abFileIndex(ix,iy,iz),2) * ionDenUsed(elementXref(2),2)*NeUsed 
            emissionHeII(i)= (gammaHeII(i) + twoPhotHeII(i) ) * grid&
                 &%elemAbun(grid%abFileIndex(ix,iy,iz),2) * ionDenUsed(elementXref(2),3)*NeUsed

            continuum(i)   = emissionHI(i) + emissionHeI(i) + &
                 & emissionHeII(i)

            ! zero out near zero contributions

            if (emissionHI(i) < 1.e-15 ) emissionHI(i) = 0.
            if (emissionHeI(i) < 1.e-15 ) emissionHeI(i) = 0.
            if (emissionHeII(i) < 1.e-15 ) emissionHeII(i) = 0.


        end do
        
        dV = getVolume(grid,ix,iy,iz)

        ! add contribution of this cell to the integrated balmer jump
        BjumpTemp = BjumpTemp + ((emissionHI(BjumpP)&
             &+emissionHeI(BjumpP)+emissionHeII(BjumpP))-&
             & (emissionHI(BjumpP-1)+emissionHeI(BjumpP-1)&
             &+emissionHeII(BjumpP-1)))* 2.2565e11*grid%Hden(grid%active(ix,iy,iz))*dV

        ! calculate the emission due to HI and HeI recombination lines
        call RecLinesEmission()

        ! calculate the emission due to the heavy elements forbidden lines
        call forLines()

        ! set the PDFs for diffuse radiation
        call setDiffusePDF()

        ! set the resonance line escape probabilities
        if (lgDust.and.convPercent>=resLinesTransfer .and. lgResLinesFirst &
             & .and. (.not.nIterateMC==1)) then
           call setResLineEscapeProb(grids,ig,ix,iy,iz)
        end if

        ! set the extra packets due to dust absorption of resonance lines
        if (lgDust.and.convPercent>=resLinesTransfer &
             & .and. (.not.nIterateMC==1)) then
           call initResLinePackets()        
        end if

        ! deallocate arrays
        if ( associated(emissionHI) ) deallocate(emissionHI)
        if ( associated(emissionHeI) ) deallocate(emissionHeI)
        if ( associated(emissionHeII) ) deallocate(emissionHeII)
        if ( associated(gammaHI) ) deallocate(gammaHI)
        if ( associated(gammaHeI) ) deallocate(gammaHeI)
        if ( associated(gammaHeII) ) deallocate(gammaHeII)
        if ( associated(gammaHeavies) ) deallocate(gammaHeavies)
        if ( associated(twoPhotHI) ) deallocate(twoPhotHI)
        if ( associated(twoPhotHeI) ) deallocate(twoPhotHeI)
        if ( associated(twoPhotHeII) ) deallocate(twoPhotHeII)
        if ( associated(ffCoeff1) ) deallocate(ffCoeff1)
        if ( associated(ffCoeff2) ) deallocate(ffCoeff2)

        grids(iG) = grid

    contains

      subroutine initResLinePackets()
        implicit none
        
        real                :: lineIntensity

        integer             :: iRes, imul
        
        ! calculate total energy deposited into grains at location icell in grid igrid
        lineIntensity=0.
        do iRes = 1, nResLines        
           do imul = 1, resLine(iRes)%nmul
              
              if (resLine(iRes)%elem==1) then
                 if ( resLine(iRes)%ion == 1 .and. resLine(iRes)%moclow(imul)==1 &
                      & .and. resLine(iRes)%mochigh(imul)==2 ) then

                    ! fits to Storey and Hummer MNRAS 272(1995)41
                    lineIntensity = lineIntensity+10**(-0.897*log10(TeUsed) + 5.05)* &
                         & grid%elemAbun(abFileUsed,1)*ionDenUsed(elementXref(1),2)*&
                         &  NeUsed*grid%Hden(grid%active(ix,iy,iz))*dV*&
                         & (1.-grid%fEscapeResPhotons(grid%active(ix,iy,iz), iRes))

                 else
                                      
                    print*, "! initResLinePackets: [warning] only dust heating from H Lyman alpha and resonance lines from heavy"
                    print*, "elements is implemented in this version. please contact author B. Ercolano -1-", ires
                                      
                 end if
              else if (resLine(iRes)%elem>2) then
                                   
                 lineIntensity = lineIntensity+&
                      &forbiddenLines(resLine(iRes)%elem,resLine(iRes)%ion,&
                      &resLine(iRes)%moclow(imul),resLine(iRes)%mochigh(imul))*&
                      & grid%Hden(grid%active(ix,iy,iz))*dV*&
                      & (1.-grid%fEscapeResPhotons(grid%active(ix,iy,iz), iRes))
                                   
              else
                                   
                 print*, "! initResLinePackets: [warning] only dust heating from H Lyman alpha and resonance lines from heavy"
                 print*, "elements is implemented in this version. please contact author B. Ercolano -2-", ires
                                   

              end if
              
           end do
        end do

        ! calculate number of dust emitted extra packets to be transfered from this location
        ! NOTE :1.e-16 from 1.e45 (from dV) * 1.e-36 (from Lstar) * 1.e-25 (from local forlines)
        grid%resLinePackets(grid%active(ix,iy,iz)) = nint(lineIntensity*1.e-16/(Lstar/Nphotons))

!print*, '! initResLinePackets: cell, extra, lineIntensity, Lstar, Nphotons'
!print*, ix, iy, iz, grid%active(ix,iy,iz), grid%resLinePackets(grid%active(ix,iy,iz)), lineIntensity, Lstar, Nphotons

      end subroutine initResLinePackets



    ! this subroutine calculates the H, He and heavy ions f-b and f-f
    ! emission coefficients
    ! - it can only be called by emission driver -
    subroutine fb_ff()
        implicit none

        ! local variables

        logical            :: lgTeOut         ! is Te out of range of
        ! available f-b coeffs

        integer            :: elem, ion       ! counters
        integer            :: g0,&            ! stat weight of (elem, nElec) ground state
&                              g1             ! stat weight of (elem, nElec-1) ground state
        integer            :: HIPnuP, &       ! pointer to the H IP in nuArray
&                              HeIPnuP,&      ! pointer to the HeI IP in nuArray
&                              HeIIPnuP,&     ! pointer to the HeII IP in nuArray
&                              highNuP,&      ! pointer to the hi en of the current ion in nuArray
&                              IPnuP          ! pointer to the IP of the current ion in nuArray
        integer            :: i, j            ! counters
        integer            :: ios             ! I/O error status
        integer            :: iup, ilow       ! counters
        integer, parameter :: nEdgesHI  = 11  ! number of series edges included in the HI data
        integer, parameter :: nEdgesHeI = 8   ! number of series edges included in the HeI data 
        integer, parameter :: nEdgesHeII= 11  ! number of series edges included in the HeII data 
        integer            :: nElec           ! number of electrons in the ion
        integer            :: nTemp           ! =1 if 5000K<Te<10000K, =2 if 10000K<Te<20000K
        integer            :: outShell        ! outer shell number (1 for k-shell)
        integer            :: xSecP           ! pointer to phot x sec of ion M in xSecArray

        integer, dimension(2*nEdgesHI)  :: HINuEdgeP   ! pointers to the HI, HeI and HeII 
        integer, dimension(2*nEdgesHeI) :: HeINuEdgeP  ! series edges in nuArray
        integer, dimension(2*nEdgesHeII):: HeIINuEdgeP !

        real, dimension(3)              :: statW      ! g0/g1 for HI, HeI and HeII

        double precision                :: aFit, &    ! general calculations factors
&                                          factor, &  ! 
&                                          expFactor  !
        double precision                :: constant   ! general calculation constant
        double precision                :: phXSecHI,& ! phot x sec for HI
&                                         phXSecHeI,& ! phot x sec for HeI
&                                         phXSecHeII,&! phot x sec for HeII 
&                                         phXSecM     ! phot x sec for ion M
        double precision, dimension(2*nEdgesHI,3)  :: gHI         ! gamma [e-40 erg/Hz] for HI, HeI and
        double precision, dimension(2*nEdgesHeI,3) :: gHeI        ! HeII at T=5,10,20 e3K
        double precision, dimension(2*nEdgesHeII,3):: gHeII       !
        double precision, dimension(2*nEdgesHI)    :: logGammaHI  ! log10(gamma)
        double precision, dimension(2*nEdgesHeI)   :: logGammaHeI ! log10(gamma)
        double precision, dimension(2*nEdgesHeII)  :: logGammaHeII! log10(gamma)
        real, dimension(2*nEdgesHI)    :: nuEdge      ! freq [Ryd] at the series edge (used
                                                      ! for HI, HeI and HeII)
     
        ! read in the f-b coefficients for T=5000, 10000 and 20000 K
        ! Ferland,PASP 92(1980)596 
        ! Osterbrock,tbl 4.8

        ! HI 11 series edges
        close(90)
        open(unit = 90, file = "data/gammaHI.dat", status = "old",&
             & position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! fb_ff: can't open file: data/gammaHI.dat"
            stop
        end if
        do i = 1, 2*nEdgesHI
            read(unit = 90, fmt = *) nuEdge(i), gHI(i,1),  gHI(i,2), &
                 & gHI(i,3)
            call locate(nuArray, nuEdge(i), HINuEdgeP(i))            
        end do
        HINuEdgeP(1) = 1       

        ! HeI 8 series edges
        close(91)
        open(unit = 91, file = "data/gammaHeI.dat", status = "old",&
             & position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! fb_ff: can't open file: data/gammaHeI.dat"
            stop
        end if
        do i = 1, 2*nEdgesHeI
            read(unit = 91, fmt = *) nuEdge(i), gHeI(i,1),  gHeI(i,2)&
                 &, gHeI(i,3)
            call locate(nuArray, nuEdge(i), HeINuEdgeP(i))
        end do
        HeINuEdgeP(1) = 1 

        ! HeII 11 series edges
        close(92)
        open(unit = 92, file = "data/gammaHeII.dat", status = "old",&
             & position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! fb_ff: can't open file: data/gammaHeII.dat"
            stop
        end if
        do i = 1, 2*nEdgesHeII
            read(unit = 92, fmt = *) nuEdge(i), gHeII(i,1),  gHeII(i&
                 &,2), gHeII(i,3)
            call locate(nuArray, nuEdge(i), HeIINuEdgeP(i))
        end do
        HeIINuEdgeP(1) = 1

        ! close files 
        close(90)
        close(91)
        close(92)

        ! interpolate the coefficients in Temperature (linear in log-log)
        
        ! check what coefficients are needed for the local electron Temperature
        lgTeOut = .false. ! (re)initialize out of range Te flag 
        if (TeUsed < 5000.) then
            lgTeOut = .true.
            do i = 1, 2*nEdgesHI
                logGammaHI(i) = log10(gHI(i, 1))
            end do
            do i = 1, 2*nEdgesHeI
                logGammaHeI(i) = log10(gHeI(i, 1))
            end do
            do i = 1, 2*nEdgesHeII
                logGammaHeII(i) = log10(gHeII(i, 1))
            end do   
        else if (TeUsed < 10000.) then
            nTemp = 1
        else if (TeUsed < 20000.) then
            nTemp = 2
        else ! Te > 20000.
            lgTeOut = .true.
            do i = 1, 2*nEdgesHI
                logGammaHI(i) = log10(gHI(i, 3))
            end do
            do i = 1, 2*nEdgesHeI
                logGammaHeI(i) = log10(gHeI(i, 3))
            end do
            do i = 1, 2*nEdgesHeII
                logGammaHeII(i) = log10(gHeII(i, 3))
            end do   
        end if

        ! carry out interpolation in temperature (if 5000K<Te<20000K)
        if (.not.lgTeOut) then    ! if Te is within the  range
            factor = log10(TeUsed/(float(nTemp)*5000.))
            do i = 1, 2*nEdgesHI            
                aFit = log10(gHI(i, nTemp+1)/gHI(i, nTemp))/ log10(2.)
                logGammaHI(i) = log10(gHI(i, nTemp)) + aFit * factor
            end do
            do i = 1, 2*nEdgesHeI
                aFit = log10(gHeI(i, nTemp+1)/gHeI(i, nTemp))/&
                     & log10(2.)
                logGammaHeI(i) = log10(gHeI(i, nTemp)) + aFit * factor
            end do
            do i = 1, 2*nEdgesHeII
                ! in the following log10(a)-log10(b) has to be used
               ! instead of 
                ! log10(a/b) to avoid division by 0 caused by
               ! underflow occuring for
                ! gHeII(22,1)=8.08e-39. This problem does not occur
               ! on a supercomputer
                aFit = ( log10(gHeII(i, nTemp+1))-log10(gHeII(i,&
                     & nTemp)) )/ log10(2.)
                logGammaHeII(i) = log10(gHeII(i, nTemp)) + aFit *&
                     & factor
            end do
        end if

        ! calculate f-b and f-f emissivity of HI, HeI and HeII 

        ! interpolate between the edges in nu, for the HI, HeI and
        ! HeII series 

        ! calculate gammaHI (nu<1.Ryd)
        do i = 1, nEdgesHI
            iup  = HINuEdgeP(2*i)
            ilow = HINuEdgeP(2*i - 1)
            if (iup == ilow) then
                print*, "! fb_ff: frequency grid too course - divide&
                     & by zero                          [nu]",&
                     & nuArray(iup)
                stop
            end if
            aFit = (logGammaHI(2*i)-logGammaHI(2*i - 1)) /&
                 & (nuArray(iup) - nuArray(ilow))
            ! carry out interpolation between edges
            do j = ilow, iup-1
                gammaHI(j) = logGammaHI(2*i)+aFit* (nuArray(j)&
                     &-nuArray(iup))
                gammaHI(j) = 10.**gammaHI(j)
            end do
        end do

        ! calculate gammaHeI (nu<1.8Ryd)
        do i = 1, nEdgesHeI
            iup  = HeINuEdgeP(2*i)
            ilow = HeINuEdgeP(2*i - 1)
            if (iup == ilow) then       
                print*, "! fb_ff: frequency grid too course - divide&
                     & by zero                          [nu]",&
                     & nuArray(iup) 
                stop
            end if 
            aFit = (logGammaHeI(2*i)-logGammaHeI(2*i - 1)) /&
                 & (nuArray(iup) - nuArray(ilow))
            ! carry out interpolation between edges
            do j = ilow, iup-1
                gammaHeI(j) = logGammaHeI(2*i)+aFit* (nuArray(j)&
                     &-nuArray(iup))
                gammaHeI(j) = 10.**gammaHeI(j)  
            end do
        end do 

        ! calculate gammaHeII (nu<4.0Ryd)
        do i = 1, nEdgesHeII
            iup  = HeIINuEdgeP(2*i)   
            ilow = HeIINuEdgeP(2*i - 1)
            if (iup == ilow) then
                print*, "! fb_ff: frequency grid too course - divide&
                     & by zero                          [nu]",&
                     & nuArray(iup)
                stop
            end if
            aFit = (logGammaHeII(2*i)-logGammaHeII(2*i - 1)) /&
                 & (nuArray(iup) - nuArray(ilow))  
            ! carry out interpolation between edges
            do j = ilow, iup-1
                gammaHeII(j) = logGammaHeII(2*i)+aFit* (nuArray(j)&
                     &-nuArray(iup))
                gammaHeII(j) = 10.**gammaHeII(j)
            end do
        end do
        
        ! calculate f-b emission for HI  (nu > 1.), for HeI (nu > 1.8 Ryd) 
        ! and for HeII (nu > 4.0 Ryd) by means of Milne relation 
        ! - note the ff contribution will then need to be added separately 
        ! so set up the free free coefficients
        call free_free()

        constant = 4.9874105e-6 ! Ryd*Ryd*(h^2/(2*Pi*Me*K))**(3/2) [cm*K^(3/2)]
        HIPnuP   = HlevNuP(1)
        HeIPnuP = HeIlevNuP(1)
        HeIIPnuP = HeIIlevNuP(1)
        statW = (/2., 0.5, 2./)
        do i = HIPnuP, nbins
            factor = (nuArray(i)*nuArray(i)*nuArray(i)) / (TeUsed*sqrTeUsed)
            expFactor = exp( (-nuArray(i) + nuArray(HIPnuP)) * hcRyd_k / TeUsed)

            phXSecHI =  xSecArray(i-HIPnuP+1+HlevXSecP(1)-1)
            gammaHI(i) = fourPi * phXSecHI * statW(1) * hcRyd * constant*&
&                          factor * expFactor * 1.e20 *1.e20
            ! add the free free contribution
            gammaHI(i) = gammaHI(i)+ffCoeff1(i)
            ! calculate gammaHeI
            if (i >= HeIPnuP) then
                expFactor = exp( (-nuArray(i) + nuArray(HeIPnuP)) * hcRyd_k / TeUsed)
                phXSecHeI = xSecArray(i-HeIPnuP+1+HeISingXSecP(1)-1)

                gammaHeI(i) = fourPi * phXSecHeI * statW(2) * hcRyd * constant*&
&                          factor * expFactor * 1.e20 *1.e20
                ! add the free free contribution
                gammaHeI(i) = gammaHeI(i) + ffCoeff1(i)
                if (i >= HeIIPnuP) then
                    expFactor  = exp( (-nuArray(i) + nuArray(HeIIPnuP)) * hcRyd_k / TeUsed)
                    phXSecHeII = xSecArray(i-HeIIPnuP+1+HeIIXSecP(1)-1)
                    gammaHeII(i) = fourPi * phXSecHeII * statW(3) * hcRyd * constant*&
&                              factor * expFactor * 1.e20 *1.e20     
                    ! add the free free contribution
                    gammaHeII(i) = gammaHeII(i) + ffCoeff2(i)
                end if               
            end if        
        end do

        ! calculate gammaHeavies
        do i = 1, nbins
           ! only elemnts up to z=19 calculated.. for now (problem with statistical
           ! weights for z>19. Contribution from elements wih Z > 19 generally small. 
           do elem = 3, 19

                do ion = 1, nstages-1
                    ! check if this element is present in the nebula
                    if (.not.lgElementOn(elem)) exit

                    if (ion > elem) exit

                    ! find the number of electron in this ion
                    nElec = elem - ion +1

                    ! find the outer shell number and stat weights
                    call getOuterShell(elem, nElec, outShell, g0, g1)

                    ! get pointer to this ion's high energy limit in nuArray
                    highNuP = elementP(elem, ion, outShell, 2)
                    if ( highNuP > nbins ) then
                        print*, "! fb_ff: high frequency limit beyond grid limit. Increase&
&                             frequency grid, or switch element off [nuMax,elem,&
&                             ion]", nuMax,elem,ion
                        stop
                    end if

                    ! get pointer to this ion's IP in nuArray
                    IPnuP = elementP(elem, ion, outShell, 1)
                    ! check IP within frequency range
                    if ( IPnuP > nbins ) then
                        print*, "! fb_ff: IP frequency beyond grid limit. Increase&
&                             frequency grid, or switch element off [nuMax,elem,&
&                             ion]", nuMax,elem,ion
                        stop
                    end if

                    ! check IP of atom against current energy
                    if ( (i >= IPnuP) .and. (i < highNuP) ) then
  
                        ! get pointer to phot xSec of this ion in xSecArray
                        xSecP = elementP(elem, ion, outShell, 3)

                        ! get phot xSec of this ion at this energy
                        phXSecM = xSecArray(i+xSecP-IPnuP+1-1) 
             
                        ! sum up the coefficients
                        expFactor = exp( (-nuArray(i)+nuArray(IPnuP)) *&
&                                  hcRyd_k / TeUsed)

                        gammaHeavies(i) = gammaHeavies(i) + &
&                                  expFactor * phXSecM * (real(g0)/real(g1)) * &
&                                  ionDenUsed(elementXref(elem), ion+1) * grid%elemAbun(grid%abFileIndex(ix,iy,iz),elem)            

                    end if
                    end do
                end do
            factor = nuArray(i)*nuArray(i)*nuArray(i) / (TeUsed*sqrTeUsed)
            gammaHeavies(i) = gammaHeavies(i) * fourPi * constant *&
&                                hcRyd * factor * 1.e20 * 1.e20

        end do

    end subroutine fb_ff

    subroutine twoPhoton()
        implicit none

        ! local variables
        integer                 :: i      ! counter    

        ! calculate the two photon emission due to HI and HeII
        ! (use the hydrogenic ions function)
        do i = 1, nbins
            twoPhotHI(i) = hydro2phot(nuArray(i), 1)
            twoPhotHeII(i) = hydro2phot(nuArray(i), 2)
        end do
        
        ! calculate the two photon emission due to HeI
        call HeI2phot()
    end subroutine twoPhoton


    ! this function calculates the two photon emission coefficient 
    ! for Hydrogenic ions. See Nussbaumer and Schmutz, A&A 138(1984)495
    function hydro2phot(nu, Z)
        implicit none

        double precision                :: hydro2phot   ! emission coeff for 2-phot cont [e-40 erg/cm^3/s/Hz]
                                            ! for the hydrogenic ion [e-40 erg*cm^3/s/Hz]
        real, intent(in)    :: nu           ! frequency [Ryd]

        integer, intent(in) :: Z            ! nuclear charge

        ! local variables
        real                :: alphaEffH22S ! effective recombination coefficient
                                            ! to the 22S state of HI [e-13 cm^3/s]
        real                :: alphaEffHeII21S ! effective recombination coefficient
                                            ! to the 21S state of HeII [e-13 cm^3/s]
        real, parameter     :: A22S12S=8.23 ! Einstein A of forbibben LyAlpha [1/s]
        real, parameter     :: A2qH = 8.2249! total hydrogen 2s->1s 2 photon
                                            ! transition probability [1/s]
        real                :: A2qZ         ! total Z-ion 2s->1s 2 photon
                                            ! transition probability [1/s]
        real                :: Ay           ! transition probability (nu-dependent)
        real                :: factor       ! general interpolation factor
        real                :: gNu          ! spec distribution of 2phot emission
                                            ! [e-27 erg/Hz]
        real                :: nu0          ! frequencyof the 2s->1s transition
        real                :: Q22S22P      ! collisional transitional rates
                                            ! for HI 22S-22P
        real                :: y            ! nu/nu0 where nu0is the frequency
                                            ! of the 2s->1s transition
        ! select the hydrogenic ion
        select case(Z)
        ! for HI
        case (1)
            ! nu [Ryd] of 2s->1s transition (=1215.7 A)
            nu0 = 0.7496
            y = nu/nu0

            ! check for nu >= nu0
            if (nu >= nu0 ) then
                hydro2phot = 0.
                return
            end if            

            ! for alphaEffH22S see Pengelly, MNRAS, 127(1964)145
            alphaEffH22S = 0.8368*((TeUsed*1.e-4)**(-0.723))

            Ay = 202.0*(y*(1-y)*(1-(4*y*(1-y))**0.8) + &
&                 0.88*((y*(1-y))**1.53)*(4*y*(1-y))**0.8)

            gNu = hPlanck*y*Ay/A2qH
            gNu = gNu*1e27

            ! calculate the collisional transitional rates for HI 22S-22P
            ! interpolate over temperature the value given in table 4.10
            ! from Osterbrock
            if (TeUsed <= 10000.) then
                Q22S22P = 5.31e-4
            else if (TeUsed >= 20000.) then
                Q22S22P = 4.71e-4
            else
                factor = log10(4.71e-4/5.31e-4) / log10(2.)
                Q22S22P = 10**( log10(5.31e-4) + factor*log10(TeUsed/10000.) )
            end if

            hydro2phot = alphaEffH22s*gNu / (1.+(NeUsed*Q22S22P)/A22S12S)
             
        case (2) ! HeII
            !  nu [Ryd] of 2s->1s transition (=303.8 A)
            nu0 = 3.00
            y = nu/nu0
                
            ! check for nu >= nu0
            if (nu >= nu0 ) then 
                hydro2phot = 0. 
                return
            end if

            ! for alphaEffHeII21S see Storey & Hummer, MNRAS, 272(1995)41
            ! the values used are for Ne = 100 cm^-3; however alphaEffHeII21S is not
            ! very sensitive to density.. only interpolate over temperature
            if (TeUsed <= 5000.) then
                alphaEffHeII21S = 6.161
            else if (TeUsed >= 30000.) then
                alphaEffHeII21S = 2.035
            else if (TeUsed > 5000. .and. TeUsed <= 10000. ) then
                factor = log10(4.091/6.161) / log10(2.)
                alphaEffHeII21S = 10**( log10(6.161) + factor*log10(TeUsed/5000.) )
            else if (TeUsed > 10000. .and. TeUsed <= 15000. ) then
                factor = log10(3.189/4.091) / log10(1.5)
                alphaEffHeII21S = 10**( log10(4.091) + factor*log10(TeUsed/10000.) )
            else if (TeUsed > 15000. .and. TeUsed < 30000. ) then
                factor = log10(2.035/3.189) / log10(2.)
                alphaEffHeII21S = 10**( log10(3.189) + factor*log10(TeUsed/15000.) )
            end if

            Ay = (Z**6)*0.9994667*202.0*(y*(1-y)*(1-(4*y*(1-y))**0.8) + &
&                 0.88*((y*(1-y))**1.53)*(4*y*(1-y))**0.8)
            
            A2qZ = 8.226*Z**6
            
            gNu = hPlanck*y*Ay/A2qZ
            gNu = gNu*1e27
            
            ! collisional de-excitation of the 22S of HeII is negligible
        
            ! caalculate the 2 photon emission
            hydro2phot = alphaEffHeII21s*gNu
        

        end select
        
    end function hydro2phot
    
    subroutine HeI2phot()
        implicit none 

        ! local variables 
        integer             :: i, j          ! counters
        integer             :: ios           ! I/O error status

        real, parameter     :: A21S11S=51.3  ! HeI 2q [1/s]
                                             ! Almog & Netzer, MNRAS 238(1989)57 
        real, parameter     :: nu0=1.514     ! 21s -> 11s frequency [Ryd] 

        real                :: alphaEff21SHeI! effective rec coeff to HeI 21S [e-14 cm^3/s]
                                             ! Almog & Netzer, MNRAS 238(1989)57 
        real                :: Ay            ! interpolated rate [1/s]
        real                :: fit1, fit2    ! interpolation coefficients 
        real                :: y             ! nu/nu0

        real, dimension(41) :: Ay_dat        ! data point in Ay (rates) [1/s]
        real, dimension(41) :: y_dat         ! data point in y (=nu/nu0) 

        ! assume all HeI singlets finally end up in the 2^1S 
        ! use total recombination cofficient to all singlets Benjamin, Skillman and SMits, ApJ, 1999, 514, 307
        alphaEff21SHeI = 6.23*((TeUsed/10000.)**(-0.827))
        
        ! read in rates from data/HeI2phot.dat
        close(93)
        open(unit = 93, file = "data/HeI2phot.dat", status = "old", position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! HeI2phot: can't open file: data/HeI2phot.dat"
            stop
        end if
        do i = 1, 41
            read(unit = 93, fmt = *) y_dat(i), Ay_dat(i)
        end do

        close(93)

        ! calculate coefficients as a function of frequency
        j = 1
        do i = 1, nbins
            ! calculate y
            y  = nuArray(i)/nu0

            if ( y > y_dat(j) ) j=j+1
            if ( j >= 41 ) exit

            fit1 = (Ay_dat(j+1)-Ay_dat(j)) / &
&                                  (y_dat(j+1)-y_dat(j))
            fit2 = Ay_dat(j) - y_dat(j)*fit1

            Ay = fit1*y + fit2

            if ( Ay < 9.2163086E-03 ) Ay = 0. 

            twoPhotHeI(i) = alphaEff21SHeI*Ay*0.66262*y/A21S11S 
        end do        
        
    end subroutine HeI2phot

    ! this subroutine sets the ff emission coefficients for z=1 (ffCoeff1)
    ! and z=2 (ffCoeff2)            
    subroutine free_free()
        implicit none
  
        ! local variables
        integer        :: i   ! counter

        do i = 1, nbins
            ! find the free-free emission coefficient for H+ and He+
            ffCoeff1(i) = ffCoeff(nuArray(i), 1, gauntFF(i))
            ! find the free-free emission coefficient for He++
            ffCoeff2(i) = ffCoeff(nuArray(i), 2, gauntFFHeII(i)) 
        end do
    end subroutine free_free
    
    ! this function calculates the ff emission coeficient for interactions between
    ! ions of nuclear charge Z and electrons. 
    ! see Allen pg 103.
    function ffCoeff(nu, Z, g)
        implicit none

        real                 :: ffCoeff     ! ff emission coefficient [e-40 erg*cm^3/s.Hz]
        real, intent(in)     :: g           ! ff gaunt 
        real, intent(in)     :: nu          ! frequency [Ryd]

        integer, intent(in)  :: Z           ! nuclear charge

        ! local variables
        real                 :: expFactor   ! exponential factor

        expFactor = exp(-nu*hcRyd_k/TeUsed)

        ffCoeff = fourPi*54.43*Z*Z*g*expFactor/sqrTeUsed

    end function ffCoeff

    subroutine RecLinesEmission()
        implicit none

        ! local variables
        integer                    :: ios         ! I/O error status
        integer                    :: i           ! counters
        integer                    :: ilow,&      ! pointer to lower level
&                                      iup        ! pointer to upper level

        real                       :: A4471, A4922! HeI reference lines 
        real                       :: aFit        ! general interpolation coeff
        real                       :: C5876, C6678! collition exc. corrections
        real                       :: Hbeta       ! Hbeta emission
        real                       :: HeII4686    ! HeII 4686 emission
        real                       :: Lalpha      ! Lalpha emission
        real                       :: T4          ! TeUsed/10000.

        real, dimension(9,3)       :: HeISingRead ! array reader for HeI sing rec lines
        real, dimension(11,3)       :: HeITripRead ! array reader for HeI trip  rec lines  

        ! do hydrogenic ions first

        ! read in HI recombination lines [e-25 ergs*cm^3/s] 
        ! (Storey and Hummer MNRAS 272(1995)41)
        close(94)
        open(unit = 94, file = "data/r1b0100.dat", status = "old", position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! RecLinesEmission: can't open file: data/r1b0100.dat"
            stop
        end if
        do iup = 15, 3, -1
            read(94, fmt=*) (HIRecLines(iup, ilow), ilow = 2, min0(8, iup-1)) 
        end do

        close(94)

        ! calculate Hbeta 
!        Hbeta = 2530./(TeUsed**0.833) ! TeUsed < 26000K CASE 
        ! fits to Storey and Hummer MNRAS 272(1995)41
        Hbeta = 10**(-0.870*log10Te + 3.57)
        Hbeta = Hbeta*NeUsed*ionDenUsed(elementXref(1),2)*grid%elemAbun(grid%abFileIndex(ix,iy,iz),1)

        ! calculate emission due to HI recombination lines [e-25 ergs/s/cm^3]
        do iup = 15, 3, -1
            do ilow = 2, min0(8, iup-1)
                HIRecLines(iup, ilow) = HIRecLines(iup, ilow)*Hbeta
            end do
        end do    

        ! add contribution of Lyman alpha 
        ! fits to Storey and Hummer MNRAS 272(1995)41
        Lalpha = 10**(-0.897*log10Te + 5.05) 
        HIRecLines(15, 8) =HIRecLines(15, 8) + grid%elemAbun(grid%abFileIndex(ix,iy,iz),1)*ionDenUsed(elementXref(1),2)*&
             & NeUsed*Lalpha 
        

        ! read in HeII recombination lines [e-25 ergs*cm^3/s]
        ! (Storey and Hummer MNRAS 272(1995)41)
        close(95)
        open(unit = 95, file = "data/r2b0100.dat", status = "old", position = "rewind", iostat=ios)
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
        HeII4686 = HeII4686*NeUsed*grid%elemAbun(grid%abFileIndex(ix,iy,iz),2)*ionDenUsed(elementXref(2),3)

        ! calculate emission due to HeI recombination lines [e-25 ergs/s/cm^3]
        do iup = 30, 3, -1
            do ilow = 2, min0(16, iup-1)
                HeIIRecLines(iup, ilow) = HeIIRecLines(iup, ilow)*HeII4686
            end do
        end do    

        ! now do HeI
        
        ! read in HeI singlet recombination lines [e-25 ergs*cm^3/s]
        ! Benjamin, Skillman and Smits ApJ514(1999)307
        close(96)
        open(unit = 96, file = "data/heIrecS.dat", status = "old", position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! RecLinesEmission: can't open file: data/heIrecS.dat"
            stop
        end if
        do i = 1, 9
            read(96, fmt=*) HeISingRead(i, :)

            ! interpolate over temperature (log10-log10 space)
            if (TeUsed <= 5000.) then
                HeIRecLinesS(i) = HeISingRead(i, 1)
            else if (TeUsed >= 20000.) then
                HeIRecLinesS(i) = HeISingRead(i, 3)
            else if ((TeUsed > 5000.) .and. (TeUsed <= 10000.)) then
                aFit = log10(HeISingRead(i, 2)/HeISingRead(i, 1)) / log10(2.)
                HeIRecLinesS(i) = 10**(log10(HeISingRead(i, 1)) + &
&                                  aFit*log10(TeUsed/5000.) ) 
            else if ((TeUsed > 10000.) .and. (TeUsed < 20000.)) then
                aFit = log10(HeISingRead(i, 3)/HeISingRead(i, 2)) / log10(2.)
                HeIRecLinesS(i) = 10**(log10(HeISingRead(i, 2)) + &
&                                  aFit*log10(TeUsed/10000.) )
            end if
        end do

        close(96)
        
        ! read in HeI triplet recombination lines [e-25 ergs*cm^3/s]
        ! Benjamin, Skillman and Smits ApJ514(1999)307
        close(97)
        open(unit = 97, file = "data/heIrecT.dat", status = "old", position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! RecLinesEmission: can't open file: data/heIrecT.dat"
            stop
        end if
        do i = 1, 11
            read(97, fmt=*) HeITripRead(i, :)

            ! interpolate over temperature (log10-log10 space)
            if (TeUsed <= 5000.) then
                HeIRecLinesT(i) = HeITripRead(i, 1)
            else if (TeUsed >= 20000) then
                HeIRecLinesT(i) = HeITripRead(i, 3)
            else if ((TeUsed > 5000.) .and. (TeUsed <= 10000.)) then
                aFit = log10(HeITripRead(i, 2)/HeITripRead(i, 1)) / log10(2.)
                HeIRecLinesT(i) = 10**(log10(HeITripRead(i, 1)) + &
&                                  aFit*log10(TeUsed/5000.) )
            else if ((TeUsed > 10000.) .and. (TeUsed < 20000.)) then
                aFit = log10(HeITripRead(i, 3)/HeITripRead(i, 2)) / log10(2.)
                HeIRecLinesT(i) = 10**(log10(HeITripRead(i, 2)) + &
&                                  aFit*log10(TeUsed/10000.) )
            end if
        end do

        close(97)

        HeIRecLinesS = HeIRecLinesS*NeUsed*grid%elemAbun(grid%abFileIndex(ix,iy,iz),2)*ionDenUsed(elementXref(2),2)
        HeIRecLinesT = HeIRecLinesT*NeUsed*grid%elemAbun(grid%abFileIndex(ix,iy,iz),2)*ionDenUsed(elementXref(2),2)

                                                   
    end subroutine RecLinesEmission


    ! this subroutine is the driver for the calculation of the emissivity
    ! from the heavy elements forbidden lines. 
    subroutine forLines()
        implicit none

        integer :: elem, ion ! counters

        ! re-initialize forbiddenLines
        forbiddenLines = 0.

        do elem = 3, nElements
           do ion = 1, min(elem+1, nstages)
              if (.not.lgElementOn(elem)) exit

              if (lgDataAvailable(elem, ion)) then

                 if (ion<nstages) then
                    call equilibrium(file_name=dataFile(elem, ion), &
                         &ionDenUp=ionDenUsed(elementXref(elem),ion+1), Te=TeUsed,&
                         &Ne=NeUsed, FlineEm=forbiddenLines(elem, ion,1:nForLevels, 1:nForLevels))
                 else
                    call equilibrium(file_name=dataFile(elem, ion), &
                         &ionDenUp=0., Te=TeUsed, Ne=NeUsed, &
                         &FlineEm=forbiddenLines(elem, ion,1:nForLevels,1:nForLevels))
                 end if

                 forbiddenLines(elem, ion, :, :) =&
                      &forbiddenLines(elem, ion, :, :)*grid%elemAbun(grid%abFileIndex(ix,iy,iz),elem)*&
                      & ionDenUsed(elementXref(elem), ion)

              end if

           end do
        end do


        ! scale the forbidden lines emissivity to give units of [10^-25 erg/s/Ngas]
        ! comment: the forbidden line emissivity is so far in units of cm^-1/s/Ngas
        !          the energy [erg] of unit wave number [cm^-1] is 1.9865e-16, hence 
        !          the right units are obtained by multiplying by 1.9865e-16*1e25 
        !          which is 1.9865*1e9 
        forbiddenLines = forbiddenLines*1.9865e9

    end subroutine forLines    
        
    subroutine setDiffusePDF()
        implicit none

        ! local variables
        integer           :: elem, ion          ! counters
        integer           :: err                ! allocation error status
        integer, parameter:: NHeIILyman = 4     ! Number of HeII Lym lines included 
        integer           :: j2TsP              ! pointer to 2s triplet state in nuArray
        integer           :: i,iup,ilow,j       ! counters       
        integer           :: ios                ! I/O error status
        integer           :: nuP                ! frequency pointer in nuArray

        real              :: A4922             ! reference line intensity for HeI Lyman series
        real              :: aFit              ! general calculations component
        real              :: alpha2tS          ! effective rec coeff to the 2s trip HeI
        real              :: alpha2tP          ! effective rec coeff to the 2p trip HeI
        real              :: bFit              ! general calculations component
        real              :: correction        ! general correction term
        real              :: normalize         ! normHI+normHeI+normHeII
        real              :: normFor           ! normalization constant for for lines
        real              :: normHI            ! normalization constant for HI
        real              :: normHeI           !                            HeI
        real              :: normHeII          !                            HeII
        real              :: normRec           !                            rec. lines
        real              :: T4                ! TeUsed/10000.

        real, dimension(NHeIILyman) &
&                          :: HeIILyman         ! HeII Lyman lines em. [e-25ergs*cm^3/s]
        real, dimension(NHeIILyman) &
&                          :: HeIILymanNu       ! HeII Lyman lines freq. [Ryd]          

        real, pointer     :: sumDiffuseHI(:)   ! summation terms for HI emission
        real, pointer     :: sumDiffuseHeI(:)  ! summation terms for HeI emission   
        real, pointer     :: sumDiffuseHeII(:) ! summation terms for HeII emission   


        logical, save     :: lgFirst =.true.   ! first time setDiffusePDF is called? 

        ! allocate space for the summation arrays
        allocate(sumDiffuseHI(1:nbins), stat = err)
        if (err /= 0) then
            print*, "! setDiffusePDF: can't allocate recPDF array memory"
            stop
        end if
        allocate(sumDiffuseHeI(1:nbins), stat = err)
        if (err /= 0) then
            print*, "! setDiffusePDF: can't allocate recPDF array memory"
            stop  
        end if
        allocate(sumDiffuseHeII(1:nbins), stat = err)
        if (err /= 0) then
            print*, "! setDiffusePDF: can't allocate recPDF array memory"
            stop  
        end if



        ! assign pointers for the He line photons
        call locate(nuArray, 1.45673, j2TsP)

        ! calculate normalization constants for the HI, HeI and HeII diffuse probs

        ! initialize constants
        normHI   = 0. 
        normHeI  = 0.
        normHeII = 0.

        sumDiffuseHI = 0.
        sumDiffuseHeI = 0.
        sumDiffuseHeII = 0.

        ! perform summations
        do i = 1, nbins

           sumDiffuseHI(i) =  cRyd*widFlx(i)*emissionHI(i)/1.e15
           
           normHI = normHI + sumDiffuseHI(i)
           
           sumDiffuseHeI(i) =  cRyd*widFlx(i)*emissionHeI(i)/1.e15
            
           normHeI = normHeI + sumDiffuseHeI(i)
   

           sumDiffuseHeII(i) = cRyd*widFlx(i)*emissionHeII(i)/1.e15

           normHeII = normHeII + sumDiffuseHeII(i)

        end do


        ! Add contribution of HeI lines able to ionize HI 
        
        ! triplet states:
        ! calculate effective recombination coefficients to triplet states[e-14 cm^3/s]
        ! Benjamin, SKillman and mits, 1999, ApJ 514, 307
        T4 = TeUsed*1.e-4
        alpha2tS = 27.2*(T4**(-0.678))       
  
        ! calculate correction for partial sums and normalization constant 
        ! triplet 2s 

        correction = alpha2tS * &
             & NeUsed*grid%elemAbun(grid%abFileIndex(ix,iy,iz),2)*&
             &ionDenUsed(elementXref(2),2)*1.45673*hcRyd*1e11

        sumDiffuseHeI(j2TsP) = sumDiffuseHeI(j2TsP) + correction
        normHeI = normHeI + correction

        ! Add contribution of HeII lines able to ionize HI and HeII

        ! read in HeII Lyman line ratios up to level n=5 [e-25 ergs*cm^3/s]
        ! (Storey and Hummer MNRAS 272(1995)41)
        close(98)
        open(unit = 98, file = "data/r2a0100.dat", status = "old", position = "rewind", iostat=ios)
        if (ios /= 0) then
            print*, "! setDiffusePDF: can't open file: data/r2a0100.dat"
            stop
        end if
        do i = 1, NHeIILyman
            read(98, fmt=*) HeIILyman(i), HeIILymanNu(i)
        end do
        
        close(98)

        ! calculate Lyman alpha first
        HeIILyman(4) = 10.**(-0.792*log10(TeUsed)+6.01)
        HeIILyman(4) = HeIILyman(4)*NeUsed*ionDenUsed(elementXref(2),3)*grid%elemAbun(grid%abFileIndex(ix,iy,iz),2)

        ! calculate emission due to HeII Lyman lines [e-25 ergs/s/cm^3]
        do i = 1, NHeIILyman-1
            HeIILyman(i) = HeIILyman(i)*HeIILyman(4)
        end do

        ! add contributions in respective bins
        do i = 1, NHeIILyman
            ! locate the frequency bin in the array
            call locate(NuArray, HeIILymanNu(i), NuP)
            ! add contribution to HeII partial sums ...
            sumDiffuseHeII(NuP) = sumDiffuseHeII(NuP) + HeIILyman(i)
            ! ... and to HeII normilazation constant
            normHeII = normHeII + HeIILyman(i)
        end do

        ! Total normalization constant
        normalize = normHI + normHeI + normHeII

        ! Sum  up energy in recombination lines
        
        normRec = 0.
        ! HI
        do iup = 3, 15
            do ilow = 2, min0(8, iup-1)
                normRec = normRec+HIRecLines(iup, ilow)
            end do
        end do
        ! HeII            
        do iup = 3, 30
            do ilow = 2, min0(16, iup-1)
                normRec = normRec+HeIIRecLines(iup, ilow)
            end do
        end do
        ! HeI singlets
        do i = 1, 9
            normRec = normRec+HeIRecLinesS(i)
        end do
        ! HeI triplets
        do i = 1, 11
            normRec = normRec+HeIRecLinesT(i)
        end do

        ! Sum  up energy in forbidden lines
            
        normFor = 0.
        do elem = 3, nElements
           do ion = 1, min(elem+1, nstages)
              do ilow = 1, nForLevels
                 do iup = 1, nForLevels
                    if (forbiddenLines(elem, ion, ilow, iup) > 1.e-35)  then 
                       normFor = normFor + forbiddenLines(elem, ion, ilow, iup)
                    end if
                 end do
              end do
           end do
        end do

        ! total energy in lines (note: this variable will later be turned into the fraction of
        ! non-ionizing line photons as in declaration)

        grid%totalLines(grid%active(ix, iy, iz)) = normRec + normFor

        do i = 1, nbins
            if (i == 1) then
                    grid%recPDF(grid%active(ix, iy, iz), i) = (sumDiffuseHI(i) + sumDiffuseHeI(i) +&
&                                          sumDiffuseHeII(i) ) 
            else

                    grid%recPDF(grid%active(ix, iy, iz), i) = grid%recPDF(grid%active(ix, iy, iz), i-1) +&
&                                         (sumDiffuseHI(i) + sumDiffuseHeI(i) +&
&                                          sumDiffuseHeII(i) )
            end if

        end do
        grid%recPDF(grid%active(ix, iy, iz), :) = grid%recPDF(grid%active(ix, iy, iz), :)/normalize
        grid%recPDF(grid%active(ix, iy, iz), nbins) = 1.0000000001

        ! calculate the PDF for recombination and forbidden lines
        if (lgDebug) then 
           i = 1
           ! HI recombination lines
           do iup = 3, 15
              do ilow = 2, min(8, iup-1)
                 if (i == 1) then

                    grid%linePDF(grid%active(ix, iy, iz), i) = HIRecLines(iup, ilow) / &
&                                                  grid%totalLines(grid%active(ix,iy,iz))



                 else
                    
                    grid%linePDF(grid%active(ix, iy, iz), i) = grid%linePDF(grid%active(ix, iy, iz), i-1) + &
&                        HIRecLines(iup, ilow) / grid%totalLines(grid%active(ix, iy, iz))



                 end if
                 i = i+1

              end do
           end do

           ! HeI singlet recombination lines
           do j = 1, 9
              grid%linePDF(grid%active(ix, iy, iz), i) = grid%linePDF(grid%active(ix, iy, iz), i-1) + &
&                           HeIRecLinesS(j) / grid%totalLines(grid%active(ix, iy, iz))
              i = i+1 
           end do

           ! HeI triplet recombination lines
           do j = 1, 11
              grid%linePDF(grid%active(ix, iy, iz), i) = grid%linePDF(grid%active(ix, iy, iz), i-1) + &
&                           HeIRecLinesT(j) / grid%totalLines(grid%active(ix, iy, iz))
              i = i+1
           end do

           ! HeII recombination lines
           do iup = 3, 30
              do ilow = 2, min(16, iup -1)
                 grid%linePDF(grid%active(ix, iy, iz), i) = grid%linePDF(grid%active(ix, iy,&
                      & iz), i-1) + HeIIRecLines(iup, ilow) / grid&
                      &%totalLines(grid%active(ix,iy,iz))
                 i = i+1
              end do
           end do

           ! heavy elements forbidden lines
           do elem = 3, nElements
              do ion = 1, min(elem+1, nstages)

                 if (.not.lgElementOn(elem)) exit
                 if (lgDataAvailable(elem, ion)) then
                    do iup = 1, nForLevels
                       do ilow = 1, nForLevels
                          grid%linePDF(grid%active(ix, iy, iz), i) = grid%linePDF(grid%active(ix, iy&
                               &, iz), i-1) + forbiddenLines(elem, ion, iup, ilow)&
                               & / grid%totalLines(grid%active(ix,iy,iz))
                       
                          if (grid%linePDF(grid%active(ix, iy, iz), i) > 1. ) grid&
                               &%linePDF(grid%active(ix, iy, iz), i) = 1. 
                          i = i+1

                       end do
                    end do
                 end if

              end do
           end do
           grid%linePDF(grid%active(ix, iy, iz), nLines) = 1.

           ! calculate the probability of Hbeta
           HbetaProb = HIRecLines(4,2) / (grid%totalLines(grid%active(ix,iy,iz))&
                &+normalize)

        end if ! debug condition

        ! now compute the fraction of non-ionizing photons

        normalize = normalize + grid%totalLines(grid%active(ix,iy,iz))

        grid%totalLines(grid%active(ix,iy,iz)) = grid%totalLines(grid%active(ix,iy,iz)) / normalize

        ! deallocate arrays
        if ( associated(sumDiffuseHI) ) deallocate(sumDiffuseHI)
        if ( associated(sumDiffuseHeI) ) deallocate(sumDiffuseHeI)
        if ( associated(sumDiffuseHeII) ) deallocate(sumDiffuseHeII)

    end subroutine setDiffusePDF  


    subroutine setDustPDF()
      implicit none

      real :: bb

      integer :: i, n, ai ! counters

      grid%dustPDF(grid%active(ix,iy,iz),:)=0.

      contShape = 'blackbody'

      do n = 1, nSpecies
         do ai = 1, nSizes
            if (grid%Tdust(n,ai,grid%active(ix,iy,iz))>0. .and. & 
                 & grid%Tdust(n,ai,grid%active(ix,iy,iz))<TdustSublime(n)) then
               bb = getFlux(nuArray(1), grid%Tdust(n,ai,grid%active(ix,iy,iz)))
               grid%dustPDF(grid%active(ix,iy,iz), 1) = xSecArray(dustAbsXsecP(n,ai))*bb*widFlx(1)
            end if
         end do
      end do

      do i = 2, nbins
         do n = 1, nSpecies
            do ai = 1, nSizes
               if (grid%Tdust(n,ai,grid%active(ix,iy,iz))>0. .and. &
                    & grid%Tdust(n,ai,grid%active(ix,iy,iz))<TdustSublime(n)) then

                  bb = getFlux(nuArray(i), grid%Tdust(n,ai,grid%active(ix,iy,iz)))
                  grid%dustPDF(grid%active(ix,iy,iz),i) = grid%dustPDF(grid%active(ix,iy,iz),i-1) + &
                       &  xSecArray(dustAbsXsecP(n,ai)+i-1)*bb*widFlx(i)

               end if
            end do
         end do
      end do

      contShape = contShapeIn

      ! normalise
      do i = 1, nbins

         if (grid%dustPDF(grid%active(ix,iy,iz),nbins)>0.) grid%dustPDF(grid%active(ix,iy,iz),i) = &
              & grid%dustPDF(grid%active(ix,iy,iz),i)/grid%dustPDF(grid%active(ix,iy,iz),nbins)

      end do

    end subroutine setDustPDF


    end subroutine emissionDriver

    subroutine equilibrium(file_name, ionDenUp, Te, Ne, fLineEm, wav)
        implicit none

        character(len = *), &
&          intent(in)  :: file_name   ! ionic data file name

        real, intent(in) &
&                      :: Te, &       ! electron temperature [K]
&                         Ne, &       ! electron density [cm^-3]
&                         ionDenUp    ! ion density of the upper ion stage 

        double precision, dimension(nForLevels,nForLevels), &
&          intent(out) :: fLineEm     ! forbidden line emissivity

        double precision, dimension(nForLevels,nForLevels), &
&          intent(out), optional :: wav         ! wavelength of transition [A]

        ! local variables 

        integer  :: gx               ! stat weight reader  
        integer  :: i, j, k, l, iT   ! counters/indeces
        integer  :: iRats            ! coll strength (iRats=0) or (coll rates)/10**iRats
        integer  :: ios              ! I/O error status
        integer  :: nLev             ! number of levels in atomic data file
        integer  :: numLines         ! number of lines for atomic data refernce

       integer  :: nTemp            ! number of temperature points in atomic data file        
        integer  :: err              ! allocation error status

        integer, parameter :: safeLim = 10000 ! loop safety limit

        integer, pointer :: g(:)              ! statistical weight array

        integer, dimension(2) :: ilow, &      ! lower index
&                                 iup          ! upper index

        character(len = 20), pointer :: &
&                                      label(:)! labels array
        character(len = 75) :: text  ! lines of text

        real                 :: alphaTotal(100)           ! maximum 100-level ion
        real                 :: a_fit, b_fit              !         
        double precision     :: ax, ex                    ! readers
        double precision     :: constant                  ! calculations constant 
        double precision     :: delTeK                    ! Boltzmann exponent
        double precision     :: expFac                    ! calculations factor
        double precision     :: Eji                       ! energy between levels j and i
        real       :: qomInt                    ! interpolated value of qom
        double precision       :: qx                        ! reader
        double precision   :: sumN                      ! normalization factor for populations
        double precision        :: sqrTe                     ! sqrt(Te)
        double precision    :: value                     ! general calculations value

        double precision, pointer          :: a(:,:)      ! transition rates array
        double precision, pointer          :: cs(:,:)     ! collisional strengths array
        double precision, pointer          :: e(:)        ! energy levels array
        real,    pointer          :: logTemp(:)  ! log10 temperature points array
        double precision, pointer          :: qeff(:,:)   ! q eff  array
        double precision,    pointer          :: qom(:,:,:)  ! qom array
        real,    pointer          :: qq(:)       ! qq array
        real,    pointer          :: qq2(:)      ! 2nd deriv qq2 array
        double precision, pointer          :: tnij(:,:)   ! tnij array
        double precision, pointer          :: x(:,:)      ! matrix arrays
        double precision, pointer          :: y(:)        ! 

       

        double precision, &
&                       pointer :: n(:), n2(:) ! level population arrays

        sqrTe   = sqrt(Te)
        log10Te = log10(Te) 

        ! open file containing atomic data
        close(11)
        open(unit=11, file = file_name, status="old", position="rewind", &
&                                                     iostat = ios)
        if (ios /= 0) then
            print*, "! equilibrium: can't open file: ", file_name
            stop
        end if

        ! read reference heading
        read(11, *) numLines
        do i = 1, numLines
            read(11, '(78A1)') text
        end do

        ! read number of levels and temperature points available
        read(11, *) nLev, nTemp

        ! allocate space for labels array
        allocate(label(nLev), stat = err)
        if (err /= 0) then
            print*, "! equilibrium: can't allocate label array memory"
            stop
        end if

        ! allocate space for tempertaure points array
        allocate(logTemp(nTemp), stat = err)
        if (err /= 0) then
            print*, "! equilibrium: can't allocate logTemp array memory"
            stop
        end if

        ! allocate space for transition probability  array
        allocate(a(nLev, nLev), stat = err)
        if (err /= 0) then 
            print*, "! equilibrium: can't allocate a array memory"
            stop
        end if

        ! allocate space for energy levels array
        allocate(e(nLev), stat = err)
        if (err /= 0) then 
            print*, "! equilibrium: can't allocate e array memory"
            stop
        end if

        ! allocate space for statistical weights array
        allocate(g(nLev), stat = err)
        if (err /= 0) then 
            print*, "! equilibrium: can't allocate g array memory"
            stop
        end if

        ! allocate space for qom array
        allocate(qom(nTemp, nLev, nLev), stat = err)
        if (err /= 0) then 
            print*, "! equilibrium: can't allocate qom array memory"
            stop
        end if

        ! allocate space for collisional strengths array
        allocate(cs(nLev, nLev), stat = err)
        if (err /= 0) then
            print*, "! equilibrium: can't allocate cm array memory"
            stop
        end if

        ! allocate space for q eff array
        allocate(qeff(nLev, nLev), stat = err)
        if (err /= 0) then
            print*, "! equilibrium: can't allocate qeff array memory"
            stop
        end if

        ! allocate space for tnij array
        allocate(tnij(nLev, nLev), stat = err)
        if (err /= 0) then
            print*, "! equilibrium: can't allocate tnij array memory"
            stop
        end if

        ! allocate space for x array
        allocate(x(nLev, nLev), stat = err)
        if (err /= 0) then
            print*, "! equilibrium: can't allocate x array memory"
            stop
        end if

        ! allocate space for y array
        allocate(y(nLev), stat = err)
        if (err /= 0) then
            print*, "! equilibrium: can't allocate y array memory"
            stop
        end if

        ! allocate space for qq array
        allocate(qq(nTemp), stat = err)
        if (err /= 0) then
            print*, "! equilibrium: can't allocate qq array memory"
            stop
        end if

        ! allocate space for qq2 array
        allocate(qq2(nTemp), stat = err)
        if (err /= 0) then
            print*, "! equilibrium: can't allocate qq2 array memory"
            stop
        end if

        ! allocate space for n array
        allocate(n(nLev), stat = err)
        if (err /= 0) then
            print*, "! equilibrium: can't allocate n array memory"
            stop
        end if

        ! allocate space for n2 array
        allocate(n2(nLev), stat = err)
        if (err /= 0) then
            print*, "! equilibrium: can't allocate n2 array memory"
            stop
        end if

        ! zero out arrays
        a = 0.
        cs = 0.
        e = 0.
        fLineEm = 0.
        g = 0
        n = 0.
        n2 = 0.
        qom = 0.
        qeff = 0.
        qq = 0.
        qq2 = 0.
        logTemp = 0.
        tnij = 0.
        x = 0.d0
        y = 0.

        ! read labels
        do i = 1, nLev
            read(11, '(A20)') label(i)

        end do

        ! read temperetaure points and get logs
        do i = 1, nTemp
            read(11, *) logTemp(i)
            logTemp(i) = log10(logTemp(i))

        end do

        ! read iRats
        read(11, *) iRats

        do i = 1, safeLim
            ! read in data
            read(11, *) ilow(2), iup(2), qx

            ! check if end of qx dat
            if (qx == 0.d0) exit
            
            ! ilow(2) always starts with a non zero value
            ! so the else condition is true at first and k initialized
            if (ilow(2) == 0) then
                ilow(2) = ilow(1)
                k = k + 1
            else 
                ilow(1) = ilow(2)
                k = 1
            end if 

            ! the same as above
            if (iup(2) == 0) then
                iup(2) = iup(1)                                                  
            else
                iup(1) = iup(2)
            end if

            qom(k, ilow(2), iup(2)) = qx 

        end do

        ! read in transition probabilities
        do k = 1, nLev-1
            do l = k+1, nLev
                read(11, *) i, j, ax

                a(j,i) = ax
            end do
        end do

        ! read statistical weights, energy levels [1/cm]
        do j = 1, nLev
            read(11, *) i, gx, ex

            g(i) = gx
            e(i) = ex
        end do

        ! read power law fit coefficients [e-13 cm^3/s]
        ! and calculate total recombination coefficient
        ! (direct + cascades)
        alphaTotal = 0.

        do j = 2, nLev
           read(unit=11,fmt=*,iostat=ios) a_fit, b_fit

           if (ios<0) then
              exit
           else
              alphaTotal(1) = 1.
           end if

           if (nLev>100) then
              print*, '! equilibrium: alphaTotal set to a &
                   &max 100 level atom, please expand'
              stop
           end if

           alphaTotal(j) = a_fit * (TeUsed/1.e4)**(b_fit) * 1.e-13

        end do

        ! close atomic data file 
        close(11)

        ! form matrices
        ! set up qeff
        do i = 2, nLev
            do j = i, nLev
                do iT = 1, nTemp
                    qq(iT) = qom(iT, i-1, j)
                end do

                if (nTemp == 1) then  
                   ! collisional strength available only for one temperature - assume constant
                   qomInt = qq(1)
                else if (nTemp == 2) then
                   ! collisional strength available only for one temperature - linear interpolation
                   qomInt = qq(1) + &
                        (qq(2)-qq(1)) / (logTemp(2)-logTemp(1)) * (log10Te - logTemp(1))
                else
                   ! set up second derivatives for spline interpolation
                   call spline(logTemp, qq, 1.e30, 1.e30, qq2)

                   ! get interpolated qom for level
                   call splint(logTemp, qq, qq2, log10Te, qomInt)
                end if

                ! set collisional strength
                cs(i-1, j) = qomInt

                ! the calculation constant here is the energy [erg] 
                ! associated to unit wavenumber [1/cm] divided by the 
                ! boltzmann constant k.
                constant = 1.4388463d0
                ! exponential factor 
                delTeK = (e(i-1)-e(j))*constant
                expFac = exp( delTeK/Te )
                qeff(i-1, j) = 8.63d-6 * cs(i-1, j) * expFac /&
&                               (g(i-1)*sqrTe)
                qeff(j, i-1) = 8.63d-6 * cs(i-1, j) / (g(j)*sqrTe)
            end do
        end do

        ! set up x
        do i= 2, nLev
            do j = 1, nLev

               x(1,:) = 1.
               y = 0.
               y(1)   = 1.

                if (j /= i) then
                    x(i, j) = x(i, j) + Ne*qeff(j, i)
                    x(i, i) = x(i, i) - Ne*qeff(i, j)
                    if (j > i) then
                        x(i, j) = x(i, j) + a(j, i)
                    else 
                        x(i, i) = x(i, i) - a(i, j)
                    end if
                end if
            end do

            ! when coefficients are available use the following:
            y(i) = -Ne*ionDenUp*alphaTotal(i)
        end do
        ! 
!        do i = 2, nLev        
!            value = 0.d0 - x(i, 1)
!            y(i-1) = value
!            do j = 2, nLev
!                value = x(i, j)
!                x(i-1,j-1) = value
!            end do
!        end do

        ! solve matrices for populations
!        call luSlv(x, y, nLev-1)
        call luSlv(x, y, nLev)

        n = y

!        do i = nLev, 2, -1
!            n(i) = y(i-1)
!        end do
!        sumN = 1.d0
!        do i = 2, nLev
!            sumN = sumN + n(i)
!        end do
!        do i = 2, nLev
!            n(i) = n(i) / sumN
!        end do
!        n(1) = 1.d0 / sumN

        sumN = 0.d0
        do i = 1, nLev
           sumN = sumN+n(i)
        end do
        do i = 1, nLev           
           n(i) = n(i)/sumN
        end do

        ! now find emissivity and wavelengths
        do i = 1, nLev-1
            do j = i+1, nLev
                if (a(j,i) /= 0.d0) then
                    Eji = (e(j)-e(i))
                    fLineEm(i,j) = a(j,i) * Eji * n(j)
                    if (Eji>0.) then
                       if (present(wav)) wav(i,j) = 1.e8/Eji
                    else
                       if (present(wav)) wav(i,j) = 0.
                    end if
                end if
             end do
        end do   

        ! deallocate arrays
        if( associated(label) ) deallocate(label)
        if( associated(logTemp) ) deallocate(logTemp)
        if( associated(a) ) deallocate(a)
        if( associated(cs) ) deallocate(cs)
        if( associated(n) ) deallocate(n)
        if( associated(n2) ) deallocate(n2)
        if( associated(qeff) ) deallocate(qeff) 
        if( associated(qq) ) deallocate(qq)
        if( associated(qq2) ) deallocate(qq2)
        if( associated(tnij) ) deallocate(tnij) 
        if( associated(x) ) deallocate(x) 
        if( associated(y) ) deallocate(y) 
        if( associated(e) ) deallocate(e)
        if( associated(g) ) deallocate(g)
        if( associated(qom) ) deallocate(qom)

    end subroutine equilibrium

    ! this procedure performs the solution of linear equations
    subroutine luSlv(a, b, n)
        implicit none

        integer, intent(in)                  :: n 

        double precision,& 
&               intent(inout), dimension(:,:) :: a
        double precision,&
&               intent(inout), dimension(:)   :: b    

        call lured(a,n)

        call reslv(a,b,n)
    end subroutine luSlv

    subroutine lured(a,n)
        implicit none

        integer, intent(in)                  :: n

        double precision,&
&            intent(inout), dimension(:,:)    :: a

        ! local variables
        integer          :: i, j, k                    ! counters

        double precision :: factor                     ! general calculation factor

        if (n == 1) return

        do i = 1, n-1
            do k = i+1, n
                factor = a(k,i)/a(i,i)
                do j = i+1, n
                   a(k, j) = a(k, j) - a(i, j) * factor
                end do
            end do
        end do
    end subroutine lured

    subroutine reslv(a,b,n)
        implicit none 

        integer, intent(in)              :: n

        double precision,&
&            intent(inout), dimension(:,:) :: a
        double precision,&        
&            intent(inout), dimension(:)   :: b

        ! local variables
        integer    :: i, j, k, l              ! counters

        if (n == 1) then
            b(n) = b(n) / a(n,n)
            return
        end if

        do i = 1, n-1
            do j = i+1, n
                b(j) = b(j) - b(i)*a(j, i)/ a(i, i)
            end do
        end do
        b(n) = b(n) / a(n,n)
        do i = 1, n-1
            k = n-i
            l = k+1
            do j = l, n
                b(k) = b(k) - b(j)*a(k, j)
            end do
            b(k) = b(k) / a(k,k)
        end do
    end subroutine reslv

    subroutine initResLines(grid)
      implicit none
      
      type(grid_type), intent(inout) :: grid(*)

      character(len=2)   :: label
      character(len=120) :: reader ! file reader

      integer :: el  ! element number
      integer :: err ! allocation error status
      integer :: i,ii,j ! counters
      integer :: iGrid ! counter
      integer :: ios ! I/O error status
      integer :: nL  ! line counter
      integer :: nmul ! number of multiplets
      integer :: nResLinesFile
      integer :: safeLimit = 1e6 !

      open(unit=19,file="data/resLines.dat", status="old", position="rewind", iostat=ios)
      if (ios /= 0) then
         print*, "! initResLines: can't open file: data/resLines"
         stop
      end if
      
      nResLines =0
      do i = 1, safeLimit
         read(unit=19,fmt='(A2,A110,I4)',iostat=ios) label,reader, nmul
         if (ios/=0) exit
         do j = 1, nmul
            read(unit=19,fmt=*, iostat=ios) 
            if (ios/=0) then
               print*, "! initResLines: error reading data/resLines.dat file"
               stop
            end if
         end do
            
         ! check what element and if it's included
         el=0
         do j = 1, nElements
            if (trim(label) == trim(elemLabel(j)) ) then
               el=j
               exit
            end if
         end do
         if (el==0) then
            print*, "! initResLines: illegal label ", label, el
            stop
         end if
         if (lgElementOn(el)) nResLines=nResLines+1
         
      end do
      if (nResLines <= 0) then
         print*, '! initResLines: number of resonant lines is zero or negative'
         stop
      end if

      if (i>=safeLimit) then
         print*, "! initResLines : safe limit exceeded in resonant line list read"
         stop
      end if
      nResLinesFile=i-1
      rewind(19)

      ! allocate space for resLineLabel and resLineNuP arrays
      allocate(resLine(1:nResLines), stat = err)
      if (err /= 0) then
         print*, "! initResLines: can't allocate array memory, resLine"
         stop
      end if
      
      nL=1
      do ii = 1, nResLinesFile
         if (nL>nResLines) exit

         read(unit=19,fmt='(A2,I2,1x,A23,11x,A7,1x,F9.4,7x,F9.6,4x,F12.6,1x,I2,1x,I2,1x,E8.2,1x,E8.2,1x,I2)', &
              &iostat=ios) resLine(nL)%species, resLine(nL)%ion,&
              & resLine(nL)%transition, resLine(nL)%multiplet, resLine(nL)%wav, resLine(nL)%Elow, resLine(nL)%Ehigh,&
              & resLine(nL)%glow, &
              & resLine(nL)%ghigh, resLine(nL)%Aik, resLine(nL)%fik, resLine(nL)%nmul  
         if (ios/=0) then
            print*, '! initResline: error reading file -1- '
            stop
         end if

         allocate(resLine(nL)%moclow(1:resLine(nL)%nmul), stat=err)
         if (err /= 0) then
            print*, "! initResLines: can't allocate array memory, resLine(nL)%moclow"
            stop
         end if
         allocate(resLine(nL)%mochigh(1:resLine(nL)%nmul), stat=err)
         if (err /= 0) then
            print*, "! initResLines: can't allocate array memory, resLine(nL)%mochigh"
            stop
         end if

         do j = 1, resLine(nL)%nmul
            read(unit=19,fmt='(A112,1x,I2,1x,I2)', iostat=ios) reader, resLine(nL)%moclow(j), resLine(nL)%mochigh(j) 
            if (ios/=0) then
               print*, "! initResLines: error reading data/resLines.dat file - 2"
               stop
            end if
         end do

         resLine(nL)%elem=0
         do j = 1, nElements

            if (trim(resLine(nL)%species) == trim(elemLabel(j)))  then
               resLine(nL)%elem=j
               exit
            end if
         end do

         if (resLine(nL)%elem==0) then
            print*, "! initResLines: illegal label -2 -"
            stop
         end if

         if (lgElementOn(resLine(nL)%elem)) then         
      
            ! convert A to Ryd
            resLine(nL)%freq = c*1.e8/(resLine(nL)%wav*fr1Ryd)
      
            ! find the mass of the ion in [g]
            resLine(nL)%m_ion = aWeight(resLine(nL)%elem)*amu

            ! locate on nuArary
            call locate(nuArray, resLine(nL)%freq, resLine(nL)%nuP)
            if (resLine(nL)%nuP<=0 .or. resLine(nL)%nuP>=nbins) then
               print*, "! initResLines: wavelength of resonant line is outside grid range", &
                    & ii, resLine(nL)%nuP, resLine(nL)%freq, nuArray(1), nuArray(nbins)
               stop
            end if
            
            nL = nL +1
         end if

      end do
      close(19)

      if (taskid == 0) then
         print*, "! initResLines: The following resonance lines have been initialised and will &
              &be included in the radiative transfe"
         do i = 1, nResLines
            do j = 1, resLine(i)%nmul
               print*, i, resLine(i)%elem, resLine(i)%ion, resLine(i)%nmul, resLine(i)%wav
               print*, resLine(i)%moclow(j), resLine(i)%mochigh(j)
            end do
         end do
      end if

      do iGrid = 1, nGrids
         ! allocate fEscapeResPhotons
         allocate(grid(iGrid)%fEscapeResPhotons(0:grid(iGrid)%nCells,1:nResLines), stat = err)
         if (err /= 0) then
            print*, "! initResLines: can't allocate array memory, fEscapeResPhotons"
            stop
         end if
         grid(iGrid)%fEscapeResPhotons=0.

         ! and the extra packets to be transferred from each location
         allocate(grid(iGrid)%resLinePackets(0:grid(iGrid)%nCells), stat = err)
         if (err /= 0) then
            print*, "! initResLines: can't allocate array memory, resLinePackets"
            stop
         end if
         grid(iGrid)%resLinePackets=0.

      end do

      allocate(dustHeatingBudget(0:nAbComponents, 0:nResLines), stat=err)
      if (err /= 0) then
         print*, "! initResLines: can't allocate array memory, dustHeatingBudget"
         stop
      end if
      dustHeatingBudget=0.

    end subroutine initResLines

    subroutine setResLineEscapeProb(grids, gPin, xPin,yPin,zPin)
      implicit none

      type(grid_type), intent(inout) :: grids(*) ! incoming grids 

      type(vector) :: uHat, rVec ! direction and position vectors
     
      real    :: a0,DeltaNuD ! line centre xSec [cm^2] and doppler width
      real    :: bFrac
      real    :: deltax    !
      real    :: k_d, k_l  ! dust and line opacities
      real    :: dSx, dSy, dSz, dS ! distances from walls
      real    :: I1, I2    ! calculation integrals
      real    :: M2_tau    ! direction-dependent escape probability
      real    :: nuFrac    ! fractional nu for profile integration
      real    :: Pline     ! prob line
      real    :: radius    ! radial distance from the centre of the grid
      real    :: tau_mu    ! optical depth in direction mu at line centre
      real    :: gasdensity,gastemperature ! properties of the cells along the line of travel
      real    :: gasdensity0,gastemperature0 ! local cell properties
      real, dimension(nbins) :: radfield,phi 

      integer, intent(in) :: xPin,yPin,zPin,gPin ! incoming cell pointers
      integer             :: xp,yp,zp ! cell pointers
      integer             :: gPmother,xPmother,yPmother,zPmother
      integer             :: icount  ! counter
      integer             :: inu, idir, jdir, kdir ! freq counter
      integer, parameter  :: safeLimit=1e5 ! loop limit
      integer             :: isub,jsub,ksub,subTot,gP ! subgrid counters
      integer             :: cellP   ! cell pointer on active grid array 
      integer             :: compoP  ! pointer to gas component index 
      integer             :: iLine   ! res line counter
      integer             :: ndirx, ndiry,ndirz ! number of directions to sample tau
      
     
      xp=xPin
      yp=yPin
      zp=zPin
      gp=gPin
      
      cellP = grids(gPin)%active(xP,yP,zP)

      if (cellP<=0) return

      ndirx = 4
      ndiry = 4
      ndirz = 4
      
      ! get local physical properties of the gas

      gasdensity0 = grids(gPin)%Hden(grids(gPin)%active(xp,yp,zp))
      gastemperature0 = grids(gPin)%Te(grids(gPin)%active(xP,yP,zP))

      radfield = 0.
      do inu = 1, nbins
         if (lgDebug) then
            radfield(inu) = radfield(inu)+grids(gPin)%Jste(grids(gPin)%active(xP,yP,zP),inu)+&
                 & grids(gPin)%Jdif(grids(gPin)%active(xP,yP,zP),inu)
         else
            radfield(inu) = radfield(inu)+grids(gPin)%Jste(grids(gPin)%active(xP,yP,zP),inu)
         end if
      end do

      do iLine = 1, nResLines
         if (.not.lgElementOn(resLine(iLine)%elem)) then
            print*, "!setResLineEscapeProb: attempted to calculate escape probability of resonant line &
                 & from an absent ion"
            stop
         end if

         xp=xPin
         yp=yPin
         zp=zPin
         gp=gPin
      
         cellP = grids(gPin)%active(xP,yP,zP)
         compoP = grids(gPin)%abFileIndex(xP,yP,zP)

         Pline=0.

         ! calculate mean line opacity for this cell
         DeltaNuD = (resLine(iLine)%freq*fr1Ryd/c)*(2.76124e-16*gastemperature0/resLine(iLine)%m_ion)**0.5
!         call meanLineOpacity(resLine(iLine), DeltaNuD, a0)
         call lineCentreXSec(resLine(iLine), DeltaNuD, a0)
         a0 = a0/1.e-14

         ! find k_d and k_l
         k_d = grids(gPin)%absOpac(cellP,resLine(iLine)%nuP)
         k_l = a0* grids(gPin)%Hden(cellP) * grids(gPin)%elemAbun(compoP,resLine(iLine)%elem)* &
              &  grids(gPin)%ionDen(cellP,elementXref(resLine(iLine)%elem),resLine(iLine)%ion)

         ! calculate line profile and calculation integrals
         I1=0.
         I2=0.
         deltax=0.

         ! profiles are too narrow for the following

         if (resLine(iLine)%nuP>=nbins) then
            print*, '! setResLineEscapeProb: res line freq pointer >= nbins', iline, &
                 & resLine(iLine)%nuP, nbins
            stop
         end if

         do inu = 1, nbins

            ! calculate profile over the bin
            nuFrac = nuArray(resLine(iLine)%nuP)+&
                 &real(inu-1)*widFlx(resLine(iLine)%nuP)/real(nbins-1)
             phi(inu) = profile("doppler", DeltaNuD, resLine(iLine), nuFrac)
            ! I1=I1+phi(inu)
            ! I2=I2+phi(inu)*radfield(inu)            

         end do
!         I1 = widFlx(resLine(iLine)%nuP)
!         I2 = I1*radfield(resLine(iLine)%nuP)

!print*, 'i1 i2, nup', i1,i2, resLine(iLine)%nuP
         ! calculate the optical depths sampling over ndirx*ndiry*ndirz directions
         uHat%x = 1./sqrt(3.) 
         uHat%y = uHat%x 
         uHat%z = uHat%x 

         do idir = 0,ndirx
            uHat = rotateX(uHat,real(idir)*Pi/2.)
            do jdir = 0,ndiry
               uHat = rotateY(uHat,real(jdir)*Pi/2.)
               do kdir = 0,ndirz
                  uHat = rotateZ(uHat,real(kdir)*Pi/2.)

                  tau_mu = 0.

                  rVec%x = grids(gPin)%xAxis(xPin)
                  rVec%y = grids(gPin)%yAxis(yPin)
                  rVec%z = grids(gPin)%zAxis(zPin)         

                  xp = xpin
                  yp = ypin
                  zp = zpin
                  gp = gpin
      

                  do icount = 1, safeLimit

                     ! check if we are in a subgrid
                     if (grids(gPin)%active(xP,yP,zP)<0) then                

                        gPmother = gP
                        xPmother = xP
                        yPmother = yP
                        zPmother = zP                        
                        
                        gP = abs(grids(gP)%active(xP,yP,zP))

                        ! where are we in the sub-grid?
                        call locate(grids(gP)%xAxis, rVec%x, xP)
                        if (xP< grids(gP)%nx) then
                           if (rVec%x >=  (grids(gP)%xAxis(xP+1)+grids(gP)%xAxis(xP))/2.) &
                                & xP = xP + 1
                        end if

                        call locate(grids(gP)%yAxis, rVec%y, yP)
                        if (yP< grids(gP)%ny) then
                           if (rVec%y >=  (grids(gP)%yAxis(yP+1)+grids(gP)%yAxis(yP))/2.) &
                                & yP = yP + 1
                        end if
                        call locate(grids(gP)%zAxis, rVec%z, zP)
                        if (zP< grids(gP)%nz) then
                           if (rVec%z >=  (grids(gP)%zAxis(zP+1)+grids(gP)%zAxis(zP))/2.) &
                                & zP = zP + 1
                        end if

                     end if


                     cellP = grids(gPin)%active(xP,yP,zP)

                     if (cellP > 0 ) then
                        compoP = grids(gPin)%abFileIndex(xP,yP,zP)
                        if (xp == xpin .and. yp == ypin .and. zp == zpin .and. gP == gpin) then
                           gasdensity = gasdensity0
                           gastemperature=gastemperature0
                        else
                           gasdensity = grids(gP)%Hden(grids(gP)%active(xp,yp,zp))
                           gastemperature = grids(gP)%Te(grids(gP)%active(xP,yP,zP))
                        end if

                     else
                        
                        gasdensity=0.
                        gastemperature=1.

                     end if
             
                     ! find distances from all walls

                     if (uHat%x>0.) then
                        if (xP<grids(gP)%nx) then
                           dSx = ( (grids(gP)%xAxis(xP+1)+grids(gP)%xAxis(xP))/2.-rVec%x)/uHat%x
                           if (abs(dSx)<1.e-10) then
                              rVec%x=(grids(gP)%xAxis(xP+1)+grids(gP)%xAxis(xP))/2.
                              xP = xP+1
                           end if
                        else
                           dSx = ( grids(gP)%xAxis(grids(gP)%nx)-rVec%x)/uHat%x
                           if (abs(dSx)<1.e-10) then
                              rVec%x=grids(gP)%xAxis(grids(gP)%nx)
                           end if
                        end if
                     else if (uHat%x<0.) then
                        if (xP>1) then
                           dSx = ( (grids(gP)%xAxis(xP)+grids(gP)%xAxis(xP-1))/2.-rVec%x)/uHat%x
                           if (abs(dSx)<1.e-10) then             
                              rVec%x=(grids(gP)%xAxis(xP)+grids(gP)%xAxis(xP-1))/2.
                              xP = xP-1
                           end if
                        else
                           dSx = (grids(gP)%xAxis(1)-rVec%x)/uHat%x
                           if (abs(dSx)<1.e-10) then             
                              rVec%x=grids(gP)%xAxis(1)
                           end if
                        end if
                     else if (uHat%x==0.) then
                        dSx = grids(gP)%xAxis(grids(gP)%nx)
                     end if
                     
                     if (.not.lg1D) then 
                        if (uHat%y>0.) then
                           if (yP<grids(gP)%ny) then
                              dSy = ( (grids(gP)%yAxis(yP+1)+grids(gP)%yAxis(yP))/2.-rVec%y)/uHat%y
                              if (abs(dSy)<1.e-10) then
                                 rVec%y=(grids(gP)%yAxis(yP+1)+grids(gP)%yAxis(yP))/2.
                                 yP = yP+1
                              end if
                           else
                              dSy = (  grids(gP)%yAxis(grids(gP)%ny)-rVec%y)/uHat%y
                              if (abs(dSy)<1.e-10) then
                                 rVec%y=grids(gP)%yAxis(grids(gP)%ny)
                              end if
                           end if
                        else if (uHat%y<0.) then
                           if (yP>1) then
                              dSy = ( (grids(gP)%yAxis(yP)+grids(gP)%yAxis(yP-1))/2.-rVec%y)/uHat%y
                              if (abs(dSy)<1.e-10) then             
                                 rVec%y=(grids(gP)%yAxis(yP)+grids(gP)%yAxis(yP-1))/2.
                                 yP = yP-1
                              end if
                           else 
                              dSy = ( grids(gP)%yAxis(1)-rVec%y)/uHat%y
                              if (abs(dSy)<1.e-10) then             
                                 rVec%y=grids(gP)%yAxis(1)
                              end if
                           end if
                        else if (uHat%y==0.) then
                           dSy = grids(gP)%yAxis(grids(gP)%ny)
                        end if

                        if (uHat%z>0.) then
                           if (zP<grids(gP)%nz) then
                              dSz = ( (grids(gP)%zAxis(zP+1)+grids(gP)%zAxis(zP))/2.-rVec%z)/uHat%z
                              if (abs(dSz)<1.e-10) then
                                 rVec%z=(grids(gP)%zAxis(zP+1)+grids(gP)%zAxis(zP))/2.
                                 zP = zP+1
                              end if
                           else
                              dSz = ( grids(gP)%zAxis(grids(gP)%nz)-rVec%z)/uHat%z
                              if (abs(dSz)<1.e-10) then
                                 rVec%z=grids(gP)%zAxis(grids(gP)%nz)
                              end if
                           end if
                        else if (uHat%z<0.) then
                           if (zP>1) then             
                              dSz = ( (grids(gP)%zAxis(zP)+grids(gP)%zAxis(zP-1))/2.-rVec%z)/uHat%z
                              if (abs(dSz)<1.e-10) then             
                                 rVec%z=(grids(gP)%zAxis(zP)+grids(gP)%zAxis(zP-1))/2.
                                 zP = zP-1
                              end if
                           else
                              dSz = ( grids(gP)%zAxis(1)-rVec%z)/uHat%z
                              if (abs(dSz)<1.e-10) then             
                                 rVec%z=grids(gP)%zAxis(1)
                              end if
                           end if
                        else if (uHat%z==0.) then
                           dSz = grids(gP)%zAxis(grids(gP)%nz)
                        end if

                     else
                       
                        print*, '! setResLineEscapeProb: oneD option not yet implemented for res line transfer routine'
                        stop

                     end if
                     

                     ! cater for cells on cell wall
                     if ( abs(dSx)<1.e-10 ) dSx = grids(gP)%xAxis(grids(gP)%nx)
                     if ( abs(dSy)<1.e-10 ) dSy = grids(gP)%yAxis(grids(gP)%ny)
                     if ( abs(dSz)<1.e-10 ) dSz = grids(gP)%zAxis(grids(gP)%nz)
                     
                     ! find the nearest wall
                     dSx = abs(dSx)
                     dSy = abs(dSy)
                     dSz = abs(dSz)

                     if (dSx<=0.) then
                        print*, '! setResLineEscapeProb: [warning] dSx <= 0.'
                        dS = amin1(dSy, dSz)
                     else if (dSy<=0.) then
                        print*, '! setResLineEscapeProb: [warning] dSy <= 0.'
                        dS = amin1(dSx, dSz)
                     else if (dSz<=0.) then
                        print*, '! setResLineEscapeProb: [warning] dSz <= 0.'
                        dS = amin1(dSx, dSy)
                     else
                        dS = amin1(dSx,dSy,dSz)
                     end if

                     ! this should now never ever happen
                     if (dS <= 0.) then
                        print*, 'setResLineEscapeProb: dS <= 0'
                        stop
                     end if

                     ! update optical depth 
                     if (cellP>0) then
                        tau_mu  = tau_mu + a0*1.e-14*dS*gasdensity*grids(gP)%elemAbun(compoP,resLine(iLine)%elem)*&
                             &grids(gP)%ionDen(cellP,elementXref(resLine(iLine)%elem),resLine(iLine)%ion)
                     end if

                     ! update position
                     rVec = rVec + dS*uHat

                     ! keep track of where you are on mother grid
                     if (gP>1) then
                        if (uHat%x>0.) then
                           if ( xPmother < grids(gPmother)%nx ) then
                              if ( rVec%x >= (grids(gPmother)%xAxis(xPmother)+& 
                                   & grids(gPmother)%xAxis(xPmother+1))/2. ) then
                                 xPmother = xPmother+1
                              end if
                           else
                              if ( rVec%x > grids(gPmother)%xAxis(xPmother)) then
                                 print*, '! pathSegment: insanity occured at mother grid transfer (x axis +)', & 
                                      & rVec%x, gP, gPmother
                                 stop
                              end if
                           end if
                        else
                           if ( xPmother > 1 ) then
                              if ( rVec%x <= (grids(gPmother)%xAxis(xPmother-1)+& 
                                   & grids(gPmother)%xAxis(xPmother))/2. ) then
                                 xPmother = xPmother-1
                              end if
                           else
                              if (rVec%x < grids(gPmother)%xAxis(1)) then
                                 print*, '! pathSegment: insanity occured at mother grids transfer (x axis -)', & 
                                      & rVec%x, gP, gPmother
                                 stop
                              end if
                           end if
                        end if
                        if (uHat%y>0.) then
                           if (  yPmother < grids(gPmother)%ny ) then
                              if ( rVec%y >= (grids(gPmother)%yAxis( yPmother)+& 
                                   & grids(gPmother)%yAxis( yPmother+1))/2. ) then
                                 yPmother =  yPmother+1
                              end if
                           else
                              if ( rVec%y > grids(gPmother)%yAxis( yPmother)) then
                                 print*, '! pathSegment: insanity occured at mother grid transfer (y axis +)', & 
                                      & rVec%y, gP, gPmother
                                 stop
                              end if
                           end if
                        else
                           if (  yPmother > 1 ) then
                              if ( rVec%y <= (grids(gPmother)%yAxis( yPmother-1)+& 
                                   & grids(gPmother)%yAxis( yPmother))/2. ) then
                                 yPmother =  yPmother-1
                              end if
                           else
                              if (rVec%y < grids(gPmother)%yAxis(1)) then
                                 print*, '! pathSegment: insanity occured at mother grid transfer (y axis -)', & 
                                      & rVec%y, gP, gPmother
                                 stop
                              end if
                           end if
                        end if
                        if (uHat%z>0.) then
                           if (  zPmother < grids(gPmother)%nz ) then
                              if ( rVec%z >= (grids(gPmother)%zAxis( zPmother)+&
                                   & grids(gPmother)%zAxis( zPmother+1))/2. ) then
                                 zPmother =  zPmother+1
                              end if
                           else
                              if ( rVec%z > grids(gPmother)%zAxis( zPmother)) then
                                 print*, '! pathSegment: insanity occured at mother grid transfer (z axis +)', &
                                      & rVec%z, gP, gPmother
                                 stop
                              end if
                           end if
                        else
                           if (  zPmother > 1 ) then
                              if ( rVec%z <= (grids(gPmother)%zAxis( zPmother-1)+&
                                   & grids(gPmother)%zAxis( zPmother))/2. ) then
                                 zPmother =  zPmother-1
                              end if
                           else
                              if (rVec%z < grids(gPmother)%zAxis(1)) then
                                 print*, '! pathSegment: insanity occured at mother grid transfer (z axis -)', &
                                      & rVec%z, gP, gPmother
                                 stop
                              end if
                           end if
                        end if
                     end if

                     ! and update current grid indeces
                     if (.not.lg1D) then
                        if ( (dS == dSx) .and. (uHat%x > 0.)  ) then
                           xP = xP+1
                        else if ( (dS == dSx) .and. (uHat%x < 0.) ) then
                           xP = xP-1
                        else if ( (dS == dSy) .and. (uHat%y > 0.) ) then
                           yP = yP+1
                        else if ( (dS == dSy) .and. (uHat%y < 0.) ) then 
                           yP = yP-1
                        else if ( (dS == dSz) .and. (uHat%z > 0.) ) then
                           zP = zP+1
                        else if ( (dS == dSz) .and. (uHat%z < 0.) ) then
                           zP = zP-1
                        end if
                     else
                        
                        print*, '! setResLineEscapeProb: oneD option not yet implemented for res line transfer routine'
                        stop
                        
                        radius = 1.e10*sqrt((rVec%x/1.e10)*(rVec%x/1.e10) + &
                             & (rVec%y/1.e10)*(rVec%y/1.e10) + &
                             & (rVec%z/1.e10)*(rVec%z/1.e10))
                        call locate(grids(gP)%xAxis, radius , xP)

                     end if                    

                     if(lgPlaneIonization) then
                        
                        if ( rVec%y <= grids(gP)%yAxis(1)-grids(gP)%geoCorrY .or. yP<1) then
                           ! exit the loop if on mother or return to mother if on sub
                           if (gP==1) then	
                              yP=1
                              exit
                           else if (gP>1) then
                              xP = xPmother
                              yP = yPmother
                              zP = zPmother
                              gP = gPmother
                           else
                              print*, '! setResLineEscapeProb: insanity occured - invalid gP', gP
                              stop
                           end if

                        else if (rVec%y >= grids(gP)%yAxis(grids(gP)%ny)+grids(gP)%geoCorrY .or. yP>grids(gP)%ny) then
                           
                           if (gP==1) then	
                              ! exit the loop if on mother or return to mother if on sub
                              yP = grids(gP)%ny
                              exit
                           else if (gP>1) then
                              xP = xPmother
                              yP = yPmother
                              zP = zPmother
                              gP = gPmother                      
                           else
                              print*, '! setResLineEscapeProb: insanity occured - invalid gP', gP
                              stop
                           end if

                        else if (rVec%x <= grids(gP)%xAxis(1) .or. xP<1) then

                           if (gP==1) then	                       
!                         rVec%x = -rVec%x 
                              xP=1
                              rVec%x = grids(gP)%xAxis(1)
                              uHat%x = -uHat%x
                           else if (gP>1) then
                              xP = xPmother
                              yP = yPmother
                              zP = zPmother
                              gP = gPmother
                           else
                              print*, '! setResLineEscapeProb: insanity occured - invalid gP', gP
                              stop
                           end if
                           
                        else if (rVec%x >=  grids(gP)%xAxis(grids(gP)%nx) .or. xP>grids(gP)%nx)then

                           if (gP==1) then
!                         rVec%x = 2.*grid(gP)%xAxis(grid(gP)%nx)- rVec%x
                              xP = grids(gP)%nx
                              rVec%x = grids(gP)%xAxis(grids(gP)%nx)
                              uHat%x = -uHat%x
                           else if (gP>1) then
                              xP = xPmother
                              yP =  yPmother
                              zP =  zPmother
                              gP = gPmother
                           else
                              print*, '! setResLineEscapeProb: insanity occured - invalid gP', gP
                              stop
                           end if
                      
                        else if (rVec%z <= grids(gP)%yAxis(1) .or.zP<1) then

                           if (gP==1) then
!                         rVec%z = -rVec%z 
                              zP=1
                              rVec%z = grids(gP)%yAxis(1)
                              uHat%z = -uHat%z
                           else if (gP>1) then
                              xP = xPmother
                              yP =  yPmother
                              zP =  zPmother
                              gP = gPmother
                           else
                              print*, '! setResLineEscapeProb: insanity occured - invalid gP', gP
                              stop
                           end if
                      
                        else if (rVec%z >=  grids(gP)%zAxis(grids(gP)%nz) .or. zP>grids(gP)%nz)then

                           if (gP==1) then
                              zP = grids(gP)%nz
!                         rVec%z = 2.*grid(gP)%zAxis(grid(gP)%nz)- rVec%z
                              rVec%z = grids(gP)%zAxis(grids(gP)%nz)
                              uHat%z = -uHat%z
                           else if (gP>1) then
                              xP = xPmother
                              yP = yPmother
                              zP = zPmother
                              gP = gPmother
                           else
                              print*, '! setResLineEscapeProb: insanity occured - invalid gP', gP
                              stop
                           end if
                      
                        end if
                    
                     end if

                     ! check if the path is still within the ionized region 
                     radius = 1.e10*sqrt((rVec%x/1.e10)*(rVec%x/1.e10) + &
                          & (rVec%y/1.e10)*(rVec%y/1.e10) + &
                          & (rVec%z/1.e10)*(rVec%z/1.e10))

                     if ( radius >= R_out .and. R_out >= 0. ) exit

                     if (.not.lgPlaneIonization) then
                        
                        if ( (abs(rVec%x) >= grids(gP)%xAxis(grids(gP)%nx)+grids(gP)%geoCorrX) .or.&
                             &(abs(rVec%y) >= grids(gP)%yAxis(grids(gP)%ny)+grids(gP)%geoCorrY) .or.&
                             &(abs(rVec%z) >= grids(gP)%zAxis(grids(gP)%nz)+grids(gP)%geoCorrZ) .or. &
                             & xP>grids(gP)%nx .or. yP>grids(gP)%ny .or. zP>grids(gP)%nz  ) then


                           if ((gP==1) .or.  (radius >= R_out .and. R_out >= 0.)) then
                              if (xP > grids(gP)%nx) xP = grids(gP)%nx
                              if (yP > grids(gP)%ny) yP = grids(gP)%ny
                              if (zP > grids(gP)%nz) zP = grids(gP)%nz
                              exit
                           else if (gP>1) then
                              xP = xPmother
                              yP = yPmother
                              zP = zPmother
                              gP = gPmother
                           else
                              print*, '! setResLineEscapeProb: insanity occured - invalid gP - ', gP
                              stop
                           end if
                        end if
                     end if

                     if (lgSymmetricXYZ .and. gP == 1) then
                        if (lgPlaneIonization) then
                           print*, '! setResLineEscapeProb: lgSymmetric and lgPlaneionization flags both raised'
                           stop
                        end if
                        if ( rVec%x <= grids(gP)%xAxis(1) .or. xP<1) then
                           if (uHat%x<0.) uHat%x = -uHat%x 
                           xP = 1
                           rVec%x = grids(gP)%xAxis(1)
                        end if
                        if ( rVec%y <= grids(gP)%yAxis(1) .or. yP<1) then
                           if (uHat%y<0.) uHat%y = -uHat%y 
                           yP = 1
                           rVec%y = grids(gP)%yAxis(1)
                        end if
                        if ( rVec%z <= grids(gP)%zAxis(1) .or. zP<1) then
                           if (uHat%z<0.) uHat%z = -uHat%z 
                           zP=1
                           rVec%z = grids(gP)%zAxis(1)
                        end if
                     end if

                     if (gP>1) then

                        if ( ( (rVec%x <= grids(gP)%xAxis(1) .or. xP<1) .and. uHat%x <=0.) .or. & 
                             & ( (rVec%y <= grids(gP)%yAxis(1) .or. yP<1) .and. uHat%y <=0.) .or. & 
                             & ( (rVec%z <= grids(gP)%zAxis(1) .or. zP<1) .and. uHat%z <=0.) ) then
                           
                           ! go back to mother grid
                           xP = xPmother
                           yP = yPmother
                           zP = zPmother
                           gP = gPmother

                        end if
                        
                     end if               
                  end do
                  if (iCount>=safeLimit) then
                     print*, '! setResLineEscapeProb: safeLimit exceeded in tau loop', icount, safelimit
                     stop
                  end if

                  if (tau_mu <= 0.) &
                       & print*, '! setResLineEscapeProb: [warning] tau_mu <= 0. ', tau_mu, uHat, xPin,yPin,zPin, cellP, xp,yp,zp


                  ! calculate direction-dependent escape probability according to Harrington, Cohen and Hess (1984) eqn A1
!                  M2_tau = getM2_tau(DeltaNuD, resLine(iLine), tau_mu, phi)
!                  Pline = M2_tau
                  
                  ! calculate direction-dependent escape probability according to Hollenbach and McKee 1978
                  if (tau_mu>=1.) then
                     bFac = sqrt(log(tau_mu))/(1.+tau_mu/1.e5)
                     Pline = 1./(tau_mu*sqrt(Pi)*(1.2+bFac))
                  else
                     Pline=(1.-exp(-2.*tau_mu))/(2.*tau_mu)
                  end if

                  if (tau_mu <= 0.) then
                     deltax=0.
                  else
                     deltax = log(tau_mu/sqrt(Pi)) 
                     if (deltax >0.) then
                        deltax=sqrt(deltax)
                     else
                        deltax=0.
                     end if
                  end if
         

print*, tau_mu, deltanud, deltax, Pline
                  ! calculate fEscapeResPhotons according to Harrington, Cohen and Hess (1984) eqn A5
                  ! for direction mu!
                  ! f = 1 + k_d Int[J_nu, {nu, 0, Inf}]/ (P_l k_l Int[J_nu phi_nu, {nu 0, Inf}])         
                  if (Pline == 0.) then
                     grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),iLine) = 0.
                  else
                     grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),iLine) = &
                          & grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),iLine) + &
                          & 1./(1. + k_d*1.e+14*2.*deltax/(k_l*Pline))
                  end if

print*, idir, jdir, kdir, iline, pline, grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),iLine)
               end do
            end do
         end do

         grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),iLine) = &
              & grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),iLine)/&
              & real((ndirx+1)*(ndiry+1)*(ndirz+1))

print*, xPin,yPin,zPin,iLine, grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),iLine), k_d,k_l
      end do


    end subroutine setResLineEscapeProb

    ! calculate direction-dependent escape probability according to Harrington, Cohen and Hess (1984) eqn A1
    function getM2_tau(width,lineIn,tauIn, profilein)
      implicit none

      type(resLine_type), intent(in) :: lineIn

      real, intent(in)              :: width, tauIn
      real, intent(in), optional    :: profileIn(nbins)
      real                          :: getM2_tau
      real                          :: dx

      integer                       :: i
      
      if (tauIn<0.) then
         print*, "! getM2_tau: negative tau", tauIn
         stop
      end if         

      if (present(profilein)) then
         getM2_tau=0.
         do i = 1, nbins
            dx = widFlx(lineIn%nuP)*fr1ryd/(real(nbins-1)*width)
            getM2_tau = getM2_tau+profileIn(i)*exp(-tauIn*profileIn(i))*dx
         end do
      else
         getM2_tau = exp(-tauIn)
      end if

    end function getM2_tau

    function profile(profileID, widthIn, lineIn, nuIn)
      implicit none
      
      type(resLine_type), intent(in) :: lineIn
      real              :: profile ! phi(nu)
      real, intent(in)  :: nuIn ! frequency [Ryd]
      real, intent(in)  :: widthIn ! 
      real              :: m_ion

      character(len=7), intent(in) :: profileID ! what type of broadening?

      select case (profileID)
         case('doppler')
            profile = (1./(sqrt(Pi)*widthIn) ) * Exp(- ((nuIn-lineIn%freq)*fr1Ryd/widthIn)**2.)
         case default
            print*, '! profile: unknown profile requested ', profileID
            stop
      end select

    end function profile

    ! use line centre xsec formula 
    ! a0 = (Pi e^2 / (m_e c)) f /( (Pi^1/2 DeltaNuD)) 
    ! DeltaNuD = Doppler width = Nu0 /c * (2 k T /m_ion)^1/2
    ! 
    subroutine lineCentreXSec(line, width, line_xSec)
      implicit none
      
      type(resLine_type), intent(inout) :: line

      real, intent(in)         :: width
      real, intent(out)        :: line_xSec
      real                     :: m_ion          ! mass of the ion
      

      if (width<=0.) then
         print*, '! lineCentreXSec: width less or equal to zero'
         stop
      end if
      line_xSec = 14.969e-3*line%fik/width
            
    end subroutine lineCentreXSec


    ! use hummer & kunasz eqn 2.3 : k_l = N_k B_kl h v_kl/(4 Pi DeltaNuP)
    ! and B_kl = e^2 f_kl / (4 E_0 m_e c h n_kl)
    ! hence k_l = N_k e^2 f_kl / (4 E_0 m_e c 4 Pi DeltaNuP)
    ! E_0 is permissivity of free space
    ! a0 = (Pi e^2 / (m_e c)) f /( (Pi^1/2 DeltaNuD)) 
    ! DeltaNuD = Doppler width = Nu0 /c * (2 k T /m_ion)^1/2
    ! 
    subroutine meanLineOpacity(line, width, line_xSec)
      implicit none
      
      type(resLine_type), intent(inout) :: line

      real, intent(in)         :: width
      real, intent(out)        :: line_xSec
      real                     :: m_ion          ! mass of the ion
      

      if (width<=0.) then
         print*, '! lineCentreXSec: width less or equal to zero'
         stop
      end if
      line_xSec = 1.681e-4*line%fik/width
            
    end subroutine meanLineOpacity
    
  end module emission_mod


