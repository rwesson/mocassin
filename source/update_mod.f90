! Copyright (C) 2005 Barbara Ercolano 
!
! Version 2.02
module update_mod
    use constants_mod
    use common_mod
    use composition_mod
    use emission_mod
    use grid_mod
    use interpolation_mod
    use xSec_mod

    contains

      subroutine updateCell(grid, xP, yP, zP)
        implicit none 

        type(grid_type), intent(inout) :: grid         ! the grid
        
        real, parameter :: Y0 = 0.5, Y1 = 0.2          ! see Baldwin et al. 1991

        real                           :: aFit, bFit   ! general fit terms
        real                           :: deltaXHI        ! delta ionDen of H0                   
        real, dimension(3)             :: ions         ! # of ionizations (1 is from
                                                       ! H0 2 is from He0 and 3 is 
                                                       ! is for He+)
        real                           :: phXSecM
        real                           :: phXSecM1
        real                           :: phXSec      ! ph xSec of ion M at freqency bin j
        real                           :: term1,term2,&! general calculation terms 
             & term3        
        real                           :: thResidual   ! thermal balance residual
        real                           :: thResidualHigh! ther bal residual-high lim
        real                           :: thResidualLow! ther bal residual-low lim
        real                           :: Thigh        ! T high limit [K]
        real                           :: Tlow         ! T low limit  [K]
        real                           :: XOldHI       ! old ionDen of H0 at this cell 

        ! dust-gas interaction heating and cooling process
        real                           :: grainEmi, grainRec ! grain emissions and recom          
        real, pointer                  :: gasDustColl_d(:,:) ! cooling/heating of the 
                                                             ! dust through collisions with grains        
        real                           :: gasDustColl_g ! cooling/heating of the 
                                                       ! gas through collisions with grains
        real,pointer                   :: photoelHeat_d(:,:) ! cooling of dust by photoelectric
                                                             ! emission 
        real                           :: photoelHeat_g ! heating of gas by dust photoelctric 
                                                       !emission 
        real,pointer                   :: grainPot(:,:)  ! [Ryd]        

        real, parameter                :: hcRyd_k = &  ! constant: h*cRyd/k (Ryd at inf used) [K]
             & 157893.94   
        real, parameter                :: hcRyd = &    ! constant: h*c*Ryd (Ryd at inf used) [erg]
             & 2.1799153e-11
        real, parameter                :: thLimit = 0.01! convergence limit T-iteration
        
        real, dimension(nElements, nstages) &
             & :: alphaTot     ! total recombination coeffs

        real, dimension(nElements, nStages) &
             & :: nPhotoSte, & ! # of stellar photoionizations
             & nPhotoDif    ! # of diffuse photoionizations



        integer, intent(in)            :: xP, yP, zP   ! cell indexes on the Cartesian axes

        ! local variables
        logical                        :: lgHit        ! has this cell been hit by a photon?

        integer ,pointer               :: grainPotP(:,:) 
        integer                        :: cellP        ! points to this cell
        integer                        :: HIPnuP       ! pointer to H IP in NuArray
        integer                        :: HeIPnuP      ! pointer to HeI IP in NuArray
        integer                        :: HeIIPnuP     ! pointer to HeII IP in NuArray
        integer                        :: highNuP      ! pointer to highest energy of shell
        integer                        :: IPnuP        ! pointer to this ion's IP in NuArray
        integer                        :: xSecP        ! pointer to an ion's xSec in xSecArray
        integer                        :: elem         ! element counter
        integer                        :: err          ! allocation error status
        integer                        :: g0, g1       ! stat weights
        integer                        :: ion          ! ionization stage counter
        integer                        :: ios          ! I/O error status
        integer                        :: i,j          ! counter
        integer                        :: ncutoff= 0   ! see clrate
        integer                        :: nElec        ! # of electrons in ion
        integer                        :: nIterateGC   ! # of GC iterations
        integer                        :: nIterateT    ! # of T iterations
        integer                        :: nIterateX    ! # of X iterations
        integer                        :: outShell     ! outer shell number (1 for k shell)
        integer                        :: nspU
        integer, parameter             ::  nTbins=300  ! number of enthalpy bins for T spike
        integer, parameter             :: maxIterateGC&! limit to number of grain charge it
             & = 100 
        integer, parameter             :: maxIterateX& ! limit to number of X-iterat ions
             & = 25 
        integer, parameter             :: maxIterateT& ! limit to number of T-iterations
             & = 25 


        ! check whether this cell is outside the nebula
        if (grid%active(xP, yP, zP)<=0) return

        cellP = grid%active(xP, yP, zP)

        if (lgMultiDustChemistry) then
           nspU = grid%dustAbunIndex(cellP)
        else
           nspU = 1
        end if

        ! initialise lgBlack
        grid%lgBlack(cellP) = 0

        ! initialise lgHit
        lgHit = .false. 

        ! find out if this cell has been hit by at least one photon
        if (.not.lgDebug) then
           do i = 1, nbins
              if ( grid%Jste(cellP,i) > 0.) then
                 lgHit = .true.
                 exit
              end if
           end do
        else 
           do i = 1, nbins
              if ( (grid%Jste(cellP,i) > 0.) .or. &
                   & (grid%Jdif(cellP,i) > 0.)) then
                 lgHit = .true.
                 exit
              end if
           end do
        end if

        ! do not update this cell if there were no hits
        if (.not.lgHit) then

            if (lgTalk) print*, "! updateCell [talk]: no photon hits, &
                 &returning...", xP,yP,zP
            
            TdustTemp(:,:,cellP)       = grid%Tdust(:,:,cellP)
            
            if (lgGas) then
               ! the grid values stay the same
               TeTemp(cellP)         = grid%Te(cellP)
               NeTemp(cellP)         = grid%Ne(cellP)
               ionDenTemp(cellP,:,:) = grid%ionDen(cellP,:,:)
            end if

            grid%noHit = grid%noHit+1.

            return
        end if        

        if (lgGas) then

           ! initialize local Te and Ne

           TeUsed     = grid%Te(cellP)
           NeUsed     = grid%Ne(cellP)
           ionDenUsed = grid%ionDen(cellP, :, :)


           ! save present value of H0 abundance in XOldHI
           XOldHI = grid%ionDen(cellP,elementXref(1),1)

           nPhotoSte = 1.e-20
           if (lgDebug) nPhotoDif = 1.e-20


           ! calculate the number of stellar and diffuse photoionizations for each species
           ! NOTE: no need to time by the frequency bin width (widflx(j)) because JSte and JDif 
           !       were calculated for each individual bin 
           do elem = 1, nElements  ! begin element loop
            do ion = 1, min(elem, nStages-1) ! begin ion loop
               if(.not.lgElementOn(elem)) exit
               
               if (elem > 2) then

                  ! find the number of electrons in this ion
                  nElec = elem - ion + 1
                  
                  ! find the outer shell number and statistical weights 
                  call getOuterShell(elem, nELec, outShell, g0, g1)
                  
                  ! get pointer to this ion's IP in NuArray
                  IPnuP = elementP(elem, ion, outShell, 1)
                  
                  ! get pointer to this cell's highest energy
                  highNuP = elementP(elem, ion, outShell, 2) - 1
                  
               else if (elem == 1) then ! HI
                   
                  IPNuP   = HlevNuP(1)
                  highNuP = nbins

               else if ( (elem == 2) .and. (ion == 1) ) then ! HeI
 
                  IPNuP   = HeIlevNuP(1)
                  highNuP = nbins 

               else if ( (elem == 2) .and. (ion == 2) ) then ! HeII

                  IPNuP   = HeIIlevNuP(1)
                  highNuP = nbins
                    
               end if

               do j = IPnuP, highNuP ! begin frequency loop
                  if ( elem > 2 ) then
                     
                     ! get pointer to phot xSec of this ion in xSecArray
                     xSecP = elementP(elem, ion, outShell, 3)
                     
                     ! get phot xSec of ion at this energy
                     phXSec = xSecArray(j+xSecP-IPnuP+1-1)
                     
                  else if ( elem == 1 ) then ! HI

                     phXSec  = xSecArray(j-IPnuP+1+HlevXSecP(1)-1)

                  else if ( (elem == 2) .and. (ion == 1) ) then ! HeI

                     phXSec  = xSecArray(j-IPnuP+1+HeISingXSecP(1)-1)
      
                  else if ( (elem == 2) .and. (ion == 2) ) then ! HeII

                     phXSec  = xSecArray(j-IPnuP+1+HeIIXSecP(1)-1)
 
                  end if

                  if ((phXSec < 1.e-35) ) then
!                     print*, "! updateCell: [warning] bad xSec value &
!                          & exiting loop (j, phXSec)",&
!                          & nuArray(j), phXSec
!                     exit
                     phXSec = 0.
                 end if

                  if ((nuArray(j) < 1.e-35) ) then
                     print*, "! updateCell: insane nu value (j, nuArray(j)", &
                          & j, nuArray(j)
                     stop
                  end if

                  ! stellar photoionizations 
                  
                  if ( grid%Jste(cellP,j) > 0. )then
                       
                     nPhotoSte(elem,ion) = nPhotoSte(elem,ion) + &
                          & grid%JSte(cellP,j)*phXSec/(hcRyd*nuArray(j))


                  end if
                    
                  ! diffuse photoionizations
                    
                  if (lgDebug) then
                     if ( grid%Jdif(cellP,j) > 0. ) then
                        nPhotoDif(elem,ion) = nPhotoDif(elem,ion) + &
                             & grid%JDif(cellP,j)*phXSec/(hcRyd*nuArray(j))

                     end if
                  end if
                  
               end do ! end frequency loop
            end do ! end ion loop
         end do ! end element loop

         ! zero out T-iteration components
         nIterateT      = 0
         thResidual     = 0.
         thResidualHigh = 0.
         thResidualLow  = 0.
         Thigh          = 0.
         Tlow           = 0.
         
         if (lgDust .and. lgGas .and. lgPhotoelectric) then
            allocate (grainPot(1:nSPecies, 1:nsizes))
            grainPot=0.
            allocate (grainPotP(1:nSPecies, 1:nsizes))
            grainPotP=0
            allocate (photoelHeat_d(1:nSPecies, 1:nsizes))
            photoelHeat_d=0.
            allocate ( gasDustColl_d(1:nSPecies, 1:nsizes))
            gasDustColl_d=0.
         end if

         ! call the recursive iteration procedure for the gas
         call iterateT()

         ! call the non recursive itaration procedure for the dust
         if (lgDust) call getDustT()

         ! this was added to help implementing MPI comunication
         if (lgDust) TdustTemp(:,:,cellP)          = grid%Tdust(:,:,cellP)

         if (lgNeInput .and. lgTalk) then
            print*, '! updateCell: [talk] cell:', xP,yP,zP, " NeInput: ", &
                 &grid%NeInput(cellP), " NeUsed: ", &
                 &grid%Ne(cellP),  " N_gas: ", &
                 &grid%Hden(cellP)
         end if

      else

         if (.not.lgDust) then
            print*, '! updateCell: no gas or dust present. the grid is empty.'
            stop
         end if

         XOldHI = grid%Tdust(0,0,cellP)

         call getDustT()

         ! determine if the model has converged at this cell
         ! NOT  variable names refer to gas phase for reasons of laziness 
         deltaXHI = (grid%Tdust(0,0,cellP) - XOldHI) / XOldHI
         if ( abs(deltaXHI) <= XHILimit ) then
            grid%lgConverged(cellP) = 1
         else
            grid%lgConverged(cellP) = 0
         end if

         if (lgTalk) &
              & print*, "updateCell: [talk] cell", xP,yP,zP, "; converged?",&
              & grid%lgConverged(cellP), "; mean dust T: ", &
              & grid%Tdust(0,0,cellP), "; &
              & mean dust T old: ", XOldHI,"; dT(dust): ", deltaXHI

         ! this was added to help implementing MPI comunication
         TdustTemp(:,:,cellP)          = grid%Tdust(:,:,cellP)
         
      end if

    contains

        recursive subroutine iterateT()
            implicit none

            ! local variables
            integer                       :: i               ! counter
            integer                       :: isp, ai         ! counters
            integer                       :: ios             ! I/O error status
            integer                       :: n               ! stopper
            integer                       :: ispbig          ! counter

            logical                       :: lgGCBConv       ! grain charge converged?
            logical                       :: lgIBConv        ! converged?

            real                          :: a, b, b2        ! calculation coefficients
            real                          :: coolInt         ! tot cooling integral [erg/s/Hden] 
            real                          :: disc            ! discriminant
            real                          :: expFact         ! general exponential factor 
            real                          :: gamma           ! He+ density
            real                          :: heatInt         ! tot heating integral [erg/s/Hden]
            real                          :: root            ! root for H and He calculations

            real, parameter               :: Xmax = 0.9999   ! max rel abundance

            ! step up T-iteration
            nIterateT = nIterateT + 1

            ! calculate the recombination coefficients for this cell
            call calcalpha()

            ! initialize X iteration counter
            nIterateX = 1

            ! calculate the ion abundances for the remaning ions
            call ionBalance(lgIBConv)

            ! calculate dust-gas interaction heating/cooling terms
            if (lgDust .and. lgGas .and. lgPhotoelectric) then
               do isp = 1, nSpeciesPart(nspU)
                  do ai = 1, nsizes
                     nIterateGC     = 0
                     grainEmi       = 0.
                     grainRec       = 0.        

                     ispbig = isp+dustComPoint(nspU)-1
                     call setGrainPotential(ispbig,ai,lgGCBConv)

                     call locate(nuArray, grainPot(isp,ai), grainPotP(isp,ai))                     
                     if (grainPotP(isp,ai) == 0) grainPotP = 1

                  end do
               end do
               call setPhotoelHeatCool()


               call setDustGasCollHeatCool()

            end if


            ! solve thermal balance
            call thermBalance(heatInt, coolInt)

            ! Now calculate new Te

            ! calculate thResidual = heatInt - coolInt
            thResidual = heatInt - coolInt
!print*, cellp, thresidual, heatint, coolint
            if ( abs(thResidual) >= thLimit*heatInt ) then ! start convergence condition
            
                if ( nIterateT < maxIterateT ) then ! start nIterateT condition
  
                    if ( thResidual < 0) then

                        thResidualLow = thResidual                       
                        Tlow          = TeUsed
                        
                        if ( Thigh /= 0. ) then
                       
                            TeUsed = Tlow-(Thigh-Tlow)*thResidualLow/&
                                 & (thResidualHigh-thResidualLow)

                        else

                            TeUsed = TeUsed/1.2
 
                        end if

                    else

                        thResidualHigh = thResidual
                        Thigh          = TeUsed
                      
                        if ( Tlow /= 0. ) then
 
                            TeUsed = Tlow-(Thigh-Tlow)*thResidualLow/&
                                 & (thResidualHigh-thResidualLow)

                        else
 
                            TeUsed = TeUsed*1.2

                        end if

                    end if


                    if (nIterateT >= 6) then
                        if (abs(Thigh-Tlow) <= (0.002*(Thigh+Tlow)) ) then

                            if ( thResidual < 0 ) then

                                TeUsed = Thigh - 100.
                                 Thigh  = 0.

                            else

                                TeUsed = Tlow+100.
                                Tlow   = 0.

                            end if
                        end if 
                    end if

                    if (TeUsed <= 0.) TeUsed = 1.

                    ! next T-iteration
                    call iterateT()
                    return                    

                else ! if nIterateT > maxIterateT
               
                    ! after maxIterateT number of iterations the temperature 
                    ! is determined as follows

                    if ( Thigh == 0. ) then
 
                        TeUsed = Tlow
  
                    else
                        
                        if (Tlow == 0.) then

                            TeUsed = Thigh

                        else 

                            TeUsed = (Thigh+Tlow)*0.5
 
                        end if

                    end if

                    
                    print*, "! iterateT: [warning] no convergence after ", &
                         & nIterateT, " steps. (cell, T)", cellP,xP,yP,zP,TeUsed

                    grid%noTeBal = grid%noTeBal+1.
                    
                    grid%lgBlack(cellP) = 1

                end if

            else ! if thResidual < thLimit 
                          
                ! the T-iteration has converged

                if (lgTalk) print*, "! iterateT: [talk] convergence achieved after ",&
                     & nIterateT, " steps. (cell, T)", xP,yP,zP,TeUsed, grid%abFileIndex(xP,yP,zP)

            end if ! end convergence condition

            if (.not.lgIBConv) then
               grid%noIonBal = grid%noIonBal+1
               grid%lgBlack(cellP) = 1
            end if

            grid%Te(cellP)         = TeUsed          

            ! determine if the model has converged at this cell
            deltaXHI = (grid%ionDen(cellP,elementXref(1),1) - XOldHI) / XOldHI
            if ( abs(deltaXHI) <= XHILimit ) then 
               grid%lgConverged(cellP) = 1
            else
               grid%lgConverged(cellP) = 0
            end if

            if (lgTalk) print*, "iterateT: [talk] cell", xP,yP,zP, "; converged?",&
                 & grid%lgConverged(cellP), "; X(H0): ", &
                 & grid%ionDen(cellP,elementXref(1),1), &
                 & "; X(H0) old: ", XOldHI,"; dX(H0): ", deltaXHI    

            ! this was added to help implementing MPI comunication
            TeTemp(cellP)          = TeUsed

        end subroutine iterateT

        recursive subroutine setGrainPotential(iSp, ai, lgGCBConv)
          implicit none 

          real,save            :: delta0, delta,delta1
          real,save            :: grainPotOld ! local copy of grai pot
          real                 :: grainEmi, grainRec ! grain emissions and recom          
          real,save            :: grainEmiOld, grainRecOld ! grain emissions and recom          
          real,parameter       :: errorLim = 0.005, dm = 0.05 ! loop convergence
          real,parameter       :: safeLim = 100 ! loop safety limit
          real                 :: threshold,fac
          real,save :: dlow,dhigh,grainpotlow,grainpothigh,dVg,slope
          
          integer, intent(in)  :: iSp, ai

          logical, intent(inout) :: lgGCBConv   ! converged?
                    
          nIterateGC = nIterateGC+1         

          if (nIterateGC==1) then 
             dVg=0.05
             grainPot(isp,ai) = 0.             
             grainPotOld = grainPot(isp,ai)
             lgGCBConv = .true.
             grainEmiOld = getGrainEmission(grainPot(isp,ai), isp, ai)
             grainRecOld = getGrainRecombination(grainPot(isp,ai), isp)
             grainPot(isp,ai) = grainPotOld+0.05             
          end if

          threshold = max(grainVn(isp)+grainPot(isp,ai),grainVn(isp))
          grainEmi = getGrainEmission(grainPot(isp,ai), isp, ai)
          grainRec = getGrainRecombination(grainPot(isp,ai), isp)
          delta = grainEmi-grainRec

          delta1 = abs(delta/(0.5*max(1.e-35, grainEmi+grainRec)))


          ! check for convergence
          if (delta1<errorLim) then
             return
          else 
             
             if (grainPot(iSp,ai) /= grainPotOld) then
                fac = (grainEmi-grainEmiOld)-(grainRec-grainRecOld)
                if (fac/=0.) slope = fac/(grainPot(iSp,ai)-grainPotOld)
             end if

             grainPotOld=grainPot(iSp,ai)
             grainRecOld = grainRec
             grainEmiOld = grainEmi

             delta1 = -delta/slope
             delta = abs(delta1)
             delta = min(delta, dm*threshold)
             delta = sign(delta, delta1)

             grainPot(iSp,ai) = grainPot(iSp,ai)+delta

             if (nIterateGC< maxIterateGC) then
                call setGrainPotential(isp,ai,lgGCBConv)
                return
             else

!                print*, '! setGrainPotential: no convergence', &
!                     & cellP,grainPot(isp,ai),grainEmi,grainRec
                lgGCBConv=.false.                 
                return
             end if
             
          end if

        end subroutine setGrainPotential

        ! calculate the grain recombination 
        ! using eqn 18 etc of Baldwin et al. (1991)
        ! cellFactor is dependant on the physical conditions of the gas, 
        function getGrainEmission(Vg, isp,ai)
          implicit none


          real :: getGrainEmission
          real, intent(in) :: Vg ! grain potential
          real :: Yn, Yhat

          real :: Qa ! grain absorption efficiency
          real :: thres, photFlux
          
          integer, intent(in) :: isp, ai
          integer :: ifreq, ip

          getGrainEmission=0.


          ! get the threshold
          thres = max(grainVn(isp)+Vg,grainVn(isp))

          call locate( nuArray, thres, ip)           
          ip = ip+1


          do ifreq = ip, nbins


             Yn = min(Y0*(1.-grainVn(isp)/nuArray(ifreq)), Y1) 

             Yhat = Yn*min(1., max(0.,1.-Vg/(nuArray(ifreq)-grainVn(isp))))

             if (.not. lgDebug) then
!                photFlux = (grid%JPEots(cellP,ifreq) + &
!                & grid%Jste(cellP,ifreq))/(hcRyd*nuArray(ifreq))
                photFlux =  grid%Jste(cellP,ifreq)/(hcRyd*nuArray(ifreq))
             else
!                photFlux = (grid%JPEots(cellP,ifreq) +&
!                &  grid%Jste(cellP,ifreq)+grid%Jdif(cellP,ifreq))/& 
!                     & (hcRyd*nuArray(ifreq))
             end if

             photFlux = photFlux*fourPi             

             Qa = XSecArray(dustAbsXsecP(isp,ai)+ifreq-1)/(1.e-8*grainRadius(ai)**2.)

             getGrainEmission = getGrainEmission+Yhat*photFlux*Qa

          end do

        end function getGrainEmission


        ! calculate the grain recombination 
        ! using eqn 23 etcof Baldwin et al. (1991)
        ! cellFactor is dependant on the physical conditions of the gas, 
        function getGrainRecombination(Vg, isp)
          implicit none

          real :: getGrainRecombination
          real, intent(in) :: Vg ! grain potential [ryd]
          
          real :: eta   ! Coulomb correction
          real :: cpDen ! colliding particle number density [cm^-3]
          real :: eightkT_pi ! 8*k * Te/Pi [erg]
          real :: mcp   ! mass of colliding particle in [g]
          real :: kT    ! k*Te [ryd]
          real :: S     ! sticking coefficient 
          real :: vmean ! colliding particle mean velocity
          real :: Z     ! colliding particle charge
          
          integer, intent(in) :: isp ! species identifier 
          
          integer :: istage, ielem
          
          getGrainRecombination = 0.
          eta = 0.
          kT =6.336e-6*TeUsed
          eightkT_pi = (1.1045e-15)*TeUsed/Pi
          
          ! add e- collisions contributions
          
          vmean = sqrt(eightkT_pi/me)

          S  =1. ! electron sticking probability


          ! eq 24 of Baldwin et al 91
          eta = -Vg/kT

          if (eta <= 0.) then
             eta = 1.-eta
          else if (eta >0.) then
             eta = exp(-eta)
          else
             print*, "! getGrainRecombination: insane eta for e-", eta
          end if
             
          getGrainRecombination = getGrainRecombination + &
               & NeUsed*vmean*S*eta

          ! add contribution from all other neutral and ionic species
          do elem = 1, nElements
             if (lgElementOn(elem)) then                
                do istage = 2, min(elem+1,nstages)
                   ! get cpDen
                   cpDen = grid%ionDen(cellP,elementXref(elem),istage)*&
                        & grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)*& 
                        & grid%Hden(cellP)
                   mcp = (aWeight(elem)*amu)
                   vmean = sqrt(eightkT_pi/mcp)

                   S = 1.
                   Z = real(istage-1)
                   
                   eta = Z*Vg/kT

                   if (eta <= 0.) then
                      eta = 1.-eta
                   else if (eta >0.) then
                      eta = exp(-eta)
                   else
                      print*, "! getGrainRecombination: insane eta", & 
                           & eta, elem, istage
                   end if
                   
                   getGrainRecombination = getGrainRecombination - &
                        & cpDen*vmean*S*eta

                end do
             end if
          end do          

        end function getGrainRecombination
          
        ! see Baldwin et al 1991; but beware that we are resolving the size distribution.
        ! so must keep size dependance, so Qa ia really Pi a^2 Qa
        subroutine setPhotoelHeatCool()
          implicit none
          
          real    :: Qa, photFlux, EY, Yhat,Yn,th

          integer :: ns,na,ifreq,thP

          ! calculate the cooling of dust by photoelectric emission 
          !  and heating of gas
          ! Baldwin et al. 1991 eqn 25-27

          photoelHeat_d=0.
          photoelHeat_g=0.          
          do ns = 1, nSpeciesPart(nspU)
             do na = 1, nSizes


                if (grainPotP(ns,na) <= 0) then
                   print*, "! setPhotoelHeatCool irregular grain potential index", & 
                        grainPotP(ns,na), grainPot(ns,na)
                   stop
                end if

                th = max(grainVn(ns)+grainPot(ns,na), grainVn(ns))
                call locate(nuArray,th,thP)
                if(thP<=0) thP=1

                do ifreq = thP, nbins                

!                do ifreq = grainPotP(ns,na), nbins

                   Yn = min(Y0*(1.-grainVn(ns)/nuArray(ifreq)), Y1) 


                   Yhat = Yn*min(1., max(0.,1.-grainPot(ns,na)/(nuArray(ifreq)-grainVn(ns))))

                   if (.not. lgDebug) then
!                      photFlux = (grid%JPEots(cellP,ifreq) + & 
!                           & grid%Jste(cellP,ifreq))/(hcRyd*nuArray(ifreq))
                      photFlux = grid%Jste(cellP,ifreq)/(hcRyd*nuArray(ifreq))
                   else
!                      photFlux = (grid%JPEots(cellP,ifreq) + grid%Jste(cellP,ifreq)+& 
!                           & grid%Jdif(cellP,ifreq))/(hcRyd*nuArray(ifreq))
                      photFlux = (grid%Jste(cellP,ifreq)+grid%Jdif(cellP,ifreq))/(hcRyd*nuArray(ifreq))

                   end if                   

                   photFlux = photFLux

                   Qa = XSecArray(dustAbsXsecP(ns,na)+ifreq-1)


                   EY = Yn*0.5* min(nuArray(ifreq)-grainVn(ns),&
                        & max(0., ((nuArray(ifreq)-grainVn(ns))**2.-grainPot(ns,na)**2.)/& 
                        & (nuArray(ifreq)-grainVn(ns))))

!                   photoelHeat_d(ns,na) = photoelHeat_d(ns,na)+Qa*photFlux*EY*hcRyd/Pi
                   photoelHeat_d(ns,na) = photoelHeat_d(ns,na)+Qa*photFlux*EY*hcRyd/Pi

                   photoelHeat_g = photoelHeat_g+Qa*photFlux*(EY-Yhat*grainPot(ns,na))*& 
                        & grainAbun(nspU,ns)*grainWeight(na)

                end do
             end do
          end do

          photoelHeat_g = photoelHeat_g*grid%Ndust(cellP)*hcRyd/grid%Hden(cellP)

        end subroutine setPhotoelHeatCool

        ! sets the cooling and heating rates of gas and dust due to collisions between the two phases        
        ! see Baldwin et al 1991
        ! only collisions with up to the 3 times ionised case are considered here
        ! (process becomes unimportant for higher ionisation cases)
        subroutine setDustGasCollHeatCool()
          implicit none         

          real :: Z , kT, eta, psi,xi,  S, IPerg
          real :: eightkT_pi ! 8*k * Te/Pi [erg]
          integer :: ss(30)
          real :: vmean, mcp ! colliding particle mean velocity and mean particle mass
          integer :: nelectrons, istage, ns, na

          ! outer shell array
          ss = (/1,1,2,2,3,3,3,3,3,3,4,4,5,5,5,5,5,5,&
               &6,6,6,6,6,6,6,6,6,6,7,7/)

          gasDustColl_g = 0.
          gasDustColl_d = 0.

          kT =6.336e-6*TeUsed ! ryd
          eightkT_pi = (1.1045e-15)*TeUsed/Pi          


          do elem = 1, nElements ! 0 for electrons
             if ( lgElementOn(elem)) then                                   
                   
                   do istage = 1, min(nstages,elem+1)
                      Z = real(istage-1)

                      mcp = aWeight(elem)*amu
                      vmean = sqrt(eightkT_pi/mcp)

                      ! number of e-'s of the ionisation stage above
                      nelectrons = elem-istage+1
                      
                      if (istage>1) then
                         if (elem == 1) then
                            IPerg = nuArray(HlevNuP(1))*ryd2erg
                         else if (elem == 2 .and. istage == 2) then
                            IPerg = nuArray(HeIlevNuP(1))*ryd2erg
                         else if (elem == 2 .and. istage == 3) then
                            IPerg = nuArray(HeIIlevNuP(1))*ryd2erg
                         else
                            IPerg = nuArray(elementP(elem,istage-1,nShells(elem,istage-1),1))*ryd2erg
                         end if
                      end if

                      do ns = 1, nSpeciesPart(nspU)
                         if (istage == 1) then
                            S = 2*mcp*MsurfAtom(ns)/(mcp+MsurfAtom(ns))**2.
                         else
                            S = 1.
                         end if
                         do na = 1, nSizes

                            psi = Z*grainPot(ns,na)/kT
                            if (psi <= 0. ) then
                               eta = 1.-psi
                               xi = 1. - psi/2.
                            else
                               eta = exp(-psi)
                               xi = (1.+psi/2.) * eta
                            end if
             
                            gasDustColl_d(ns,na) = gasDustColl_d(ns,na)+& 
                                 & ionDenUsed(elementXref(elem),istage)*& 
                                 & grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)*&
                                 & grid%Hden(cellP)*& 
                                 & Pi*grainRadius(na)*grainRadius(na)*1.e-8*&
                                 & S*vmean*(2.*kT*Ryd2erg*xi-eta*&
                                 & (Z*grainPot(ns,na)*Ryd2erg-IPerg+& 
                                 & 2.*kBoltzmann*grid%Tdust(ns,na,cellP)))     

                            gasDustColl_g = gasDustColl_g + & 
                                 & ionDenUsed(elementXref(elem),istage)*& 
                                 & grainWeight(na)*grainAbun(nspU,ns)*grid%Ndust(cellP)*&
                                 & grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)*&
                                 & Pi*grainRadius(na)*grainRadius(na)*1.e-8*&
                                 & S*vmean*(2.*kT*Ryd2erg*xi-eta*2.*kBoltzmann*& 
                                 & grid%Tdust(ns,na,cellP))
                         end do
                      end do
                   end do
                end if
             end do
                   
             ! add contribution of e- collisions
             vmean = sqrt(eightkT_pi/me)
             do ns = 1, nSpeciesPart(nspU)
                S = 1.
                Z=-1.
                do na = 1, nSizes


                   psi = Z*grainPot(ns,na)/kT
                   if (psi <= 0. ) then
                      eta = 1.-psi
                      xi = 1. - psi/2.
                   else
                      eta = exp(-psi)
                      xi = (1.+psi/2.) * eta
                   end if
                   
                   gasDustColl_d(ns,na) = gasDustColl_d(ns,na)+& 
                        & NeUsed* &
                        & Pi*grainRadius(na)*grainRadius(na)*1.e-8*&
                        & S*vmean*(2.*kT*Ryd2erg*xi-eta*&
                        & (Z*grainPot(ns,na)*ryd2erg))
                   
                   gasDustColl_g = gasDustColl_g + & 
                        & NeUsed* &
                        & grainWeight(na)*grainAbun(nspU, ns)*grid%Ndust(cellP)*&
                        & Pi*grainRadius(na)*grainRadius(na)*1.e-8*&
                        & S*vmean*(2.*kT*Ryd2erg*xi)/&
                        & grid%Hden(cellP)
                end do
             end do

             ! factor of fourpi to make up for the lack of fourpi 
             ! in the balance eqns for J
             gasDustColl_d(:,:)= gasDustColl_d(:,:)/fourpi
             gasDustColl_g= gasDustColl_g/fourpi
             
           end subroutine setDustGasCollHeatCool


        subroutine thermBalance(heatInt, coolInt)
            implicit none

            real, intent(out)      :: heatInt, &    ! total heating and
                 & coolInt       ! cooling integrals
            ! local variables
            
            integer                :: i,j,k         ! counters
            integer                :: elem, ion     ! counters
 
            real                   :: betaFF        ! energy loss coeff due to ff rad
            real                   :: betaRec       ! energy loss coeff due to recomb
            real                   :: ch12, ch13, & ! collision excitation of 
                 & ex12, ex13,&                     ! hydrogen data
                 & th12, th13                       ! (Mathis, Ly alpha, beta)
            real                   :: coolFF        ! cool due to FF radiation [erg/s/Hden]
            real                   :: coolColl      ! cool due to coll excit   [erg/s/Hden]
            real                   :: coolRec       ! cool due to recombination[erg/s/Hden]
            real                   :: fcool
            real                   :: heatIonSte    ! heat due to this ion stellar phot
            real                   :: heatIonDif    ! heat due to this ion diffuse phot
            real                   :: heatSte       ! tot heat gain due to stellar phot
            real                   :: heatDif       ! tot heat gain due to diffuse phot
            real                   :: hTerm1        ! heat calculation term
            real                   :: hTerm2        ! heat calculation term
            real                   :: log10Te       ! log10(Te)
            real                   :: log10TeScaled ! log10(Te/Z^2)
            real                   :: Np            ! proton density
            real                   :: Te4           ! Te/10000.


            ! find log10Te and Te4 at this cell
            log10Te = log10(TeUsed)
            Te4     = TeUsed / 1.e4

            ! open ff, rec beta files for H+ and He2+ (Hummer, MNRAS 268(1994) 109, Table 1.  

            ! find the N(H+)
            Np = grid%ionDen(cellP,elementXref(1),2)*grid%elemAbun(grid%abFileIndex(xP,yP,zP),1)

            ! cooling of gas due to FF radiation from H+ 
            ! fits to Hummer, MNRAS 268(1994) 109, Table 1. or  least square fitting to m=4

            betaFF = 1.0108464E-11 + 9.7930778E-13*log10Te - &
                 & 6.6433144E-13*log10Te*log10Te + 2.4793747E-13*log10Te*log10Te*log10Te -&
                 & 2.3938215E-14*log10Te*log10Te*log10Te*log10Te

            coolFF = Np*NeUsed*betaFF*kBoltzmann*TeUsed/sqrt(TeUsed)

            if (lgTraceHeating.and.taskid==0) then
               write(57,*) 'Cell: ', xp,yp,zp
               write(57,*) 'FF H+: ', coolFF
            end if

            ! cooling of gas due to recombination of H+
            ! fits to Hummer, MNRAS 268(1994) 109, Table 1.  
            ! least square fitting to m=4
            betaRec = 9.4255985E-11 -4.04794384E-12*log10Te &
                 & -1.0055237E-11*log10Te*log10Te +  1.99266862E-12*log10Te*log10Te*log10Te&
                 & -1.06681387E-13*log10Te*log10Te*log10Te*log10Te

            coolRec = Np*NeUsed*betaRec*kBoltzmann*TeUsed/sqrt(TeUsed)

            if (lgTraceHeating.and.taskid==0) then
               if (nIterateT > 1) then
                  do i = 1, 12
                     backspace 57
                  end do
               end if

               write(57,*) 'Rec H+: ', coolrec
            end if

            ! cooling of gas due to FF radiation from He++
            ! fits to Hummer, MNRAS 268(1994) 109, Table 1. least square fitting to m=4 and 
            ! scaled to Z=2

            ! find N(He++)
            Np =  grid%ionDen(cellP,elementXref(2),3)*grid%elemAbun(grid%abFileIndex(xP, yP,zP),2)

            log10TeScaled = log10(TeUsed/4.)

            betaFF = 2.*(1.0108464E-11 + 9.7930778E-13*log10TeScaled - &
                 & 6.6433144E-13*log10TeScaled*log10TeScaled + &
                 & 2.4793747E-13*log10TeScaled*log10TeScaled*log10TeScaled -&
                 & 2.3938215E-14*log10TeScaled*log10TeScaled*log10TeScaled*log10TeScaled)


            coolFF = coolFF + Np*NeUsed*betaFF*kBoltzmann*TeUsed/sqrt(TeUsed/4.)

            if (lgTraceHeating.and.taskid==0) then
               write(57,*) 'FF He++: ',Np*NeUsed*betaFF*kBoltzmann*TeUsed/sqrt(TeUsed/4.)
            end if


            ! cooling of gas due to recombination of He++
            ! fits to Hummer, MNRAS 268(1994) 109, Table 1.  least square fitting to m=4 
            ! and scaled to Z=2
            betaRec = 2.*(9.4255985E-11 -4.04794384E-12*log10TeScaled &
                 & -1.0055237E-11*log10TeScaled*log10TeScaled  &
                 & +1.99266862E-12*log10TeScaled*log10TeScaled*log10TeScaled&
                 & -1.06681387E-13*log10TeScaled*log10TeScaled*log10TeScaled*log10TeScaled)

            coolRec = coolRec + Np*NeUsed*betaRec*kBoltzmann*TeUsed/sqrt(TeUsed/4.)

            if (lgTraceHeating.and.taskid==0) then
               write(57,*) 'Rec He++: ',Np*NeUsed*betaRec*kBoltzmann*TeUsed/sqrt(TeUsed/4.)
            end if

            ! cooling of gas due to FF radiation from He+
            ! fits to Hummer and Storey, MNRAS 297(1998) 1073, Table 6. least square fitting to m=4 

            ! find N(He+)
            Np =  grid%ionDen(cellP,elementXref(2),2)*grid%elemAbun(grid%abFileIndex(xp, yP, zP),2)

            betaFF =  1.070073e-11    -2.5730207e-13*log10Te + &
                 & 2.109134e-13*log10Te*log10Te 

            coolFF = coolFF + Np*NeUsed*betaFF*kBoltzmann*TeUsed/sqrt(TeUsed)

            if (lgTraceHeating.and.taskid==0) then
               write(57,*) 'FF He+: ',Np*NeUsed*betaFF*kBoltzmann*TeUsed/sqrt(TeUsed)
            end if


            ! cooling of gas due to recombination of He+
            ! fits to Hummer and Storey, MNRAS 297(1998) 1073, Table 6. least square fitting to m=4
            betaRec =    9.4255985E-11 -4.04794384E-12*log10Te &
                 & -1.0055237E-11*log10Te*log10Te  &
                 & +1.99266862E-12*log10Te*log10Te*log10Te &
                 & -1.06681387E-13*log10Te*log10Te*log10Te*log10Te 

            coolRec = coolRec + Np*NeUsed*betaRec*kBoltzmann*TeUsed/sqrt(TeUsed)

            if (lgTraceHeating.and.taskid==0) then
               write(57,*) 'Rec He+: ',Np*NeUsed*betaRec*kBoltzmann*TeUsed/sqrt(TeUsed)
            end if

            ! collisional excitation of Hydrogen
            ! Mathis, Ly alpha, beta
            ch12 = 2.47e-8
            ch13 = 1.32e-8
            ex12 = -0.228
            ex13 = -0.460
            th12 = 118338.
            th13 = 140252.

            if (TeUsed > 5000.) then 
             
                coolColl = (ch12*exp(-th12/TeUsed)*Te4**ex12 + &
                     & ch13*exp(-th13/TeUsed)*Te4**ex13) * &
                     & hcRyd*grid%ionDen(cellP,elementXref(1),1)*NeUsed

             else

                coolColl = 0.
     
             end if

            if (lgTraceHeating.and.taskid==0) then
               write(57,*) 'Coll exc H: ',coolColl
               fcool = 0.
            end if


             ! collisional excitation of Heavies

             ! get the emissivities of the forb lines
             call forLines()

             ! sum all contributions from the heavies to coolColl 
             do elem = 3, nElements
                do ion = 1, min(elem+1, nstages)
                   do j = 1, 10
                      do k = 1, 10
                         coolColl = coolColl + forbiddenLines(elem,ion,j,k)                         
                         if (lgTraceHeating.and.taskid==0) then
                            fcool = fcool + forbiddenLines(elem,ion,j,k)                         
                         end if
                      end do
                   end do
                end do
             end do

            if (lgTraceHeating.and.taskid==0) then
               write(57,*) 'CELs cool: ',fcool
               fcool = 0.
            end if

             ! heating due to photoionization

             ! re-initialize heatSte and heatDif
             heatSte = 0.
             if (lgDebug) heatDif = 0.

             do elem = 1, nElements  ! begin element loop
                 do ion = 1, min(elem, nStages-1) ! begin ion loop
                     if(.not.lgElementOn(elem)) exit
                         
                     if (elem > 2) then

                         ! find the number of electrons in this ion
                         nElec = elem - ion + 1

                         ! find the outer shell number and statistical weights 
                         call getOuterShell(elem, nELec, outShell, g0, g1)
                         
                         ! get pointer to this ion's IP in NuArray
                         IPnuP = elementP(elem, ion, outShell, 1)

                         ! get pointer to this cell's highest energy
                         highNuP = elementP(elem, ion, outShell, 2) 

                     else if (elem == 1) then ! HI
            
                         IPNuP   = HlevNuP(1)
                         highNuP = nbins

                     else if ( (elem == 2) .and. (ion == 1) ) then ! HeI
 
                         IPNuP   = HeIlevNuP(1)
                         highNuP = nbins 

                     else if ( (elem == 2) .and. (ion == 2) ) then ! HeII

                         IPNuP   = HeIIlevNuP(1)
                         highNuP = nbins

                     end if

                     ! re-initialize heatIonSte and heatIonDif
                     heatIonSte = 0.
                     if (lgDebug) heatIonDif = 0.
                     do j = IPnuP, highNuP ! begin frequency loop
                             if ( elem > 2 ) then

                                 ! get pointer to phot xSec of this ion in xSecArray
                                 xSecP = elementP(elem, ion, outShell, 3)

                                 ! get phot xSec of ion at this energy
                                 phXSec = xSecArray(j+xSecP-IPnuP+1-1)

                             else if ( elem == 1 ) then ! HI

                                 phXSec  = xSecArray(j-IPnuP+1+HlevXSecP(1)-1)

                             else if ( (elem == 2) .and. (ion == 1) ) then ! HeI


                                 phXSec  = xSecArray(j-IPnuP+1+HeISingXSecP(1)-1)
      

                             else if ( (elem == 2) .and. (ion == 2) ) then ! HeII

                                 phXSec  = xSecArray(j-IPnuP+1+HeIIXSecP(1)-1)
 
                             end if

                             if ((phXSec < 1.e-35) ) then
                                 exit
                             end if

                             if ((nuArray(j) < 1.e-35) ) then
                                 print*, "! thermBalance: insane nu value (j, nuArray(j)", &
                                      & j, nuArray(j)
                                 stop
                             end if    

                             ! stellar photoionizations 
                             if ( grid%Jste(cellP,j) > 0. )then

                                heatIonSte = heatIonSte + phXSec*grid%Jste(cellP,j)*&
                                     & (nuArray(j)-nuArray(IPNuP)) / (nuArray(j))
 
                             end if 

                             ! diffuse photoionizations

                             if (lgDebug) then
                                if ( grid%Jdif(cellP,j) > 0. ) then

                                   heatIonDif = heatIonDif + phXSec*grid%JDif(cellP,j)* &
                                        & (nuArray(j)-nuArray(IPNuP)) / (nuArray(j))

                                end if
                             end if
                                 
                    end do ! end frequency loop

                    ! get total heat gain from stellar and diffuse photoionizations
 
                    heatIonSte = heatIonSte*ionDenUsed(elementXref(elem),ion)*& 
                         & grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)
                    if (lgDebug) &
                         & heatIonDif = heatIonDif*ionDenUsed(elementXref(elem),ion)*& 
                         &grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)

                    heatSte    = heatSte + heatIonSte
                    if (lgDebug) &
                         & heatDIf    = heatDif + heatIonDif

                end do ! end ion loop   

            end do ! end element loop

            ! calculate the total heating and cooling integrals
 
            if (lgDebug) then
               heatInt = heatSte + heatDIf
            else
               heatInt = heatSte
            end if
           
            if (lgTraceHeating.and.taskid==0) then
               write(57,*) 'Dust gas coll cool: ',gasDustColl_g
               write(57,*) 'Heat photionization: ', heatInt
               write(57,*) 'Heat photoelectric: ', photoelHeat_g
            end if

            coolInt = coolFF + coolRec + coolColl

            if (lgDust .and. lgPhotoelectric) then
               coolInt = coolInt+gasDustColl_g
               heatInt = heatInt+photoelHeat_g
            end if
           
            
        end subroutine thermBalance


        ! this subroutine is the driver for the calculation of the emissivity
        ! from the heavy elements forbidden lines. 
        subroutine forLines()
          implicit none

          integer        :: elem, ion ! counters

          ! re-initialize forbiddenLines
          forbiddenLines = 0.

          do elem = 3, nElements
             do ion = 1, min(elem+1, nstages)
                if (.not.lgElementOn(elem)) exit

                if (lgDataAvailable(elem, ion)) then
                  
                   if (elem == 26 .and. ion == 2) then
                      if (nstages > 2) then
                         call equilibrium(file_name=dataFile(elem, ion), &
                              &ionDenUp=ionDenUsed(elementXref(elem),ion+1)/&
                              &ionDenUsed(elementXref(elem),ion), Te=TeUsed,&
                              &Ne=NeUsed, FlineEm=forbiddenLinesLarge(&
                              &1:nForLevelsLarge, 1:nForLevelsLarge))
                      else
                         call equilibrium(file_name=dataFile(elem, ion), &
                              &ionDenUp=0., Te=TeUsed,&
                              &Ne=NeUsed, FlineEm=forbiddenLinesLarge(&
                              &1:nForLevelsLarge, 1:nForLevelsLarge))
                      end if
                      forbiddenLinesLarge(:, :) = forbiddenLinesLarge(:, :)*& 
                           & grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)*&
                           & ionDenUsed(elementXref(elem), ion)                      
                   else
                      if (ion<nstages) then
                         call equilibrium(dataFile(elem, ion), &
                              &ionDenUsed(elementXref(elem), ion+1)/&
                              &ionDenUsed(elementXref(elem),ion), &
                              & TeUsed, NeUsed, forbiddenLines(elem, ion,:,:))
                      else
                         call equilibrium(dataFile(elem, ion), 0., &
                            & TeUsed, NeUsed, forbiddenLines(elem, ion,:,:))
                      end if
                      forbiddenLines(elem, ion, :, :) = forbiddenLines(elem, ion, :, :)*& 
                           & grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)*&
                           & ionDenUsed(elementXref(elem), ion)                      
                   end if
                end if
             end do
          end do

          ! scale the forbidden lines emissivity to give units of [erg/s/Ngas]
          ! comment: the forbidden line emissivity is so far in units of cm^-1/s/Ngas
          !          the energy [erg] of unit wave number [cm^-1] is 1.9865e-16, hence 
          !          the right units are obtained by multiplying by 1.9865e-16 
          forbiddenLines = forbiddenLines*1.9865e-16   
          if (lgElementOn(26) .and. ion>2) forbiddenLinesLarge = forbiddenLinesLarge*1.9865e-16

        end subroutine forLines

        recursive subroutine ionBalance(lgConv)
            implicit none

            logical, intent(out)   :: lgConv        ! did ion bal converge? 

            ! local variables
            real                   :: collIon       ! collisional ionization of H 
            real                   :: correction    ! used in lgNeInput = .t.
            real                   :: deltaHI       ! delta(X(H0))
            real                   :: deltaHeI      ! delta(X(He0))
            real                   :: deltaHeII     ! delta(X(HeII))
            real                   :: expFact       ! exponential factor
            real                   :: fac0,fac1     ! calculation factors
            real                   :: t4            ! TeUsed/10000.
            real, save             :: HIOld         ! X(H0) from last iteration
            real, save             :: HeIOld        ! X(He0) from last iteration
            real, save             :: HeIIOld       ! X(He+) from last iteration
            real, parameter        :: limit = 0.01  ! convergence limit

            integer                :: elem, ion     ! element and ion counters
            integer                :: g0,g1         ! stat weights
            integer                :: i             ! counter
            integer                :: nElec         ! of  e's in the ion
            integer                :: outShell      ! byproduct of proc to get stat weights


            real                   :: revRate       ! reverse charge exchange rate

            double precision, dimension(nELements) :: &
                 & denominator   ! denominator of final ion abundance


            real, dimension(nElements, nstages) :: &
                 & deltaE_k      ! deltaE/k [K]

            double precision, dimension(nElements, nstages) :: &
                 & ionRatio, &   ! X(i+1)/X(+i)
                 & ionProd       ! ionRatio products
            real, dimension(nElements, nstages,4) :: &
                 & chex         ! ch exchange coeff in cm^3/s


            ! initialize variables
            correction  = 0.
            ionRatio    = 0.
            ionProd     = 1.
            denominator = 1.
            lgConv      = .true.

            ! take into account collisional ionization of H
            ! Drake & Ulrich, ApJS42(1980)351
            expFact = 157893.94/TeUsed
 
            if (expFact > 75.) then
                ! prevents underflow of exponential factor in collIon
                expFact = 75.
            end if
         
            collIon = 2.75E-16*TeUsed*sqrt(TeUsed)*&
                 & (157893.94/TeUsed+2.)*exp(-expFact)                           


            ! get HIOld, HeIOld, HeIIOld from last iteration
            HIOld   = ionDenUsed(elementXref(1),1)
            HeIOld  = ionDenUSed(elementXref(2),1)
            HeIIOld = ionDenUsed(elementXref(2),2)

            ! the set of charge exchange coeffs is not complete; the following might need
            ! to be changed when a more complete set is available



            chex          = 0.
            deltaE_k      = 0.

            chex(2,1,:)  = (/7.47e-6, 2.06, 9.93,-3.89/)! He0
            chex(2,2,:)  = (/1.e-5  , 0.  , 0.  , 0./)  ! He+
            chex(3,2,:)  = (/1.26   , 0.96,3.02 ,-0.65/)! Li+
            chex(3,3,:)  = (/1.e-5  , 0.  , 0.  , 0. /) ! Li+2
            chex(4,2,:)  = (/1.e-5  , 0.  , 0.  , 0. /) ! Be+
            chex(4,3,:)  = (/1.e-5  , 0.	 , 0.  , 0. /) ! Be+2
            chex(4,4,:)  = (/5.17   , 0.82, -.69, -1.12 /)! Be+3
            chex(5,2,:)  = (/2.e-2  , 0.  , 0.  , 0. /) ! B+
            chex(5,3,:)  = (/1.e-5  , 0.  , 0.  , 0. /) ! B+2
            chex(5,4,:)  = (/5.27e-1, 0.76,-0.63,-1.17/)! B+3
            chex(6,1,:)  = (/1.76e-9, 8.33, 4278.78, -6.41/)! C0
            chex(6,2,:)  = (/1.67e-4, 2.79, 304.72, -4.07/)! C+
            chex(6,3,:)  = (/3.25   , 0.21, 0.19, -3.29/)! C+2
            chex(6,4,:)  = (/332.46 ,-0.11,-0.995,-1.58e-3/)! C+3
            chex(7,1,:)  = (/1.01e-3,-0.29,-0.92, -8.38/)! N0
            chex(7,2,:)  = (/3.05e-1, 0.60, 2.65, -0.93/)! N+
            chex(7,3,:)  = (/4.54   , 0.57,-0.65, -0.89/)! N2+
            chex(7,4,:)  = (/3.28   , 0.52,-0.52, -0.19/)! N3+
            chex(8,1,:)  = (/1.04   , 3.15e-2, -0.61, -9.73/)! O0
            chex(8,2,:)  = (/1.04   , 0.27, 2.02, -5.92/)! O+
            chex(8,3,:)  = (/3.98   , 0.26, 0.56, -2.62/)! O2+
            chex(8,4,:)  = (/2.52e-1, 0.63, 2.08, -4.16/)! O3+
            chex(9,2,:)  = (/1.e-5  , 0.	 , 0.  , 0./) ! F+
            chex(9,3,:)  = (/9.86   , 0.29,-0.21,-1.15/) ! F+2
            chex(9,4,:)  = (/7.15e-1, 1.21,-0.70,-0.85/) ! F3+
            chex(10,2,:) = (/1.e-5  , 0.  , 0.  , 0.  /) ! Ne+
            chex(10,3,:) = (/14.73  , 4.52e-2, -0.84, -0.31 /) ! Ne+2
            chex(10,4,:) = (/6.47   , 0.54 , 3.59 , -5.22 /) ! Ne+3
            chex(11,2,:) = (/1.e-5  , 0.   , 0.   , 0. /) ! Na+
            chex(11,3,:) = (/1.33   , 1.15 , 1.20 , -0.32 /)! Na+2
            chex(11,4,:) = (/1.01e-1, 1.34 , 10.05, -6.41 /)! Na+3
            chex(12,2,:) = (/8.58e-5, 2.49e-3, 2.93e-2, -4.33 /)! Mg+
            chex(12,3,:) = (/6.49   , 0.53 , 2.82, -7.63 /) ! Mg+2
            chex(12,4,:) = (/6.36   , 0.55 , 3.86, -5.19 /) ! Mg+3
            chex(13,2,:) = (/1.e-5  , 0.   , 0.  , 0./) ! Al+
            chex(13,3,:) = (/7.11e-5, 4.12 , 1.72e4, -22.24/)! Al+2
            chex(13,4,:) = (/7.52e-1, 0.77 , 6.24, -5.67/) ! Al+3
            chex(14,2,:) = (/1.23   , 0.24 , 3.17, 4.18e-3/) ! Si+
            chex(14,3,:) = (/4.900e-1, -8.74e-2, -0.36, -0.79/)! Si+2
            chex(14,4,:) = (/7.58   , 0.37 , 1.06, -4.09/)! Si+3
            chex(16,1,:) = (/3.82e-7, 11.10, 2.57e4, -8.22/)! S0
            chex(16,2,:) = (/1.e-5  , 0.   , 0.   ,0. /)! S+
            chex(16,3,:) = (/2.29   , 4.02e-2, 1.59, -6.06/)! S+2
            chex(16,4,:) = (/6.44   , 0.13 , 2.69 , -5.69/)! S+3
            chex(18,2,:) = (/1.e-5  , 0.   , 0.    , 0.	/) ! Ar+
            chex(18,3,:) = (/4.57   , 0.27 , -0.18 , -1.57/)! Ar+2
            chex(18,4,:) = (/6.37   , 2.12 , 10.21 , -6.22/)! Ar+3
            chex(18,3,:) = (/3.17e-2, 2.12 , 12.06 , -0.40/)! Ca+2
            chex(18,4,:) = (/2.68   , 0.69 , -0.68 , -4.47/)! Ca+3
            chex(26,2,:) = (/1.26   , 7.72e-2, -0.41, -7.31/)! Fe+
            chex(26,3,:) = (/3.42   , 0.51 , -2.06 , -8.99/)! Fe+2.	


            deltaE_k(7,1) = 10863.                    
            deltaE_k(8,1) = 2205.

            chex(:,:,1) = chex(:,:,1)*1.e-9


            t4 = TeUsed/10000.
            
            
            ! calculate the X(i+1)/X(i) ratio            
            do elem = 1, nElements
                do ion = 1, min(elem, nstages-1)
                    if (.not.lgElementOn(elem)) exit

                    chex(elem,ion,1) = chex(elem,ion,1)*(t4**chex(elem,ion,2))*& 
                         & (1.+chex(elem,ion,3)*exp(chex(elem,ion,4)*t4))

                    if (chex(elem,ion,1) < 0. ) chex(elem,ion,1) = 0.
                    if (TeUsed < 6000. .or. TeUsed>5.e4) chex(elem,ion,1) = 0.

                    ! find the number of electron in this ion
                    nElec = elem - ion +1

                    ! find the stat weights
                    call getOuterShell(elem, nElec, outShell, g0, g1)                        

                    ! calculate the reverse charge exchange rate (only if deltaE_k > 1.)
                    if ( deltaE_k(elem,ion) > 1.) then

                       revRate = chex(elem,ion,1) * &
                            & (1.-grid%ionDen(cellP, &
                            & elementXref(1),1))*grid%Hden(cellP)*&
                            & 2. * exp(-deltaE_k(elem,ion)/TeUsed)/(real(g0)/real(g1))
                    else
                       revRate = 0.
                    end if

                    if ( (elem>1) .or. (ion>1) ) collIon = 0.

                    ! calculate the X(i+1)/X(i) ratio  
                    if (lgDebug) then
                       ionRatio(elem,ion) = (nPhotoSte(elem,ion)+nPhotoDif(elem,ion)+&
                            & collIon+revRate)/&
                            & (NeUsed*alphaTot(elem,ion)& 
                            & +chex(elem,ion,1)*grid%Hden(grid%active(xP,yp,zP)) * &
                            & grid%ionDen(cellP,elementXref(1),1))
                    else
                       ionRatio(elem,ion) = (nPhotoSte(elem,ion)+&
                            & collIon+revRate)/&
                            & (NeUsed*alphaTot(elem,ion)& 
                            & +chex(elem,ion,1)*grid%Hden(grid%active(xP,yp,zP)) * &
                            & grid%ionDen(cellP,elementXref(1),1))

                    end if
                 end do
              end do
           
            ! calculate the products of ionRatio
            do elem = 1, nElements
                do ion = 1, min(elem, nstages-1)
                   if (.not.lgElementOn(elem)) exit
 
                   ! generate the product
                   do i = 1, ion
                       ionProd(elem, ion) = ionProd(elem,ion)*ionRatio(elem, i)

                   end do

                end do
            end do
 
            ! calculate denominators for final ion abundances
            do elem = 1, nElements
                do ion = 1, min(elem, nstages-1)
                    if (.not.lgElementOn(elem)) exit

                    denominator(elem) = denominator(elem) + ionProd(elem, ion)

               end do
            end do


            ! calculate the new abundances for all ions
            do elem = 1, nElements
                do ion = 1, min(elem+1, nstages)

                    if (.not.lgElementOn(elem)) exit


                    if (ion == 1) then 
                        grid%ionDen(cellP,elementXref(elem),ion) = 1./denominator(elem)
                    else

                        grid%ionDen(cellP,elementXref(elem),ion) = &
                             & ionProd(elem, ion-1)/denominator(elem)
 
                   end if               

                   ! take back within the limit
                   if (grid%ionDen(cellP,elementXref(elem),ion) > xMax) &
                        & grid%ionDen(cellP,elementXref(elem),ion) = xMax
                   if (grid%ionDen(cellP,elementXref(elem),ion) < 1.e-20) &
                        & grid%ionDen(cellP,elementXref(elem),ion) = 1.e-20


                   ionDenUsed(elementXref(elem), ion) = &
                        & grid%ionDen(cellP,elementXref(elem),ion)
                   ! this was added to help MPI communication
                   ionDenTemp(cellP,elementXref(elem),ion) = &
                        & grid%ionDen(cellP,elementXref(elem),ion)


                end do
            end do            

            ! calculate new Ne 
          NeUsed = 0.
            do elem = 1, nElements
                do ion = 2, min(elem+1, nstages)
                    if (lgElementOn(elem)) then
                      if( ionDenUsed(elementXref(elem),ion) >= 1.e-10) &
                           & NeUsed = NeUsed + (ion-1)*& 
                           &grid%elemAbun(grid%abFileIndex(xP,yP,zP),elem)*&
                           &ionDenUsed(elementXref(elem), ion)
                    end if

                end do
            end do

            NeUsed = NeUsed * grid%Hden(cellP)

            if (NeUsed==0.) print*, '! ionBalance [warning]: cell ', xP,yP,zP, &
                 &'; NeUsed = ',  NeUsed

            if (NeUsed == 0.) then
               NeUsed = 1.
            end if

            if (LgNeInput) then
               correction = NeUsed/grid%NeInput(cellP)
               grid%Hden(cellP) = grid%Hden(cellP)/correction
            end if

            ! this was added to help MPI implementation
            NeTemp(cellP) = NeUsed

            grid%Ne(cellP) = NeUsed


            ! calculate the residuals
            deltaHI   = (ionDenUsed(elementXref(1),1) - HIOld)   / HIOld
            deltaHeI  = (ionDenUsed(elementXref(2),1) - HeIOld)  / HeIOld 
            deltaHeII = (ionDenUsed(elementXref(2),2) - HeIIOld) / HeIIOld

            ! check for convergence
            if ( ( (abs(deltaHI)>limit) .or. (abs(deltaHeI)>limit) .or. &
                 &(abs(deltaHeII)>limit) ) .and. (nIterateX<maxIterateX) ) then
                
                ! stepa up X iteration
                nIterateX = nIterateX + 1
                call ionBalance(lgConv)
                return

            else if (( (abs(deltaHI)>limit) .or. (abs(deltaHeI)>limit) .or. &
                 &(abs(deltaHeII)>limit) ) .and. (nIterateX == maxIterateX) ) then
            
               if (lgTalk)  print*, "! ionBalance: [warning] convergence not reached after ", &
                    &maxIterateX, " steps. Finishing up..."

                lgConv = .false.
            end if
         
            ! this was added to help MPI implementation
            NeTemp(cellP) = NeUsed

        end subroutine ionBalance

        subroutine calcAlpha()
            implicit none

            real, dimension(nElements, nstages) &
                 &:: diRec  ! total recombination coeffs
           

            ! zero out alphaTot and contributors
            alphaTot = 0.
            diRec    = 0.

            ! calculate radiative recombination part first
            do elem = 1, nElements
                do ion = 1, min(nstages-1, elem)
                    alphaTot(elem, ion) = radRecFit(elem, elem-ion+1)    
                end do
            end do

            call dielectronic(diRec)
            
            ! calculate dielectronic recombination part
            do elem = 3, nElements
                do ion = 1, min(nstages-1, elem)
                   if (diRec(elem,ion) == 0.) diRec(elem,ion) = diRec(8,ion)
                    alphaTot(elem, ion) = alphaTot(elem, ion) + diRec(elem, ion)
                end do
            end do

            

       end subroutine calcAlpha
            

       ! This subroutine calculates rates of radiative recombination for all ions
       ! of all elements from H through Zn by use of the following fits:
       ! H-like, He-like, Li-like, Na-like - Verner & Ferland, 1996, ApJS, 103, 467
       ! Other ions of C, N, O, Ne - Pequignot et al. 1991, A&A, 251, 680,
       ! refitted by Verner & Ferland formula to ensure correct asymptotes
       ! Fe XVII-XXIII - Arnaud & Raymond, 1992, ApJ, 398, 394
       ! Fe I-XV - refitted by Verner & Ferland formula to ensure correct asymptotes
       ! Other ions of Mg, Si, S, Ar, Ca, Fe, Ni -
       !                     - Shull & Van Steenberg, 1982, ApJS, 48, 95
       ! Other ions of Na, Al - Landini & Monsignori Fossi, 1990, A&AS, 82, 229
       ! Other ions of F, P, Cl, K, Ti, Cr, Mn, Co (excluding Ti I-II, Cr I-IV,
       ! Mn I-V, Co I)        - Landini & Monsignori Fossi, 1991, A&AS, 91, 183
       ! All other species    - interpolations of the power-law fits
       !  Input parameters: z - atomic number
       !                    n - number of electrons from 1 to z
       ! Output parameter:  radRecFit  - rate coefficient, cm^3 s^(-1)

       function radRecFit(z, n)
            implicit none

            integer, intent(in)             :: z     ! atomic weight of the element
            integer, intent(in)             :: n     ! number of electrons (from 1 to z)
                

            real                            :: radRecFit    ! rad rate coeff [cm^3/s]

            ! local variables
            integer                              :: elem    ! element counter
            integer                              :: i       ! counter
            integer                              :: ion     ! ion stage counter
            integer                              :: ios     ! I/O error status

 
            logical, save                        :: lgFirst = .true.
                                                            ! first time this is evaluated?

            real                                 :: tt      ! temp dep fact in interpolation

            real, dimension(2, nElements, nElements),&      ! coefficients for the  
                 & save  :: rrec    ! calculation of the 
            real, dimension(4, nElements, nElements),&      ! radiative rates
                 & save  :: rnew    ! 
            real, dimension(3, 4:13), save         :: fe    ! 

            ! if this is the first time this procedure is called 
            ! read in radiative recombination coefficient file 
            ! and the dielectronic recombination data
            if (lgFirst) then
                close(17)
                open (unit=17, file='data/radrec.dat', status='old',position='rewind', &
                     & iostat = ios, action="read")
   
                do ion = 4, 30
                    if (ion /= 11) then
                        do elem = ion, 30
                            if ( (elem /= 26) .or. (ion >= 14) ) then
                                read(unit=17, fmt=*, iostat=ios) (rrec(i,elem,ion), i=1,2)
                            end if
                        end do
                    end if
                end do

                do ion = 1, 3
                    do elem = ion, 30
                        read(unit=17, fmt=*, iostat=ios) (rnew(i,elem,ion), i=1,4)
                    end do
                end do

                do elem = 11, 30
                   read(unit=17, fmt=*, iostat=ios) (rnew(i,elem,11), i=1,4) 
                end do

                do ion = 4, 10
                    do elem = max(6,ion), 8
                        read(unit=17, fmt=*, iostat=ios) (rnew(i,elem,ion), i=1,4)
                    end do
                    read(unit=17, fmt=*, iostat=ios) (rnew(i,10,ion), i=1,4)
                end do

                do ion = 12, 26
                    read(unit=17, fmt=*, iostat=ios) (rnew(i,26,ion), i=1,4)
                end do

                ! close the files
                close(17)

                ! set data for Fe
                fe(1,:) = (/4.33e-10,3.91e-10,3.49e-10,3.16e-10,&
                     &2.96e-10,2.59e-10,2.24e-10,1.91e-10,1.68e-10,1.46e-10/)
                fe(2,:) = (/0.531,0.523,0.521,0.534, &
                     &0.557,0.567,0.579,0.601,0.602,0.597/)
                fe(3,:) = (/5.77e-02,6.15e-02,6.22e-02,6.02e-02,&
                     &5.79e-02,5.65e-02,5.49e-02,5.10e-02,5.07e-02,5.22e-02/)
                 
                ! set lgFirst to .false.
                lgFirst = .false.
            end if

            ! check right element and number of electron reference
            if ( (z<1) .or. (z>30) ) then
                print*, "! radRecFit: insane atomic number", z
                stop
            end if
            if ( (n<1) .or. (n>z) ) then
                print*, "! radRecFit: insane number of electrons", n
                stop
            end if

            ! calculate the rates
            if ( (n<=3) .or. (n==11) .or. ((z>5) .and. (z<9)) .or. (z==10) .or.&
                 &((z==26) .and. (n>11)) ) then

                tt = sqrt(TeUsed/rnew(3,z,n))
                radRecFit = rnew(1,z,n) / (tt*(tt+1.)**(1.-rnew(2,z,n))*&
                     &(1.+sqrt(TeUsed/rnew(4,z,n)))**(1.+rnew(2,z,n)))
            else    
               tt = TeUsed*1.e-4
               if ( (z==26) .and. (n<=13) ) then
                   radRecFit = fe(1,n)/tt**(fe(2,n)+fe(3,n)*log10(tt))
               else 
                   radRecFit = rrec(1,z,n)/tt**rrec(2,z,n)
               end if
            end if
                           
        end function radRecFit
                 
        subroutine dielectronic(diRec)
            implicit none

            real, intent(out), dimension(nElements, nstages) &
                 &:: diRec   ! total recombination coeffs

            ! local variables

            integer              :: elem      ! atomic number
            integer              :: ion       ! ionization stage       
            integer              :: n         ! number of electrons
            integer              :: ios       ! I/O error status
            integer              :: g         ! temperature flag
            
            real, dimension(nElements, nstages) :: &
                                    aldroPequi! high T dielec rec coeff by A&P73

            real                 :: a,b,c,d,f ! fitting coefficients
            real                 :: t,t0,t1   ! t = TeUsed/10000., t0,t1 are fitting par
            real                 :: alpha

            t = TeUsed/10000.

            alpha = 0.
            aldroPequi=0.

            close(18)
            open (unit=18, file='data/dielectronic.dat', status='old',position='rewind', &
                 &iostat = ios, action="read")
            do i = 1, 10000
               read(unit=18, fmt=*, iostat=ios) elem, n, a, b, c, d, f, g
               if (ios < 0) exit ! end of file reached
               ion = elem + 1 - n

               if (ion <= nstages) then
                  
                  if (g == 0) then
                     diRec(elem, ion) = (10.**(-12))*(a/t+b+c*t+d*t**2.)*t**(-3./2.)*exp(-f/t)
                  else if (g == 1) then
                     if (TeUsed < 20000.) then
                        diRec(elem, ion) = (10.**(-12))*(a/t+b+c*t+d*t**2.)*t**(-3./2.)*exp(-f/t)
                     end if
                  else if (g == 2) then
                     if (TeUsed >= 20000.) then
                        diRec(elem, ion) = (10.**(-12))*(a/t+b+c*t+d*t**2.)*t**(-3./2.)*exp(-f/t)
                     end if
                  end if
               end if
            end do
            close(18)
               
            ! calculate the high temperatures dielectronic recombination coeficients of 
            ! Aldrovandi and Pequignot 1973
            t = TeUsed

            close(17)
            open (unit=17, file='data/aldrovandi.dat', status='old',position='rewind', iostat = ios, action="read")
            if (ios /= 0) then
               print*, "! dielectronic: can't open file data/alrovandi.dat"
               stop
            end if

            do i = 1, 100000
               read(unit=17, fmt=*, iostat=ios) elem, n, a, b, t0, t1
               if (ios < 0) exit ! end of file reached
               ion = elem + 1 - n
     
               alpha = a*t**(-3./2.)*exp(-t0/t)*(1.+b*exp(-t1/t)) 
               
               ion = elem-n+1
               
               if (ion <= nstages) aldroPequi(elem, ion) = alpha

            end do
            
            close(17)

            
            if (TeUsed>60000.) then               
               
               diRec = aldroPequi
               
            else
               do elem = 1, nElements
                  do ion = 1, nstages
                     if (diRec(elem,ion) == 0.) diRec(elem,ion) = aldroPequi(elem,ion) 
                  end do
               end do
            end if
            

          end subroutine dielectronic
                  
         subroutine getDustT()
            implicit none

            real                   :: Tspike(nTbins),Pspike(nTbins)
            real                   :: dustAbsIntegral   ! dust absorption integral
            real                   :: dabs
            real                   :: resLineHeat       ! resonance line heating
            real, dimension(nbins) :: radField,yint     ! radiation field

            integer :: nS, i, ai ! counters
            integer :: iT        ! pointer to dust temp in dust temp array

            ! radiation field at this location
            if (lgDebug) then
               radField = ((grid%Jste(cellP,:) + grid%Jdif(cellP,:)))/Pi
            else
               radField = grid%Jste(cellP,:)/Pi
            end if

            ! zero out dust temperature arrays
            grid%Tdust(:,:,cellP) = 0.

            
            ! calculate absorption integrals for each species
            do nS = 1, nSpeciesPart(nspU)

               do ai = 1, nSizes

                  dustAbsIntegral=0.
                  if (lgTraceHeating.and.taskid==0) dabs=0.

                  do i = 1, nbins
                     dustAbsIntegral = dustAbsIntegral+xSecArray(dustAbsXsecP(nS,ai)+i-1)*radField(i)
                     if (lgTraceHeating.and.taskid==0) then
                        dabs = dabs+xSecArray(dustAbsXsecP(nS,ai)+i-1)*radField(i)                        
                     end if
                  end do
                 
                  if (lgGas .and. convPercent>=resLinesTransfer .and. (.not.lgResLinesFirst) .and. &
                       & (.not.nIterateMC==1) ) then
                     dustHeatingBudget(grid%abFileIndex(xp,yp,zp),0) = dustHeatingBudget(grid%abFileIndex(xp,yp,zp),0)+&  
                          & dustAbsIntegral*grainWeight(ai)*grainAbun(nspU,nS)*grid%Ndust(cellP)
                     dustHeatingBudget(0,0) = dustHeatingBudget(0,0)+&                                                    
                          & dustAbsIntegral*grainWeight(ai)*grainAbun(nspU,nS)*grid%Ndust(cellP)
                     resLineHeat = resLineHeating(ai,ns,nspU)
                     dustAbsIntegral = dustAbsIntegral+resLineHeat
                  end if

                  if (lgGas .and. lgPhotoelectric) then
                     dustAbsIntegral = dustAbsIntegral+gasDustColl_d(nS,ai)-& 
                          & photoelHeat_d(nS,ai)
                  end if

                  if (lgTraceHeating.and.taskid==0) then
                     write(57,*) 'Dust Species: ', ns, ai
                     write(57,*) 'Abs of cont rad field  ', dabs
                     write(57,*) 'Res Lines Heating: ', reslineheat
                     if (lgGas .and. lgPhotoelectric) then

                        write(57,*) 'Grain potential ', grainPot(ns,ai)
                        write(57,*) 'Gas-dust collision heat', gasDustColl_d(nS,ai)
                        write(57,*) 'Photoelctric cooling: ', photoelHeat_d(nS,ai)
                     end if
                  end if

                  call locate(dustEmIntegral(nS,ai,:), dustAbsIntegral, iT)
                  
                  if (iT<=0) then
                     print*, "getDustT: [warning] temperature of grain = 0. K!!!!"
                     print*, cellP
                     print*, nS, dustAbsIntegral
                     print*, dustEmIntegral(nS,ai,1)
                     !stop
                     iT=1
                     grid%Tdust(nS,ai,cellP) = 1.
                  else if (iT>=nTemps) then
                     grid%Tdust(nS,ai,cellP) = real(nTemps)
                  else
                     grid%Tdust(nS,ai,cellP) = real(iT) + &
                          & (dustAbsIntegral-dustEmIntegral(nS,ai,iT))*&
                          & (real(iT+1)-real(iT))/(dustEmIntegral(nS,ai,iT+1)-&
                          & dustEmIntegral(nS,ai,iT))
                  end if

                  if (lgTraceHeating.and.taskid==0) then
                     write(57,*) 'Radiative cooling: ', dustEmIntegral(ns,ai,iT)
                     write(57,*) 'Grain temperature: ', grid%Tdust(nS,ai,cellP), &
                          & " species ", grainLabel(nS), " size:", grainRadius(ai)

                  end if

                  if (lgTalk) &
                       & print*, "! getDustT: [talk] cell ", xP,yP,zP,"; Grain temperature: "&
                       &, grid%Tdust(nS,ai,cellP), " species ", grainLabel(nS), " size:", grainRadius(ai)

                  ! find weighted mean
                  grid%Tdust(nS,0,cellP) = grid%Tdust(nS,0,cellP)+&
                       & grid%Tdust(nS,ai,cellP)*grainWeight(ai)

               end do
               
               grid%Tdust(0,0,cellP) = grid%Tdust(0,0,cellP)+&
                    & grid%Tdust(nS,0,cellP)*grainAbun(nspU,nS)

            end do
            

          end subroutine getDustT

          function resLineHeating(sizeP,speciesP, icompP)
            implicit none

            real                :: resLineHeating ! dust heating due to res lines

            real                :: Gline   ! energy is the line [erg sec^-1]
            real                :: heat    ! sub calculation var
            
            integer, intent(in) :: sizeP   ! size pointer
            integer, intent(in) :: speciesP! species pointer
            integer, intent(in) :: icompP! dust componenent pointer

            integer             :: iL      ! line counter
            integer             :: imul    ! multiplet counter 

            resLineHeating = 0.
            do iL = 1, nResLines

               Gline = 0.
               do imul = 1, resLine(iL)%nmul
                  if (resLine(iL)%elem==1) then

                     if ( resLine(iL)%ion == 1 .and. resLine(iL)%moclow(imul)==1 .and. resLine(iL)%mochigh(imul)==2 ) then

                        ! fits to Storey and Hummer MNRAS 272(1995)41
                        Gline = Gline + 10**(-0.897*log10(TeUsed) + 5.05)*grid%elemAbun(grid%abFileIndex(xP,yP,zP),1)&
                             & *1.e-25*ionDenUsed(elementXref(1),2)*&
                             & NeUsed

                     else
                        
                        print*, "! resLineHeating: [warning] only dust &
                             &heating from H Lyman alpha and resonance lines from heavy"
                        print*, "elements is implemented in this version. please contact author B. Ercolano -1-", &
                             &iL, resLine(iL)%elem, &
                             & resLine(iL)%ion

                     end if
                     
                  else if (resLine(iL)%elem>2) then
                     
                     Gline = Gline+forbiddenLines(resLine(iL)%elem,resLine(iL)%ion, &
                          &resline(iL)%moclow(imul), resline(iL)%mochigh(imul))

                  else
                  
                     print*, "! resLineHeating: [warning] only dust heating from H ",&
                          &"Lyman alpha and resonance lines from heavy"
                     print*, "elements is implemented in this version. please contact author B. Ercolano -2-", &
                          &iL, resLine(iL)%elem, &
                             & resLine(iL)%ion

                  end if


               end do

               Gline=Gline/Pi

               heat = Gline*grid%Hden(cellP)* &
                    & (1.-grid%fEscapeResPhotons(cellP, iL))*&
                    & xSecArray(dustAbsXSecP(speciesP,sizeP)+resLine(iL)%nuP-1)/&
                    & (grid%Ndust(cellP)*&
                    & absOpacSpecies(speciesP,resLine(iL)%nuP))

               ! Harrington Monk and Clegg 1988 (section 3.2)
               resLineHeating = resLineHeating + heat

               dustHeatingBudget(grid%abFileIndex(xP,yP,zP),iL) = &
                    & dustHeatingBudget(grid%abFileIndex(xP,yP,zP),iL) + &
                    & heat*grainWeight(sizeP)*grainAbun(icompP,speciesP)*grid%Ndust(cellP)
               dustHeatingBudget(0,iL) = dustHeatingBudget(0,iL) + &
                    & heat*grainWeight(sizeP)*grainAbun(icompP,speciesP)*grid%Ndust(cellP)


            end do
            
          end function resLineHeating


        end subroutine updateCell

end module update_mod

