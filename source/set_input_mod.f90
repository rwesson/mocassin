! Copyright (C) 2005 Barbara Ercolano 
!
! Version 2.02
module set_input_mod
    use common_mod
    implicit none

    character(len=50)   :: shapeDiffuse     ! diffuse spectrum shape

    contains

    subroutine readInput()
        implicit none

        ! local variables 

        integer             :: err              ! allocation error status
        integer             :: i, j             ! counters
        integer             :: ios              ! I/O error status
        integer             :: iValue           ! integer value corresponding to keyword
        integer             :: maxLim = 1000    ! max limit for reading loop

        real                :: rValue           ! real value corresponding to keyword

        character(len=50)   :: cValue           ! character value corresponding to keyword
        character(len=50)   :: in_file          ! input file
        character(len=50)   :: keyword          ! input parameter keyword
        character(len=50)   :: multiPhotoSources ! stars file
        
        
        ! set default values and set non oprional values to 0 or 0. or "zero"

        lgIsotropic   = .false.       
        lgRecombination = .false.
        lgAutoPackets = .false.
        lgMultiDustChemistry = .false.
        lgMultiChemistry = .false.
        lgTalk        = .false.
        lgHdenConstant= .false.
        lgDfile       = .false.
        lgDebug       = .false.
        lgDlaw        = .false.
        lgDust        = .false.
        lgDustConstant= .false.
        lgGas         = .true.
        lgMdMg        = .false.
        lgMdMh        = .false.
        lgNeInput     = .false.
        lgNeutral     = .true.
        lgOutput      = .false. 
        lgPlaneIonization = .false.
        lg1D          = .false.
        lgDustScattering = .true.
        lgSymmetricXYZ= .false.
        lgEcho        = .false.
        lgNosource    = .false.
        lgforceTDust = .false.

        nPhotonsDiffuse = 0        
        nStars        = 0
        nxIn(:)       = 0
        nyIn(:)       = 0
        nzIn(:)       = 0
        nxIn(1)       = 30
        nyIn(1)       = 30
        nzIn(1)       = 30        
        maxIterateMC  = 30
        maxPhotons    = 0
        nbins         = 600        
        MdMgValue     = 0.
        MdMgFile      = "none"               
        nAngleBins    = 0        
        nGrids        = 1
        NdustValue     = 0.
        NdustFile      = "none"
        nstages        = 7

        ! qheat
        Tmax           =700.
        nTbins         =300
        lgWritePss     =.false.
        minaQHeat     = 1.e-3
        minConvQHeat  = 99.

        fillingFactor = 1.
        contCube      = -1.
        convIncPercent= 0.
        convWriteGrid = 0.
        nPhotIncrease = 0.

        densityLaw    = 0.

        multiPhotoSources = "none"
        densityFile   = "none"
        dustFile      = "none"
        gridList      = "none"

        nu0           = 0.  
        Ldiffuse      = 0. 
        Tdiffuse      = 0.
        shapeDiffuse  = "none"
        Hdensity      = 0.
        H0Start       = 3.e-5
        LPhot         = 0.
        minConvergence= 95.
        NeStart       = 0.
        nuMax         = 15.
        nuMin         = 1.001e-5
        nuStepSize    = 0.075
        resLinesTransfer = 101.
        Rnx           = -1.
        Rny           = -1.
        Rnz           = -1.
        R_in          = -1.
        R_out         = 0.
        TeStart       = 10000.
        XHILimit      = 0.05

        ! ask the user to input file name and check for errors (three chances)
        do i=1,3
             in_file = 'input/input.in'

             close(10)
             open (unit = 10, file = in_file, status = "old", &
                  &position="rewind",  iostat = ios, action="read")

             if (ios==0) exit
 
             print *, "! readInput: can't open file - please try again", ios

        end do

        if (ios /= 0) then
            print*, "!out_in: can't open file - terminating"
            stop
        end if

        print*, taskid, "input file: ", in_file

        do i = 1, maxLim

            ! check if the end of file has been reached
            read(unit=10, fmt=*, iostat=ios) keyword
            if (ios < 0) exit ! end of file reached

            select case (keyword)
            case ("home")
               backspace 10
               read(unit=10,fmt=*, iostat=ios) keyword, home
            case ("TDust")           ! BEKS2010
               backspace 10
               read(unit=10,fmt=*, iostat=ios) keyword, forceTDust
               lgforceTDust = .true.
               print*,keyword,lgforceTDust,forceTDust
            case ("echo")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, echot1, echot2, echoTemp
               print "(1X,A4,1X,3es12.4)",keyword,echot1,echot2,echoTemp
               print*, '! readInput [Reminder]: Light echo option - is the inclination correctly specified? [see Readme]' 
               print*, '! readInput [Warning]: Light echo option - '
               print*, '!   The actual volume of material illuminated by the echo is probably'
               print*, '!   much smaller than that illuminated in this grid.  Be sure to perform'
               print*, '!   volume correction. See echo.out file for volume illuminated in this grid.'
!               echot1 = echot1 * 2.59020684e15 !ct
!               echot2 = echot2 * 2.59020684e15
               lgEcho = .true.                 
            case("NoSourceSED")
               lgNosource = .true.
               print*, keyword, lgNosource
            case("isotropicScattering")
               lgIsotropic = .true.
               print*, keyword, lgIsotropic
            case("2D")
               lg2D = .true.
               print*, keyword, lg2D
            case("nstages")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, nstages
               if (nstages>10) then
                  print*, " !out_in: Max 10 ionisation stages allowed. &
                       &If more are required (and you have the relevant atomic data), &
                       & please update the data/ directory (and let B. Ercolano know ;-) &
                       & and edit the data/fileNames.dat file"  
                  stop
               end if
               print*, keyword, nstages
            case ("getEquivalentTau")
               lgEquivalentTau = .true.               
               print*, keyword, lgEquivalentTau
            case ("diffuseSource") 
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, Ldiffuse, Tdiffuse, shapeDiffuse, nPhotonsDiffuse, emittingGrid
               print*, keyword, Ldiffuse, Tdiffuse, shapeDiffuse, nPhotonsDiffuse, emittingGrid
            case ("traceHeating")
               lgTraceHeating = .true.
               print*, keyword, lgTraceHeating
            case("quantumHeatGrain")
               lgQheat = .true.
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword,minaQHeat,minConvQHeat
               print*, keyword, minaQHeat,minConvQHeat
            case("quantumHeatGrainParameters") 
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword,Tmax,nTbins,lgWritePss
               print*, keyword,Tmax,nTbins,lgWritePss              
            case ("noPhotoelectric")
               print*, keyword
               lgPhotoelectric = .false.
            case ("multiPhotoSources") 
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword,multiPhotoSources
               print*, keyword,multiPhotoSources
               lgMultiStars=.true.
            case("continuumCube")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, contCube(1), contCube(2)
               print*, keyword, contCube
            case("noScattering")
               print*, keyword
               lgDustScattering = .false.
            case("resLinesTransfer")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, resLinesTransfer
               print*, keyword, resLinesTransfer
            case("multiGrids")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, nGrids, gridList
               print*, keyword, nGrids, gridList
               if (nGrids > maxGrids) then
                  print*, "! readInput: please increse the maxGrids parameter in constants_mod.f90", &
                       & nGrids, maxGrids
                  stop
               end if
               call readGridList(gridList)
            case("recombinationLines")
               lgRecombination=.true.
               print*, keyword
            case ("fillingFactor") 
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, fillingFactor
               print*, keyword, fillingFactor
            case ("slit")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, dxSlit, dySlit
               print*, keyword, dxSlit, dySlit
            case ("inputNe")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword
               lgNeInput = .true.
               print*, keyword
            case ("edges")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, Rnx, Rny, Rnz
               print*, keyword,  Rnx, Rny, Rnz
            case ("output")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword
               lgOutput = .true.
               print*, keyword
            case ("multiChemistry")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, nAbComponents
               lgMultiChemistry = .true.
               backspace(10)
               allocate(abundanceFile(1:nAbComponents))
               read(unit=10, fmt=*, iostat=ios) keyword, nAbComponents, &
                    & (abundanceFile(j), j=1,nAbComponents)
               print*, keyword, nAbComponents, (abundanceFile(j), j=1,nAbComponents)
            case ("multiDustChemistry")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, nDustComponents
               lgMultiDustChemistry = .true.
               lgDust = .true.
               backspace(10)
               allocate(dustSpeciesFile(1:nDustComponents))
               allocate(nSpeciesPart(1:nDustComponents))
               read(unit=10, fmt=*, iostat=ios) keyword, nDustComponents, (dustSpeciesFile(j), j=1,nDustComponents), dustFile(2)
               print*, keyword, nDustComponents, (dustSpeciesFile(j), j=1,nDustComponents)
               print*, dustFile(2)
            case ("planeIonization")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, meanField, nu0
               lgPlaneIonization = .true.
               print*, keyword, meanField, nu0
               allocate(Lstar(1))
               Lstar=0.
            case ("autoPackets")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, convIncPercent, nPhotIncrease, maxPhotons
               lgAutoPackets = .true.
               print*, keyword, convIncPercent, nPhotIncrease, maxPhotons
            case ("oneD")
               backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword
                lg1D = .true.
                print*, keyword        
                print*, 'ERROR: oneD option is not currently avalable!!!!'
                stop
             case ("densityFile")
                 backspace 10
                 read(unit=10, fmt=*, iostat=ios) keyword, densityFile
                 lgDfile = .true.
                 print*, keyword, densityFile
                   open(unit=11,file=densityFile,status='old')
                   read(unit=11, fmt=*, iostat=ios) keyword
                   if (keyword=="#") then
                      backspace 11
                      read(unit=11,fmt=*, iostat=ios) keyword,nxin(1),nyin(1),nzin(1)
                      print*,'grid axes read in from file:',trim(densityFile),nxin(1),nyin(1),nzin(1)
                   endif
                   close(11)
             case ("dustFile")
                 backspace 10
                 read(unit=10, fmt=*, iostat=ios) keyword, dustFile(1), dustFile(2)
                 print*, keyword, trim(dustFile(1))," | ",trim(dustFile(2))," |"
                 allocate(dustSpeciesFile(1:1))
                 allocate(nspeciesPart(1:1))
                 nDustComponents = 1
                 dustSpeciesFile(1) = dustFile(1)
             case ("debug")
                 backspace 10
                 read(unit=10, fmt=*, iostat=ios) keyword
                 lgDebug = .true.
                 print*, keyword
            case ("nbins")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, nbins
                print*, keyword, nbins                
            case ("talk")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword
                lgTalk = .true.
                print*, keyword
            case ("symmetricXYZ")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword
                lgSymmetricXYZ = .true.
                print*, keyword
            case ("contShape")
                backspace 10
                nstars=1
                allocate(contShape(0:1))
                allocate(contShapeIn(0:1))
                allocate(spID(0:1))
                allocate(tstep(0:1))
                contShape = 'none'
                read(unit=10, fmt=*, iostat=ios) keyword, contShape(1)
                if (contShape(1) == 'powerlaw') then
                   backspace 10
                   read(unit=10, fmt=*, iostat=ios) keyword, contShape(1), pwlIndex                   
                   print*, keyword, contShape(1), pwlIndex                   
                else
                   print*, keyword, contShape(1)
                end if
                if (contShape(1) == 'blackbody' .or. contShape(1) == 'powerlaw') then
                   spID(1) = contShape(1)
                else
                   spID(1) = 'rauch'
                end if
                tstep = 0.                
            case ("nebComposition")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, cValue
                if (cValue=="noGas") lgGas=.false.
                print*, keyword, cValue
                allocate(abundanceFile(1:nAbComponents))
                abundanceFile = cValue
            case ("Hdensity") 
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, Hdensity
                lgHdenConstant = .true.
                print*, keyword, Hdensity
            case ("densityLaw")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, (densityLaw(j),j=1,3)
                lgDlaw = .true.
                print*, keyword, densityLaw
            case ("maxIterateMC")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, maxIterateMC, minConvergence
                if (taskid==0) print*, keyword, maxIterateMC, minConvergence
            case ("nPhotons")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, nPhotonsTot
                print*, keyword, nPhotonsTot
            case ("nx")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, nxin(1)
            case ("ny")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, nyin(1)
            case ("nz")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, nzin(1)
             case ("starPosition")
                backspace 10
                allocate(starPosition(1))
                read(unit=10, fmt=*, iostat=ios) keyword, starPosition(1)%x, starPosition(1)%y, starPosition(1)%z
                print*,  keyword, starPosition(1)%x, starPosition(1)%y, starPosition(1)%z
             case ("TeStart")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, TeStart
                print*, keyword, TeStart
            case ("NeStart")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, NeStart
                print*, keyword, NeStart
            case ("TStellar")
                backspace 10
                allocate(TStellar(0:1))
                read(unit=10, fmt=*, iostat=ios) keyword, TStellar(1)
                print*, keyword, TStellar(1)
            case ("convLimit")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword,  XHILimit
                print*, keyword, XHILimit
            case ("H0Start")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, H0Start
                print*, keyword, H0Start
            case ("LStar")
                backspace 10
                allocate(Lstar(0:1))
                read(unit=10, fmt=*, iostat=ios) keyword, LStar(1)
                print*, keyword, LStar(1)
            case ("LPhot")
                backspace 10
                if (.not. associated(LStar)) allocate(Lstar(0:1))                
                read(unit=10, fmt=*, iostat=ios) keyword, LPhot
                print*, keyword, LPhot
            case ("nuMax")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, nuMax
                print*, keyword, nuMax
            case ("nuMin")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, nuMin
                print*, keyword, nuMin
            case ("Rin")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, R_in
                print*, keyword, R_in
            case ("Rout")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, R_out
                print*, keyword, R_out
            case ("MdMg")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, keyword
                if (keyword=="constant") then
                   backspace 10
                   read(unit=10, fmt=*, iostat=ios) keyword, keyword, MdMgValue
                   print*, 'MdMg constant ', MdMgValue
                   lgDustConstant = .true.
                   lgDust = .true.
                else if (keyword=="file") then
                   backspace 10
                   read(unit=10, fmt=*, iostat=ios) keyword, keyword, MdMgFile
                   print*, 'MdMg file ', MdMgFile
                   lgDust = .true.
                   open(unit=11,file=MdMgFile,status='old')
                   read(unit=11, fmt=*, iostat=ios) keyword
                   if (keyword=="#") then
                      backspace 11
                      read(unit=11,fmt=*, iostat=ios) keyword,nxin(1),nyin(1),nzin(1)
                      print*,'grid axes read in from file:',trim(MdMgFile),nxin(1),nyin(1),nzin(1)
                   endif
                   close(11)
                else
                   print*, "! readInput: invalid keyword in MdMg field ", keyword
                   stop
                end if
                lgMdMg = .true.
            case ("MdMh")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, keyword
               if (keyword=="constant") then
                   backspace 10
                   read(unit=10, fmt=*, iostat=ios) keyword, keyword, MdMgValue
                   print*, 'MdMh constant ', MdMgValue
                   lgDustConstant = .true.
                   lgDust = .true.
                else if (keyword=="file") then
                   backspace 10
                   read(unit=10, fmt=*, iostat=ios) keyword, keyword, MdMgFile
                   print*, 'MdMh file ', MdMgFile
                   lgDust = .true.
                   open(unit=11,file=MdMgFile,status='old')
                   read(unit=11, fmt=*, iostat=ios) keyword
                   if (keyword=="#") then
                      backspace 11
                      read(unit=11,fmt=*, iostat=ios) keyword,nxin(1),nyin(1),nzin(1)
                      print*,'grid axes read in from file:',trim(MdMgFile),nxin(1),nyin(1),nzin(1)
                   endif
                   close(11)
                else
                   print*, "! readInput: invalid keyword in MdMh field ", keyword
                   stop
                end if
                lgMdMh = .true.
            case ("Ndust")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, keyword
                if (keyword=="constant") then
                   backspace 10
                   read(unit=10, fmt=*, iostat=ios) keyword, keyword, NdustValue
                   print*, 'Ndust constant ', NdustValue
                   lgDustConstant = .true.
                   lgDust = .true.
                else if (keyword=="file") then
                   backspace 10
                   read(unit=10, fmt=*, iostat=ios) keyword, keyword, NdustFile
                   print*, 'Ndust file ', NdustFile
                   open(unit=11,file=NdustFile,status='old')
                   read(unit=11, fmt=*, iostat=ios) keyword
                   if (keyword=="#") then
                      backspace 11
                      read(unit=11,fmt=*, iostat=ios) keyword,nxin(1),nyin(1),nzin(1)
                      print*,'grid axes read in from file:',trim(NdustFile),nxin(1),nyin(1),nzin(1)
                   endif
                   close(11)
                   
                   lgDust = .true.
                else
                   print*, "! readInput: invalid keyword in Ndust field ", keyword
                   stop
                end if
            case ("writeGrid")
                backspace 10
                read(unit=10, fmt=*, iostat=ios) keyword, convWriteGrid
                if (taskid==0) print*, keyword, convWriteGrid
            case ("inclination")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, nAngleBins 
               if (nAngleBins>2) then
                  print*, '! readInput: only two inclination anges are allowed per simulation. &
                       & If more inclinations are required please run mocassinWarm for the other&
                       & viewing angles.'
                  stop
               end if
               allocate(viewPointTheta(0:nAngleBins), stat=err)
               if (err /= 0) then
                  print*, '! readInput: allocation error for viewPointTheta pointer'
                  stop
               end if
               viewPointTheta=0.
               allocate(viewPointPhi(0:nAngleBins), stat=err)
               if (err /= 0) then
                  print*, '! readInput: allocation error for viewPointPhi pointer'
                  stop
               end if
               viewPointPhi=0.
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, nAngleBins, (viewPointTheta(j), viewPointPhi(j), j=1,nAngleBins)
               if (taskid==0) print*, keyword, nAngleBins, (viewPointTheta(j), viewPointPhi(j), j = 1, nAngleBins)
            case ("dustMass")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, inputDustMass
               lginputDustMass = .true.  
            case ("gasMass")
               backspace 10
               read(unit=10, fmt=*, iostat=ios) keyword, inputGasMass
               lginputGasMass = .true.  
               print*,"Using gasMass keyword. Be sure to specify gas density file or Hdensity 1.0"
            case default
                print*, "! readInput: unrecognised keyword in model parameter input file", &
&                        in_file, keyword
            end select

         end do

        close(10)

        allocate(lgDataAvailable(3:nElements, 1:nstages), stat=err)
        if (err /= 0) then
           print*, '! readInput: allocation error for lgDataAvailable pointer'
           stop
        end if

        if (lginputDustMass) then
          inputDustMass = inputDustMass / 5.028e11
          if (lgSymmetricXYZ) then
            inputDustMass = inputDustMass / 8
          end if
        end if

        if (lginputGasMass) then
          inputGasMass = inputGasMass / 5.028e11
          if (lgSymmetricXYZ) then
            inputGasMass = inputGasMass / 8
          end if
        end if

        if (lginputDustMass .and. .not.lgDust) then
           print*,"! readInput: you have asked to scale dustMass but there is no dust"
           stop
        endif
        if (lginputGasMass .and. .not.lgGas) then
           print*,"! readInput: you have asked to scale gasMass but there is no gas"
           stop
        endif
        
        if (lginputDustMass .and. lginputGasMass) then
           if (lgDust .and. lgGas) then
              if (lgMdMg .or. lgMdMh) then
                 print*,'! readInput: dust mass is defined by your gas-to-dust ratio'
                 print*,'! readInput: please remove dustMass keyword from input.in'
                 stop
              else
                 print*,'! readInput: Simulation will scale gas to gasMass and dust to dustMass.'
                 print*,'! readInput: This will change your gas-to-dust ratio so BE SURE THIS IS WHAT YOU WANT...'
              endif
           endif
        endif

        if (lgMultiStars) then
           call setMultiPhotoSources(multiPhotoSources)
        else if (lgMultiStars .and. nStars==1) then
           print*, '! readInput: multiPhotoSources keyword and Lstar, Tstellar, ContShape are cannot be specified together'
           stop
        else
           allocate(nPhotons(1))
           allocate(deltaE(0:1))
           nPhotons(1)=nPhotonsTot

           if (Ldiffuse>0. .and. nPhotonsDiffuse>0) then
              if (nStars<1) then
                 if(.not.associated(contShape)) allocate(contShape(0:0))
                 if(.not.associated(spID)) allocate(spID(0:0))
                 if(.not.associated(LStar)) allocate(LStar(0:0))
              end if
              deltaE(0) = Ldiffuse/nPhotonsDiffuse
              contShape(0) = shapeDiffuse
              if (contShape(0) == 'blackbody' .or. contShape(0) == 'powerlaw') then
                 spID(0) = contShape(0)
              else
                 spID(0) = 'rauch'
              end if
           end if

           if (.not.associated(starPosition)) then
              allocate (starPosition(1))
              starPosition(1)%x=0.
              starPosition(1)%y=0.
              starPosition(1)%z=0.
              print* , '! readInput: [talk] ionising source positioned at 0.,0.,0.'
           end if
        end if

        if (.not.associated(contShape)) then
           print*, "! readInput: no information about the continuum shape - contShape -"
           stop
        end if
        if (.not.associated(TStellar)) then
           print*, "! readInput: no information about the stellar effective temperature - Tstellar -"
           stop
        end if
        if (.not.associated(nPhotons)) then
           print*, "! readInput: no information about the number of energy packets to be used - nPhotons -"
           stop
        end if
        if (.not.associated(Lstar) .and. Lphot<=0.) then
           print*, "! readInput: no information about the stellar luminosity - Lstar or Lphot -"
           stop
        end if

        contShapeIn = contShape

        print*, "mother nx, ny, nz" , nxin(1), nyin(1), nzin(1)

        ! check for missing or invalid values in the model parameters input file        
        if (lg1D) then
           nyIn = 1
           nzIn = 1
           lgSymmetricXYZ = .true.
        else if (lgPlaneIonization .and. lg2d) then
           print*, "! readInput: planeIonization and 2D options are not compatible"
           stop
        else if (resLinesTransfer >= minConvergence .and. resLinesTransfer /= 101. .and. lgDust) then
           print*, "! readInput: the min convergence level assigned to the calculation of &
                &the resonant lines transfer is higher or equal to that assigned by maxIterateMC"
           stop
        else if (resLinesTransfer >= 100.5 .and. lgDust .and. lgGas) then
           print*, "! readInput: [warning] both dust and gas processes have been activated, but the &
                &resonanceLinesTransfer keyword is absent. The dust temperatures may be underestimated!!!!"
        else if (resLinesTransfer <= 100 .and. (.not.lgDust .or. .not.lgGas)) then
           print*, "! readInput: [warning] both dust and gas processes must be activated, to use the&
                &resonanceLinesTransfer keyword. Reset to 101."
           resLinesTransfer = 101.
        else if (.not.associated(abundanceFile) .and. lgGas) then
           print*, "! readInput: no elemental abundance file has been specified"
           stop
        else if ((lgMdMh .and. lgMdMg) .and. (NdustFile /= "none" .or. NdustValue /=0) )then
           print*, "! readInput: MdMg and Ndust cannot be specified together"
           stop
        else if ( (.not.lgDust) .and. (.not.lgGas) )then
           print*, "! readInput: the grid is completely empty. no gas or dust present.",  &
                lgGas,lgDust
           stop                      
        else if ( lgDebug .and. .not.(lgGas) ) then
           print*, "! readInput: debugging mode can only be used when a gas component is present. &
                & debugging mode will be turned off."
           lgDebug = .false.
        else if (fillingFactor > 1.) then
           print*, "! readInput: filling factor, epsilon, specified is greater than unity.", &
                & fillingFactor
           stop
        else if ((.not.lgDfile) .and. (Rnx<0. .or. Rny<0. .or. Rnz<0.) .and. (NdustFile == "none")) then
           print*,  "! readInput: Grid edges unspecified or non-valid grid &
                &edges", Rnx,Rny,Rnz
           stop
        else if (lgDfile .and. (Rnx>=0. .or. Rny>=0. .or. Rnz>=0.) ) then
           print*,  "! readInput: [warning] when the density distribution is&
                & specified via an external density file, no grid edges need &
                & to be specified. The specified values will be ignored."
        else if (lgPlaneIonization .and. .not.lgDfile) then
           print*, "! readInput: plane ionizing field specified - must use user defined &
                &density distribution via the densityFile option"
           stop
        else if (lgPlaneIonization .and. lgSymmetricXYZ) then
           print*, "! readInput [warning]: plane ionizing field specified - cannot use&
                & symmetricXYZ - removed."
           lgSymmetricXYZ = .false.
        else if (TStellar(1) == 0. .and. Tdiffuse==0. .and. contShape(1) /= 'powerlaw') then
            print*, "! readInput: TStellar missing from model parameter input file" 
            stop
        else if (R_in < 0.) then
            print*, "! readInput: Invalid Rin parameter in the input file", R_in
            stop
       else if (LPhot == 0. .and. LStar(1) == 0. .and. .not.lgPlaneIonization & 
            & .and. Ldiffuse<=0.) then
            print*, "! readInput: LPhot and LStar missing from model parameter input file"
            stop
        else if (LPhot /= 0. .and. LStar(1) /= 0.) then 
            print*, "! readInput: [warning] LPhot and LStar both entered. "
        else if (maxIterateMC < 1) then
            print*, "! readInput: invalid maxIterateMC in model parameter input &
                 & file", in_file, maxIterateMC
            stop
         else if (minConvergence < 0. .or. minConvergence > 100.) then
            print*, "! readInput: invalid minConvergence input (must be between 0. and 100.) ", & 
                 &  in_file, minConvergence
            stop
        else if (nPhotons(1) < 0) then
            print*, "! readInput: invalid nPhotons in model parameter input &
                 & file", in_file, nPhotons
            stop
        else if (lgDust) then
            if ((dustFile(1) == "none" .or. dustFile(2)=="none") .and. &
                 &.not.lgMultiDustChemistry) then
               print*, "! readInput: dust files were not specified; please include dustFile field &
                    & in the input files and specify name and location of dust files"
               stop
            end if
        else if (nxin(1) < 1 .or. nyin(1) < 1 .or. nzin(1) < 1) then
            print*, "! readInput: invalid grid boundaries in model parameter input & 
                 & file", in_file, nxin(1), nyin(1), nzin(1)
            stop
        else if (TStellar(1) < 0.) then
            print*, "! readInput: invalid TStellar  in model parameter input &
                 & file", in_file, TStellar
            stop
        else if (XHILimit < 0.) then
            print*, "! readInput: invalid XHILimit in model parameter input &
                 & file", in_file, XHILimit
            stop
        else if (LPhot < 0.) then
            print*, "! readInput: invalid LPhot  in model parameter input &
                 & file", in_file, LPhot
            stop
        else if (LStar(1) < 0.) then
            print*, "! readInput: invalid LStar  in model parameter input &
                 & file", in_file, LStar
            stop
        else if (nuMax < 0.) then
            print*, "! readInput: invalid nuMax  in model parameter input &
                 & file", in_file, nuMax
            stop
        else if (nuMin < 0.) then
            print*, "! readInput: invalid nuMin  in model parameter input &
                 & file", in_file, nuMin
            stop
        else if (nuMin > 1.4e-6) then
            print*, "! readInput [waring]:  nuMin  in model parameter input &
                 & file is too large to include 5GHz flux calculations", in_file,nuMin            
        else if (R_out < 0.) then
            print*, "! readInput: invalid R_out  in model parameter input &
                 & file", in_file, R_out
            stop
        else if (R_in <0.) then
            print*, "! readInput: invalid R_in  in model parameter input &
                 & file", in_file, R_in
            stop
        else if ((R_in >= R_out) .and. (R_out > 0.) ) then
            print*, "! readInput: the inner radius must be smaller than &
                 &the outer radius ", in_file, R_in, R_out
            stop
        else if (lgMultiChemistry .and. .not.lgDfile) then
           print*, "! readInput: the MultiChemistry option has been specified &
                &therefore an user density file of the type x,y,z,den,abIndex must also &
                &be specified."
           stop
        end if

        TStellar(0) = Tdiffuse
        
        contains

          subroutine readGridList(filename)
            implicit none
            
            real                          :: skipR

            character(len=50), intent(in) :: filename
            character(len=50)             ::  skipC            

            integer                       :: i,iG
            integer                       :: err, ios
            integer                       :: skip
            
            close(17)
            open (unit = 17, file = filename, status = "old", &
                 &position="rewind",  iostat = ios, action="read")
            if (ios/=0) then
               print*, "! readGridList: can't open file for reading ", filename
               stop
            end if
            
            do i = 2, nGrids
              read(17, *) skip, nxIn(i), nyIn(i), nzIn(i), skipC, skipR
              read(17, *) skipR, skipR, skipR, skipR, skipR, skipR

            end do
            close(17)

         end subroutine readGridList

       end subroutine readInput
       
       subroutine setMultiPhotoSources(infile)
        implicit none
        
        character(len=50), intent(in) :: infile

        integer                       :: ios, i, iloop,nSafeLimit=10000

        real                          :: E0


        close(13)
        open(file=infile, unit=13, iostat=ios, action="read")
        if (ios /= 0) then
           print*, "! setMultiPhotoSources: can't open ionising sources file", infile
           stop
        end if
        
        read(13,*) nStars
        print*, "Multiple Ionising sources", nStars
        print*, "(i, Tstellar, Lstar, ContShape, Nphotons, position)"

        allocate(deltaE(0:nStars))
        allocate(spID(0:nStars))        
        allocate(tStep(nStars))        
        allocate(TStellar(0:nStars))        
        allocate(Lstar(nStars))
        allocate(ContShape(0:nStars))        
        allocate(ContShapeIn(0:nStars))        
        allocate(nPhotons(nStars))        
        allocate(starPosition(nStars))       

        TStellar=0.
        Lstar=0.
        ContShape='none'
        ContShapeIn='none'
        spID='none'
        tStep=0.
        nPhotons=0.

        do i = 1, nStars
           read(13,*) TStellar(i), Lstar(i), ContShape(i), starPosition(i)%x, starPosition(i)%y, starPosition(i)%z, &
                & spID(i), tStep(i)
           contShapeIn(i) = contShape(i)
        end do
        
        
        print*, 'i, Tstellar(i), Lstar(i), contShape(i), Nphotons(i), starPosition(i), deltaE(i)'
        do i = 1, nStars
           nPhotons(i) = NphotonsTot/nStars
           deltaE(i) = Lstar(i)/nPhotons(i)
           print*, i, Tstellar(i), Lstar(i), contShape(i), Nphotons(i), starPosition(i), deltaE(i)
        end do                

        if (Ldiffuse>0. .and. nPhotonsDiffuse>0) then
           deltaE(0) = Ldiffuse/nPhotonsDiffuse
           contShape(0) = shapeDiffuse
           if (contShape(0) == 'blackbody' .or. contShape(0) == 'powerlaw') then
              spID(0) = contShape(0)
           else
              spID(0) = 'rauch'
           end if
        end if

      end subroutine setMultiPhotoSources
        
end module set_input_mod

