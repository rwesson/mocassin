! Copyright (C) 2005 Barbara Ercolano 
!  
! Version 2.02
module common_mod
    use constants_mod
    use vector_mod
    implicit none

    ! fluorescence variables
    integer          :: nFluo ! # of fluorescent lines to be transfered
    integer          :: nVP   ! # of fluorescent viewing angles

    character(len=4), pointer :: fluorescenceLabel(:) ! ID labels of fluorescent transitions
    real(kind=dp), pointer    :: fluorescenceVP(:,:) ! viewing angles for fluorescence cube &
                                            ! (for the first dimension 1=theta;2=phi)

    ! MPI variables
    integer, pointer :: lgConvergedTemp(:)  ! temporary converged? flag
    integer, pointer :: lgBlackTemp(:)      ! temporary converged? flag

    logical         :: lg2D=.false.             ! 2D?
    logical         :: lgIsotropic

    logical             :: lgMultiStars=.false.
    logical         :: lgEquivalentTau          ! calculate equivalent tau?
    logical         :: lgWarm=.false.           ! warm started?
    logical         :: lgFluorescence=.false.   ! fluorescence run?
    logical         :: lgNeInput=.false.        ! Ne distribution entered
    logical         :: lgResLinesFirst = .true. ! first time the res lines transfer proc is called? 
    logical         :: lgPhotoelectric = .true. ! photoelectric effect on
    logical         :: lgQheat = .false.        ! temperature spiking on? 
    logical         :: lgWritePss=.false.       ! write Pss for qHeat file?
    logical         :: lgTraceHeating = .false. ! trace thermal balance?
    logical         :: lgGrainSpiking = .true.  ! temperature spiking for small grains
    logical         :: lgEcho = .false.         ! light time travel included?
    logical         :: lgforceTDust = .false.   !
    logical         :: lgNosource = .false.     ! exclude sources from SED?
    logical         :: lginputDustMass = .false.! user sets input dust mass?
    logical         :: lginputGasMass = .false. ! user sets input gas mass?

    real(kind=dp), pointer :: gSca(:)                 ! gSca(freq)

    real(kind=dp), pointer   :: TdustTemp(:,:,:)     ! temporary dust temperature array (species,size,cell)
    real(kind=dp), pointer   :: ionDenTemp(:,:,:)    ! temporary ion density array    
    real(kind=dp), pointer   :: NeTemp(:)            ! temporary electron density array
    real(kind=dp), pointer   :: TeTemp(:)            ! temporary electron Temperature array    
    real(kind=dp)            :: HeIrecLineCoeff(34,3,4)
    real(kind=dp), pointer   :: KNsigmaArray(:,:)    ! Klein Nishina calc PDFs (nbins,180)
    real(kind=dp), pointer   :: KNsigmaT(:)          ! Klein Nishina cross-sections integrated over the solid angle
    real(kind=dp), pointer   :: PcompArray (:,:) 
    real(kind=dp), pointer   :: twoDscaleJ(:)

    integer         :: radio4p9GHzP   ! 
      
    integer         :: ierr           ! MPI error status
    integer         :: numtasks       ! total # of processes
    integer         :: taskid         ! process identification #
    real(kind=dp)            :: starttime      ! start time [sec]

    real(kind=dp)            :: endtime        ! end time [sec]    
    real(kind=dp)            :: absInt         ! total number of absorption events
    real(kind=dp)            :: scaInt         ! total number of scattering events
    real(kind=dp)            :: echot1=0., echot2=0.   ! light time travel parameters
    real(kind=dp)            :: echoTemp=0.       ! temperature outside the illuminated parabula
    real(kind=dp)            :: dTheta                  ! 
    real(kind=dp)            :: dPhi                    ! 
    real(kind=dp)            :: nu0                     ! 
    real(kind=dp)            :: nu0Add                  !
    real(kind=dp)            :: totalDustMass 
    real(kind=dp)            :: totalGasMass 
    real(kind=dp)            :: inputDustMass           ! for when user sets desired dust mass in input
    real(kind=dp)            :: inputGasMass           ! for when user sets desired gas mass in input
    real(kind=dp)            :: forceTDust = 0.         ! for TDust feature
    real(kind=dp)            :: convPercent=0.          ! total convergence percentage
    real(kind=dp)            :: pwlIndex = 0.           ! power law input spectrum index
    real(kind=dp)            :: pwlMin=0.               ! power law lower cutoff [Ryd]
    real(kind=dp)            :: pwlMax=0.               ! power law higher cutoff [Ryd]
    real(kind=dp)            :: totPercent=0.           ! 
    real(kind=dp)            :: minaQHeat = 0.          ! min radius for Q Heat in [um]
    real(kind=dp)            :: NeUsed, TeUsed, HdenUsed ! local properties of the gas
    real(kind=dp)            :: Tmax                    ! max temp for qheat
    real(kind=dp)            :: auger(30,30,10,10)      ! auger yields



    integer, pointer  :: HINuEdgeP(:)     ! pointers to the HI, HeI and HeII 
    integer, pointer  :: HeINuEdgeP(:)    ! series edges in nuArray
    integer, pointer  :: HeIINuEdgeP(:)   !


    integer :: nTbins                     ! # of temp bins in qheat
    integer :: nlimGammaHI,nlimGammaHeI, nlimGammaHeII, ntkGamma

    integer :: nIterateMC = 1                  ! current MC interation # 

    integer :: nflux = 1                       ! number of energy bins to the high energy end 
                                               ! of the incident spectrum
    integer :: nLines = 0                      ! total number of lines 

    integer :: iOrigin, jOrigin, kOrigin       ! pointers to the origin of the axes

    ! stack pointers

    integer :: KshellLimitP                    ! KshellLimitP
    integer :: FeKaColdP                       ! pointer to 6.4keV in nuArray
    integer :: secIonP
    integer :: cRecoilP                        ! pointer to compton recoil (194.Ryd) in nuArray
    integer :: xrayP                           ! pointer to carbon k-shell ionization in nuArray
    integer :: cRecoilHP, cRecoilHeP           ! recoil pointers

    integer :: BjumpP                          ! pointer to balmer jump

    integer :: nAuger(30,30,10)                ! max number of auger electrons freed
    
    integer :: noGpLoc = -1
    integer, dimension(3) :: noCellLoc=-1

    integer, dimension(nElements, nElements, 7, 3) &
         &:: elementP                          ! first dim is atomic number of element,
                                               ! second dimension is ion stage, 1 for atom
                                               ! third is shell, 1 for k shell, up to 7
                                               ! last is 1 for pointer to energy of threshold
                                               ! 2 is highest energy for shell
                                               ! 3 is opacity offset


    integer,pointer :: dustScaXsecP(:,:)       ! pointer to dust scatterring x-sec in xSecArray (species,size)
    integer,pointer :: dustAbsXsecP(:,:)       ! pointer to dust absorption x-sec in xSecArray  (species,size)



    integer, dimension(nHlevel+1) &
         & :: HlevNuP = 0      ! pointer to the nth H level in nuArray

    integer, dimension(nHeIlevel+1) &
         & :: HeIlevNuP = 0    ! pointer to the nth HeI level in nuArray


    integer, dimension(nHeIIlevel+1) &
         & :: HeIIlevNuP = 0   ! pointer to the nth HeII level in nuArray

   
    logical, pointer :: &
         & lgDataAvailable(:,:)! is the atomic data available for this species?

    character(len=50), dimension(3:nElements, 1:10) :: &
         & dataFile       ! name of the file containing the atomic data for this species

    integer, dimension(nElements, nElements)&
         & :: nShells          ! number of shells for each element and ionization stage
                                               ! first dimension is atomic weight, second is ionization stage     

    character(len=2)        :: elemLabel(30) 
    

    ! energies temperatures etc..
    real(kind=dp) :: Bjump = 0.                        ! Balmer jump


    ! slit stuff
    real(kind=dp) :: dxSlit = 0.
    real(kind=dp) :: dySlit = 0.

    ! photoionization data 
    real(kind=dp), dimension(5,30,30) :: CF=0.
    real(kind=dp), dimension(30,30)   :: cionTherm = 0.
    real(kind=dp), dimension(6, 30, 30, 7) :: ph1
    real(kind=dp), dimension(7, 30, 30) :: ph2
    real(kind=dp), dimension(nElements) :: aWeight
    real(kind=dp), dimension(450) :: ionEdge 

    real(kind=dp), pointer :: logGammaHI(:,:),logGammaHeI(:,:),logGammaHeII(:,:), tkGamma(:), nuGammaHI(:),&
         & nuGammaHeI(:),nuGammaHeII(:)

    ! dust arrays
    real(kind=dp), pointer :: TdustSublime(:)        ! sublimation T (nspecies)
    real(kind=dp), pointer :: grainVn(:)             ! potential of neutral grain
    real(kind=dp), pointer :: MsurfAtom(:)           ! mass of surf atom of a grain [amu]
    real(kind=dp), pointer :: rho(:)                 ! intrinsic density (nspecies)    
    real(kind=dp), pointer :: dustHeatingBudget(:,:) ! heating budget of grains (nAbComponents, nResLines+1)       
    real(kind=dp), pointer :: SEDnoExt(:)            ! SED no extinction
    real(kind=dp), pointer :: equivalentTau(:)       ! SED no extinction

    real(kind=dp), save,pointer::&
          & forbiddenLines(:,:,:,:)                                 ! emissivity from heavies  rec lines
    real(kind=dp), save,pointer::&
          & forbiddenLinesLarge(:,:) 

    ! linear increments
    real(kind=dp), pointer :: dl(:)

    ! calculate recombination lines?
    logical :: lgRecombination=.false.

    ! composition switches
    logical, dimension(nElements) :: lgElementOn ! is this element present?

    ! dust arrays
    real(kind=dp), pointer :: dustEmIntegral(:,:,:)   ! dust emission integral (species,size,temperature)
    
    ! allocatable arrays
    real(kind=dp), pointer :: absOpacSpecies(:,:)     ! absOpacSpecies (nspecies, freq)

    real(kind=dp), pointer :: continuum(:)            ! continuum array (freq)

    real(kind=dp), pointer :: gauntFF(:)              ! gaunt factor

    real(kind=dp), pointer :: gauntFFHeII(:)          ! gaunt ff for HeII

    real(kind=dp), pointer :: H2phot(:)               ! H 2-photon emission
 
    real(kind=dp), pointer :: HeI2phot(:)             ! HeI 2-photon emission

    real(kind=dp), pointer :: HeII2phot(:)            ! HeII 2-photon emission

    real(kind=dp), pointer :: ionDenUsed(:,:)         ! local ionDen array

    real(kind=dp), pointer :: nuArray(:)              ! energy array [Ryd]

    real(kind=dp), pointer :: xSecArray(:)            ! x Section array calculated in xSec_mod

    real(kind=dp), pointer :: xSecRecoil(:)           ! compton recoil xSec

    real(kind=dp), pointer :: comXSecH(:)             ! compton Xsec Heat array

    real(kind=dp), pointer :: comXSecC(:)             ! compton Xsec Cool array

    real(kind=dp), pointer :: widFlx(:)               ! widFlx array

    real(kind=dp), pointer :: photoRateFluo(:,:)      ! photoionisation rate (nfluo, nstages) [sec^-1]


    integer, pointer :: planeIonDistribution(:,:) ! initial distribution of ionising photons 



    ! derived types

    type grid_type             ! derived grid type
        sequence

        character(len=50) :: composition            ! chemical composition

        integer :: nx                               ! size of cartesian grid in x
        integer :: ny                               ! size of cartesian grid in y
        integer :: nz                               ! size of cartesian grid in z
        integer :: nCells                           ! number of active cells
        integer :: motherP                          ! pointer to the motherGrid

        integer, pointer :: dustAbunIndex(:)        ! dust chemistry file index (cell) 
        integer, pointer :: resLinePackets(:)       ! extraPackets to be generated cell (cell)
        integer, pointer :: abFileIndex(:,:,:)      ! abundance file index axxay
        integer, pointer :: lgConverged(:)          ! has the model converged at this grid cell?
        integer, pointer :: lgBlack(:)              ! is this to remain a black cell?

        integer, pointer :: active(:,:,:)           ! point to active cell in 1d array; 
                                                    ! returns 0 for inactive cells        

        real(kind=dp)    :: geoCorrX,geoCorrY,geoCorrZ       ! geometric correction
        real(kind=dp)    :: noHit                            ! cells not sampled by rad field
        real(kind=dp)    :: noIonBal                         ! Ion Balance not raeched 
        real(kind=dp)    :: noTeBal                          ! Te Balance not reached


        real(kind=dp), pointer :: absOpac(:,:)               ! dust absorption opacity 
        real(kind=dp), pointer :: scaOpac(:,:)               ! dust scattering opacity 
        real(kind=dp), pointer :: dustPDF(:,:)               ! dust emission PDF 
        real(kind=dp), pointer :: fEscapeResPhotons(:,:)     ! frac of res line photons escaping (cell,resline)
        real(kind=dp), pointer :: fluorescenceCube(:,:,:)    ! escaped fluorescence packets (ncell,nFluo,nVP)
        real(kind=dp), pointer :: Hden(:)                    ! density of H [1/cm^3]
        real(kind=dp), pointer :: Ne(:)                      ! electron density [1/cm^3]
        real(kind=dp), pointer :: Te(:)                      ! electron temperature [K]
        real(kind=dp), pointer :: Tdust(:,:,:)               ! dust temperature (species,size,cell) [K] 
        real(kind=dp), pointer :: ionDen(:,:,:)              ! fractional ion density (x,y,z,elem,stage)
        real(kind=dp), pointer :: Jste(:,:)                  ! MC estimator of stellar J (cell,nu) 
        real(kind=dp), pointer :: Jdif(:,:)                  ! MC estimator of diffuse J (cell,nu)
        real(kind=dp), pointer :: JPEots(:,:)                ! OTS line contribution to photoelectric emission
        real(kind=dp), pointer :: escapedPackets(:,:,:)      ! escaped packets (cell,nu, angle)
        real(kind=dp), pointer :: linePackets(:,:)           ! line packets (x,y,z,n)
        real(kind=dp), pointer :: LFluorescence(:,:)         ! local fluorescence luminosity of cell [e36 erg/sec]
        real(kind=dp), pointer :: LdiffuseLoc(:)             ! loc luminosity of diffuse source [e36 erg/sec]
        real(kind=dp), pointer :: Ndust(:)                   ! number density for dust

        real(kind=dp), pointer :: xAxis(:)                   ! x-axis
        real(kind=dp), pointer :: yAxis(:)                   ! y-axiss
        real(kind=dp), pointer :: zAxis(:)                   ! z-axis

        real(kind=dp), pointer :: NeInput(:)                 ! only used if lgNeInput=.t.


        real(kind=dp), pointer :: opacity(:,:)               ! total  opacity (cell,nu) [1/cm]
        real(kind=dp), pointer :: recPDF(:,:)                ! rec prob distribution function
        real(kind=dp), pointer :: linePDF(:,:)               ! line prob distribution function
        real(kind=dp), pointer :: totalLines(:)              ! fraction of non-ionizing line phots
    
        real(kind=dp), pointer :: elemAbun(:,:)              ! elemental abundance (specified
                                                    ! by number relative to total 
                                                    ! hydrogen density)
        real(kind=dp), pointer :: echoVol(:,:,:)  ! BEKS 2010. Contains volume of one
        ! grid cell enclosed by echo.  Only used if lgEcho=.true.
    end type grid_type

    ! photon packet type
    type photon_packet

        integer       :: nuP       ! pointer
        integer, dimension(2) & 
             & :: xP,yP,zP         ! grids position indeces 1= mother 2=sub
        integer, dimension(2) :: origin ! 1=ig, 2=icell
        integer       :: iG        ! grid index 

        real(kind=dp)          :: nu
        logical       :: lgStellar
        logical       :: lgLine
    
        type(vector)  :: position
        type(vector)  :: direction
  

    end type photon_packet

    ! plot type
    type plot_type
       real(kind=dp), pointer         :: intensity(:,:,:)
       real(kind=dp), pointer         :: xAxis(:), yAxis(:)

       logical               :: lgFilter
       logical               :: lgMap
       logical, pointer      :: lgLine(:)

       integer, pointer      :: lineNumber(:)
       integer               :: nPlots
       integer, pointer      :: nuP(:,:)
       
    end type plot_type

    type resLine_type
       character(len=2)      :: species
       character(len=30)     :: transition, multiplet! label

       integer               :: nuP  ! pointer to freq of line in nuArray
       integer               :: elem ! elemnt
       integer               :: ion  ! ionization stage
       integer               :: glow, ghigh ! lower and upper statistical weights
       integer               :: nmul ! # of components in the multiplet
       integer,pointer       :: moclow(:), mochigh(:) ! mocassin's lower and higher levels code

       real(kind=dp)                  :: m_ion ! mass of the ion
       real(kind=dp)                  :: wav  ! wavelength [A]
       real(kind=dp)                  :: freq ! frequency [Ryd]
       real(kind=dp)                  :: Elow,Ehigh ! lower and higher level energies
       real(kind=dp)                  :: Aik,fik ! trans prob and oscillator strength
    end type resLine_type

    !  model parameters

    type(vector), pointer :: starPosition(:) ! ionising source(s) position
    type(vector)       :: nullUnitVector 

    logical            :: lgCompton        ! evaluate Compton energy exchange?
    logical            :: lgAutoPackets    ! automatic increase of packets on the fly? 
    logical            :: lgTalk           ! talk on?
    logical            :: lgDfile          ! use an external density file?
    logical            :: lgDebug          ! debigging mode? (memory consuming)
    logical            :: lgDlaw           ! use the density law routine?
    logical            :: lgDust           ! is dust included in this model?
    logical            :: lgDustScattering ! is dust scattering included in this model?
    logical            :: lgDustConstant   ! is the md/Mg constant throughtout the grid?
    logical            :: lgGas            ! is there any gas in teh grid? 
    logical            :: lgHdenConstant   ! use the constant density parameters?
    logical            :: lgMdMg           ! dust to gas mass ratio used>? 
    logical            :: lgMdMh           ! dust to hydrogen mass ratio used>? 
    logical            :: lgOutput         ! output line fluxes,temp & ion struct files?
    logical            :: lg1D             ! 1-D switch
    logical            :: lgMultiDustChemistry ! do we have a inhomogeneous dust species distribution? 
    logical            :: lgMultiChemistry ! do we have a chemically inhomogeneous gas? 
    logical            :: lgNeutral        ! starting from neutral gas?
    logical            :: lgPlaneIonization! plane parallel ionization?    
    logical            :: lgSymmetricXYZ   ! symmetric in x, y, and z?

    character(len=50)  :: gridList ! grid list file name
    character(len=50), pointer  :: abundanceFIle(:) ! abundance file names
    character(len=50), pointer  :: dustSpeciesFIle(:) ! abundance file names
    character(len=50)           :: contDiffuse      ! shape of the diffuse ionising spectrum
    character(len=50), pointer  :: contShape(:)     ! continuumShape
    character(len=50), pointer  :: contShapeIn(:)   ! continuumShape
    character(len=50), pointer  :: spID(:)          ! input spectrum generator ID
    character(len=50)  :: densityFile      ! density file
    character(len=50)  :: dustFile(2)      ! dust files
    character(len=50)  :: MdMgFile         ! name of MdMg file
    character(len=50)  :: NdustFile        ! name of Ndust file    
    character(len=50)  :: Qfile            ! name of Qfile 
    character(len=30),pointer       :: grainLabel(:)    ! name of this species

    integer,pointer    :: viewPointPtheta(:), viewPointPphi(:)       ! viewing angles

    integer            :: nAngleBins=0     ! number of viewing angles for SED
    integer            :: TotAngleBinsTheta=10 ! total # of theta angle bins for SED
    integer            :: TotAngleBinsPhi=20 ! total # of phi angle bins for SED
    integer            :: nGrids           ! total number of grids to be used in the simulation
    integer            :: maxIterateMC     ! limit on number of MC iterations
    integer            :: maxPhotons       ! limit to packets to be use
    integer            :: nAbComponents=1  ! number of abundance components            
    integer            :: nDustComponents=1 ! number of dust components            
    integer, pointer   :: dustComPoint(:)  !  
    integer            :: nBins            ! number of energy bins 
    integer            :: nElementsUsed    ! actual number of elements used
    integer, pointer   :: nPhotons(:)      ! # of packets to be used in the sim
    integer            :: nPhotonsTot      ! # of packets to be used in the sim
    integer            :: nPhotonsDiffuse  ! # of packets to be used by the diffuse ionisation source
    integer            :: nPhotonsDiffuseLoc! # of packets to be used by the diffuse ionisation source
    integer            :: nPacketsFluo     ! # of packets to be used for each fluorescent transition 
    integer            :: nPacketsFluoLoc  ! # of packets to be used for each fluorescent trans in one cell
    integer, dimension(maxGrids) &
                      &:: nxIn,nyIn, nzIn  ! x, y and z dimensions of the grids
    integer            :: nStars           ! number of ionising sources
    integer            :: elementXref(nElements) ! x reference index array for elements actually used
    integer            :: emittingGrid     ! grid emiting illuminating radiation [0 for all]
    integer, pointer   :: starIndeces(:,:) ! (nstars, 3) 1=x; 2=y; 3=z, 4=gp
    integer            :: nstages          ! # of ionisation stages to be included

    real(kind=dp)               :: fluoCubeMineV, fluoCubeMaxeV ! limits of the fluorescence cube band
    real(kind=dp)               :: fillingFactor    ! filling factor epsilon
    real(kind=dp)               :: contCube(2)      ! continuum cube
    real(kind=dp)               :: convIncPercent   ! percentage by  which conv must increase
    real(kind=dp)               :: convWriteGrid    ! min conv level before starting to write the grid
    real(kind=dp)               :: nPhotIncrease    ! nPhoton increase factor
    real(kind=dp)               :: densityLaw(3)    ! density law parameters (R1,n,f,N0)
    real(kind=dp),pointer       :: grainAbun(:,:)     ! abundance of this species
    real(kind=dp),pointer       :: grainRadius(:)   ! grain radius [um]
    real(kind=dp),pointer       :: grainWeight(:)   ! grain weight normalised to 1
    real(kind=dp),pointer       :: viewPointTheta(:),viewPointPhi(:)     ! viewing angles
    real(kind=dp)               :: Hdensity         ! constant H density values (cm-^3)
    real(kind=dp)               :: H0Start          ! initial guess at X(H0) for regions I and II
    real(kind=dp)               :: Lphot            ! L of ionizing source [e36 phot/sec]
    real(kind=dp)               :: Ldiffuse         ! total luminosity of diffuse source [e36 erg/sec]
    real(kind=dp),pointer       :: LStar(:)         ! L of ionizing source [e36 erg/sec]
    real(kind=dp)               :: meanField        ! mean ionizing field to be used with plane parallel 
                                           ! geometry [erg/sec/cm^2]
    real(kind=dp)               :: meanFieldin
    real(kind=dp)               :: minConvergence   ! stop when this level of convergence has been reached  
    real(kind=dp)               :: minConvQHeat     ! min convergence level for Qheat routines to run
    real(kind=dp)               :: NeStart          ! initial guess at Ne
    real(kind=dp)               :: MdMgValue        ! dust-to-gas ratio by mass
    real(kind=dp)               :: NdustValue       ! dust number density [cm^-3]
    real(kind=dp)               :: nuMax            ! max nu limit [Ryd]
    real(kind=dp)               :: nuMin            ! min nu limit [Ryd]
    real(kind=dp)               :: nuStepSize       ! 
    real(kind=dp)               :: resLinesTransfer ! after what conv percent should the resonant line transfer start?    
    real(kind=dp)               :: Rnx,Rny,Rnz      ! edges [cm]
    real(kind=dp)               :: R_in             ! inner radius [cm]
    real(kind=dp)               :: R_out            ! outer radius [cm]
    real(kind=dp)               :: SEDfreq(2)       ! 
    real(kind=dp)               :: TeStart          ! initial guess at Te for regions I and II
    real(kind=dp), pointer      :: deltaE(:)        ! energy carried by a single photon
    real(kind=dp)               :: Tdiffuse         ! energy of diffuse ionising source [K]
    real(kind=dp), pointer      :: Tstellar(:)      ! T of ionizing source [K]    
    real(kind=dp), pointer      :: tStep(:)         ! time step for sb99 inputs [yrs]
    real(kind=dp)               :: XHILimit         ! convergence limit on X(HI)
    
    integer             :: lymanP                  ! frequency pointer to Lyman limit
    ! dust parameters
    integer            :: nSizes           ! number of grain sizes
    integer,pointer    :: nSpeciesPart(:)  ! number of grain specie
    integer            :: nSpecies         ! number of grain specie
    integer            :: nSpeciesMax      ! max number of grain species
    integer            :: nResLines=0      ! number of resonant lines to be transfered

end module common_mod
   

    

