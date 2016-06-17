! Copyright (C) 2005 Barbara Ercolano 
!  
! Version 2.02
module common_mod
    use constants_mod
    use vector_mod
    implicit none
    save

    ! MPI variables
    integer, pointer :: lgConvergedTemp(:)  ! temporary converged? flag
    integer, pointer :: lgBlackTemp(:)      ! temporary converged? flag

    logical         :: lgWarm=.false.           ! warm started?
    logical         :: lgNeInput=.false.        ! Ne distribution entered
    logical         :: lgResLinesFirst = .true. ! first time the res lines transfer proc is called? 
    logical         :: lgPhotoelectric = .true. ! photoelectric effect on
    logical         :: lgQheat = .false.        ! temperature spiking on? 
    logical         :: lgWritePss=.false.       ! write Pss for qHeat file?
    logical         :: lgTraceHeating = .false. ! trace thermal balance?
    logical         :: lgGrainSpiking = .true.  ! temperature spiking for small grains


    real, pointer   :: TdustTemp(:,:,:)     ! temporary dust temperature array (species,size,cell)
    real, pointer   :: ionDenTemp(:,:,:)    ! temporary ion density array    
    real, pointer   :: NeTemp(:)            ! temporary electron density array
    real, pointer   :: TeTemp(:)            ! temporary electron Temperature array    
   

    integer         :: radio4p9GHzP   ! 
      
    integer         :: ierr           ! MPI error status
    integer         :: numtasks       ! total # of processes
    integer         :: taskid         ! process identification #
    real            :: starttime      ! start time [sec]

    real            :: endtime        ! end time [sec]    

    real            :: dTheta                  ! 
    real            :: dPhi                    ! 
    real            :: totalDustMass 
    real            :: convPercent=0.          ! total convergence percentage
    real            :: totPercent=0.           ! 
    real            :: minaQHeat = 0.          ! min radius for Q Heat in [um]
    real            :: NeUsed, TeUsed, HdenUsed ! local properties of the gas
    real            :: Tmax                    ! max temp for qheat

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

    integer :: xrayP                           ! pointer to carbon k-shell ionization in nuArray

    integer :: BjumpP                          ! pointer to balmer jump

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
         & :: HlevNuP = 0.     ! pointer to the nth H level in nuArray

    integer, dimension(nHeIlevel+1) &
         & :: HeIlevNuP = 0.   ! pointer to the nth HeI level in nuArray


    integer, dimension(nHeIIlevel+1) &
         & :: HeIIlevNuP = 0.  ! pointer to the nth HeII level in nuArray

   
    logical, dimension(3:nElements, 1:nstages) :: &
         & lgDataAvailable! is the atomic data available for this species?

    character(len=20), dimension(3:nElements, 1:10) :: &
         & dataFile       ! name of the file containing the atomic data for this species

    integer, dimension(nElements, nElements)&
         & :: nShells          ! number of shells for each element and ionization stage
                                               ! first dimension is atomic weight, second is ionization stage     

    character(len=2)        :: elemLabel(30) 
    

    ! energies temperatures etc..
    real :: Bjump = 0.                        ! Balmer jump


    ! slit stuff
    real :: dxSlit = 0.
    real :: dySlit = 0.

    ! photoionization data 
    real, dimension(6, 30, 30, 7) :: ph1
    real, dimension(7, 30, 30) :: ph2
    real, dimension(nElements) :: aWeight
    real, dimension(450) :: ionEdge 

    real, pointer :: logGammaHI(:,:),logGammaHeI(:,:),logGammaHeII(:,:), tkGamma(:), nuGammaHI(:),&
         & nuGammaHeI(:),nuGammaHeII(:)

    ! dust arrays
    real, pointer :: TdustSublime(:)        ! sublimation T (nspecies)
    real, pointer :: grainVn(:)             ! potential of neutral grain
    real, pointer :: MsurfAtom(:)           ! mass of surf atom of a grain [amu]
    real, pointer :: rho(:)                 ! intrinsic density (nspecies)    
    real, pointer :: dustHeatingBudget(:,:) ! heating budget of grains (nAbComponents, nResLines+1)       

    ! linear increments
    real, pointer :: dl(:)

    ! calculate recombination lines?
    logical :: lgRecombination=.false.

    ! composition switches
    logical, dimension(nElements) :: lgElementOn ! is this element present?

    ! dust arrays
    real, pointer :: dustEmIntegral(:,:,:)   ! dust emission integral (species,size,temperature)
    
    ! allocatable arrays
    real, pointer :: absOpacSpecies(:,:)     ! absOpacSpecies (nspecies, freq)

    real, pointer :: continuum(:)            ! continuum array (freq)

    real, pointer :: gauntFF(:)              ! gaunt factor

    real, pointer :: gauntFFHeII(:)          ! gaunt ff for HeII

    real, pointer :: H2phot(:)               ! H 2-photon emission
 
    real, pointer :: HeI2phot(:)             ! HeI 2-photon emission

    real, pointer :: HeII2phot(:)            ! HeII 2-photon emission

    real, pointer :: ionDenUsed(:,:)         ! local ionDen array

    real, pointer :: nuArray(:)              ! energy array [Ryd]

    real, pointer :: xSecArray(:)            ! x Section array calculated in xSec_mod

    real, pointer :: widFlx(:)               ! widFlx array

    integer, pointer :: planeIonDistribution(:,:) ! initial distribution of ionising photons 



    ! derived types

    type grid_type             ! derived grid type
        sequence

        character(len=20) :: composition            ! chemical composition

        integer :: nx                               ! size of cartesian grid in x
        integer :: ny                               ! size of cartesian grid in y
        integer :: nz                               ! size of cartesian grid in z
        integer :: nCells                           ! number of active cells
        integer :: motherP                          ! pointer to the motherGrid

        integer, pointer :: resLinePackets(:)       ! extraPackets to be generated cell (cell)
        integer, pointer :: abFileIndex(:,:,:)      ! abundance file index axxay
        integer, pointer :: lgConverged(:)          ! has the model converged at this grid cell?
        integer, pointer :: lgBlack(:)              ! is this to remain a black cell?
        

        real    :: geoCorrX,geoCorrY,geoCorrZ       ! geometric correction
        real    :: noHit                            ! cells not sampled by rad field
        real    :: noIonBal                         ! Ion Balance not raeched 
        real    :: noTeBal                          ! Te Balance not reached


        integer, pointer :: active(:,:,:)           ! point to active cell in 1d array; 
                                                    ! returns 0 for inactive cells

        real, pointer :: absOpac(:,:)               ! dust absorption opacity 
        real, pointer :: scaOpac(:,:)               ! dust scattering opacity 
        real, pointer :: dustPDF(:,:)               ! dust emission PDF 
        real, pointer :: fEscapeResPhotons(:,:)     ! frac of res line photons escaping (cell,resline)
        real, pointer :: Hden(:)                    ! density of H [1/cm^3]
        real, pointer :: Ne(:)                      ! electron density [1/cm^3]
        real, pointer :: Te(:)                      ! electron temperature [K]
        real, pointer :: Tdust(:,:,:)               ! dust temperature (species,size,cell) [K] 
        real, pointer :: ionDen(:,:,:)              ! fractional ion density (x,y,z,elem,stage)
        real, pointer :: Jste(:,:)                  ! MC estimator of stellar J (cell,nu) 
        real, pointer :: Jdif(:,:)                  ! MC estimator of diffuse J (cell,nu)
        real, pointer :: JPEots(:,:)                ! OTS line contribution to photoelectric emission
        real, pointer :: escapedPackets(:,:,:)      ! escaped packets (cell,nu, angle)
        real, pointer :: linePackets(:,:)           ! line packets (x,y,z,n)
        real, pointer :: Ndust(:)                   ! number density for dust

        real, pointer :: xAxis(:)                   ! x-axis
        real, pointer :: yAxis(:)                   ! y-axiss
        real, pointer :: zAxis(:)                   ! z-axis

        real, pointer :: NeInput(:)                 ! only used if lgNeInput=.t.


        real, pointer :: opacity(:,:)               ! total  opacity (cell,nu) [1/cm]
        real, pointer :: recPDF(:,:)                ! rec prob distribution function
        real, pointer :: linePDF(:,:)               ! line prob distribution function
        real, pointer :: totalLines(:)              ! fraction of non-ionizing line phots
    
        real, pointer :: elemAbun(:,:)              ! elemental abundance (specified
                                                    ! by number relative to total 
                                                    ! hydrogen density)
    end type grid_type

    ! photon packet type
    type photon_packet

        integer       :: nuP       ! pointer
        integer, dimension(maxGrids) & 
             & :: xP,yP,zP         ! grids position indeces
        integer, dimension(2) :: origin ! 1=ig, 2=icell
        integer       :: iG        ! grid index 

        real          :: nu

        logical       :: lgStellar
        logical       :: lgLine
    
        type(vector)  :: position
        type(vector)  :: direction
  
    end type photon_packet

    ! plot type
    type plot_type
       real, pointer         :: intensity(:,:,:)
       real, pointer         :: xAxis(:), yAxis(:)

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

       real                  :: m_ion ! mass of the ion
       real                  :: wav  ! wavelength [A]
       real                  :: freq ! frequency [Ryd]
       real                  :: Elow,Ehigh ! lower and higher level energies
       real                  :: Aik,fik ! trans prob and oscillator strength
    end type resLine_type

    !  model parameters

    type(vector), pointer :: starPosition(:) ! ionising source(s) position

    logical            :: lg3DextinctionMap! is 3D extinction to be accounted for?
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
    logical            :: lgMultiChemistry ! do we have a chemically inhomogeneous gas? 
    logical            :: lgNeutral        ! starting from neutral gas?
    logical            :: lgPlaneIonization! plane parallel ionization?    
    logical            :: lgSymmetricXYZ   ! symmetric in x, y, and z?

    character(len=50)  :: gridList ! grid list file name
    character(len=50), pointer  :: abundanceFIle(:) ! abundance file names
    character(len=50), pointer  :: contShape(:)     ! continuumShape
    character(len=50), pointer  :: contShapeIn(:)   ! continuumShape
    character(len=50), pointer  :: spID(:)          ! input spectrum generator ID
    character(len=50)  :: densityFile      ! density file
    character(len=50)  :: dustFile(2)      ! dust files
    character(len=50)  :: extMapFile       ! name of 3D extinction map file
    character(len=50)  :: MdMgFile         ! name of MdMg file
    character(len=50)  :: NdustFile        ! name of Ndust file    
    character(len=50)  :: Qfile            ! name of Qfile 
    character(len=30),pointer       :: grainLabel(:)    ! name of this species

    integer,pointer    :: viewPointPtheta(:), viewPointPphi(:)       ! viewing angles
    integer            :: nAngleBins=0     ! number of viewing angles for SED
    integer            :: TotAngleBinsTheta=180 ! total # of theta angle bins for SED
    integer            :: TotAngleBinsPhi=360 ! total # of phi angle bins for SED
    integer            :: nGrids           ! total number of grids to be used in the simulation
    integer            :: maxIterateMC     ! limit on number of MC iterations
    integer            :: maxPhotons       ! limit to packets to be used
    integer            :: nAbComponents=1  ! number of abundance components            
    integer            :: nBins            ! number of energy bins 
    integer            :: nElementsUsed    ! actual number of elements used
    integer, pointer   :: nPhotons(:)      ! # of packets to be used in the sim
    integer            :: nPhotonsTot      ! # of packets to be used in the sim
    integer, dimension(maxGrids) &
                      &:: nxIn,nyIn, nzIn  ! x, y and z dimensions of the grids
    integer            :: nStars           ! number of ionising sources
    integer            :: elementXref(nElements) ! x reference index array for elements actually used
    integer, pointer   :: starIndeces(:,:) ! (nstars, 3) 1=x; 2=y; 3=z

    real               :: fillingFactor    ! filling factor epsilon
    real               :: contCube(2)      ! continuum cube
    real               :: convIncPercent   ! percentage by  which conv must increase
    real               :: convWriteGrid    ! min conv level before starting to write the grid
    real               :: nPhotIncrease    ! nPhoton increase factor
    real               :: densityLaw(3)    ! density law parameters (R1,n,f,N0)
    real,pointer       :: grainAbun(:)     ! abundance of this species
    real,pointer       :: grainRadius(:)   ! grain radius [um]
    real,pointer       :: grainWeight(:)   ! grain weight normalised to 1
    real,pointer       :: viewPointTheta(:),viewPointPhi(:)     ! viewing angles
    real               :: Hdensity         ! constant H density values (cm-^3)
    real               :: H0Start          ! initial guess at X(H0) for regions I and II
    real               :: Lphot            ! L of ionizing source [e36 phot/sec]
    real,pointer       :: LStar(:)         ! L of ionizing source [e36 erg/sec]
    real               :: meanField        ! mean ionizing field to be used with plane parallel 
                                           ! geometry [erg/sec/cm^2]
    real               :: minConvergence   ! stop when this level of convergence has been reached  
    real               :: minConvQHeat     ! min convergence level for Qheat routines to run
    real               :: NeStart          ! initial guess at Ne
    real               :: MdMgValue        ! dust-to-gas ratio by mass
    real               :: NdustValue       ! dust number density [cm^-3]
    real               :: nuMax            ! max nu limit [Ryd]
    real               :: nuMin            ! min nu limit [Ryd]
    real               :: nuStepSize       ! 
    real               :: resLinesTransfer ! after what conv percent should the resonant line transfer start?    
    real               :: Rnx,Rny,Rnz      ! edges [cm]
    real               :: R_in             ! inner radius [cm]
    real               :: R_out            ! outer radius [cm]
    real               :: SEDfreq(2)       ! 
    real               :: TeStart          ! initial guess at Te for regions I and II
    real, pointer      :: deltaE(:)        ! energy carried by a single photon
    real, pointer      :: Tstellar(:)      ! T of ionizing source [K]
    real, pointer      :: tStep(:)         ! time step for sb99 inputs [yrs]
    real               :: XHILimit         ! convergence limit on X(HI)
    

    ! dust parameters
    integer            :: nSizes           ! number of grain sizes
    integer            :: nSpecies         ! number of grain species
    integer            :: nResLines=0      ! number of resonant lines to be transfered

end module common_mod
   

    

