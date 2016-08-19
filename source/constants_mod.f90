! Copyright (C) 2005 Barbara Ercolano 
! changed nStages from 6 to 9 for the modelling of higher energy ions (nw2006)
!
! Version 2.02
module constants_mod

    public
 
    ! physical constant
    real, parameter :: A0 = 6.30e-18               ! H0 phot x sec at the edge [cm^2]
    real, parameter :: amu = 1.66053e-24           ! [g]           
    real, parameter :: c = 2.9979250e10            ! speed of light [cm/s]
    real, parameter :: fr1Ryd = 3.28984e15         ! frequency a 1 Ryd [Hz]
    real, parameter :: hPlanck = 6.6262e-27        ! Planck constant [erg s]
    real, parameter :: HIonPot = 0.99946           ! ionization potential of H [Ryd]
    real, parameter :: hnu3c2 = 0.524166
    real, parameter :: me = 9.109382e-28           ! mass of e- [g]
    real, parameter :: kBoltzmann = 1.380658e-16     ! Boltzmann constant
    real, parameter :: radio4p9GHz = 1.489434e-6     ! 4.9 GHz in Ryd
    real, parameter :: Ryd2erg = 2.1799e-11        ! converts ryd to erg
    real, parameter :: RydInf = 109737.315         ! Ryd for infinitme mass nuclei [1/cm]
    real, parameter :: sigma = 5.66956e-5          ! Stefan-Boltzmann constant [erg*cm^-2deg^-4s-1]
    real, parameter :: Te1Ryd = 1.578866e5         ! Te at 1 Ryd [K]
    real, parameter :: Y1 = 0.2                    ! see Baldwin et al. 1991, ApJ 374, 580
    real, parameter :: FeKaCold = 470.388          ![ryd] =6.4keV -> energy of Fe K alpha for cold iron

    ! numerical constants
    real, parameter :: rootThree = 1.732050808
    real, parameter :: kShellLimit = 7.35e4        ! high energy limit to code 
    
    ! pi 
    real, parameter :: pi = 3.141592654
    real, parameter :: twoPi = 2.*pi
    real, parameter :: fourPi = 4.*pi
    real, parameter :: piByFour = pi/4.
    real, parameter :: oneOnFourPi = 1./(4.*pi)
    real, parameter :: oneOnTwoPi = 1./(2.*pi)

    ! conversion factors
    real, parameter :: degToRad = pi/180.
    real, parameter :: RydToeV = 13.6056981
    real, parameter :: radToDeg = 180./pi

    ! hard limits
    integer, parameter :: maxGrids = 50000           ! limit to the number of grids to be used 
    integer, parameter :: maxTau = 10000000         ! limit to the optical depth arrays
    integer, parameter :: nElements = 30            ! number of elements
    integer, parameter :: nForLevels = 17           ! number of levels
    integer, parameter :: nForLevelsLarge = 142     ! number of levels for FeII 
    integer, parameter :: nHlevel =  10             ! number of levels in the H atom 
    integer, parameter :: nHeIlevel = 9             ! number of levels in the HeI singlets 
    integer, parameter :: nHeIIlevel = 9            ! number of levels in the HeII singlets    
    integer, parameter :: NnuMax = 3000             ! limit to the possible number of energy bins (nbins) 
    integer, parameter :: xSecMax = 1000000         ! max limit to xSecArray 
    integer, parameter :: nTemps=3000
    integer, parameter :: recursionLimit=5000

    real, parameter :: xMax = 0.99999               ! max limit for a relative ionic abundance

end module constants_mod



