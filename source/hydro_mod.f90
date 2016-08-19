! Copyright (C) 2005 Barbara Ercolano 
!
! Version 2.02
module elements_mod
    use common_mod
    use interpolation_mod
    implicit none

    real :: HlevEn(nHlevel+1), HeIlevEn(nHeIlevel+1), HeIIlevEn(nHeIIlevel+1)

    contains

    subroutine makeCollIonData()
      implicit none
      
      integer :: elem,ion, i

      close(12)
      open(unit = 12, file = trim(home)//'data/cion.dat', status='old', position='rewind')
      
      do elem = 1, nElements
         do ion = elem, 1, -1
            read(12,*) (CF(i,elem,ion),i=1,5)
         end do
      end do
      read(12,*)
      read(12,*)
      do elem = 3, nElements
         read(12,*) (cionTherm(elem,ion) , ion=1, elem)
      end do
      close(12)

    end subroutine makeCollIonData

    subroutine makeAugerData()
      implicit none
      
      integer :: elem,ion,shell,nelec,imax,i,j,ios

      close(12)
      open(unit = 12, file = trim(home)//'data/auger.dat', status='old', position='rewind')

      nAuger=0
      auger =0.

      do i = 1, 1696
         read(12,*) elem,ion,shell,nelec
         nauger(elem,ion,shell)=nelec
         read(12,*) (auger(elem,ion,shell,j), j = 1, 10)
      end do
      close(12)

    end subroutine makeAugerData

    ! this subroutine makes the data for Hydrogen and Helium
    subroutine makeHydro()

        ! local variables
        integer :: iup, ilow           ! counters
        
        ! define excitation and ionization temperatures for hydrogen
        do iup = 1, nHlevel
            HlevEn(iup) = HIonPot / float(iup*iup)
            ! fill in stat weight for H
        end do
        HlevEn(nHlevel+1) = HlevEn(2)


        ! define excitation and ionization temperatures for HeI
        do iup = 1, nHeIlevel
            HeIlevEn(iup) = 1. / float(iup*iup)
            ! fill in stat weight for HeI
        end do
        ! ground state and the n=2 are not hydrogenic
        HeIlevEn(1) = 1.80802
        HeIlevEn(2) = 0.2478
        HeIlevEn(nHeIlevel+1) = HeIlevEn(2)

        ! define excitation and ionization temperatures for HeII   
        do iup = 1, nHeIIlevel
            HeIIlevEn(iup) = 4. / float(iup*iup)
        end do
        HeIIlevEn(nHeIIlevel+1) = HeIIlevEn(2)

     end subroutine makeHydro
     
    ! this subroutine determines the outer shell and statistical weights
     ! LS coupling
    subroutine getOuterShell(z, nElec, outShell, g0, g1)
        implicit none
        
        integer, intent(in)  :: z        ! atomic number from 1 to 30
        integer, intent(in)  :: nElec    ! number of electrons from 1 to z
        integer, intent(out) :: outShell ! number of the outer shell
        integer, intent(out) :: g0       ! statistical weight of (z, nelec) ground state
        integer, intent(out) :: g1       ! statistical weight of (z, nelec-1) ground state
            
        ! local variables
        integer, dimension(30) :: ss, gl ! data for outer shells and statistical weights
        integer, dimension(19:30) :: glhigh
        
 
        ! assign the data to the local arrays
        ss = (/1,1,2,2,3,3,3,3,3,3,4,4, &
              &  5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7/)

        gl = (/2, 1, 2, 1, 6, 9, 4, 9, 6, 1, 2, 1, 6, 9, 4, 9, 6, 1, 2, 1, &
             & 10, 21, 28, 7, 6, 25, 28, 21, 2, 1/) 

        glhigh = (/10,21,28,25,6,25,28,21,10,1,2,1/)
        

        if ( .not.lgElementOn(z) ) then
            print*, "! getOuterShell: element is switched off [elem]", z
            stop
        end if

        if ( z > 30 ) then
            print*, "! getOuterShell: atomic number out of range"
            stop
         end if
        if ( (nElec < 1) .or. (nElec > z) ) then
            print*, "! getOuterShell: number of electrons out of range"
            stop
        end if
        
        ! number of the outer shell and statistical weight
        outShell = ss(nElec)
        if (z==nElec .and.z>18) outShell = 7
        if (z==nElec+1 .and. &
             &(z==20.or.z==21.or.z==22.or.z==25.or.z==26)) outShell = 7

        g0 = gl(nElec)
        if (nElec == 1) then
            g1 = 1
        else
            g1 = gl(nElec-1)
        end if        
        
        if ( (z>20) .and. nElec >= 19 .and. (z-nElec)>=1) then
           g0 = glhigh(nElec) 
           if (nElec == 19) then
              g1 = gl(nElec-1)
           else
              g1 = glHigh(nElec-1)
           end if
        end if
         
        if ( z == nElec ) then 
           if (z>20 .and. nElec >= 21) then
              g1 = glHigh(nElec-1)
           end if

           select case (z)
           case (21)
              g1 = 15
           case (22)
              g1 = 28
           case (25) 
              g1 = 7
           case (26)
              g1 = 30
           end select

        end if


        if ( z-nElec == 1 ) then
           select case (z)
           case (21)
              g0 = 15
           case (22)
              g0 = 28
           case (25) 
              g0 = 7
           case (26)
              g0 = 30
           end select
        end if     


     end subroutine getOuterShell

      subroutine readHeIRecLines()
       implicit none
       
       ! iden = 1 is 100cm^-2
       ! iden = 2 is 10000cm^-2
       ! iden = 3 is 1000000cm^-2
       integer :: iline, iden, j

       close(13)
       open(file=trim(home)//'data/HeIrecLines.dat', unit = 13, status='old')
       
       do iden = 1, 3
          do iline = 1, 34
             read(13,*) (HeIrecLineCoeff(iline,iden,j), j = 1, 4)
             HeIrecLineCoeff(iline,iden,1) = HeIrecLineCoeff(iline,iden,1)*1.e25 
          end do
       end do


       close(13)
     end subroutine readHeIRecLines
     
     ! this subroutine assignes continuum energy pointers 
     ! to shells for all atoms     
     subroutine setShells(nElem)
         implicit none
 
         integer, intent(in) :: nElem
 
         ! local variables
         integer            :: g0             ! statistical weight of (z, n) ground state
         integer            :: g1             ! statistical weight of (z, n-1) ground state
         integer            :: ion            ! ionic stage (=1 for atom up to totElem for hydrogenic)
         integer            :: nElec          ! number of boundn electrons
         integer            :: shell          ! shell number
         integer, parameter :: totElem = 30
         integer            :: outShell       ! outer shell number


         real               :: thres ! threshold [Ryd]

         if ( .not.lgElementOn(nElem) ) then
             print*, "! setShells: element is switched off [elem]", nElem
             stop
         end if

         ! check nElem is within the array bounds
         if ( (nElem > totElem) .or. (nElem < 1) ) then
             print*, "! setShells: nElement out of range [1, 30]"         
             stop
         end if

         do ion = 1, nElem

             ! find number of bound electrons
             nElec = nElem - ion + 1
             call getOuterShell(nElem, nElec, outShell, g0, g1)
             nShells(nELem, ion) = outShell

             ! loop on all inner shells, valence shell
             do shell = 1, outShell
                 thres = ph1(1, nElem, nElec, shell) / RydToeV

                 ! negative IP shell doesn't exist, so set 
                 ! the upper limit lower then the 
                 ! lower limit so this never loop upon. 
                 ! used as flags by limitShell to check
                 ! whether this is a real shell
!print*, nelem, ion, shell,  thres, nuArray(1)
                 if ( thres <= 0.1 )  then
                     elementP(nElem, ion, shell, 1) = 2
                     elementP(nElem, ion, shell, 2) = 1
                 else
                     ! this is the lower limit to the range
                     call locate(nuArray, thres, elementP(nElem, ion, shell, 1))
                     if ( elementP(nElem, ion, shell, 1) < nbins ) then
                         elementP(nElem, ion, shell, 1) = elementP(nElem, ion, shell, 1)
                     end if

                     ! this is the higher limit to the range.
                     ! limitShell returns the pointer to the threshold
                     ! of the next major shell. Fo the k-shell returns
                     ! kShellLimit (= nuMax)
                     elementP(nElem, ion, shell, 2) = limitShell(ion, shell, nElem)
                 end if
             end do
             ! this is the valence pointer
             call locate(nuArray, thres, elementP(nElem, ion, outShell, 1))
             if ( elementP(nElem, ion, outShell, 1) < nbins ) then
                 elementP(nElem, ion, outShell, 1) = elementP(nElem, ion, outShell, 1)
             end if

          end do


    end subroutine setShells

    ! this function returns the high energy limit to the energy range     
    function limitShell(ion, shell, nElem)
        implicit none

        integer, intent(in) :: ion      ! ionic stage (=1 for atom up to totElem for hydrogenic)
        integer, intent(in) :: nElem    ! element number
        integer, intent(in) :: shell    ! shell number

        integer             :: limitShell
        
        ! local variables

        if (.not.lgElementOn(nElem)) then
            print*, "! limitShell: element is not switched on [elem]", nElem
            stop
        end if 


        select case (shell)
        ! high energy limit to code kShellLimit = nuMax
        case (1)
            limitShell = kShellLimitP
        ! 2s shell, upper limit set to high energy limit
        case (2)
            limitShell = kShellLimitP
        ! 2p shell,  upper limit set to high energy limit 
        case (3)
            limitShell = kShellLimitP
        ! 3s shell, upper limit set to k-shell edge
        case (4)
            limitShell = elementP(nElem, ion, 1, 1) -1
        ! 3p shell,  upper limit set to k-shell edge
        case (5)
            limitShell = elementP(nElem, ion, 1, 1) -1 
        ! 3d sgell, upper limit set to k-shell edge
        case (6)
            limitShell = elementP(nElem, ion, 1, 1) -1 
        ! 4s shell, upper limit set to 3d
        case (7)
            ! if the shell 6 is empty (3d) then set it to 5 (3p)
            if (elementP(nElem, ion, 6, 1)<3) then
                limitShell = elementP(nElem, ion, 5, 1) -1
            else
                limitShell = elementP(nElem, ion, 6, 1) -1
            end if
        case default
            print*, "! limitShell: shell number out of range [1, 7]"
            stop
        end select
        
      end function limitShell


     ! this subroutine sets up the pointers for the lines and 
     ! the cotinua
     subroutine setPointers()
     
         ! local variables
         integer :: i, j          ! counters

         ! set pointer for K-shell limit 
         call locate(nuArray, KshellLimit, KshellLimitP)

         ! set pointer for secondary ionisation
         call locate(nuArray,7.353,secIonP)

         ! set pointer for Compton recoil
         call locate(nuArray,194., cRecoilHP)
         call locate(nuArray,260., cRecoilHeP)
         
         ! set xRay pointer
         call locate(nuArray,20.6, xrayP)

         ! set the pointer to the Balmer jump in the nuArray
         call locate(nuArray, 0.25, BjumpP)

         ! set data for hydrogen and helium atoms
         call makeHydro()  
         
         ! read in data from benj skil smit 99
         call readHeIRecLines()

         ! set HlevNuP (pointer to the nth H level in nuArray)           
         do i = 1, nHlevel
            call locate(nuArray, HlevEn(i), HlevNuP(i))     
         end do
         HlevNuP(nHlevel+1) = HlevNuP(2)

         ! set HeIlevNuP (pointer to the nth HeI level in nuArray)
         do i = 1, nHeIlevel
            call locate(nuArray, HeIlevEn(i), HeIlevNuP(i))
         end do
         HeIlevNuP(nHeIlevel+1) = HeIlevNuP(2)

         ! set HeIIlevNuP (pointer to the nth HeII level in nuArray)
         do i = 1, nHeIIlevel
            call locate(nuArray, HeIIlevEn(i), HeIIlevNuP(i))
            HeIIlevNuP(i) = HeIIlevNuP(i)         
         end do
         HeIIlevNuP(nHeIIlevel+1) = HeIIlevNuP(2)
 
         ! set up shell pointers
         ! the folloowing order of elements is more or less in decreasing
         ! abundance

         ! H and He
         nShells(1,1) = 1
         nShells(2,1) = 1
         nShells(2,2) = 1

         ! carbon
         if (lgElementOn(6)) call setShells(6)

         ! oxygen
         if (lgElementOn(8)) call setShells(8)
   
         ! nitrogen -after O so that O gets more accurate pointer
         if (lgElementOn(7)) call setShells(7)

         ! neon
         if (lgElementOn(10)) call setShells(10)

         ! sodium
         if (lgElementOn(11)) call setShells(11) 

         ! magnesium
         if (lgElementOn(12)) call setShells(12) 

         ! alluminium
         if (lgElementOn(13)) call setShells(13) 

         ! silicon
         if (lgElementOn(14)) call setShells(14)

         ! phosphorus
         if (lgElementOn(15)) call setShells(15)

         ! sulphur
         if (lgElementOn(16)) call setShells(16)

         ! chlorine
         if (lgElementOn(17)) call setShells(17)

         ! iron
         if (lgElementOn(26)) call setShells(26)

         ! argon
         if (lgElementOn(18)) call setShells(18)

         ! potassium
         if (lgElementOn(19)) call setShells(19)

         ! calcium
         if (lgElementOn(20)) call setShells(20) 

         ! scandium
         if (lgElementOn(21)) call setShells(21)

         ! titanium
         if (lgElementOn(22)) call setShells(22) 
 
         ! vanadium
         if (lgElementOn(23)) call setShells(23)

         ! chromium
         if (lgElementOn(24)) call setShells(24)

         ! manganese
         if (lgElementOn(25)) call setShells(25)

         ! fluorine
         if (lgElementOn(9)) call setShells(9)

         ! lithium
         if (lgElementOn(3)) call setShells(3)

         ! beryllium
         if (lgElementOn(4)) call setShells(4)

         ! boron
         if (lgElementOn(5)) call setShells(5)

         ! cobalt
         if (lgElementOn(27)) call setShells(27)

         ! nickel
         if (lgElementOn(28)) call setShells(28)

         ! copper
         if (lgElementOn(29)) call setShells(29)

         ! zinc
         if (lgElementOn(30)) call setShells(30)  

     end subroutine setPointers

     subroutine makeElements()
       implicit none

       integer           :: elem, ion  ! counters
       integer           :: i, iup, ilow ! counters
       integer           :: ios        ! I/O error status
       
       close(17)
       open(file=trim(home)//'data/fileNames.dat', action="read", unit=17, status='old', &
            & position='rewind', iostat=ios)
       
       if (ios/=0) then
          
          print*, "! makeElements: can't open ",trim(home),"data/fileNames.dat"
          stop
          
       end if


       ! read in the file names for the atomica data
       ! (even the non existing ones.. never know might add them later)
       do elem = 3, nElements
          do ion = 1, min(elem+1, 10)

             read(17, '(A20)') dataFile(elem, ion)

          end do
       end do
        
       close(17)
       

       ! find out the number of lines and what files are available
       nLines = 0

       ! HI recombination lines
       do iup = 3, 15
          do ilow = 2, min(8,iup-1)
             nLines = nLines + 1
          end do
       end do

       ! HeI singlet recombination lines
       do i = 1, 9
          nLines = nLines + 1
       end do

       ! HeI triplet recombination lines
       do i = 1, 11
          nLines = nLines + 1
       end do

       ! HeII recombination lines
       do iup = 3, 30
          do ilow = 2, min(16, iup-1)
             nLines = nLines+1
          end do
       end do

       ! initialise lgDataAvailable
       lgDataAvailable = .false.

       ! heavies collisional lines
       do elem = 3, nELements
          do ion = 1, min(elem+1, nstages)

             if(.not.lgElementOn(elem)) exit

             close(18)
             open(file=trim(home)//dataFile(elem,ion), unit=18, action="read", status='old', &
                  & position='rewind', iostat=ios)

             if (ios == 0) then
                
                close(18)
             
                if (elem == 26 .and. ion ==2) then
                   nLines = nLines + nForLevelsLarge*nForLevelsLarge! nForlLevels
                   lgDataAvailable(elem, ion) = .true.
                else
                   nLines = nLines + nForLevels*nForLevels ! nForlLevels
                   lgDataAvailable(elem, ion) = .true.
                end if
             else

                lgDataAvailable(elem, ion) = .false.

             end if

          end do
       end do

     end subroutine makeElements

end module elements_mod




