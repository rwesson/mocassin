!
! this program needs the input file makeSubGrids.in
! it will create a number of input grid files to be 
! used by the multigrid routines of MOCASSIN 2.0
!
! Barbara Ercolano December 2004
program makeSubGrids
  implicit none

  logical           :: newFile

  character(len=50) :: gridList
  character(len=50), pointer :: dFileS(:)    

  integer           :: ios, nGrids, iG, jG, motherP
  integer           :: nx,ny,nz,abIndex
  integer           :: i,j,k

  real              :: x1,xN,y1,yN,z1,zN
  real              :: hDen,x,y,z

  gridList = 'subGrids.in'

  open(unit=10, file='makeSubGrids.in',  position='rewind',status='unknown', iostat=ios)
  if (ios/=0) then
     print*, 'error opening makeSubGrids'
     stop
  end if
  open(unit=11, file=gridList, position='rewind',status='unknown', iostat=ios)
  if (ios/=0) then
     print*, 'error opening gridList'
     stop
  end if

  read(10,*) nGrids

  allocate(dFileS(1:nGrids), stat=ios)
  if (ios/=0) then
     print*, 'error allocating memory to dFileS'
     stop
  end if

  abIndex = 0

  do iG = 1, nGrids

     read(10,*) motherP, nx,ny,nz,dFileS(iG)
     read(10,*) x1,xN,y1,yN,z1,zN
     read(10,*) hDen, abIndex

     newFile = .true.

     write(11,*) motherP, nx,ny,nz,trim(dFileS(iG))
     write(11,*) x1,xN,y1,yN,z1,zN
    
     do jG = 1, iG-1
        if (dFileS(iG) == dFileS(jG) )then
           newFile = .false.
           exit
        end if
     end do

     if (newFile) then

        open(unit=12, file=dFileS(iG), position='rewind',status='unknown', iostat = ios)
        if (ios/=0) then
           print*, 'error opening ', dFileS(iG)
           stop
        end if

        do i = 1, nx
           do j = 1, ny
              do k = 1, nz
                 
                 x = (real(i)-1)/(real(nx)-1)
                 y = (real(j)-1)/(real(ny)-1)
                 z = (real(k)-1)/(real(nz)-1)
             
                 if (abIndex>0) then                 
                    write(12,*) x,y,z,hDen,abIndex                 
                 else
                    write(12,*) x,y,z,hDen                 
                 end if
              end do
           end do
        end do
        close(11)
           
     end if
     
  end do
        
end program makeSubGrids
