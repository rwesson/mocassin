! 
! this program calculates T0 and t2 according 
! to Peimbert (1967) formulation
! compatible with mocassin.2.02.16
! expects to read in grid0.out grid1.out grid2.out
!
program getT0t2
  implicit none 
  
  real ::  dx,dy,dz
  real, pointer :: xaxis(:),yaxis(:),zaxis(:),Te(:,:,:),Ne(:,:,:),& 
       & ionDen(:,:,:,:,:), Hden(:,:,:), dV(:,:,:),T0(:,:), t2(:,:), den(:,:)

  integer :: i,j,k,nx,ny,nz,iread, ncells,nel=11,l,li,m,nions=6

  open(unit=10,file='grid0.out')
  
  read(10,*) iread
  read(10,*) nx, ny, nz, ncells
print*, nx,ny,nz

  allocate(xaxis(nx))
  allocate(yaxis(ny))
  allocate(zaxis(nz))
  allocate(ionDen(nx,ny,nz,nel,nions))
  allocate(Hden(nx,ny,nz))
  allocate(Te(nx,ny,nz))
  allocate(Ne(nx,ny,nz))
  allocate(T0(nel,nions))
  allocate(t2(nel,nions))
  allocate(dV(nx,ny,nz))
  allocate(den(nel,nions))
  dV=0.

  do i = 1, nx
     read(10,*) xaxis(i)
  end do
  do i = 1, ny
     read(10,*) yaxis(i)
  end do
  do i = 1, nz
     read(10,*) zaxis(i)
  end do
  xaxis = xaxis/1.e15
  yaxis = yaxis/1.e15
  zaxis = zaxis/1.e15

  close(10)

  open(unit=11,file='grid1.out')
  do i= 1, nx
     do j= 1, ny
        do k= 1, nz
           read(11,*) Te(i,j,k), Ne(i,j,k), Hden(i,j,k)
           if (i<nx.and.i>1) then
              dx = (xaxis(i+1)-xaxis(i-1))/2.
           else if (i==nx) then
              dx = xaxis(i)-xaxis(i-1)
           else if (i==1) then
              dx = xaxis(i+1)-xaxis(i)
           end if
           if (j<ny.and.j>1) then
              dy = (yaxis(j+1)-yaxis(j-1))/2.
           else if (j==ny) then
              dy = yaxis(j)-yaxis(j-1)
           else if (i==1) then
              dy = yaxis(j+1)-yaxis(j)
           end if
           if (k<nz.and.k>1) then
              dz = (zaxis(k+1)-zaxis(k-1))/2.
           else if (k==nz) then
              dz = zaxis(k)-zaxis(k-1)
           else if (k==1) then
              dz = zaxis(k+1)-zaxis(k)
           end if
           dV(i,j,k) = dx*dy*dz
        end do
     end do
  end do

  close(11)

  open(unit=12,file='grid2.out') 

  do i= 1, nx
     do j= 1, ny
        do k= 1, nz
           do l = 1, nel
              if (l==1) then
                 li = 2
              else if (l==2) then
                 li = 3
              else
                 li = 6
              end if

              read(12,*) (ionDen(i,j,k,l,m), m=1, li)

           end do
        end do
     end do
  end do

  close(12)

  den=0.
  T0=0.
  t2=0.
  ! calculate T0 and t2
  do i = 1,  nx
     do j= 1, ny
        do k= 1, nz
           do l = 1, nel
              if (l==1) then
                 li = 2
              else if (l==2) then
                 li = 3
              else
                 li = 6
              end if

              do m = 1, li
                 
                 den(l,m) = den(l,m)+&
                      & Ne(i,j,k)*ionDen(i,j,k,l,m)*Hden(i,j,k)*dV(i,j,k)
                 
                 T0(l,m) = T0(l,m)+&
                      & Te(i,j,k)*Ne(i,j,k)*ionDen(i,j,k,l,m)*Hden(i,j,k)*dV(i,j,k)

              end do
           end do
        end do
     end do
  end do

  do l = 1, nel
     if (l==1) then
        nions = 2
     else if (l==2) then
        nions = 3
     else
        nions = 6
     end if

     do li = 1, nions
        T0(l,li) = T0(l,li)/den(l,li)
        do i = 1, nx
           do j = 1, ny
              do k = 1, nz
                 t2(l,li) = t2(l,li) + ((Te(i,j,k)-T0(l,li))**2.)*&
                      & Ne(i,j,k)*ionDen(i,j,k,l,li)*Hden(i,j,k)*dV(i,j,k)
              end do
           end do
        end do
        t2(l,li) = t2(l,li)/(den(l,li)*T0(l,li)**2.)
     end do
  end do

  print*, '           Element             Ion   T0              t^2'
  do l = 1, nel
     if (l==1) then
        nions = 2
     else if (l==2) then
        nions = 3
     else
        nions = 6
     end if

     do li = 1, nions
        print*, l, li, T0(l,li),t2(l,li)
     end do
  end do

end program getT0t2

  
  
 
  
