program makeHomCube
  implicit none
  integer :: i,j,k,nx,ny,nz
  real::Rx,Ry,Rz,den,dx,dy,dz
  real,pointer::x(:),y(:),z(:)

  nx=6
  ny=60
  nz=6

  Rx=1.e17
  Ry=1.e18
  Rz=1.e17
  den=10000.

  allocate(x(nx))
  allocate(y(ny))
  allocate(z(nz))

  open(unit=10, file='cuboid.dat', status='unknown')

  dx = Rx/(nx-1)
  dy = Ry/(ny-1)
  dz = Rz/(nz-1)
  
  x(1)=0.
  do i = 2, nx
     x(i)=x(i-1)+dx
  end do
  y(1)=0.
  do i = 2, ny
     y(i)=y(i-1)+dy
  end do
  z(1)=0.
  do i = 2, nz
     z(i)=z(i-1)+dz
  end do
  
  do i = 1, nx
     do j = 1, ny
        do k = 1, nz
           
           write(10,*) x(i), ' ', y(j), ' ', z(k), ' ', den

        end do
     end do
  end do
  
  close(10)
end program makeHomCube
