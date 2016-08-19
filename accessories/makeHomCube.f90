program makeHomCube
  implicit none
  integer :: i,j,k,nx,ny,nz
  real::R

  nx=50
  ny=50
  nz=50
  R=1.e17
  

  open(unit=10, file='cube.dat', status='unknown')

  do i = 1, nx
     do j = 1, nx
        do k = 1, nx
           
           write(10,*) real(i-1)*R/(nx-1),' ',real(j-1)*R/(ny-1), ' ', real(k-1)*R/(nz-1), 1.e-7

        end do
     end do
  end do
  
  close(10)
end program makeHomCube
