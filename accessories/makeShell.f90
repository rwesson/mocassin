program makeShell
  implicit none


  integer :: nx, i,j,k, nA

  real   :: Rin, Rout, dr, Nh, x, y, z, inr, radius, V

  open(unit=10, file='inputShell.in', status='old')
  open(unit=11, file='shell.out', status='unknown')

  read(10,*) nx
  read(10,*) Rin, Rout
  read(10,*) Nh

  dr = 1./(real(nx-1))
  inr = Rin/Rout

  V = dr*dr*dr

  nA = 0
  do i = 1, nx
     do j = 1, nx
        do k = 1, nx
           
           x = (i-1)*dr
           y = (j-1)*dr
           z = (k-1)*dr

           radius = sqrt(x*x+y*y+z*z)
  
           if (radius<=1. .and. radius >= inr) then
           
              write(11,*) x*Rout, y*Rout, z*Rout, Nh
              nA = nA + 1
      
           else

              write(11,*) x*Rout, y*Rout, z*Rout, 0.

           end if
  
        end do
     end do
  end do

  print*, 'Total number of active cells: ', nA 
  print*, 'Total volume                : ', V*nA*Rout
  

  close(10)
  close(11)

end program makeShell
  
  
  
