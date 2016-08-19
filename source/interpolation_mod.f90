! Copyright (C) 2005 Barbara Ercolano 
!
! Version 2.02
module interpolation_mod
    contains

    subroutine sortUp(arr)
      implicit none
      
      real, dimension(:), intent(inout) :: arr
      real, pointer :: tmp(:)

      real       :: min
      

      integer    :: i, j
      integer    :: n ! size 

      n = size(arr)

      allocate(tmp(n))
      
      min = 1.e30
      do i = 1, n
         if (arr(i) < min) min = arr(i)
      end do
      tmp(1) = min
      do j = 2, n
         min = 1.e30
         do i = 1, n
            if (arr(i) < min .and. arr(i) > tmp(j-1)) min = arr(i)
         end do
         tmp(j) = min
      end do
      
      arr = tmp
      
      deallocate(tmp)

    end subroutine sortUp            
      
 
    ! subroutine locate uses the method of bisection   
    ! given an array xa of length n, and given a value x,
    ! it returns a value ns, such that x is located between 
    ! xa(ns) and xa(ns+1) 
    ! ns = 0 or ns=n is returned to indicate that x is out
    ! of range. 
    subroutine locate(xa,  x, ns)
        implicit none

        integer, intent(out) :: ns
        
        real, dimension(:), intent(in) :: xa  ! input array
        real, intent(in) :: x                 ! variable to be located
  
        ! local variables

        integer :: i                          ! counter
        integer :: jlow, jup   ! lower and upper limits
        integer :: jm                       ! mid point
        integer :: imax = 10000   ! max # of operations
        integer :: n                 ! size of array xa
 
        n = size(xa)


        ! first check if x is out of range
        if ( x > xa(n) ) then
            ns = n
            return
        end if

        if ( x < xa(1) ) then
            ns = 0
            return
        end if        


        ! initialize lower and upper limits
        jlow = 0                 
        jup = n+1

        do i = 1, imax
            ! if we are not done yet
            if ( (jup-jlow) < 1) exit 
            ! compute another mid point and either
            jm = (jup+jlow)/2             

            ! the following if statement caters for the case when x
            ! falls exactly on one of the array members
            if ( (xa(jm) == x) ) then
                ns = jm
                return
            end if

            if ( (xa(n) > xa(1)) .eqv. (x > xa(jm)) ) then
                jlow = jm ! reaplace the lower limit, or
            else
                jup = jm  ! replace the upper limit
            end if
        end do

        ns = jlow ! set the output
        
    end subroutine locate


    ! this routine will map y(x) onto x_new and return y_new(x_new)
    ! mapping is carried out by means of linear interpolation
    subroutine linearMap(y, x, nx, y_new, x_new, nx_new)
      implicit none

      real, intent(out)   :: y_new(*)
      real, intent(in)    :: y(*), x(*), x_new(*)
      
      integer, intent(in) :: nx, nx_new
      integer             :: i, ii 

      do i = 1, nx_new
         call locate(x(1:nx), x_new(i), ii)
         if (ii==0) then
            y_new(i) = y(1)
         else if (ii==nx) then
            y_new(i) = y(nx)
         else 
            y_new(i) = y(ii) +(y(ii+1)-y(ii))*(x_new(i)-x(ii))/(x(ii+1)-x(ii))
         end if
      end do

    end subroutine linearMap      

    ! subroutine polint carries out polynomial interpoation 
    ! or extrapolation. given arrays xa nd ya, each of length 
    ! n, and given a value x, it returns a value y and an 
    ! error estimate dy.
    subroutine polint(xa, ya,  x, y, dy)
        implicit none

        real, dimension(:), intent(in) :: xa, ya
        real, intent(in) :: x
        real, intent(out) :: y, dy

        ! local variables

        integer :: i, m                       ! counters
        integer :: n                          ! size of the arrays
        integer :: ns = 1            

        real :: den, ho, hp, w
        real, dimension(size(xa)) :: c, d 

        ! locate the index ns closest to the table entry
        call locate(xa, x, ns) 

        c = ya
        d = ya
        y = ya(ns)
        ns = ns-1

        n = size(xa)

        ! for each entry in the table, loop over the current
        ! c's and d's and update them
        do m=1, n-1
            do i = 1, n-m
                ho = xa(i) - x
                hp = xa(i+m) - x
                w = c(i+1) - d(i)
                den = ho - hp

                if (den == 0) then
                    print*, "! polint: two input xa's &
&                             identical (to within roundoff)"
                    stop
                end if

                den = w/den
                ! update c and d
                d(i) = hp*den
                c(i) = ho*den
            end do
    
            ! after each column in the table, decide the 
            ! correction c or d to be added to the 
            ! accumulating value of y. 
            if ( (2*ns) < (n-m) ) then
                dy = c(ns+1)
            else
                dy = d(ns)
                ns = ns-1
            end if

            y = y+dy
        end do

    end subroutine polint        

    ! subroutine spline given arrays wa and ya each of length
    ! n and given yp1 and ypn for the first derivatives of
    ! the interpolating function at points 1 and n
    ! respectively, returns the array y2a of length n
    ! containing the second derivatives of the interpolating
    ! funtion at the tabulated points. to use the condition
    ! of natural spline (2nd derivatives = 0 at the end
    ! points)  set yp1 and ypn to a value > or = to 1.e30.
    subroutine spline(xa, ya, yp1,ypn, y2a)
    implicit none

        real, dimension(:), intent(inout) :: xa, ya
        real, intent(in) :: yp1, ypn
        real, dimension(:), intent(out) :: y2a
                
        ! local variables
            
        integer :: i, k               ! counters
        integer :: n                  ! size of the arrays
        integer, parameter :: nmax=500! safety limit

        real                  :: p, qn, sig, un
        real, dimension(nmax) :: u
        n = size(xa)

        if (yp1 > .99e30) then
            y2a(1) = 0.
            u(1) = 0.
        else
            y2a(1) = -.5
            if ((xa(2)-xa(1))==0.) then
                print*, "! spline: bad xa input or x out of range [2,1] "
                stop
            end if
            u(1) = (3./(xa(2)-xa(1)))* &
                 & ((ya(2)-ya(1))/(xa(2)-xa(1))-yp1)
        end if
        do i = 2, n-1
            if ((xa(i+1)-xa(i-1))==0.) then
                print*, "! spline: bad xa input or x out of range [i+1, i-1] ", (i+1), (i-1)
                stop   
            end if               
            sig = (xa(i)-xa(i-1))/(xa(i+1)-xa(i-1))
            p=sig*y2a(i-1)+2.
            y2a(i)=(sig-1.)/p
            if ((xa(i+1)-xa(i))==0.) then    
                print*, "! spline: bad xa input or x out of range [i+1, i] ", (i+1), (i)    
                stop
            end if
            if ((xa(i)-xa(i-1))==0.) then    
                print*, "! spline: bad xa input or x out of range [i, i-1] ", (i), (i-1)    
                stop
            end if

            
            u(i)= ya(i+1)
            u(i)=u(i) - ya(i) 
            u(i) = u(i) / (xa(i+1)-xa(i)) 
            u(i)= u(i) - (ya(i)-ya(i-1))/(xa(i)-xa(i-1)) 
            u(i)=(6.*u(i)/(xa(i+1)-xa(i-1)) - sig*u(i-1))/p
       end do
        if (ypn > .99e30) then
            qn=0.
            un=0.
        else
            qn=0.5
            if ((xa(n)-xa(n-1))==0.) then
                print*, "! spline: bad xa input or x out of range [n, n-1] ", n, (n-1)              
                stop
            end if
            un=(3./(xa(n)-xa(n-1))) * &
                 & (ypn-(ya(n)-ya(n-1))/(xa(n)-xa(n-1)))
        end if
        y2a(n)=(un-qn*u(n-1))/(qn*y2a(n-1)+1.)
        
        do k = (n-1), 1, -1
            y2a(k) = y2a(k)*y2a(k+1)+u(k)
        end do

    end subroutine spline

    ! subroutine splint, given the arrays xa and ya each of
    ! length n and given the array y2a, which is the output
    ! from the subroutine spline, it returns a cubic spline
    ! interpolated value y for the given x 
    ! NOTE: if x is out of range splint will return the limit
    ! values
    subroutine splint(xa, ya, y2a, x, y)

       real, dimension(:), intent(in) :: xa, ya, y2a
       real, intent(in) :: x
       real, intent(out) :: y

       ! local variables

       integer :: klo, khi, i, n
       integer, parameter :: imax = 1000

       real :: a, b, h

       n = size(xa)

       ! check if x is within range
       if ( x <= xa(1) ) then
           y = ya(1)
           return
       else if ( x >= xa(n) ) then
           y = ya(n)
           return
       end if
      
       klo = 1
       khi = n
   
       do i = 1, imax
           if ( (khi-klo) <= 1) exit
           k = (khi+klo)/2  
           if (xa(k) > x) then
               khi=k
           else
               klo=k
           end if
       end do

       h=xa(khi)-xa(klo)
       
       if (h==0.) then
           print*, "! splint: bad xa input or x out of range [khi, klow] ", khi, klo
           stop
       end if

       a = (xa(khi) - x)/h        
       b = (x - xa(klo))/h
       y = a*ya(klo)+b*ya(khi) + &
            & ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) * &
            & (h**2)/6.

    end subroutine splint

    ! this function performs a linear interpolation of a scalar quantity in a 
    ! frequency dependent 3D grid. the interpolation coefficients t1, t2, t3 are
    ! are left to be calculated in the calling program, so that only a small 
    ! number or arguments need to be passed to the procedure
    function interpGrid(s3Dgrid, s4DGrid, s5Dgrid, xP, yP, zP, freqP, elem, ion, t1, t2, t3)
        implicit none

        integer, intent(in)                  :: xP, yP, zP      ! x, y and z axes  indeces
        integer, intent(in), optional        :: freqP           ! frequency index
        integer, intent(in), optional        :: elem            ! element index
        integer, intent(in), optional        :: ion             ! ion index             

        real, intent(in), dimension(:,:,:),&
&                          optional           :: s3Dgrid          ! 3D scalar grid
        real, intent(in), dimension(:,:,:,:),&  
&                          optional           :: s4Dgrid          ! 4D scalar grid
        real, intent(in), dimension(:,:,:,:,:),&  
&                          optional           :: s5Dgrid          ! 5D scalar grid
        real, intent(in)                     :: t1, t2, t3      ! interpolation coefficients
        real                                 :: interpGrid      ! interpolated value

        if ( present(freqP) .and. (.not.present(elem)) .and. (.not.present(ion)) ) then 
            if(.not.present(s4Dgrid)) then
                print*, "! interpGrid: insanity occurred - arguments incompatible"
                stop
            end if
            interpGrid = &  
&                   ((1.-t1)  * (1.-t2) * (1.-t3))* s4Dgrid(xP  , yP   , zP   , freqP) + &
&                   ((t1   )  * (1.-t2) * (1.-t3))* s4Dgrid(xP+1, yP   , zP   , freqP) + &
&                   ((t1   )  * (t2   ) * (1.-t3))* s4Dgrid(xP+1, yP+1 , zP   , freqP) + &
&                   ((1.-t1)  * (t2   ) * (t3   ))* s4Dgrid(xP  , yP+1 , zP+1 , freqP) + &
&                   ((1.-t1)  * (t2   ) * (1.-t3))* s4Dgrid(xP  , yP+1 , zP   , freqP) + &   
&                   ((t1   )  * (1.-t2) * (t3   ))* s4Dgrid(xP+1, yP   , zP+1 , freqP) + &   
&                   ((1.-t1)  * (1.-t2) * (t3   ))* s4Dgrid(xP  , yP   , zP+1 , freqP) + &   
&                   ((t1   )  * (t2   ) * (t3   ))* s4Dgrid(xP+1, yP+1 , zP+1 , freqP)
        else if ((.not.present(freqP)) .and. present(elem) .and. present(ion) ) then
           if(.not.present(s5Dgrid)) then
                print*, "! interpGrid: insanity occurred - arguments incompatible"
                stop
            end if
            interpGrid = &
&                   ((1.-t1)  * (1.-t2) * (1.-t3))* s5Dgrid(xP  , yP   , zP   , elem, ion) + &
&                   ((t1   )  * (1.-t2) * (1.-t3))* s5Dgrid(xP+1, yP   , zP   , elem, ion) + &
&                   ((t1   )  * (t2   ) * (1.-t3))* s5Dgrid(xP+1, yP+1 , zP   , elem, ion) + &
&                   ((1.-t1)  * (t2   ) * (t3   ))* s5Dgrid(xP  , yP+1 , zP+1 , elem, ion) + &
&                   ((1.-t1)  * (t2   ) * (1.-t3))* s5Dgrid(xP  , yP+1 , zP   , elem, ion) + &
&                   ((t1   )  * (1.-t2) * (t3   ))* s5Dgrid(xP+1, yP   , zP+1 , elem, ion) + &
&                   ((1.-t1)  * (1.-t2) * (t3   ))* s5Dgrid(xP  , yP   , zP+1 , elem, ion) + &
&                   ((t1   )  * (t2   ) * (t3   ))* s5Dgrid(xP+1, yP+1 , zP+1 , elem, ion)
        else if ((.not.present(freqP)) .and. (.not.present(elem)) .and. (.not.present(ion)) ) then
           if(.not.present(s3Dgrid)) then
                print*, "! interpGrid: insanity occurred - arguments incompatible"
                stop
            end if
            interpGrid = &
&                   ((1.-t1)  * (1.-t2) * (1.-t3))* s3Dgrid(xP  , yP   , zP  ) + &
&                   ((t1   )  * (1.-t2) * (1.-t3))* s3Dgrid(xP+1, yP   , zP  ) + &
&                   ((t1   )  * (t2   ) * (1.-t3))* s3Dgrid(xP+1, yP+1 , zP  ) + &
&                   ((1.-t1)  * (t2   ) * (t3   ))* s3Dgrid(xP  , yP+1 , zP+1) + &
&                   ((1.-t1)  * (t2   ) * (1.-t3))* s3Dgrid(xP  , yP+1 , zP  ) + &
&                   ((t1   )  * (1.-t2) * (t3   ))* s3Dgrid(xP+1, yP   , zP+1) + &
&                   ((1.-t1)  * (1.-t2) * (t3   ))* s3Dgrid(xP  , yP   , zP+1) + &
&                   ((t1   )  * (t2   ) * (t3   ))* s3Dgrid(xP+1, yP+1 , zP+1)
         end if

    end  function interpGrid

end module interpolation_mod

