! Copyright (C) 2005 Barbara Ercolano 
!
! Version 2.02
module vector_mod


  use constants_mod

  implicit none

  public

  ! The definition of the vector type

  type vector
     real :: x
     real :: y
     real :: z
  end type vector

  ! Define multiply

  interface operator(*)
     module procedure mult
  end interface

  ! and divide

  interface operator(/)
     module procedure divideVec
  end interface

  ! add

  interface operator(+)
     module procedure add
  end interface

  ! subtract

  interface operator(-)
     module procedure subtract
  end interface

  ! equivalence

  interface operator(==)
      module procedure equivalence
  end interface 

  ! notEquivalence(/=)

  interface operator(/=)
      module procedure notEquivalence
  end interface

  ! note that the following two operators are costumised 
  ! for use in the mocassinPlot routines, these are only 
  ! meaningful in the context of this program

  ! >=

  interface operator(>=)
     module procedure greaterEqual
  end interface

  ! <= 
  
  interface operator(<=)
     module procedure lessEqual
  end interface

  ! dot product

  interface operator(.dot.)
     module procedure dotProd
  end interface

  ! cross product

  interface operator(.cross.)
     module procedure crossProd
  end interface

contains

  ! the dot product function

  real function dotProd(a , b)
    type(vector), intent(in) :: a
    type(vector), intent(in) :: b

    dotProd = a%x*b%x + a%y*b%y + a%z*b%z

  end function dotProd

  ! the cross product function

  type(vector) function crossProd(a ,b)
    type(vector), intent(in) :: a
    type(vector), intent(in) :: b

    crossProd%x =  (a%y*b%z - a%z*b%y)
    crossProd%y = -(a%x*b%z - a%z*b%x)
    crossProd%z =  (a%x*b%y - a%y*b%x)
  end function crossProd

  ! normalization subroutine - checks for zero vector

  subroutine normalize(a)
    type(vector) :: a
    real :: m

    m = modulus(a)

    if (m == 0.) then
       write(*,'(a)') "! Attempt to normalize the zero vector"
       STOP
    endif

    a%x = a%x / m
    a%y = a%y / m
    a%z = a%z / m

  end subroutine normalize

  ! find the modulus of a vector

  real function modulus(a)
    type(vector) :: a

    modulus = a%x*a%x + a%y*a%y + a%z*a%z
    modulus = sqrt(modulus)

  end function modulus

  ! multiply function

  type(vector) function mult(a,b)
    real, intent(in) :: a
    type(vector), intent(in) :: b

    mult%x = a * b%x
    mult%y = a * b%y
    mult%z = a * b%z

  end function mult


  ! divide vector by a scalar

  type(vector) function divideVec(a,b)
    type(vector), intent(in) :: a
    real, intent(in) :: b

    divideVec%x = a%x / b
    divideVec%y = a%y / b
    divideVec%z = a%z / b

  end function divideVec

  ! add two vectors
  
  type(vector) function add(a,b)
    type(vector), intent(in) :: a
    type(vector), intent(in) :: b

    add%x = a%x + b%x
    add%y = a%y + b%y
    add%z = a%z + b%z

  end function add

  ! subtract two vectors

  type(vector) function subtract(a,b)
    type(vector), intent(in) :: a
    type(vector), intent(in) :: b

    subtract%x = a%x - b%x
    subtract%y = a%y - b%y
    subtract%z = a%z - b%z

  end function subtract

  logical function equivalence(a,b)
    type(vector), intent(in) :: a
    type(vector), intent(in) :: b
    
    if( (a%x == b%x) .and.&
         & (a%y == b%y) .and.&
         & (a%z == b%z) ) then
      equivalence = .true.
    else
      equivalence = .false.
    end if
  end function equivalence

  logical function greaterEqual(a,b)
    type(vector), intent(in) :: a 
    type(vector), intent(in) :: b

    if( (a%z <= b%z) .and.&
         & ( a%x/(sqrt(1.-a%z*a%z)) >= &
         & b%x/(sqrt(1.-b%z*b%z))) ) then
      greaterEqual = .true.
    else
      greaterEqual = .false.
    end if
  end function greaterEqual

  logical function lessEqual(a,b)
    type(vector), intent(in) :: a 
    type(vector), intent(in) :: b

    if( (a%z >= b%z) .and.&
         & ( a%x/(sqrt(1.-a%z*a%z)) <= &
         & b%x/(sqrt(1.-b%z*b%z))) ) then
      lessEqual = .true.
    else
      lessEqual = .false.
    end if
  end function lessEqual


  logical function notEquivalence(a,b)
    type(vector), intent(in) :: a
    type(vector), intent(in) :: b

    if( (a%x == b%x) .and.&
         & (a%y == b%y) .and.&
         & (a%z == b%z) ) then
      notEquivalence = .false.   
    else
      notEquivalence = .true.
    end if
  end function notEquivalence


  ! get polar form of a cartesian vector
  
  subroutine getPolar(vec, r, theta, phi)

    implicit none
    type(vector) :: vec
    real :: r, theta, phi, cosTheta 

    r = modulus(vec)
    if ((vec%y == 0.) .and. (vec%x == 0)) then
       phi = 0.
    else
       phi = atan2(vec%y, vec%x)
    endif
    if (phi < 0.) phi = phi + twoPi
    cosTheta = vec%z/r
    theta = acos(cosTheta)
  end subroutine getPolar

  ! rotate a vector "a" about the z-axis by angle b

  type(vector) function rotateZ(a,b)
    type(vector), intent(in) :: a
    real, intent(in) :: b   ! angle in radians
    real :: cosb, sinb

    cosb = cos(b)
    sinb = sin(b)

    rotateZ%x = cosb * a%x + sinb * a%y
    rotateZ%y =-sinb * a%x + cosb * a%y
    rotateZ%z = a%z

  end function rotateZ

  type(vector) function rotateX(a,b)
    type(vector), intent(in) :: a
    real, intent(in) :: b   ! angle in radians
    real :: cosb, sinb

    cosb = cos(b)
    sinb = sin(b)

    rotateX%x = a%x 
    rotateX%y = cosb * a%y + sinb * a%z
    rotateX%z =-sinb * a%y + cosb * a%z

  end function rotateX

  type(vector) function rotateY(a,b)
    type(vector), intent(in) :: a
    real, intent(in) :: b   ! angle in radians
    real :: cosb, sinb

    cosb = cos(b)
    sinb = sin(b)

    rotateY%x = cosb * a%x + sinb * a%z
    rotateY%y = a%y
    rotateY%z =-sinb * a%x + cosb * a%z

  end function rotateY

    type(vector) function randomUnitVector()
    real :: r1, r2, u, v, w, t, ang
    call random_number(r1)
    w = 2.*r1 - 1.
    t = sqrt(1.-w*w)
    call random_number(r2)
    ang = 3.141592654*(2.*r2-1.)
    u = t*cos(ang)
    v = t*sin(ang)

    randomUnitVector = vector(u,v,w)
  end function randomUnitVector


end module vector_mod

