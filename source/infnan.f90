module Inf_NaN_Detection

!!     Inf_NaN_Detection module 
!!     Copyright(c) 2003, Lahey Computer Systems, Inc.
!!     Copies of this source code, or standalone compiled files 
!!     derived from this source may not be sold without permission
!!     from Lahey Computers Systems. All or part of this module may be 
!!     freely incorporated into executable programs which are offered
!!     for sale. Otherwise, distribution of all or part of this file is
!!     permitted, provided this copyright notice and header are included.

!!     This module exposes four elemental functions:
!!
!!     isnan(x)    - test for a "not a number" value
!!
!!     isinf(x)    - test for either a positive or negative "infinite" value
!!
!!     isposinf(x) - test for a positive "infinite" value
!!
!!     isneginf(x) - test for a negative "infinite" value
!!
!!     Each function accepts a single or double precision real argument, and
!!     returns a true or false value to indicate the presence of the value 
!!     being tested for. If the argument is array valued, the function returns
!!     a conformable logical array, suitable for use with the ANY function, or
!!     as a logical mask.
!!
!!     Each function operates by transferring the bit pattern from a real 
!!     variable to an integer container. Unless testing for + or - infinity,
!!     the sign bit is cleared to zero. The value is exclusive ORed with
!!     the value being tested for. The integer result of the IEOR function is
!!     converted to a logical result by comparing it to zero.
!!

    implicit none

    private

    public :: isnan, isinf, isposinf, isneginf

    ! Kind numbers for single and double precision integer containers
    integer, parameter :: Single = selected_int_kind(precision(1.e0))
    integer, parameter :: Double = selected_int_kind(precision(1.d0))

    ! Single precision IEEE values
    integer(Single), parameter :: sNaN    = Z"7FC00000"
    integer(Single), parameter :: sPosInf = Z"7F800000"
    integer(Single), parameter :: sNegInf = Z"FF800000"

    ! Double precision IEEE values
    integer(Double), parameter :: dNaN    = Z"7FF8000000000000"
    integer(Double), parameter :: dPosInf = Z"7FF0000000000000"
    integer(Double), parameter :: dNegInf = Z"FFF0000000000000"

    ! Locatation of single and double precision sign bit (Intel)
    ! Subtract one because bit numbering starts at zero
    integer, parameter :: SPSB = bit_size(sNaN) - 1
    integer, parameter :: DPSB = bit_size(dNaN) - 1
    
   interface isnan
      module procedure sisnan
      module procedure disnan
   end interface   

   interface isinf
      module procedure sisinf
      module procedure disinf
   end interface   
   
   interface isposinf
      module procedure sisposinf
      module procedure disposinf
   end interface   
   
   interface isneginf
      module procedure sisneginf
      module procedure disneginf
   end interface   
   
contains    

  ! Single precision test for NaN
  elemental function sisnan(x) result(res)
    real(kind(1.e0)), intent(in) :: x
    logical :: res
    res = ieor(ibclr(transfer(x,sNan),SPSB), sNaN) == 0
  end function  

  ! Double precision test for NaN
  elemental function disnan(d) result(res)
    real(kind(1.d0)), intent(in) :: d
    logical :: res
    res = ieor(ibclr(transfer(d,dNaN),DPSB), dNaN) == 0
  end function  
  
  ! Single precision test for Inf
  elemental function sisinf(x) result(res)
    real(kind(1.e0)), intent(in) :: x
    logical :: res
    res = ieor(ibclr(transfer(x,sPosInf),SPSB), sPosInf) == 0
  end function  

  ! Double precision test for Inf
  elemental function disinf(d) result(res)
    real(kind(1.d0)), intent(in) :: d
    logical :: res
    res = ieor(ibclr(transfer(d,dPosInf),DPSB), dPosInf) == 0
  end function  
  
  ! Single precision test for +Inf
  elemental function sisposinf(x) result(res)
    real(kind(1.e0)), intent(in) :: x
    logical :: res
    res = ieor(transfer(x,sPosInf), sPosInf) == 0
  end function  

  ! Double precision test for +Inf
  elemental function disposinf(d) result(res)
    real(kind(1.d0)), intent(in) :: d
    logical :: res
    res = ieor(transfer(d,dPosInf), dPosInf) == 0
  end function  
  
  ! Single precision test for -Inf
  elemental function sisneginf(x) result(res)
    real(kind(1.e0)), intent(in) :: x
    logical :: res
    res = ieor(transfer(x,sNegInf), sNegInf) == 0
  end function  

  ! Double precision test for -Inf
  elemental function disneginf(d) result(res)
    real(kind(1.d0)), intent(in) :: d
    logical :: res
    res = ieor(transfer(d,dNegInf), dNegInf) == 0
  end function  
  
end module


