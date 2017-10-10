!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

      module symfuncshort 

      implicit none

! this does not work because num_functions needs to be known for declaring the arrays
      integer num_functions 
      integer function_type(num_functions)

      real*8 funccutoff(num_functions)
      real*8 lambda(num_functions)
      real*8 zeta(num_functions)
      real*8 eta(num_functions)
      real*8 rshift(num_functions)
      real*8 maxcutoff

      end module symfuncshort 

