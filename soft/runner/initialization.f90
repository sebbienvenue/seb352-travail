!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - main.f90
!!
      subroutine initialization(ielem,lelement)
!!
      use fileunits
      use nnflags
      use globaloptions
!!
      implicit none
!!
      integer ielem                           ! out number of elements in input.data
!!
      logical lelement(102)                   ! out
!!
!! initializations
      nelem         = 0
      ielem         = 0
      max_num_pairs = 0
!!
!! CHANGE ANDI: GFORTRAN: initialize "lelement" already here, not inside of checkstructures(...).
!!              Otherwise readinput(...) uses uninitialized array in mode 2.
!!
      lelement(:)=.false.
!! END CHANGE
!!
      call checkfiles()
!!
      call getdimensions()
!!
!!     get dimensions for structure-related arrays
      call structurecount()
!!
!!     get dimensions for pair-related arrays (and within it will call to get triplets)
      call paircount()
!!
!!
!! check structures in input.data for inconsistencies in modes 1 and 3
      if((mode.eq.1).or.(mode.eq.3))then
        call checkstructures(ielem,lelement)
      endif
!!
      return
      end
