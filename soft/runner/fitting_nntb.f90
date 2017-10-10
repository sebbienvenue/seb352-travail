!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - mode2.f90
!!
      subroutine fitting_nntb(iseed)
!!
      use fileunits
      use nnflags 
      use globaloptions
      use timings
!!
      implicit none
!!
      integer iseed                ! in
!!
!!=====================================================
!! fit the overlap-dependent part of the Hamiltonian (analogous to short range pair case)
!!=====================================================
      if(nntb_flag(1))then
        call fit_density(iseed)

!!=====================================================
!! fit the hexton potential part of the Hamiltonian 
!!=====================================================
      elseif(nntb_flag(2))then
        call fit_hexton(iseed)

!!=====================================================
!! fit the hextoff potential part of the Hamiltonian
!!=====================================================
      elseif(nntb_flag(3))then
        call fit_hextoff(iseed)
!!=====================================================
!! fit the Overlap matrix 
!!=====================================================
      elseif(nntb_flag(4))then
        call fit_overlap(iseed)



      else
        write(ounit,*)'ERROR: unexpected nn_type_nntb in fitting_nntb'
        stop
      endif
!!
      return
      end
