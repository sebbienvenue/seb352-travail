!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - main.f90
!!
      subroutine mode2(iseed)
!!
      use fileunits
      use timings
      use fittingoptions
      use nnflags
      use globaloptions
!!
      implicit none
!!
      integer iseed
!!
      call zerotime(daymode2,timemode2start,timemode2end) 
      call abstime(timemode2start,daymode2)
!!
!!
      if((lshort.and.(lelec.and.nn_type_elec.eq.1)).or.(lshort.and.lnntb).or.(lelec.and.lnntb))then
        write(ounit,*)'ERROR: simultaneous fitting of short and elec or hamiltonian  does not make sense'
        stop !'
      endif
!!
!!======================================================================
!! fitting short range part
!!======================================================================
      if(lshort)then
        if(nn_type_short.eq.1)then
          if(fitmode.eq.1)then
            call fitting_short_atomic(iseed)
!! TODO:
!! fit atomic short range part and possibly charges as second output nodes
!!            call fitting_atomic(iseed)
          elseif(fitmode.eq.2)then
            write(ounit,*)'ERROR: batch fitting is not well tested and not up to date'
            stop
          endif
        elseif(nn_type_short.eq.2)then
          call fitting_short_pair(iseed)
        endif
!!
!!======================================================================
!! fitting electrostatic NN 
!!======================================================================
      elseif(lelec.and.(nn_type_elec.eq.1))then
        call fitting_electrostatic(iseed)
!!
!!======================================================================
!! fitting Hamiltonian NN 
!!======================================================================
      elseif(lnntb)then
        call fitting_nntb(iseed)
      endif
!!
      call abstime(timemode2end,daymode2)
!!
      return
      end
