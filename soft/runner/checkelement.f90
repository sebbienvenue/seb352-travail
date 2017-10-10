!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - readinput.f90
!!
      subroutine checkelement(elementtemp)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1            ! internal
!!
      character*2 elementtemp    ! in
!!
      logical lfound             ! internal
!!
      lfound=.false.
!!
      do i1=1,nelem
        if(elementtemp.eq.element(i1))then
          lfound=.true.
        endif
      enddo
!!
      if(.not.lfound) then
        write(ounit,*)'Error: illegal element specified in some keyword in input.nn: ',elementtemp
        stop !'
      endif
!!
      return
      end
