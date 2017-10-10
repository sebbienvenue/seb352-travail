!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!!
!! This subroutine checks if all relevant files are present
!!
!! called by:
!! - initialization.f90
!!
      subroutine checkfiles()
!!
      implicit none
!!
      logical lexist !
!!
      inquire(file='input.nn',exist=lexist)
      if(.not.lexist) then
        write(*,*)'Error: file input.nn not found'
        stop
      endif
!!
      return
      end
