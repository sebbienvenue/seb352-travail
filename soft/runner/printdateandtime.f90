!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!!
      subroutine printdateandtime(countepoch)
!!
      use fileunits
      use timings
!!
      implicit none
!!
      integer countepoch ! in
!!
      call date_and_time(fulldate,fulltime,zone,timevalues)

      write(ounit,'(a8,i5,6x,i4,x,i2,x,i2,a7,i2,a3,i2,a4,i2,a3)')' DATE   ',&
        countepoch,&
        timevalues(1),timevalues(2),timevalues(3),&
        ' TIME: ',timevalues(5),' h ',&
        timevalues(6),' m ',timevalues(7),' s '

      end

