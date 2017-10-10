!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - prediction.f90
!! - predictionpair.f90
!! - checkonestructure.f90
!!
      subroutine getvolume(lattice,volume)
!!
      use fileunits
!!
      implicit none
!!
      real*8 lattice(3,3)              ! in
      real*8 volume                    ! out
      real*8 tempvec(3)                ! internal
!!
      volume=0.0d0
!!
      tempvec(1)=lattice(1,2)*lattice(2,3)-lattice(1,3)*lattice(2,2)
      tempvec(2)=lattice(1,3)*lattice(2,1)-lattice(1,1)*lattice(2,3)
      tempvec(3)=lattice(1,1)*lattice(2,2)-lattice(1,2)*lattice(2,1)
      volume=tempvec(1)*lattice(3,1)&
            +tempvec(2)*lattice(3,2)&
            +tempvec(3)*lattice(3,3)
      volume=abs(volume)
!!
      if(volume.lt.0.0d0)then
        write(ounit,*)'Error: volume<0 ',volume
        stop
      endif
!!
      return
      end
