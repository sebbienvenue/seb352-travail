!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - getkalmanmatrices.f90
!!
      subroutine initialcorrmatrix(num_weights,corrmatrix)
!!
      implicit none
!!
      integer num_weights                                               ! in
      integer i1,i2                                                     ! internal
      integer icount                                                    ! internal
!!
      real*8 corrmatrix(num_weights*(num_weights+1)/2)  ! out
!!
      icount=1
!!
      do i1=1,num_weights
        do i2=i1,num_weights
          if(i1.eq.i2)then
            corrmatrix(icount)=1.0d0 ! diagonal element
          else
            corrmatrix(icount)=0.0d0 ! off-diagonal element
          endif
          icount=icount+1
        enddo ! i2
      enddo ! i1
!!
      return
      end
