!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fitting.f90
!!
      subroutine getidx(npoints,iseed,idx)
!!
      use fittingoptions
      use globaloptions
!!
      implicit none
!!
      integer npoints        ! in
      integer idx(nblock)    ! out
      integer iseed          ! in/out
      integer itemp          ! internal     
      integer i1             ! internal
!!
      real*8 ran0            ! internal
      real*8 z(nblock)       ! internal
      real*8 ztemp           ! internal
!!
!! initialization
      do i1=1,nblock
        idx(i1)=i1
      enddo
!!
!! assign random order
      if(lrandomtrain)then
        do i1=1,npoints
          z(i1)=ran0(iseed)
!! debugging:
!!          write(*,'(i4,f14.8)')i1,z(i1)
        enddo
 10     continue
        do i1=1,npoints-1
          if(z(i1).gt.z(i1+1))then
            ztemp=z(i1)
            itemp=idx(i1)
            z(i1)=z(i1+1)
            z(i1+1)=ztemp
            idx(i1)=idx(i1+1)
            idx(i1+1)=itemp
            goto 10
          endif
        enddo
!! debugging
!!        do i1=1,npoints      
!!          write(*,'(i4,i4,f14.8)')i1,idx(i1),z(i1)
!!        enddo
      endif ! lrandomtrain
!!
      return
      end
