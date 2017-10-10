!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! called by:
!! - fittingpair.f90 (4x)
!!
      subroutine writeweightspair(npairs,&
              maxnum_weightspair,num_weightspair,&
              maxnum_layerspair,num_layerspair,&
              nodes_pair,weights_pair,filenamewp)
!!
      use fileunits
!!
      implicit none
!!
      integer npairs                  ! in
      integer i1,i2,i3,i4
      integer icount,jcount
      integer maxnum_weightspair
      integer maxnum_layerspair
      integer num_weightspair(npairs)
      integer num_layerspair(npairs)
      integer nodes_pair(0:maxnum_layerspair,npairs)
!!
      real*8 weights_pair(maxnum_weightspair,npairs)
      character*23 filenamewp(npairs)
!!
!!
      do i1=1,npairs 
        icount=0
        open(woutunit,file=filenamewp(i1),form='formatted',status='replace')
          do i2=1,num_layerspair(i1) ! '
!! connecting weights
            do i3=1,nodes_pair(i2-1,i1)
              do i4=1,nodes_pair(i2,i1)
                icount=icount+1
                write(woutunit,'(f18.10,x,a1,i10,4i6)') weights_pair(icount,i1),'a',icount,i2-1,i3,i2,i4
              enddo ! i4
            enddo ! i3
!! bias weights
            do i3=1,nodes_pair(i2,i1)
              icount=icount+1
              write(woutunit,'(f18.10,x,a1,i10,2i6)')weights_pair(icount,i1),'b',icount,i2,i3
            enddo ! i3
          enddo ! i2
        close(woutunit)
      enddo
!!
!!
      return
      end
