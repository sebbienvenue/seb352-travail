!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!! - optimize_ewald.f90
!! - optimize_short_combined.f90
!! - optimize_short_combinedpair.f90
!!
      subroutine debugweights(ndim,iswitch,&
        countepoch,point,ielem,&
        maxnum_weights_local,&
        maxnum_layers_local,num_layers_local,&
        nodes_local,weights_local)
!!
      use fileunits
!!
      implicit none
!!
      integer ndim                   ! in
      integer maxnum_weights_local          ! in
      integer maxnum_layers_local
      integer num_layers_local(ndim)
      integer i2,i3,i4
      integer nodes_local(0:maxnum_layers_local,ndim)
      integer icount
      integer iswitch
      integer point
      integer countepoch
      integer ielem
!!
      real*8 weights_local(maxnum_weights_local,ndim)      
!!
      character*15 dummy 
!!
      if(iswitch.eq.0)then
        dummy=' weights_short '
      elseif(iswitch.eq.1)then
        dummy=' weights_ewald '
      elseif(iswitch.eq.2)then
        dummy=' weights_pair '
      else
        write(ounit,*)'ERROR: unknown iswitch in debugweights ',iswitch
        stop
      endif
!!
        icount=0
          do i2=1,num_layers_local(ielem)
!! connecting weights
            do i3=1,nodes_local(i2-1,ielem)
              do i4=1,nodes_local(i2,ielem)
                icount=icount+1
                write(debugunit,'(3i6,a,f18.10,x,a1,i10,4i6)')&
                  countepoch,point,ielem,dummy,weights_local(icount,ielem),'a',icount,i2-1,i3,i2,i4
              enddo ! i4
            enddo ! i3
!! bias weights
            do i3=1,nodes_local(i2,ielem)
              icount=icount+1
              write(debugunit,'(3i6,a,f18.10,x,a1,i10,2i6)')&
                countepoch,point,ielem,dummy,weights_local(icount,ielem),'b',icount,i2,i3
            enddo ! i3
          enddo ! i2
!!
      return
      end
