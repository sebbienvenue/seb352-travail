!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!! - fitting.f90 
!! - fitting_batch.f90 
!! - fittingpair.f90 
!!
      subroutine writeweights(iswitch,ndim,&
        maxnum_weights_local,&
        maxnum_layers_local,num_layers_local,&
        nodes_local,weights_local)
!!
      use fileunits
      use fittingoptions
      use globaloptions
!!
      implicit none
!!
      integer iswitch                                    ! in
      integer ndim                                       ! in
      integer maxnum_weights_local                       ! in
      integer maxnum_layers_local
      integer num_layers_local(ndim)
      integer i1,i2,i3,i4
      integer nodes_local(0:maxnum_layers_local,ndim)
      integer icount
!!
      real*8 weights_local(maxnum_weights_local,ndim)      
!!
!!
!!
      do i1=1,ndim
        icount=0
        if(iswitch.eq.0)then
          if(lwriteunformatted)then
            open(woutunit,file=filenamews(i1),form='unformatted',status='replace')
          else !'
            open(woutunit,file=filenamews(i1),form='formatted',status='replace')
          endif !'
        elseif(iswitch.eq.1)then
          if(lwriteunformatted)then
            open(woutunit,file=filenamewe(i1),form='unformatted',status='replace')
          else !'
            open(woutunit,file=filenamewe(i1),form='formatted',status='replace')
          endif !'
        elseif(iswitch.eq.2)then
          if(lwriteunformatted)then
            open(woutunit,file=filenamewp(i1),form='unformatted',status='replace')
          else !'
            open(woutunit,file=filenamewp(i1),form='formatted',status='replace')
          endif !'
        elseif(iswitch.eq.3)then
          if(lwriteunformatted)then
            open(woutunit,file=filenamewhextoff,form='unformatted',status='replace')
          else !'
            open(woutunit,file=filenamewhextoff,form='formatted',status='replace')
          endif !'
        endif
!! write weights
        do i2=1,num_layers_local(i1) !'
!! connecting weights_local
          do i3=1,nodes_local(i2-1,i1)
            do i4=1,nodes_local(i2,i1)
              icount=icount+1
              if(lwriteunformatted)then
                write(woutunit)weights_local(icount,i1)
              else
                write(woutunit,'(f18.10,x,a1,i10,4i6)')&
                  weights_local(icount,i1),'a',icount,i2-1,i3,i2,i4
              endif
            enddo ! i4
          enddo ! i3
!! bias weights_local
          do i3=1,nodes_local(i2,i1)
            icount=icount+1
            if(lwriteunformatted)then
              write(woutunit)weights_local(icount,i1)
            else
              write(woutunit,'(f18.10,x,a1,i10,2i6)')&
                weights_local(icount,i1),'b',icount,i2,i3
            endif
          enddo ! i3
        enddo ! i2
        close(woutunit)
      enddo
!!
!!
      return
      end
