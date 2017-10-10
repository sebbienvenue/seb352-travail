!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!!
      subroutine getkaldims_short(ndim,&
        kaldim,corrdim,corrfdim,&
        maxkaldim,maxcorrdim,maxcorrfdim,&
        num_weights_short_free)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use nnflags
      use globaloptions
!!
      integer num_weights_short_free(ndim) ! in
      integer ndim              ! in
      integer n_start           ! internal, just dummy here
      integer n_end             ! internal, just dummy here
!! Kalman matrix dimensions:
      integer corrdim(ndim)     ! out
      integer corrfdim(ndim)    ! out
      integer kaldim(ndim)      ! out
      integer i1
      integer isum              ! internal
      integer maxkaldim         ! out
      integer maxcorrdim        ! out
      integer maxcorrfdim        ! out
!!
!!
!! initializations
      isum=0
!!
!! non-parallel case
      if((mpisize.eq.1).or.lompmkl)then
        do i1=1,ndim 
          if(lshort)then
            if((optmodee.eq.1).or.(optmodef.eq.1))then ! corrdim and corrfdim must have the same size in case Kalman filter is used only for E or F
              corrdim(i1)=num_weights_short_free(i1)*(num_weights_short_free(i1)+1)/2
              corrfdim(i1)=num_weights_short_free(i1)*(num_weights_short_free(i1)+1)/2
            else
              corrdim(i1)=1
              corrfdim(i1)=1
            endif
            if((optmodee.ne.1).and.(optmodef.ne.1))then
              kaldim(i1)=1
            else
              kaldim(i1)=num_weights_short_free(i1)
            endif
          endif ! lshort
        enddo ! i1
!!
      else ! parallel case
!!
        do i1=1,ndim 
          if(lshort)then
            if(optmodee.eq.1)then
              call mpifitdistribution(num_weights_short_free(i1),kaldim(i1),n_start,n_end)
              corrdim(i1)=kaldim(i1)*num_weights_short_free(i1)
            else
              corrdim(i1)=1
            endif
            if(optmodef.eq.1)then
              call mpifitdistribution(num_weights_short_free(i1),kaldim(i1),n_start,n_end)
              corrfdim(i1)=kaldim(i1)*num_weights_short_free(i1)
            else
              corrfdim(i1)=1
            endif
            if((optmodee.ne.1).and.(optmodef.ne.1))then
              kaldim(i1) =1
            endif
          endif ! lshort
        enddo ! i1
!!
      endif ! mpisize
!!
      maxkaldim=0
      maxcorrdim=0
      maxcorrfdim=0
      do i1=1,ndim
        maxkaldim=max(maxkaldim,kaldim(i1))
        maxcorrdim=max(maxcorrdim,corrdim(i1))
        maxcorrfdim=max(maxcorrfdim,corrfdim(i1))
      enddo
!!
      return
      end      
