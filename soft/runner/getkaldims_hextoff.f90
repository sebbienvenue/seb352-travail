!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################

!! multipurpose subroutine

!! called by:
!!
      subroutine getkaldims_hextoff(ndim,&
        kaldim,corrdim,&
        maxkaldim,maxcorrdim,&
        num_weights_hextoff_free)
!!
      use mpi_mod
      use fileunits
      use fittingoptions
      use nnflags
      use globaloptions
!!
      integer num_weights_hextoff_free(1) ! in
      integer ndim              ! in
      integer n_start           ! internal, just dummy here
      integer n_end             ! internal, just dummy here
!! Kalman matrix dimensions:
      integer corrdim           ! out
      integer kaldim            ! out
      integer i1
      integer isum              ! internal
      integer maxkaldim         ! out
      integer maxcorrdim        ! out
!!
!!
!! initializations
      isum=0
!!
!! non-parallel case
      if((mpisize.eq.1).or.lompmkl)then
            if(optmodehextoff.eq.1)then
              corrdim=num_weights_hextoff_free(1)*(num_weights_hextoff_free(1)+1)/2
            else
              corrdim=1
            endif
            if(optmodehextoff.ne.1)then
              kaldim=1
            else
              kaldim=num_weights_hextoff_free(1)
            endif
!!
      else ! parallel case
!!
            if(optmodehextoff.eq.1)then
              call mpifitdistribution(num_weights_hextoff_free,kaldim,n_start,n_end)
              corrdim=kaldim*num_weights_hextoff_free(1)
            else
              corrdim=1
            endif
            if(optmodehextoff.ne.1)then
              kaldim =1
            else
            endif
!!
      endif ! mpisize
!!
      maxkaldim=0
      maxcorrdim=0
      maxcorrfdim=0
        maxkaldim=max(maxkaldim,kaldim)
        maxcorrdim=max(maxcorrdim,corrdim)
!!
      return
      end      
