!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################
!!
!! called by:
!!
      subroutine getrmse_hextoff(matrixsize,ndone,&
        rmse_hextoff,mad_hextoff)
!!
      use fileunits
      use nnflags
      use globaloptions
!!
      implicit none
!!
      integer matrixsize                        ! in
      integer ndone
!!
      real*8 rmse_hextoff                       ! in/out
      real*8 mad_hextoff                        ! in/out
!!
!!
!! some security checks if settings for maxenergy or maxforce are inappropriate:
!!      if(lshort)then
!!        if(nenergies.eq.0)then
!!          write(ounit,*)'ERROR in getrmse: nenergies is zero '
!!          write(ounit,*)'Probably maxenergy and/or energy_threshold is set too low'
!!          stop  !'
!!        endif
!!      endif
!!
!! hextoff 
      if(lnntb.and.nntb_flag(3))then
        rmse_hextoff        =rmse_hextoff/dble(ndone*matrixsize)
        rmse_hextoff        =dsqrt(rmse_hextoff)
        mad_hextoff         =mad_hextoff/dble(ndone*matrixsize)
      endif
!!
      return
      end

