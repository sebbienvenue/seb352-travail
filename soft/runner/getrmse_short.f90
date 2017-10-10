!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################
!!
!! called by:
!!
      subroutine getrmse_short(nenergies,&
        nforces,rmse_short,rmse_force_s,mad_short,mad_force_s)
!!
      use fileunits
      use nnflags
      use globaloptions
!!
      implicit none
!!
      integer nenergies,nforces                 ! in
!!
      real*8 rmse_short                         ! in/out
      real*8 rmse_force_s                       ! in/out
      real*8 mad_short                          ! in/out
      real*8 mad_force_s                        ! in/out
!!
!!
!! some security checks if settings for maxenergy or maxforce are inappropriate:
      if(lshort)then
        if(nenergies.eq.0)then
          write(ounit,*)'ERROR in getrmse: nenergies is zero '
          write(ounit,*)'Probably maxenergy and/or energy_threshold is set too low'
          stop  !'
        endif
      endif
!!
!! short range
      if(lshort)then
        rmse_short          =rmse_short/dble(nenergies)
        rmse_short          =dsqrt(rmse_short)
        rmse_force_s        =rmse_force_s/dble(nforces)
        rmse_force_s        =sqrt(rmse_force_s)
        mad_short           =mad_short/dble(nenergies)
        mad_force_s         =mad_force_s/dble(nforces)
      endif
!!
      return
      end

