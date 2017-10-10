!######################################################################
! This routine is part of
! RuNNer - Ruhr-University Neural Network Energy Representation
! (c) Dr. Joerg Behler 2008
!######################################################################
!!
!! called by:
!!
      subroutine getrmse_elec(ntrain,ncharges,&
        rmse_charge,rmse_totalcharge,rmse_elec,&
        mad_charge,mad_totalcharge,mad_elec)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer ntrain,ncharges                   ! in
!!
      real*8 rmse_elec                         ! in/out
      real*8 rmse_charge                        ! in/out
      real*8 rmse_totalcharge                   ! in/out
      real*8 mad_elec                          ! in/out
      real*8 mad_charge                         ! in/out
      real*8 mad_totalcharge                    ! in/out
!!
!! electrostatic
      rmse_charge         =rmse_charge/dble(ncharges)
      rmse_charge         =dsqrt(rmse_charge)
      rmse_totalcharge    =rmse_totalcharge/dble(ntrain)
      rmse_totalcharge    =dsqrt(rmse_totalcharge)
      rmse_elec          =rmse_elec/dble(ntrain)
      rmse_elec          =dsqrt(rmse_elec)
      mad_charge          =mad_charge/dble(ncharges)
      mad_totalcharge     =mad_totalcharge/dble(ntrain)
      mad_elec           =mad_elec/dble(ntrain)
!!
      return
      end

